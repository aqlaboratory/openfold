# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from functools import partial
import math
from typing import Optional, List

import torch
import torch.nn as nn

from openfold.model.primitives import Linear, LayerNorm, Attention
from openfold.model.dropout import (
    DropoutRowwise,
    DropoutColumnwise,
)
from openfold.model.pair_transition import PairTransition
from openfold.model.triangular_attention import (
    TriangleAttentionStartingNode,
    TriangleAttentionEndingNode,
)
from openfold.model.triangular_multiplicative_update import (
    TriangleMultiplicationOutgoing,
    TriangleMultiplicationIncoming,
)
from openfold.utils.checkpointing import checkpoint_blocks
from openfold.utils.chunk_utils import (
    chunk_layer,
    ChunkSizeTuner,
)
from openfold.utils.feats import (
    build_template_angle_feat,
    build_template_pair_feat,
)
from openfold.utils.tensor_utils import (
    add,
    permute_final_dims,
    flatten_final_dims,
    tensor_tree_map,
)


class TemplatePointwiseAttention(nn.Module):
    """
    Implements Algorithm 17.
    """
    def __init__(self, c_t, c_z, c_hidden, no_heads, inf, **kwargs):
        """
        Args:
            c_t:
                Template embedding channel dimension
            c_z:
                Pair embedding channel dimension
            c_hidden:
                Hidden channel dimension
        """
        super(TemplatePointwiseAttention, self).__init__()

        self.c_t = c_t
        self.c_z = c_z
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.inf = inf

        self.mha = Attention(
            self.c_z,
            self.c_t,
            self.c_t,
            self.c_hidden,
            self.no_heads,
            gating=False,
        )

    def _chunk(self,
        z: torch.Tensor,
        t: torch.Tensor,
        biases: List[torch.Tensor],
        chunk_size: int,
        use_lma: bool = False,
    ) -> torch.Tensor:
        mha_inputs = {
            "q_x": z,
            "kv_x": t,
            "biases": biases,
        }
        return chunk_layer(
            partial(self.mha, use_lma=use_lma),
            mha_inputs,
            chunk_size=chunk_size,
            no_batch_dims=len(z.shape[:-2]),
        )


    def forward(self, 
        t: torch.Tensor, 
        z: torch.Tensor, 
        template_mask: Optional[torch.Tensor] = None,
        # This module suffers greatly from a small chunk size
        chunk_size: Optional[int] = 256,
        use_lma: bool = False,
    ) -> torch.Tensor:
        """
        Args:
            t:
                [*, N_templ, N_res, N_res, C_t] template embedding
            z:
                [*, N_res, N_res, C_t] pair embedding
            template_mask:
                [*, N_templ] template mask
        Returns:
            [*, N_res, N_res, C_z] pair embedding update
        """
        if template_mask is None:
            template_mask = t.new_ones(t.shape[:-3])

        bias = self.inf * (template_mask[..., None, None, None, None, :] - 1)

        # [*, N_res, N_res, 1, C_z]
        z = z.unsqueeze(-2)

        # [*, N_res, N_res, N_temp, C_t]
        t = permute_final_dims(t, (1, 2, 0, 3))

        # [*, N_res, N_res, 1, C_z]
        biases = [bias]
        if chunk_size is not None and not self.training:
            z = self._chunk(z, t, biases, chunk_size, use_lma=use_lma)
        else:
            z = self.mha(q_x=z, kv_x=t, biases=biases, use_lma=use_lma)

        # [*, N_res, N_res, C_z]
        z = z.squeeze(-2)

        return z


class TemplatePairStackBlock(nn.Module):
    def __init__(
        self,
        c_t: int,
        c_hidden_tri_att: int,
        c_hidden_tri_mul: int,
        no_heads: int,
        pair_transition_n: int,
        dropout_rate: float,
        inf: float,
        **kwargs,
    ):
        super(TemplatePairStackBlock, self).__init__()

        self.c_t = c_t
        self.c_hidden_tri_att = c_hidden_tri_att
        self.c_hidden_tri_mul = c_hidden_tri_mul
        self.no_heads = no_heads
        self.pair_transition_n = pair_transition_n
        self.dropout_rate = dropout_rate
        self.inf = inf

        self.dropout_row = DropoutRowwise(self.dropout_rate)
        self.dropout_col = DropoutColumnwise(self.dropout_rate)

        self.tri_att_start = TriangleAttentionStartingNode(
            self.c_t,
            self.c_hidden_tri_att,
            self.no_heads,
            inf=inf,
        )
        self.tri_att_end = TriangleAttentionEndingNode(
            self.c_t,
            self.c_hidden_tri_att,
            self.no_heads,
            inf=inf,
        )

        self.tri_mul_out = TriangleMultiplicationOutgoing(
            self.c_t,
            self.c_hidden_tri_mul,
        )
        self.tri_mul_in = TriangleMultiplicationIncoming(
            self.c_t,
            self.c_hidden_tri_mul,
        )

        self.pair_transition = PairTransition(
            self.c_t,
            self.pair_transition_n,
        )

    def forward(self, 
        z: torch.Tensor, 
        mask: torch.Tensor, 
        chunk_size: Optional[int] = None, 
        use_lma: bool = False,
        inplace_safe: bool = False,
        _mask_trans: bool = True,
        _attn_chunk_size: Optional[int] = None,
    ):
        if(_attn_chunk_size is None):
            _attn_chunk_size = chunk_size

        single_templates = [
            t.unsqueeze(-4) for t in torch.unbind(z, dim=-4)
        ]
        single_templates_masks = [
            m.unsqueeze(-3) for m in torch.unbind(mask, dim=-3)
        ]

        for i in range(len(single_templates)):
            single = single_templates[i]
            single_mask = single_templates_masks[i]
            
            single = add(single,
                self.dropout_row(
                    self.tri_att_start(
                        single,
                        chunk_size=_attn_chunk_size,
                        mask=single_mask,
                        use_lma=use_lma,
                        inplace_safe=inplace_safe,
                    )
                ),
                inplace_safe,
            )

            single = add(single,
                self.dropout_col(
                    self.tri_att_end(
                        single,
                        chunk_size=_attn_chunk_size,
                        mask=single_mask,
                        use_lma=use_lma,
                        inplace_safe=inplace_safe,
                    )
                ),
                inplace_safe,
            )

            tmu_update = self.tri_mul_out(
                single,
                mask=single_mask,
                inplace_safe=inplace_safe,
                _add_with_inplace=True,
            )
            if(not inplace_safe):
                single = single + self.dropout_row(tmu_update)
            else:
                single = tmu_update
            
            del tmu_update

            tmu_update = self.tri_mul_in(
                single,
                mask=single_mask,
                inplace_safe=inplace_safe,
                _add_with_inplace=True,
            )
            if(not inplace_safe):
                single = single + self.dropout_row(tmu_update)
            else:
                single = tmu_update
            
            del tmu_update
      
            single = add(single,
                self.pair_transition(
                    single,
                    mask=single_mask if _mask_trans else None,
                    chunk_size=chunk_size,
                ),
                inplace_safe,
            )

            if(not inplace_safe):
                single_templates[i] = single

        if(not inplace_safe):
            z = torch.cat(single_templates, dim=-4)

        return z


class TemplatePairStack(nn.Module):
    """
    Implements Algorithm 16.
    """
    def __init__(
        self,
        c_t,
        c_hidden_tri_att,
        c_hidden_tri_mul,
        no_blocks,
        no_heads,
        pair_transition_n,
        dropout_rate,
        blocks_per_ckpt,
        tune_chunk_size: bool = False,
        inf=1e9,
        **kwargs,
    ):
        """
        Args:
            c_t:
                Template embedding channel dimension
            c_hidden_tri_att:
                Per-head hidden dimension for triangular attention
            c_hidden_tri_att:
                Hidden dimension for triangular multiplication
            no_blocks:
                Number of blocks in the stack
            pair_transition_n:
                Scale of pair transition (Alg. 15) hidden dimension
            dropout_rate:
                Dropout rate used throughout the stack
            blocks_per_ckpt:
                Number of blocks per activation checkpoint. None disables
                activation checkpointing
        """
        super(TemplatePairStack, self).__init__()

        self.blocks_per_ckpt = blocks_per_ckpt

        self.blocks = nn.ModuleList()
        for _ in range(no_blocks):
            block = TemplatePairStackBlock(
                c_t=c_t,
                c_hidden_tri_att=c_hidden_tri_att,
                c_hidden_tri_mul=c_hidden_tri_mul,
                no_heads=no_heads,
                pair_transition_n=pair_transition_n,
                dropout_rate=dropout_rate,
                inf=inf,
            )
            self.blocks.append(block)

        self.layer_norm = LayerNorm(c_t)

        self.tune_chunk_size = tune_chunk_size
        self.chunk_size_tuner = None
        if(tune_chunk_size):
            self.chunk_size_tuner = ChunkSizeTuner()

    def forward(
        self,
        t: torch.tensor,
        mask: torch.tensor,
        chunk_size: int,
        use_lma: bool = False,
        inplace_safe: bool = False,
        _mask_trans: bool = True,
    ):
        """
        Args:
            t:
                [*, N_templ, N_res, N_res, C_t] template embedding
            mask:
                [*, N_templ, N_res, N_res] mask
        Returns:
            [*, N_templ, N_res, N_res, C_t] template embedding update
        """
        if(mask.shape[-3] == 1):
            expand_idx = list(mask.shape)
            expand_idx[-3] = t.shape[-4]
            mask = mask.expand(*expand_idx)

        blocks = [
            partial(
                b,
                mask=mask,
                chunk_size=chunk_size,
                use_lma=use_lma,
                inplace_safe=inplace_safe,
                _mask_trans=_mask_trans,
            )
            for b in self.blocks
        ]

        if(chunk_size is not None and self.chunk_size_tuner is not None):
            assert(not self.training)
            tuned_chunk_size = self.chunk_size_tuner.tune_chunk_size(
                representative_fn=blocks[0],
                args=(t.clone(),),
                min_chunk_size=chunk_size,
            )
            blocks = [
                partial(b, 
                    chunk_size=tuned_chunk_size,
                    _attn_chunk_size=max(chunk_size, tuned_chunk_size // 4),
                ) for b in blocks
            ]

        t, = checkpoint_blocks(
            blocks=blocks,
            args=(t,),
            blocks_per_ckpt=self.blocks_per_ckpt if self.training else None,
        )

        t = self.layer_norm(t)

        return t


def embed_templates_offload(
    model, 
    batch, 
    z, 
    pair_mask, 
    templ_dim, 
    template_chunk_size=256,
    inplace_safe=False,
):
    """
    Args:
        model: 
            An AlphaFold model object
        batch: 
            An AlphaFold input batch. See documentation of AlphaFold.
        z: 
            A [*, N, N, C_z] pair embedding
        pair_mask: 
            A [*, N, N] pair mask
        templ_dim: 
            The template dimension of the template tensors in batch
        template_chunk_size: 
            Integer value controlling how quickly the offloaded pair embedding
            tensor is brought back into GPU memory. In dire straits, can be
            lowered to reduce memory consumption of this function even more.
    Returns:
        A dictionary of template pair and angle embeddings.
    
    A version of the "embed_templates" method of the AlphaFold class that
    offloads the large template pair tensor to CPU. Slower but more frugal 
    with GPU memory than the original. Useful for long-sequence inference.
    """
    # Embed the templates one at a time (with a poor man's vmap)
    pair_embeds_cpu = []
    n = z.shape[-2]
    n_templ = batch["template_aatype"].shape[templ_dim]
    for i in range(n_templ):
        idx = batch["template_aatype"].new_tensor(i)
        single_template_feats = tensor_tree_map(
            lambda t: torch.index_select(t, templ_dim, idx).squeeze(templ_dim),
            batch,
        )

        # [*, N, N, C_t]
        t = build_template_pair_feat(
            single_template_feats,
            use_unit_vector=model.config.template.use_unit_vector,
            inf=model.config.template.inf,
            eps=model.config.template.eps,
            **model.config.template.distogram,
        ).to(z.dtype)
        t = model.template_pair_embedder(t)

        # [*, 1, N, N, C_z]
        t = model.template_pair_stack(
            t,
            pair_mask.unsqueeze(-3).to(dtype=z.dtype), 
            chunk_size=model.globals.chunk_size,
            use_lma=model.globals.use_lma,
            _mask_trans=model.config._mask_trans,
        )

        pair_embeds_cpu.append(t.cpu())

        del t

    # Preallocate the output tensor
    t = z.new_zeros(z.shape)

    for i in range(0, n, template_chunk_size):
        pair_chunks = [
            p[..., i: i + template_chunk_size, :, :] for p in pair_embeds_cpu
        ]
        pair_chunk = torch.cat(pair_chunks, dim=templ_dim).to(device=z.device)
        z_chunk = z[..., i: i + template_chunk_size, :, :]
        att_chunk = model.template_pointwise_att(
            pair_chunk,
            z_chunk,
            template_mask=batch["template_mask"].to(dtype=z.dtype),
            use_lma=model.globals.use_lma,
        )

        t[..., i: i + template_chunk_size, :, :] = att_chunk
    
    del pair_chunks

    if(inplace_safe):
        t = t * (torch.sum(batch["template_mask"], dim=-1) > 0)
    else:
        t *= (torch.sum(batch["template_mask"], dim=-1) > 0)

    ret = {}
    if model.config.template.embed_angles:
        template_angle_feat = build_template_angle_feat(
            batch,
        )

        # [*, N, C_m]
        a = model.template_angle_embedder(template_angle_feat)
 
        ret["template_angle_embedding"] = a 

    ret.update({"template_pair_embedding": t})

    return ret


def embed_templates_average(
    model, 
    batch, 
    z, 
    pair_mask, 
    templ_dim,
    templ_group_size=2,
    inplace_safe=False,
):
    """
    Args:
        model: 
            An AlphaFold model object
        batch: 
            An AlphaFold input batch. See documentation of AlphaFold.
        z: 
            A [*, N, N, C_z] pair embedding
        pair_mask: 
            A [*, N, N] pair mask
        templ_dim: 
            The template dimension of the template tensors in batch
        templ_group_size: 
            Granularity of the approximation. Larger values trade memory for 
            greater proximity to the original function
    Returns:
        A dictionary of template pair and angle embeddings.

    A memory-efficient approximation of the "embed_templates" method of the 
    AlphaFold class. Instead of running pointwise attention over pair 
    embeddings for all of the templates at the same time, it splits templates 
    into groups of size templ_group_size, computes embeddings for each group 
    normally, and then averages the group embeddings. In our experiments, this 
    approximation has a minimal effect on the quality of the resulting 
    embedding, while its low memory footprint allows the number of templates 
    to scale almost indefinitely.
    """
    # Embed the templates one at a time (with a poor man's vmap)
    n = z.shape[-2]
    n_templ = batch["template_aatype"].shape[templ_dim]
    out_tensor = z.new_zeros(z.shape)
    for i in range(0, n_templ, templ_group_size): 
        def slice_template_tensor(t):
            s = [slice(None) for _ in t.shape]
            s[templ_dim] = slice(i, i + templ_group_size)
            return t[s]
        
        template_feats = tensor_tree_map(
            slice_template_tensor,
            batch,
        )

        # [*, N, N, C_t]
        t = build_template_pair_feat(
            template_feats,
            use_unit_vector=model.config.template.use_unit_vector,
            inf=model.config.template.inf,
            eps=model.config.template.eps,
            **model.config.template.distogram,
        ).to(z.dtype)

        # [*, S_t, N, N, C_z]
        t = model.template_pair_embedder(t)
        t = model.template_pair_stack(
            t, 
            pair_mask.unsqueeze(-3).to(dtype=z.dtype), 
            chunk_size=model.globals.chunk_size,
            use_lma=model.globals.use_lma,
            _mask_trans=model.config._mask_trans,
        )

        t = model.template_pointwise_att(
            t,
            z,
            template_mask=template_feats["template_mask"].to(dtype=z.dtype),
            use_lma=model.globals.use_lma,
        )

        denom = math.ceil(n_templ / templ_group_size)
        if(inplace_safe):
            t /= denom
        else:
            t = t / denom

        if(inplace_safe):
            out_tensor += t
        else:
            out_tensor = out_tensor + t

        del t

    if(inplace_safe):
        out_tensor *= (torch.sum(batch["template_mask"], dim=-1) > 0)
    else:
        out_tensor = out_tensor * (torch.sum(batch["template_mask"], dim=-1) > 0)

    ret = {}
    if model.config.template.embed_angles:
        template_angle_feat = build_template_angle_feat(
            batch,
        )

        # [*, N, C_m]
        a = model.template_angle_embedder(template_angle_feat)
 
        ret["template_angle_embedding"] = a 

    ret.update({"template_pair_embedding": out_tensor})

    return ret
