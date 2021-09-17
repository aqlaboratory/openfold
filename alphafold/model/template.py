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
import torch
import torch.nn as nn

from alphafold.model.primitives import Linear, Attention
from alphafold.utils.deepspeed import checkpoint_blocks
from alphafold.model.dropout import  (
    DropoutRowwise,
    DropoutColumnwise,
)
from alphafold.model.pair_transition import PairTransition
from alphafold.model.triangular_attention import (
    TriangleAttentionStartingNode,
    TriangleAttentionEndingNode,
)
from alphafold.model.triangular_multiplicative_update import (
    TriangleMultiplicationOutgoing,
    TriangleMultiplicationIncoming,
)
from alphafold.utils.tensor_utils import (
    chunk_layer,
    permute_final_dims, 
    flatten_final_dims,
)


class TemplatePointwiseAttention(nn.Module):
    """
        Implements Algorithm 17.
    """
    def __init__(self, 
        c_t, 
        c_z, 
        c_hidden, 
        no_heads, 
        chunk_size,
        **kwargs
    ):
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
        self.chunk_size = chunk_size

        self.mha = Attention(
            self.c_z, self.c_t, self.c_t, 
            self.c_hidden, self.no_heads,
            gating=False,
        )

    def forward(self, t, z, template_mask=None):
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
        if(template_mask is None):
            # NOTE: This is not the "template_mask" from the supplement, but a
            # [*, N_templ] mask from the code. I'm pretty sure it's always just 1,
            # but not sure enough to remove it. It's nice to have, I guess.
            template_mask = torch.ones(t.shape[:-3], device=t.device)
        
        bias = (1e9 * (template_mask[..., None, None, None, None, :] - 1))

        # [*, N_res, N_res, 1, C_z]
        z = z.unsqueeze(-2)

        # [*, N_res, N_res, N_temp, C_t]
        t = permute_final_dims(t, 1, 2, 0, 3)

        # [*, N_res, N_res, 1, C_z]
        mha_inputs = {
            "q_x": z,
            "k_x": t,
            "v_x": t,
            "biases": [bias],
        }
        if(not self.training and self.chunk_size is not None):
            z = chunk_layer(
                self.mha,
                mha_inputs,
                chunk_size=self.chunk_size,
                no_batch_dims=len(z.shape[:-2])
            )
        else:
            z = self.mha(**mha_inputs)

        # [*, N_res, N_res, C_z]
        z = z.squeeze(-2)
         
        return z


class TemplatePairStackBlock(nn.Module):
    def __init__(self, 
        c_t, 
        c_hidden_tri_att,
        c_hidden_tri_mul,
        no_heads, 
        pair_transition_n, 
        dropout_rate,
        chunk_size,
    ):
        super(TemplatePairStackBlock, self).__init__()
        
        self.c_t = c_t
        self.c_hidden_tri_att = c_hidden_tri_att
        self.c_hidden_tri_mul = c_hidden_tri_mul
        self.no_heads = no_heads
        self.pair_transition_n = pair_transition_n
        self.dropout_rate = dropout_rate
        self.chunk_size = chunk_size

        self.dropout_row = DropoutRowwise(self.dropout_rate)
        self.dropout_col = DropoutColumnwise(self.dropout_rate)
        
        self.tri_att_start = TriangleAttentionStartingNode(
            self.c_t, 
            self.c_hidden_tri_att, 
            self.no_heads, 
            chunk_size=chunk_size,
        )
        self.tri_att_end = TriangleAttentionEndingNode(
            self.c_t,
            self.c_hidden_tri_att,
            self.no_heads,
            chunk_size=chunk_size,
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
            chunk_size=chunk_size,
        )

    def forward(self, z, mask, _mask_trans=True):
        z = z + self.dropout_row(self.tri_att_start(z, mask=mask))
        z = z + self.dropout_col(self.tri_att_end(z, mask=mask))
        z = z + self.dropout_row(self.tri_mul_out(z, mask=mask))
        z = z + self.dropout_row(self.tri_mul_in(z, mask=mask))
        z = z + self.pair_transition(z, mask=mask if _mask_trans else None)
        
        return z



class TemplatePairStack(nn.Module):
    """
        Implements Algorithm 16.
    """
    def __init__(self, 
        c_t, 
        c_hidden_tri_att,
        c_hidden_tri_mul,
        no_blocks, 
        no_heads, 
        pair_transition_n, 
        dropout_rate,
        blocks_per_ckpt,
        chunk_size,
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
                chunk_size:
                    Size of subbatches. A higher value increases throughput at
                    the cost of memory
        """  
        super(TemplatePairStack, self).__init__()

        self.blocks_per_ckpt = blocks_per_ckpt

        self.blocks = nn.ModuleList()
        for i in range(no_blocks):
            block = TemplatePairStackBlock(
                c_t=c_t,
                c_hidden_tri_att=c_hidden_tri_att,
                c_hidden_tri_mul=c_hidden_tri_mul,
                no_heads=no_heads,
                pair_transition_n=pair_transition_n,
                dropout_rate=dropout_rate,
                chunk_size=chunk_size,
            )
            self.blocks.append(block)

        self.layer_norm = nn.LayerNorm(c_t)

    def forward(self, 
        t: torch.tensor,
        mask: torch.tensor,
        _mask_trans: bool = True,
    ):
        """
            Args:
                t:
                    [*, N_res, N_res, C_t] template embedding
                mask:
                    [*, N_res, N_res] mask
            Returns:
                [*, N_res, N_res, C_t] template embedding update
        """
        t, = checkpoint_blocks(
            blocks=[
                partial(
                    b, 
                    mask=mask, 
                    _mask_trans=_mask_trans,
                ) for b in self.blocks
            ], 
            args=(t),
            blocks_per_ckpt=self.blocks_per_ckpt,
        )

        t = self.layer_norm(t)

        return t


if __name__ == "__main__":
    template_angle_dim = 51
    c_m = 256
    batch_size = 4
    n_templ = 4
    n_res = 256

    tae = TemplateAngleEmbedder(
        template_angle_dim,
        c_m,
    )
 
    x = torch.rand((batch_size, n_templ, n_res, template_angle_dim))
    x.shape_before = x.shape
    x = tae(x)
    x.shape_after = x.shape

    assert(shape_before == shape_after)

    batch_size = 2
    s_t = 4
    c_t = 64
    c_z = 128
    c = 32
    no_heads = 3
    n = 100
    
    tpa = TemplatePointwiseAttention(c_t, c_z, c)

    t = torch.rand((batch_size, s_t, n, n, c_t))
    z = torch.rand((batch_size, n, n, c_z))

    z_update = tpa(t, z)

    assert(z_update.shape == z.shape)
