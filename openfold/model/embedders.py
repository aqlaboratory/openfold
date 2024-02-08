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

import torch
import torch.nn as nn
from typing import Tuple, Optional

from openfold.utils import all_atom_multimer
from openfold.utils.feats import (
    pseudo_beta_fn,
    dgram_from_positions,
    build_template_angle_feat,
    build_template_pair_feat,
)
from openfold.model.primitives import Linear, LayerNorm
from openfold.model.template import (
    TemplatePairStack,
    TemplatePointwiseAttention,
)
from openfold.utils import geometry
from openfold.utils.tensor_utils import add, one_hot, tensor_tree_map, dict_multimap


class InputEmbedder(nn.Module):
    """
    Embeds a subset of the input features.

    Implements Algorithms 3 (InputEmbedder) and 4 (relpos).
    """

    def __init__(
        self,
        tf_dim: int,
        msa_dim: int,
        c_z: int,
        c_m: int,
        relpos_k: int,
        **kwargs,
    ):
        """
        Args:
            tf_dim:
                Final dimension of the target features
            msa_dim:
                Final dimension of the MSA features
            c_z:
                Pair embedding dimension
            c_m:
                MSA embedding dimension
            relpos_k:
                Window size used in relative positional encoding
        """
        super(InputEmbedder, self).__init__()

        self.tf_dim = tf_dim
        self.msa_dim = msa_dim

        self.c_z = c_z
        self.c_m = c_m

        self.linear_tf_z_i = Linear(tf_dim, c_z)
        self.linear_tf_z_j = Linear(tf_dim, c_z)
        self.linear_tf_m = Linear(tf_dim, c_m)
        self.linear_msa_m = Linear(msa_dim, c_m)

        # RPE stuff
        self.relpos_k = relpos_k
        self.no_bins = 2 * relpos_k + 1
        self.linear_relpos = Linear(self.no_bins, c_z)

    def relpos(self, ri: torch.Tensor):
        """
        Computes relative positional encodings

        Implements Algorithm 4.

        Args:
            ri:
                "residue_index" features of shape [*, N]
        """
        d = ri[..., None] - ri[..., None, :]
        boundaries = torch.arange(
            start=-self.relpos_k, end=self.relpos_k + 1, device=d.device
        ) 
        reshaped_bins = boundaries.view(((1,) * len(d.shape)) + (len(boundaries),))
        d = d[..., None] - reshaped_bins
        d = torch.abs(d)
        d = torch.argmin(d, dim=-1)
        d = nn.functional.one_hot(d, num_classes=len(boundaries)).float()
        d = d.to(ri.dtype)
        return self.linear_relpos(d)

    def forward(
        self,
        tf: torch.Tensor,
        ri: torch.Tensor,
        msa: torch.Tensor,
        inplace_safe: bool = False,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            batch: Dict containing
                "target_feat":
                    Features of shape [*, N_res, tf_dim]
                "residue_index":
                    Features of shape [*, N_res]
                "msa_feat":
                    Features of shape [*, N_clust, N_res, msa_dim]
        Returns:
            msa_emb:
                [*, N_clust, N_res, C_m] MSA embedding
            pair_emb:
                [*, N_res, N_res, C_z] pair embedding

        """
        # [*, N_res, c_z]
        tf_emb_i = self.linear_tf_z_i(tf)
        tf_emb_j = self.linear_tf_z_j(tf)

        # [*, N_res, N_res, c_z]
        pair_emb = self.relpos(ri.type(tf_emb_i.dtype))
        pair_emb = add(pair_emb, 
            tf_emb_i[..., None, :], 
            inplace=inplace_safe
        )
        pair_emb = add(pair_emb, 
            tf_emb_j[..., None, :, :], 
            inplace=inplace_safe
        )

        # [*, N_clust, N_res, c_m]
        n_clust = msa.shape[-3]
        tf_m = (
            self.linear_tf_m(tf)
            .unsqueeze(-3)
            .expand(((-1,) * len(tf.shape[:-2]) + (n_clust, -1, -1)))
        )
        msa_emb = self.linear_msa_m(msa) + tf_m

        return msa_emb, pair_emb


class InputEmbedderMultimer(nn.Module):
    """
    Embeds a subset of the input features.

    Implements Algorithms 3 (InputEmbedder) and 4 (relpos).
    """

    def __init__(
        self,
        tf_dim: int,
        msa_dim: int,
        c_z: int,
        c_m: int,
        max_relative_idx: int,
        use_chain_relative: bool,
        max_relative_chain: int,
        **kwargs,
    ):
        """
        Args:
            tf_dim:
                Final dimension of the target features
            msa_dim:
                Final dimension of the MSA features
            c_z:
                Pair embedding dimension
            c_m:
                MSA embedding dimension
            relpos_k:
                Window size used in relative positional encoding
        """
        super(InputEmbedderMultimer, self).__init__()

        self.tf_dim = tf_dim
        self.msa_dim = msa_dim

        self.c_z = c_z
        self.c_m = c_m

        self.linear_tf_z_i = Linear(tf_dim, c_z)
        self.linear_tf_z_j = Linear(tf_dim, c_z)
        self.linear_tf_m = Linear(tf_dim, c_m)
        self.linear_msa_m = Linear(msa_dim, c_m)

        # RPE stuff
        self.max_relative_idx = max_relative_idx
        self.use_chain_relative = use_chain_relative
        self.max_relative_chain = max_relative_chain
        if(self.use_chain_relative):
            self.no_bins = (
                2 * max_relative_idx + 2 +
                1 +
                2 * max_relative_chain + 2
            )
        else:
            self.no_bins = 2 * max_relative_idx + 1
        self.linear_relpos = Linear(self.no_bins, c_z)

    def relpos(self, batch):
        pos = batch["residue_index"]
        asym_id = batch["asym_id"]
        asym_id_same = (asym_id[..., None] == asym_id[..., None, :])
        offset = pos[..., None] - pos[..., None, :]

        clipped_offset = torch.clamp(
            offset + self.max_relative_idx, 0, 2 * self.max_relative_idx
        )

        rel_feats = []
        if(self.use_chain_relative):
            final_offset = torch.where(
                asym_id_same, 
                clipped_offset,
                (2 * self.max_relative_idx + 1) * 
                torch.ones_like(clipped_offset)
            )
            boundaries = torch.arange(
                start=0, end=2 * self.max_relative_idx + 2, device=final_offset.device
            )
            rel_pos = one_hot(
                final_offset,
                boundaries,
            )

            rel_feats.append(rel_pos)

            entity_id = batch["entity_id"]
            entity_id_same = (entity_id[..., None] == entity_id[..., None, :])
            rel_feats.append(entity_id_same[..., None].to(dtype=rel_pos.dtype))

            sym_id = batch["sym_id"]
            rel_sym_id = sym_id[..., None] - sym_id[..., None, :]

            max_rel_chain = self.max_relative_chain
            clipped_rel_chain = torch.clamp(
                rel_sym_id + max_rel_chain,
                0,
                2 * max_rel_chain,
            )

            final_rel_chain = torch.where(
                entity_id_same,
                clipped_rel_chain,
                (2 * max_rel_chain + 1) *
                torch.ones_like(clipped_rel_chain)
            )

            boundaries = torch.arange(
                start=0, end=2 * max_rel_chain + 2, device=final_rel_chain.device
            )
            rel_chain = one_hot(
                final_rel_chain,
                boundaries,
            )

            rel_feats.append(rel_chain)
        else:
            boundaries = torch.arange(
                start=0, end=2 * self.max_relative_idx + 1, device=clipped_offset.device
            )
            rel_pos = one_hot(
                clipped_offset, boundaries,
            )
            rel_feats.append(rel_pos)

        rel_feat = torch.cat(rel_feats, dim=-1).to(
            self.linear_relpos.weight.dtype
        )

        return self.linear_relpos(rel_feat)

    def forward(self, batch) -> Tuple[torch.Tensor, torch.Tensor]:
        tf = batch["target_feat"] 
        msa = batch["msa_feat"]

        # [*, N_res, c_z]
        tf_emb_i = self.linear_tf_z_i(tf)
        tf_emb_j = self.linear_tf_z_j(tf)

        # [*, N_res, N_res, c_z]
        pair_emb = tf_emb_i[..., None, :] + tf_emb_j[..., None, :, :]
        pair_emb = pair_emb + self.relpos(batch)

        # [*, N_clust, N_res, c_m]
        n_clust = msa.shape[-3]
        tf_m = (
            self.linear_tf_m(tf)
            .unsqueeze(-3)
            .expand(((-1,) * len(tf.shape[:-2]) + (n_clust, -1, -1)))
        )
        msa_emb = self.linear_msa_m(msa) + tf_m

        return msa_emb, pair_emb


class PreembeddingEmbedder(nn.Module):
    """
    Embeds the sequence pre-embedding passed to the model and the target_feat features.
    """

    def __init__(
        self,
        tf_dim: int,
        preembedding_dim: int,
        c_z: int,
        c_m: int,
        relpos_k: int,
        **kwargs,
    ):
        """
        Args:
            tf_dim:
                End channel dimension of the incoming target features
            preembedding_dim:
                End channel dimension of the incoming embeddings
            c_z:
                Pair embedding dimension
            c_m:
                Single-Seq embedding dimension
            relpos_k:
                Window size used in relative position encoding
        """
        super(PreembeddingEmbedder, self).__init__()

        self.tf_dim = tf_dim
        self.preembedding_dim = preembedding_dim

        self.c_z = c_z
        self.c_m = c_m

        self.linear_tf_m = Linear(tf_dim, c_m)
        self.linear_preemb_m = Linear(self.preembedding_dim, c_m)
        self.linear_preemb_z_i = Linear(self.preembedding_dim, c_z)
        self.linear_preemb_z_j = Linear(self.preembedding_dim, c_z)

        # Relative Positional Encoding
        self.relpos_k = relpos_k
        self.no_bins = 2 * relpos_k + 1
        self.linear_relpos = Linear(self.no_bins, c_z)

    def relpos(self, ri: torch.Tensor):
        """
        Computes relative positional encodings
        Args:
            ri:
                "residue_index" feature of shape [*, N]
        Returns:
                Relative positional encoding of protein using the
                residue_index feature
        """
        d = ri[..., None] - ri[..., None, :]
        boundaries = torch.arange(
            start=-self.relpos_k, end=self.relpos_k + 1, device=d.device
        )
        reshaped_bins = boundaries.view(((1,) * len(d.shape)) + (len(boundaries),))
        d = d[..., None] - reshaped_bins
        d = torch.abs(d)
        d = torch.argmin(d, dim=-1)
        d = nn.functional.one_hot(d, num_classes=len(boundaries)).float()
        d = d.to(ri.dtype)
        return self.linear_relpos(d)

    def forward(
        self,
        tf: torch.Tensor,
        ri: torch.Tensor,
        preemb: torch.Tensor,
        inplace_safe: bool = False,
    ) -> Tuple[torch.Tensor, torch.Tensor]:

        tf_m = (
            self.linear_tf_m(tf)
            .unsqueeze(-3)
        )
        preemb_emb = self.linear_preemb_m(preemb[..., None, :, :]) + tf_m
        preemb_emb_i = self.linear_preemb_z_i(preemb)
        preemb_emb_j = self.linear_preemb_z_j(preemb)

        pair_emb = self.relpos(ri.type(preemb_emb_i.dtype))
        pair_emb = add(pair_emb,
                       preemb_emb_i[..., None, :],
                       inplace=inplace_safe)
        pair_emb = add(pair_emb,
                       preemb_emb_j[..., None, :, :],
                       inplace=inplace_safe)

        return preemb_emb, pair_emb


class RecyclingEmbedder(nn.Module):
    """
    Embeds the output of an iteration of the model for recycling.

    Implements Algorithm 32.
    """
    def __init__(
        self,
        c_m: int,
        c_z: int,
        min_bin: float,
        max_bin: float,
        no_bins: int,
        inf: float = 1e8,
        **kwargs,
    ):
        """
        Args:
            c_m:
                MSA channel dimension
            c_z:
                Pair embedding channel dimension
            min_bin:
                Smallest distogram bin (Angstroms)
            max_bin:
                Largest distogram bin (Angstroms)
            no_bins:
                Number of distogram bins
        """
        super(RecyclingEmbedder, self).__init__()

        self.c_m = c_m
        self.c_z = c_z
        self.min_bin = min_bin
        self.max_bin = max_bin
        self.no_bins = no_bins
        self.inf = inf

        self.linear = Linear(self.no_bins, self.c_z)
        self.layer_norm_m = LayerNorm(self.c_m)
        self.layer_norm_z = LayerNorm(self.c_z)

    def forward(
        self,
        m: torch.Tensor,
        z: torch.Tensor,
        x: torch.Tensor,
        inplace_safe: bool = False,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            m:
                First row of the MSA embedding. [*, N_res, C_m]
            z:
                [*, N_res, N_res, C_z] pair embedding
            x:
                [*, N_res, 3] predicted C_beta coordinates
        Returns:
            m:
                [*, N_res, C_m] MSA embedding update
            z:
                [*, N_res, N_res, C_z] pair embedding update
        """
        # [*, N, C_m]
        m_update = self.layer_norm_m(m)
        if(inplace_safe):
            m.copy_(m_update)
            m_update = m

        # [*, N, N, C_z]
        z_update = self.layer_norm_z(z)
        if(inplace_safe):
            z.copy_(z_update)
            z_update = z

        # This squared method might become problematic in FP16 mode.
        bins = torch.linspace(
            self.min_bin,
            self.max_bin,
            self.no_bins,
            dtype=x.dtype,
            device=x.device,
            requires_grad=False,
        )
        squared_bins = bins ** 2
        upper = torch.cat(
            [squared_bins[1:], squared_bins.new_tensor([self.inf])], dim=-1
        )
        d = torch.sum(
            (x[..., None, :] - x[..., None, :, :]) ** 2, dim=-1, keepdims=True
        )

        # [*, N, N, no_bins]
        d = ((d > squared_bins) * (d < upper)).type(x.dtype)

        # [*, N, N, C_z]
        d = self.linear(d)
        z_update = add(z_update, d, inplace_safe)

        return m_update, z_update


class TemplateSingleEmbedder(nn.Module):
    """
    Embeds the "template_angle_feat" feature.

    Implements Algorithm 2, line 7.
    """

    def __init__(
        self,
        c_in: int,
        c_out: int,
        **kwargs,
    ):
        """
        Args:
            c_in:
                Final dimension of "template_angle_feat"
            c_out:
                Output channel dimension
        """
        super(TemplateSingleEmbedder, self).__init__()

        self.c_out = c_out
        self.c_in = c_in

        self.linear_1 = Linear(self.c_in, self.c_out, init="relu")
        self.relu = nn.ReLU()
        self.linear_2 = Linear(self.c_out, self.c_out, init="relu")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: [*, N_templ, N_res, c_in] "template_angle_feat" features
        Returns:
            x: [*, N_templ, N_res, C_out] embedding
        """
        x = self.linear_1(x)
        x = self.relu(x)
        x = self.linear_2(x)

        return x


class TemplatePairEmbedder(nn.Module):
    """
    Embeds "template_pair_feat" features.

    Implements Algorithm 2, line 9.
    """

    def __init__(
        self,
        c_in: int,
        c_out: int,
        **kwargs,
    ):
        """
        Args:
            c_in:

            c_out:
                Output channel dimension
        """
        super(TemplatePairEmbedder, self).__init__()

        self.c_in = c_in
        self.c_out = c_out

        # Despite there being no relu nearby, the source uses that initializer
        self.linear = Linear(self.c_in, self.c_out, init="relu")

    def forward(
        self,
        x: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, C_in] input tensor
        Returns:
            [*, C_out] output tensor
        """
        x = self.linear(x)

        return x


class ExtraMSAEmbedder(nn.Module):
    """
    Embeds unclustered MSA sequences.

    Implements Algorithm 2, line 15
    """
    def __init__(
        self,
        c_in: int,
        c_out: int,
        **kwargs,
    ):
        """
        Args:
            c_in:
                Input channel dimension
            c_out:
                Output channel dimension
        """
        super(ExtraMSAEmbedder, self).__init__()

        self.c_in = c_in
        self.c_out = c_out

        self.linear = Linear(self.c_in, self.c_out)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x:
                [*, N_extra_seq, N_res, C_in] "extra_msa_feat" features
        Returns:
            [*, N_extra_seq, N_res, C_out] embedding
        """
        x = self.linear(x)

        return x


class TemplateEmbedder(nn.Module):
    def __init__(self, config):
        super(TemplateEmbedder, self).__init__()
        
        self.config = config
        self.template_single_embedder = TemplateSingleEmbedder(
            **config["template_single_embedder"],
        )
        self.template_pair_embedder = TemplatePairEmbedder(
            **config["template_pair_embedder"],
        )
        self.template_pair_stack = TemplatePairStack(
            **config["template_pair_stack"],
        )
        self.template_pointwise_att = TemplatePointwiseAttention(
            **config["template_pointwise_attention"],
        )
    
    def forward(
        self,
        batch,
        z, 
        pair_mask,
        templ_dim,
        chunk_size,
        _mask_trans=True,
        use_deepspeed_evo_attention=False,
        use_lma=False,
        inplace_safe=False
    ):
        # Embed the templates one at a time (with a poor man's vmap)
        pair_embeds = []
        n = z.shape[-2]
        n_templ = batch["template_aatype"].shape[templ_dim]

        if (inplace_safe):
            # We'll preallocate the full pair tensor now to avoid manifesting
            # a second copy during the stack later on
            t_pair = z.new_zeros(
                z.shape[:-3] +
                (n_templ, n, n, self.config.template_pair_embedder.c_out)
            )

        for i in range(n_templ):
            idx = batch["template_aatype"].new_tensor(i)
            single_template_feats = tensor_tree_map(
                lambda t: torch.index_select(t, templ_dim, idx).squeeze(templ_dim),
                batch,
            )

            # [*, N, N, C_t]
            t = build_template_pair_feat(
                single_template_feats,
                use_unit_vector=self.config.use_unit_vector,
                inf=self.config.inf,
                eps=self.config.eps,
                **self.config.distogram,
            ).to(z.dtype)
            t = self.template_pair_embedder(t)

            if (inplace_safe):
                t_pair[..., i, :, :, :] = t
            else:
                pair_embeds.append(t)

            del t

        if (not inplace_safe):
            t_pair = torch.stack(pair_embeds, dim=templ_dim)

        del pair_embeds

        # [*, S_t, N, N, C_z]
        t = self.template_pair_stack(
            t_pair,
            pair_mask.unsqueeze(-3).to(dtype=z.dtype),
            chunk_size=chunk_size,
            use_deepspeed_evo_attention=use_deepspeed_evo_attention,
            use_lma=use_lma,
            inplace_safe=inplace_safe,
            _mask_trans=_mask_trans,
        )
        del t_pair

        # [*, N, N, C_z]
        t = self.template_pointwise_att(
            t,
            z,
            template_mask=batch["template_mask"].to(dtype=z.dtype),
            use_lma=use_lma,
        )

        t_mask = torch.sum(batch["template_mask"], dim=-1) > 0
        # Append singletons
        t_mask = t_mask.reshape(
            *t_mask.shape, *([1] * (len(t.shape) - len(t_mask.shape)))
        )

        if (inplace_safe):
            t *= t_mask
        else:
            t = t * t_mask

        ret = {}

        ret.update({"template_pair_embedding": t})

        del t

        if self.config.embed_angles:
            template_angle_feat = build_template_angle_feat(
                batch
            )

            # [*, S_t, N, C_m]
            a = self.template_single_embedder(template_angle_feat)

            ret["template_single_embedding"] = a

        return ret


class TemplatePairEmbedderMultimer(nn.Module):
    def __init__(self,
        c_in: int,
        c_out: int,
        c_dgram: int,
        c_aatype: int,
    ):
        super(TemplatePairEmbedderMultimer, self).__init__()

        self.dgram_linear = Linear(c_dgram, c_out, init='relu')
        self.aatype_linear_1 = Linear(c_aatype, c_out, init='relu')
        self.aatype_linear_2 = Linear(c_aatype, c_out, init='relu')
        self.query_embedding_layer_norm = LayerNorm(c_in)
        self.query_embedding_linear = Linear(c_in, c_out, init='relu')
        
        self.pseudo_beta_mask_linear = Linear(1, c_out, init='relu')
        self.x_linear = Linear(1, c_out, init='relu')
        self.y_linear = Linear(1, c_out, init='relu')
        self.z_linear = Linear(1, c_out, init='relu')
        self.backbone_mask_linear = Linear(1, c_out, init='relu')

    def forward(self,
        template_dgram: torch.Tensor,
        aatype_one_hot: torch.Tensor,
        query_embedding: torch.Tensor,
        pseudo_beta_mask: torch.Tensor,
        backbone_mask: torch.Tensor,
        multichain_mask_2d: torch.Tensor,
        unit_vector: geometry.Vec3Array,
    ) -> torch.Tensor:
        act = 0.

        pseudo_beta_mask_2d = (
            pseudo_beta_mask[..., None] * pseudo_beta_mask[..., None, :]
        )
        pseudo_beta_mask_2d *= multichain_mask_2d
        template_dgram *= pseudo_beta_mask_2d[..., None]
        act += self.dgram_linear(template_dgram)
        act += self.pseudo_beta_mask_linear(pseudo_beta_mask_2d[..., None])
       
        aatype_one_hot = aatype_one_hot.to(template_dgram.dtype)
        act += self.aatype_linear_1(aatype_one_hot[..., None, :, :])
        act += self.aatype_linear_2(aatype_one_hot[..., None, :])

        backbone_mask_2d = (
            backbone_mask[..., None] * backbone_mask[..., None, :]
        )
        backbone_mask_2d *= multichain_mask_2d
        x, y, z = [(coord * backbone_mask_2d).to(dtype=query_embedding.dtype) for coord in unit_vector]
        act += self.x_linear(x[..., None])
        act += self.y_linear(y[..., None])
        act += self.z_linear(z[..., None])
       
        act += self.backbone_mask_linear(backbone_mask_2d[..., None].to(dtype=query_embedding.dtype))

        query_embedding = self.query_embedding_layer_norm(query_embedding)
        act += self.query_embedding_linear(query_embedding)

        return act


class TemplateSingleEmbedderMultimer(nn.Module):
    def __init__(self,
        c_in: int,
        c_out: int,
    ):
        super(TemplateSingleEmbedderMultimer, self).__init__()
        self.template_single_embedder = Linear(c_in, c_out)
        self.template_projector = Linear(c_out, c_out)
    
    def forward(self,
        batch,
        atom_pos,
        aatype_one_hot,
    ):
        out = {}

        dtype = batch["template_all_atom_positions"].dtype

        template_chi_angles, template_chi_mask = (
            all_atom_multimer.compute_chi_angles(
                atom_pos,
                batch["template_all_atom_mask"],
                batch["template_aatype"],
            )
        )

        template_features = torch.cat(
            [
                aatype_one_hot,
                torch.sin(template_chi_angles) * template_chi_mask,
                torch.cos(template_chi_angles) * template_chi_mask,
                template_chi_mask,
            ],
            dim=-1,
        ).to(dtype=dtype)

        template_mask = template_chi_mask[..., 0].to(dtype=dtype)

        template_activations = self.template_single_embedder(
            template_features
        )
        template_activations = torch.nn.functional.relu(
            template_activations
        )
        template_activations = self.template_projector(
            template_activations,
        )

        out["template_single_embedding"] = (
            template_activations
        )
        out["template_mask"] = template_mask

        return out 

 
class TemplateEmbedderMultimer(nn.Module):
    def __init__(self, config):
        super(TemplateEmbedderMultimer, self).__init__()
        
        self.config = config
        self.template_pair_embedder = TemplatePairEmbedderMultimer(
            **config["template_pair_embedder"],
        )
        self.template_single_embedder = TemplateSingleEmbedderMultimer(
            **config["template_single_embedder"],
        )
        self.template_pair_stack = TemplatePairStack(
            **config["template_pair_stack"],
        )

        self.linear_t = Linear(config.c_t, config.c_z)
    
    def forward(self, 
        batch, 
        z, 
        padding_mask_2d, 
        templ_dim,
        chunk_size,
        multichain_mask_2d,
        _mask_trans=True,
        use_deepspeed_evo_attention=False,
        use_lma=False,
        inplace_safe=False
    ):
        template_embeds = []
        n_templ = batch["template_aatype"].shape[templ_dim]
        for i in range(n_templ):
            idx = batch["template_aatype"].new_tensor(i)
            single_template_feats = tensor_tree_map(
                lambda t: torch.index_select(t, templ_dim, idx),
                batch,
            )

            single_template_embeds = {}
            act = 0.

            template_positions, pseudo_beta_mask = pseudo_beta_fn(
                single_template_feats["template_aatype"],
                single_template_feats["template_all_atom_positions"],
                single_template_feats["template_all_atom_mask"])

            template_dgram = dgram_from_positions(
                template_positions,
                inf=self.config.inf,
                **self.config.distogram,
            )

            aatype_one_hot = torch.nn.functional.one_hot(
                single_template_feats["template_aatype"], 22,
            )
            
            raw_atom_pos = single_template_feats["template_all_atom_positions"]

            # Vec3Arrays are required to be float32
            atom_pos = geometry.Vec3Array.from_array(raw_atom_pos.to(dtype=torch.float32))

            rigid, backbone_mask = all_atom_multimer.make_backbone_affine(
                atom_pos,
                single_template_feats["template_all_atom_mask"],
                single_template_feats["template_aatype"],
            )
            points = rigid.translation
            rigid_vec = rigid[..., None].inverse().apply_to_point(points)
            unit_vector = rigid_vec.normalized()

            pair_act = self.template_pair_embedder(
                template_dgram,
                aatype_one_hot,
                z,
                pseudo_beta_mask,
                backbone_mask,
                multichain_mask_2d,
                unit_vector,
            )

            single_template_embeds["template_pair_embedding"] = pair_act
            single_template_embeds.update(
                self.template_single_embedder(
                    single_template_feats,
                    atom_pos,
                    aatype_one_hot,
                )
            )
            template_embeds.append(single_template_embeds)

        template_embeds = dict_multimap(
            partial(torch.cat, dim=templ_dim),
            template_embeds,
        )

        # [*, S_t, N, N, C_z]
        t = self.template_pair_stack(
            template_embeds["template_pair_embedding"], 
            padding_mask_2d.unsqueeze(-3).to(dtype=z.dtype), 
            chunk_size=chunk_size,
            use_deepspeed_evo_attention=use_deepspeed_evo_attention,
            use_lma=use_lma,
            inplace_safe=inplace_safe,
            _mask_trans=_mask_trans,
        )
        # [*, N, N, C_z]
        t = torch.sum(t, dim=-4) / n_templ
        t = torch.nn.functional.relu(t)
        t = self.linear_t(t)
        template_embeds["template_pair_embedding"] = t

        return template_embeds
