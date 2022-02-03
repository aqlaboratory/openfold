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

from openfold.utils.feats import (
    pseudo_beta_fn,
    build_extra_msa_feat,
    build_template_angle_feat,
    build_template_pair_feat,
    atom14_to_atom37,
)
from openfold.model.embedders import (
    InputEmbedder,
    RecyclingEmbedder,
    TemplateAngleEmbedder,
    TemplatePairEmbedder,
    ExtraMSAEmbedder,
)
from openfold.model.evoformer import EvoformerStack, ExtraMSAStack
from openfold.model.heads import AuxiliaryHeads
import openfold.np.residue_constants as residue_constants
from openfold.model.structure_module import StructureModule
from openfold.model.template import (
    TemplatePairStack,
    TemplatePointwiseAttention,
)
from openfold.utils.loss import (
    compute_plddt,
)
from openfold.utils.tensor_utils import (
    dict_multimap,
    tensor_tree_map,
)


class AlphaFold(nn.Module):
    """
    Alphafold 2.

    Implements Algorithm 2 (but with training).
    """

    def __init__(self, config):
        """
        Args:
            config:
                A dict-like config object (like the one in config.py)
        """
        super(AlphaFold, self).__init__()

        self.globals = config.globals
        config = config.model
        template_config = config.template
        extra_msa_config = config.extra_msa

        # Main trunk + structure module
        self.input_embedder = InputEmbedder(
            **config["input_embedder"],
        )
        self.recycling_embedder = RecyclingEmbedder(
            **config["recycling_embedder"],
        )
        self.template_angle_embedder = TemplateAngleEmbedder(
            **template_config["template_angle_embedder"],
        )
        self.template_pair_embedder = TemplatePairEmbedder(
            **template_config["template_pair_embedder"],
        )
        self.template_pair_stack = TemplatePairStack(
            **template_config["template_pair_stack"],
        )
        self.template_pointwise_att = TemplatePointwiseAttention(
            **template_config["template_pointwise_attention"],
        )
        self.extra_msa_embedder = ExtraMSAEmbedder(
            **extra_msa_config["extra_msa_embedder"],
        )
        self.extra_msa_stack = ExtraMSAStack(
            **extra_msa_config["extra_msa_stack"],
        )
        self.evoformer = EvoformerStack(
            **config["evoformer_stack"],
        )
        self.structure_module = StructureModule(
            **config["structure_module"],
        )

        self.aux_heads = AuxiliaryHeads(
            config["heads"],
        )

        self.config = config

    def embed_templates(self, batch, z, pair_mask, templ_dim): 
        # Embed the templates one at a time (with a poor man's vmap)
        template_embeds = []
        n_templ = batch["template_aatype"].shape[templ_dim]
        for i in range(n_templ):
            idx = batch["template_aatype"].new_tensor(i)
            single_template_feats = tensor_tree_map(
                lambda t: torch.index_select(t, templ_dim, idx),
                batch,
            )

            single_template_embeds = {}
            if self.config.template.embed_angles:
                template_angle_feat = build_template_angle_feat(
                    single_template_feats,
                )

                # [*, S_t, N, C_m]
                a = self.template_angle_embedder(template_angle_feat)

                single_template_embeds["angle"] = a

            # [*, S_t, N, N, C_t]
            t = build_template_pair_feat(
                single_template_feats,
                inf=self.config.template.inf,
                eps=self.config.template.eps,
                **self.config.template.distogram,
            ).to(z.dtype)
            t = self.template_pair_embedder(t)

            single_template_embeds.update({"pair": t})

            template_embeds.append(single_template_embeds)

        template_embeds = dict_multimap(
            partial(torch.cat, dim=templ_dim),
            template_embeds,
        )

        # [*, S_t, N, N, C_z]
        t = self.template_pair_stack(
            template_embeds["pair"], 
            pair_mask.unsqueeze(-3).to(dtype=z.dtype), 
            chunk_size=self.globals.chunk_size,
            _mask_trans=self.config._mask_trans,
        )

        # [*, N, N, C_z]
        t = self.template_pointwise_att(
            t, 
            z, 
            template_mask=batch["template_mask"].to(dtype=z.dtype),
            chunk_size=self.globals.chunk_size,
        )
        t = t * (torch.sum(batch["template_mask"]) > 0)

        ret = {}
        if self.config.template.embed_angles:
            ret["template_angle_embedding"] = template_embeds["angle"]

        ret.update({"template_pair_embedding": t})

        return ret

    def iteration(self, feats, m_1_prev, z_prev, x_prev, _recycle=True):
        # Primary output dictionary
        outputs = {}

        # This needs to be done manually for DeepSpeed's sake
        dtype = next(self.parameters()).dtype
        for k in feats:
            if(feats[k].dtype == torch.float32):
                feats[k] = feats[k].to(dtype=dtype)

        # Grab some data about the input
        batch_dims = feats["target_feat"].shape[:-2]
        no_batch_dims = len(batch_dims)
        n = feats["target_feat"].shape[-2]
        n_seq = feats["msa_feat"].shape[-3]
        device = feats["target_feat"].device

        # Prep some features
        seq_mask = feats["seq_mask"]
        pair_mask = seq_mask[..., None] * seq_mask[..., None, :]
        msa_mask = feats["msa_mask"]

        # Initialize the MSA and pair representations

        # m: [*, S_c, N, C_m]
        # z: [*, N, N, C_z]
        m, z = self.input_embedder(
            feats["target_feat"],
            feats["residue_index"],
            feats["msa_feat"],
        )

        # Initialize the recycling embeddings, if needs be
        if None in [m_1_prev, z_prev, x_prev]:
            # [*, N, C_m]
            m_1_prev = m.new_zeros(
                (*batch_dims, n, self.config.input_embedder.c_m),
                requires_grad=False,
            )

            # [*, N, N, C_z]
            z_prev = z.new_zeros(
                (*batch_dims, n, n, self.config.input_embedder.c_z),
                requires_grad=False,
            )

            # [*, N, 3]
            x_prev = z.new_zeros(
                (*batch_dims, n, residue_constants.atom_type_num, 3),
                requires_grad=False,
            )

        x_prev = pseudo_beta_fn(
            feats["aatype"], x_prev, None
        ).to(dtype=z.dtype)

        # m_1_prev_emb: [*, N, C_m]
        # z_prev_emb: [*, N, N, C_z]
        m_1_prev_emb, z_prev_emb = self.recycling_embedder(
            m_1_prev,
            z_prev,
            x_prev,
        )

        # If the number of recycling iterations is 0, skip recycling
        # altogether. We zero them this way instead of computing them
        # conditionally to avoid leaving parameters unused, which has annoying
        # implications for DDP training.
        if(not _recycle):
            m_1_prev_emb *= 0
            z_prev_emb *= 0

        # [*, S_c, N, C_m]
        m[..., 0, :, :] += m_1_prev_emb

        # [*, N, N, C_z]
        z += z_prev_emb

        # Possibly prevents memory fragmentation
        del m_1_prev, z_prev, x_prev, m_1_prev_emb, z_prev_emb

        # Embed the templates + merge with MSA/pair embeddings
        if self.config.template.enabled:
            template_feats = {
                k: v for k, v in feats.items() if k.startswith("template_")
            }
            template_embeds = self.embed_templates(
                template_feats,
                z,
                pair_mask.to(dtype=z.dtype),
                no_batch_dims,
            )

            # [*, N, N, C_z]
            z = z + template_embeds["template_pair_embedding"]

            if self.config.template.embed_angles:
                # [*, S = S_c + S_t, N, C_m]
                m = torch.cat(
                    [m, template_embeds["template_angle_embedding"]], 
                    dim=-3
                )

                # [*, S, N]
                torsion_angles_mask = feats["template_torsion_angles_mask"]
                msa_mask = torch.cat(
                    [feats["msa_mask"], torsion_angles_mask[..., 2]], 
                    dim=-2
                )

        # Embed extra MSA features + merge with pairwise embeddings
        if self.config.extra_msa.enabled:
            # [*, S_e, N, C_e]
            a = self.extra_msa_embedder(build_extra_msa_feat(feats))

            # [*, N, N, C_z]
            z = self.extra_msa_stack(
                a,
                z,
                msa_mask=feats["extra_msa_mask"].to(dtype=a.dtype),
                chunk_size=self.globals.chunk_size,
                pair_mask=pair_mask.to(dtype=z.dtype),
                _mask_trans=self.config._mask_trans,
            )

        # Run MSA + pair embeddings through the trunk of the network
        # m: [*, S, N, C_m]
        # z: [*, N, N, C_z]
        # s: [*, N, C_s]
        m, z, s = self.evoformer(
            m,
            z,
            msa_mask=msa_mask.to(dtype=m.dtype),
            pair_mask=pair_mask.to(dtype=z.dtype),
            chunk_size=self.globals.chunk_size,
            _mask_trans=self.config._mask_trans,
        )

        outputs["msa"] = m[..., :n_seq, :, :]
        outputs["pair"] = z
        outputs["single"] = s

        # Predict 3D structure
        outputs["sm"] = self.structure_module(
            s,
            z,
            feats["aatype"],
            mask=feats["seq_mask"].to(dtype=s.dtype),
        )
        outputs["final_atom_positions"] = atom14_to_atom37(
            outputs["sm"]["positions"][-1], feats
        )
        outputs["final_atom_mask"] = feats["atom37_atom_exists"]
        outputs["final_affine_tensor"] = outputs["sm"]["frames"][-1]

        # Save embeddings for use during the next recycling iteration

        # [*, N, C_m]
        m_1_prev = m[..., 0, :, :]

        # [*, N, N, C_z]
        z_prev = z

        # [*, N, 3]
        x_prev = outputs["final_atom_positions"]

        return outputs, m_1_prev, z_prev, x_prev

    def _disable_activation_checkpointing(self):
        self.template_pair_stack.blocks_per_ckpt = None
        self.evoformer.blocks_per_ckpt = None

        for b in self.extra_msa_stack.blocks:
            b.ckpt = False

    def _enable_activation_checkpointing(self):
        self.template_pair_stack.blocks_per_ckpt = (
            self.config.template.template_pair_stack.blocks_per_ckpt
        )
        self.evoformer.blocks_per_ckpt = (
            self.config.evoformer_stack.blocks_per_ckpt
        )

        for b in self.extra_msa_stack.blocks:
            b.ckpt = self.config.extra_msa.extra_msa_stack.ckpt

    def forward(self, batch):
        """
        Args:
            batch:
                Dictionary of arguments outlined in Algorithm 2. Keys must
                include the official names of the features in the
                supplement subsection 1.2.9.

                The final dimension of each input must have length equal to
                the number of recycling iterations.

                Features (without the recycling dimension):

                    "aatype" ([*, N_res]):
                        Contrary to the supplement, this tensor of residue
                        indices is not one-hot.
                    "target_feat" ([*, N_res, C_tf])
                        One-hot encoding of the target sequence. C_tf is
                        config.model.input_embedder.tf_dim.
                    "residue_index" ([*, N_res])
                        Tensor whose final dimension consists of
                        consecutive indices from 0 to N_res.
                    "msa_feat" ([*, N_seq, N_res, C_msa])
                        MSA features, constructed as in the supplement.
                        C_msa is config.model.input_embedder.msa_dim.
                    "seq_mask" ([*, N_res])
                        1-D sequence mask
                    "msa_mask" ([*, N_seq, N_res])
                        MSA mask
                    "pair_mask" ([*, N_res, N_res])
                        2-D pair mask
                    "extra_msa_mask" ([*, N_extra, N_res])
                        Extra MSA mask
                    "template_mask" ([*, N_templ])
                        Template mask (on the level of templates, not
                        residues)
                    "template_aatype" ([*, N_templ, N_res])
                        Tensor of template residue indices (indices greater
                        than 19 are clamped to 20 (Unknown))
                    "template_all_atom_positions"
                        ([*, N_templ, N_res, 37, 3])
                        Template atom coordinates in atom37 format
                    "template_all_atom_mask" ([*, N_templ, N_res, 37])
                        Template atom coordinate mask
                    "template_pseudo_beta" ([*, N_templ, N_res, 3])
                        Positions of template carbon "pseudo-beta" atoms
                        (i.e. C_beta for all residues but glycine, for
                        for which C_alpha is used instead)
                    "template_pseudo_beta_mask" ([*, N_templ, N_res])
                        Pseudo-beta mask
        """
        # Initialize recycling embeddings
        m_1_prev, z_prev, x_prev = None, None, None

        # Disable activation checkpointing for the first few recycling iters
        is_grad_enabled = torch.is_grad_enabled()
        self._disable_activation_checkpointing()

        # Main recycling loop
        num_iters = batch["aatype"].shape[-1]
        for cycle_no in range(num_iters):
            # Select the features for the current recycling cycle
            fetch_cur_batch = lambda t: t[..., cycle_no]
            feats = tensor_tree_map(fetch_cur_batch, batch)

            # Enable grad iff we're training and it's the final recycling layer
            is_final_iter = cycle_no == (num_iters - 1)
            with torch.set_grad_enabled(is_grad_enabled and is_final_iter):
                if is_final_iter:
                    self._enable_activation_checkpointing()
                    # Sidestep AMP bug (PyTorch issue #65766)
                    if torch.is_autocast_enabled():
                        torch.clear_autocast_cache()

                # Run the next iteration of the model
                outputs, m_1_prev, z_prev, x_prev = self.iteration(
                    feats,
                    m_1_prev,
                    z_prev,
                    x_prev,
                    _recycle=(num_iters > 1)
                )

        # Run auxiliary heads
        outputs.update(self.aux_heads(outputs))

        return outputs
