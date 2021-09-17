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

import torch
import torch.nn as nn

from alphafold.utils.feats import (
    pseudo_beta_fn,
    atom37_to_torsion_angles,
    build_extra_msa_feat,
    build_template_angle_feat,
    build_template_pair_feat,
    atom14_to_atom37,
)
from alphafold.model.embedders import (
    InputEmbedder, 
    RecyclingEmbedder,
    TemplateAngleEmbedder,
    TemplatePairEmbedder,
    ExtraMSAEmbedder,
)
from alphafold.model.evoformer import EvoformerStack, ExtraMSAStack
from alphafold.model.heads import AuxiliaryHeads
import alphafold.np.residue_constants as residue_constants
from alphafold.model.structure_module import StructureModule
from alphafold.model.template import (
    TemplatePairStack, 
    TemplatePointwiseAttention,
)
from alphafold.utils.loss import (
    compute_plddt,
)

from alphafold.utils.tensor_utils import (
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

    def embed_templates(self, batch, z, pair_mask):
        # Build template angle feats
        angle_feats = atom37_to_torsion_angles(
            batch["template_aatype"], 
            batch["template_all_atom_positions"], 
            batch["template_all_atom_masks"], 
            eps=1e-8
        )

        # Stow this away for later
        batch["torsion_angles_mask"] = angle_feats["torsion_angles_mask"]

        template_angle_feat = build_template_angle_feat(
            angle_feats,
            batch["template_aatype"],
        )
 
        # [*, S_t, N, C_m]
        a = self.template_angle_embedder(template_angle_feat)

        # [*, S_t, N, N, C_t]
        t = build_template_pair_feat(
            batch,
            eps=self.config.template.eps,
            **self.config.template.distogram
        )
        t = self.template_pair_embedder(t)
        t = self.template_pair_stack(
            t, 
            pair_mask.unsqueeze(-3),
            _mask_trans=self.config._mask_trans
        )

        # [*, N, N, C_z]
        t = self.template_pointwise_att(
            t, 
            z, 
            template_mask=batch["template_mask"]
        )
        t *= torch.sum(batch["template_mask"]) > 0
    
        return a, t

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
                        "template_all_atom_pos" ([*, N_templ, N_res, 37, 3])
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
        # Recycling embeddings
        m_1_prev, z_prev, x_prev = None, None, None

        # Primary output dictionary
        outputs = {}

        # Main recycling loop
        for cycle_no in range(self.config.no_cycles):
            # Select the features for the current recycling cycle
            fetch_cur_batch = lambda t: t[..., cycle_no] 
            feats = tensor_tree_map(fetch_cur_batch, batch)

            # Grab some data about the input
            batch_dims = feats["target_feat"].shape[:-2]
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

            # Inject information from previous recycling iterations
            if(self.config.no_cycles > 1):
                # Initialize the recycling embeddings, if needs be
                if(None in [m_1_prev, z_prev, x_prev]):
                    # [*, N, C_m]
                    m_1_prev = torch.zeros(
                        (*batch_dims, n, self.config.c_m), 
                        requires_grad=False, 
                        device=device,
                    )

                    # [*, N, N, C_z]
                    z_prev = torch.zeros(
                        (*batch_dims, n, n, self.config.c_z),
                        requires_grad=False,
                        device=device,
                    )

                    # [*, N, 3]
                    x_prev = torch.zeros(
                        (*batch_dims, n, residue_constants.atom_type_num, 3),
                        requires_grad=False,
                        device=device,
                    )

                x_prev = pseudo_beta_fn(
                    feats["aatype"],
                    x_prev,
                    None  # TODO: figure this part out
                )

                # m_1_prev_emb: [*, N, C_m]
                # z_prev_emb: [*, N, N, C_z]
                m_1_prev_emb, z_prev_emb = self.recycling_embedder(
                    m_1_prev, 
                    z_prev, 
                    x_prev,
                )

                # [*, S_c, N, C_m]
                m[..., 0, :, :] += m_1_prev_emb

                # [*, N, N, C_z]
                z += z_prev_emb

            # Embed the templates + merge with MSA/pair embeddings
            if(self.config.template.enabled):
                a, t = self.embed_templates(feats, z, pair_mask)

                # [*, N, N, C_z]
                z += t
 
                if(self.config.template.embed_angles):
                    # [*, S = S_c + S_t, N, C_m]
                    m = torch.cat([m, a], dim=-3)

                    # [*, S, N]
                    torsion_angles_mask = feats["torsion_angles_mask"]
                    msa_mask = torch.cat(
                        [feats["msa_mask"], torsion_angles_mask[..., 2]], axis=-2
                    )

            # Embed extra MSA features + merge with pairwise embeddings 
            if(self.config.extra_msa.enabled):
                # [*, S_e, N, C_e]
                a = self.extra_msa_embedder(build_extra_msa_feat(feats))
           
                # [*, N, N, C_z]
                z = self.extra_msa_stack(
                    a, 
                    z, 
                    msa_mask=feats["extra_msa_mask"],
                    pair_mask=pair_mask,
                    _mask_trans=self.config._mask_trans,
                )

            # Run MSA + pair embeddings through the trunk of the network
            # m: [*, S, N, C_m]
            # z: [*, N, N, C_z]
            # s: [*, N, C_s]
            m, z, s = self.evoformer(
                m, 
                z, 
                msa_mask=msa_mask, 
                pair_mask=pair_mask,
                _mask_trans=self.config._mask_trans
            )

            outputs["msa"] = m[..., :n_seq, :, :]
            outputs["pair"] = z
            outputs["single"] = s

            # Predict 3D structure
            outputs["sm"] = self.structure_module(
                s, z, feats["aatype"], mask=feats["seq_mask"],
            )        
            outputs["final_atom_positions"] = atom14_to_atom37(
                outputs["sm"]["positions"][-1], feats
            )
            outputs["final_atom_mask"] = feats["atom37_atom_exists"]

            # Save embeddings for use during the next recycling iteration 

            # [*, N, C_m]
            m_1_prev = m[..., 0, :, :]

            # [* N, N, C_z]
            z_prev = z

            # [*, N, 3]
            x_prev = outputs["final_atom_positions"]

        outputs.update(self.aux_heads(outputs))

        return outputs
