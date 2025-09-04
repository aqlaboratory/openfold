# Copyright 2021 AlQuraishi Laboratory
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

"""
Unit tests to compare components of OpenFold run with the cuEquivariance memory-efficient
attention kernel vs. a stock PyTorch attention implementation.
"""

import unittest
import numpy as np
import pickle
import torch
from torch.nn import functional as F

from openfold.data import data_transforms
from openfold.model.primitives import (
    lecun_normal_init_,
    Attention
)
from openfold.utils.tensor_utils import tensor_tree_map

from tests.config import consts
import tests.compare_utils as compare_utils
from tests.data_utils import random_template_feats, random_attention_inputs



class TestCuEquivarianceKernel(unittest.TestCase):

    def test_compare_template_stack(self):
        """
        Compare Template Stack output with and without using DeepSpeed Evoformer attention kernel.
        Kernel can be used for Triangle Attention in the Template Pair Stack.
        """
        n_templ = consts.n_templ
        n_res = 20
        eps = 2e-2

        batch = random_template_feats(n_templ, n_res)
        batch["template_all_atom_masks"] = batch["template_all_atom_mask"]
        if consts.is_multimer:
            batch["asym_id"] = batch['asym_id'][0]

        pair_act = np.random.rand(n_res, n_res, consts.c_z).astype(np.float32)
        pair_mask = np.random.randint(0, 2, (n_res, n_res)).astype(np.float32)

        batch = {k: torch.as_tensor(v).cuda() for k, v in batch.items()}
        template_feats = {
            k: v for k, v in batch.items() if k.startswith("template_")
        }

        with torch.no_grad():
            model = compare_utils.get_global_pretrained_openfold()
            model.globals.use_deepspeed_evo_attention = False
            out_repro = model.embed_templates(
                template_feats,
                batch,
                torch.as_tensor(pair_act).cuda(),
                torch.as_tensor(pair_mask).cuda(),
                templ_dim=0,
                inplace_safe=False
            )
            out_repro = out_repro["template_pair_embedding"].cpu()

            model.globals.use_cuequivariance_attention = True
            model.globals.use_cuequivariance_multiplicative_update = True
            
            out_repro_ds = model.embed_templates(
                template_feats,
                batch,
                torch.as_tensor(pair_act).cuda(),
                torch.as_tensor(pair_mask).cuda(),
                templ_dim=0,
                inplace_safe=False
            )
            out_repro_ds = out_repro_ds["template_pair_embedding"].cpu()

            compare_utils.assert_max_abs_diff_small(out_repro, out_repro_ds, eps)

    def test_compare_model(self):
        """
        Run full model with and without using CuEquivariance Evoformer attention kernel
        and compare output coordinates.
        """
        eps = 0.2
        with open("tests/test_data/sample_feats.pickle", "rb") as fp:
            batch = pickle.load(fp)

        # atom37_to_atom14 doesn't like batches
        batch["residx_atom14_to_atom37"] = batch["residx_atom14_to_atom37"][0]
        batch["atom14_atom_exists"] = batch["atom14_atom_exists"][0]

        batch["no_recycling_iters"] = np.array([3., 3., 3., 3., ])

        if consts.is_multimer:
            n_res = batch['aatype'].shape[1]
            n_extra_seq = batch['extra_msa'].shape[1]
            batch["asym_id"] = np.ones((4, n_res))
            batch["entity_id"] = np.ones((4, n_res))
            batch["sym_id"] = np.ones((4, n_res))
            batch["extra_deletion_matrix"] = np.random.randint(0, 2, size=(4, n_extra_seq, n_res))
        
        batch = {k: torch.as_tensor(v).cuda() for k, v in batch.items()}

        batch["aatype"] = batch["aatype"].long()
        batch["template_aatype"] = batch["template_aatype"].long()
        batch["extra_msa"] = batch["extra_msa"].long()
        batch["residx_atom37_to_atom14"] = batch[
            "residx_atom37_to_atom14"
        ].long()
        batch["target_feat"] = torch.nn.functional.one_hot(batch["aatype"], consts.msa_logits - 1).to(torch.float32)
        batch["template_all_atom_mask"] = batch["template_all_atom_masks"]
        batch.update(
            data_transforms.atom37_to_torsion_angles("template_")(batch)
        )

        # Move the recycling dimension to the end
        move_dim = lambda t: t.permute(*range(len(t.shape))[1:], 0)
        batch = tensor_tree_map(move_dim, batch)
        # Restrict this test to use only torch.float32 precision due to instability with torch.bfloat16
        # https://github.com/aqlaboratory/openfold/issues/532
        with torch.no_grad(), torch.cuda.amp.autocast(dtype=torch.float32):
                model = compare_utils.get_global_pretrained_openfold()
                model.globals.use_cuequivariance_attention = False
                model.globals.use_cuequivariance_multiplicative_update = False
                out_repro = model(batch)
                out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)
                out_repro = out_repro["sm"]["positions"][-1].squeeze(0)

                # Enable attention
                model.globals.use_cuequivariance_attention = True
                out_repro_attn = model(batch)
                out_repro_attn = tensor_tree_map(lambda t: t.cpu(), out_repro_attn)
                out_repro_attn = out_repro_attn["sm"]["positions"][-1].squeeze(0)

                compare_utils.assert_mean_abs_diff_small(out_repro, out_repro_attn, eps)

                # Enable multiplication
                model.globals.use_cuequivariance_attention = True
                model.globals.use_cuequivariance_multiplicative_update = True
                out_repro_mul = model(batch)
                out_repro_mul = tensor_tree_map(lambda t: t.cpu(), out_repro_mul)
                out_repro_mul = out_repro_mul["sm"]["positions"][-1].squeeze(0)

                compare_utils.assert_mean_abs_diff_small(out_repro_attn, out_repro_mul, eps)


if __name__ == "__main__":
    unittest.main()
