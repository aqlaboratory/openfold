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
Unit tests to compare components of OpenFold run with the DeepSpeed memory-efficient
attention kernel, DS4Sci_EvoformerAttention vs. a stock PyTorch attention implementation.
"""

import torch
import unittest
import numpy as np
import pickle

from openfold.model.primitives import (
    Attention,
)
from tests.config import consts
import tests.compare_utils as compare_utils
from tests.data_utils import (
    random_template_feats,
    random_extra_msa_feats,
)

from openfold.config import model_config
from openfold.data import data_transforms
from openfold.model.model import AlphaFold
from openfold.utils.tensor_utils import tensor_tree_map


class TestDeepSpeedKernel(unittest.TestCase):
    def test_ds_kernel_vs_attention(self):
        """Compare attention with and without using DeepSpeed Evoformer kernel."""
        batch_size = consts.batch_size
        c_hidden = 32
        n = 2 ** 12
        n_seq = 12
        no_heads = 4
        eps = 2e-2

        q = torch.rand(batch_size, n_seq, n, c_hidden).cuda()
        kv = torch.rand(batch_size, n_seq, n, c_hidden).cuda()

        bias = [torch.rand(batch_size, n_seq, 1, 1, n), torch.rand(batch_size, 1, no_heads, n, n)]
        bias = [b.cuda() for b in bias]

        a = Attention(
            c_hidden, c_hidden, c_hidden, c_hidden, no_heads
        ).cuda()

        with torch.no_grad():
            l = a(q, kv, biases=bias, use_deepspeed_evo_attention=True)
            real = a(q, kv, biases=bias)

        self.assertTrue(torch.max(torch.abs(l - real)) < eps)

    def compare_evoformer(self, dtype):
        """
        Compare Evoformer output with and without using DeepSpeed Evoformer attention kernel.
        Set dtype to confirm the kernel can be used during both training (BF16) and inference (FP32),
        since the kernel itself can run with either BF16 or FP16 precision.
        """
        n_res = consts.n_res
        n_seq = consts.n_seq
        eps = 2e-2

        activations = {
            "msa": torch.rand(n_seq, n_res, consts.c_m, device='cuda', dtype=dtype),
            "pair": torch.rand(n_res, n_res, consts.c_z, device='cuda', dtype=dtype)
        }

        masks = {
            "msa": torch.randint(0, 2, (n_seq, n_res), device='cuda', dtype=dtype),
            "pair": torch.randint(0, 2, (n_res, n_res), device='cuda', dtype=dtype),
        }

        with torch.cuda.amp.autocast(dtype=dtype):
            model = compare_utils.get_global_pretrained_openfold()
            out_repro_msa, out_repro_pair = model.evoformer.blocks[0](
                activations["msa"],
                activations["pair"],
                masks["msa"],
                masks["pair"],
                use_deepspeed_evo_attention=False,
                chunk_size=4,
                _mask_trans=False,
                inplace_safe=False,
            )

            out_repro_msa = out_repro_msa.cpu()
            out_repro_pair = out_repro_pair.cpu()

            out_repro_msa_ds, out_repro_pair_ds = model.evoformer.blocks[0](
                activations["msa"],
                activations["pair"],
                masks["msa"],
                masks["pair"],
                use_deepspeed_evo_attention=True,
                chunk_size=4,
                _mask_trans=False,
                inplace_safe=False,
            )
            out_repro_msa_ds = out_repro_msa_ds.cpu()
            out_repro_pair_ds = out_repro_pair_ds.cpu()

            self.assertTrue(torch.allclose(torch.abs(out_repro_msa), torch.abs(out_repro_msa_ds), atol=eps))
            self.assertTrue(torch.allclose(torch.abs(out_repro_pair), torch.abs(out_repro_pair_ds), atol=eps))

    def test_compare_evoformer_bf16(self):
        """Run evoformer comparison test with BF16 precision."""
        self.compare_evoformer(torch.bfloat16)

    def test_compare_evoformer_fp32(self):
        """Run evoformer comparison test with FP32 precision."""
        self.compare_evoformer(torch.float32)

    def test_compare_model(self):
        """
        Run full model with and without using DeepSpeed Evoformer attention kernel
        and compare output coordinates
        """
        eps = 2e-2
        with open("tests/test_data/sample_feats.pickle", "rb") as fp:
            batch = pickle.load(fp)

        # atom37_to_atom14 doesn't like batches
        batch["residx_atom14_to_atom37"] = batch["residx_atom14_to_atom37"][0]
        batch["atom14_atom_exists"] = batch["atom14_atom_exists"][0]

        batch["no_recycling_iters"] = np.array([3., 3., 3., 3., ])
        batch = {k: torch.as_tensor(v).cuda() for k, v in batch.items()}

        batch["aatype"] = batch["aatype"].long()
        batch["template_aatype"] = batch["template_aatype"].long()
        batch["extra_msa"] = batch["extra_msa"].long()
        batch["residx_atom37_to_atom14"] = batch[
            "residx_atom37_to_atom14"
        ].long()
        batch["template_all_atom_mask"] = batch["template_all_atom_masks"]
        batch.update(
            data_transforms.atom37_to_torsion_angles("template_")(batch)
        )

        # Move the recycling dimension to the end
        move_dim = lambda t: t.permute(*range(len(t.shape))[1:], 0)
        batch = tensor_tree_map(move_dim, batch)

        with torch.no_grad():
            model = compare_utils.get_global_pretrained_openfold()
            out_repro = model(batch)

            # Enable kernel
            model.globals.use_deepspeed_evo_attention = False
            out_repro_ds = model(batch)

            out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)
            out_repro_ds = tensor_tree_map(lambda t: t.cpu(), out_repro_ds)

            out_repro = out_repro["sm"]["positions"][-1].squeeze(0)
            out_repro_ds = out_repro_ds["sm"]["positions"][-1].squeeze(0)

            self.assertTrue(torch.max(torch.abs(out_repro - out_repro_ds)) < eps)


if __name__ == "__main__":
    unittest.main()
