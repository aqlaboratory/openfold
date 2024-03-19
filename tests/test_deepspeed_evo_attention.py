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


@compare_utils.skip_unless_ds4s_installed()
class TestDeepSpeedKernel(unittest.TestCase):
    def compare_attention_types(self, use_flash=False):
        """Compare attention with and without using DeepSpeed Evoformer kernel."""
        batch_size = consts.batch_size
        n_seq = 18
        n_res = 20
        c_hidden = 32
        no_heads = 4
        eps = 2e-2

        q, kv, mask, biases = random_attention_inputs(batch_size=batch_size,
                                                      n_seq=n_seq,
                                                      n=n_res,
                                                      no_heads=no_heads,
                                                      c_hidden=c_hidden)

        a = Attention(
            c_hidden, c_hidden, c_hidden, c_hidden, no_heads
        ).cuda()

        # Change output params init for testing since they are initialized with 'final' init (zeros)
        # Otherwise both will just return zero.
        with torch.no_grad():
            lecun_normal_init_(a.linear_g.weight)
            lecun_normal_init_(a.linear_o.weight)

            if use_flash:
                biases = [biases[0]]
                flash_mask = mask.reshape(batch_size * n_seq, n_res)
                real_out = a(q, kv, use_flash=True, flash_mask=flash_mask).cpu()
            else:
                real_out = a(q, kv, biases=biases).cpu()

            ds_out = a(q, kv, biases=biases, use_deepspeed_evo_attention=True).cpu()

        err = torch.max(torch.abs(ds_out - real_out))
        self.assertTrue(err < eps, f'Error: {err}')

    def test_ds_kernel_vs_attention_forward(self):
        """Compare regular attention vs. DeepSpeed Evoformer kernel."""
        self.compare_attention_types(use_flash=False)

    @compare_utils.skip_unless_flash_attn_installed()
    def test_ds_kernel_vs_flash_attn_forward(self):
        """Compare Flash Attention vs. DeepSpeed Evoformer kernel."""
        self.compare_attention_types(use_flash=True)

    def test_ds_kernel_vs_attention_backward(self):
        """Compare backward pass for regular attention vs. DeepSpeed Evoformer kernel."""
        batch_size = consts.batch_size
        n_seq = 18
        n_res = 20
        c_hidden = 32
        no_heads = 4
        eps = consts.eps

        q, kv, mask, biases = random_attention_inputs(batch_size=batch_size,
                                                      n_seq=n_seq,
                                                      n=n_res,
                                                      no_heads=no_heads,
                                                      c_hidden=c_hidden,
                                                      requires_grad=True)

        attn = Attention(
            c_hidden, c_hidden, c_hidden, c_hidden, no_heads
        ).cuda()

        with torch.no_grad():
            lecun_normal_init_(attn.linear_g.weight)
            lecun_normal_init_(attn.linear_o.weight)

        def clone(t):
            # Create new params, clone values
            t = t.clone()
            if t.requires_grad:
                t.retain_grad()
            return t

        def init_attn():
            # Create new attention object with same initial weights
            a_clone = Attention(
                c_hidden, c_hidden, c_hidden, c_hidden, no_heads
            ).cuda()

            a_clone.load_state_dict(attn.state_dict())
            return a_clone

        # Clone param values and run attention with DS kernel
        q_repro = clone(q)
        kv_repro = clone(kv)
        biases_repro = [clone(b) for b in biases]

        a_repro = init_attn()
        out_repro = a_repro(q_repro, kv_repro, biases=biases_repro, use_deepspeed_evo_attention=True)
        loss_repro = torch.mean(out_repro)
        loss_repro.backward()

        q_gt = clone(q)
        kv_gt = clone(kv)
        biases_gt = [clone(b) for b in biases]

        # Clone param values and run attention without DS kernel
        a_gt = init_attn()
        out_gt = a_gt(q_gt, kv_gt, biases=biases_gt)

        loss_gt = torch.mean(out_gt)
        loss_gt.backward()

        # Compare the grads of attention inputs
        pairs = zip([q_repro, kv_repro, biases_repro[1]],
                    [q_gt, kv_gt, biases_gt[1]])
        for i, item in enumerate(pairs):
            t_repro, t_gt = item
            err = torch.max(torch.abs(t_repro.grad.cpu() - t_gt.grad.cpu()))
            self.assertTrue(err < eps, f'Error item #{i}: {err}')

        # Compare the grads of model weights
        a_repro_params = dict(a_repro.named_parameters())
        a_gt_params = dict(a_gt.named_parameters())
        for name in a_gt_params.keys():
            t_repro = a_repro_params[name]
            t_gt = a_gt_params[name]
            err = torch.max(torch.abs(t_repro.grad.cpu() - t_gt.grad.cpu()))
            self.assertTrue(err < eps, f'Error item {name}: {err}')

    def compare_evoformer(self, dtype, eps):
        """
        Compare Evoformer output with and without using DeepSpeed Evoformer attention kernel.
        Set dtype to confirm the kernel can be used during both training (BF16) and inference (FP32),
        since the kernel itself can run with either BF16 or FP16 precision.
        """
        n_res = 20
        n_seq = 18
        c_m_shape = (consts.c_m,)
        c_z_shape = (consts.c_z,)

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

            # In practice, layer norms applied later in the network make any
            # kernel rounding errors negligible
            out_repro_msa = F.layer_norm(out_repro_msa, c_m_shape).cpu()
            out_repro_pair = F.layer_norm(out_repro_pair, c_z_shape).cpu()

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
            out_repro_msa_ds = F.layer_norm(out_repro_msa_ds, c_m_shape).cpu()
            out_repro_pair_ds = F.layer_norm(out_repro_pair_ds, c_z_shape).cpu()

            err = torch.mean(torch.abs(out_repro_msa - out_repro_msa_ds))
            self.assertTrue(err < eps, f'MSA Error: {err}')

            err = torch.mean(torch.abs(out_repro_pair - out_repro_pair_ds))
            self.assertTrue(err < eps, f'Pair Error {err}')

    def test_compare_evoformer_bf16(self):
        """Run evoformer comparison test with BF16 precision."""
        self.compare_evoformer(dtype=torch.bfloat16, eps=4e-2)

    def test_compare_evoformer_fp32(self):
        """Run evoformer comparison test with FP32 precision."""
        self.compare_evoformer(dtype=torch.float32, eps=2e-2)

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

            model.globals.use_deepspeed_evo_attention = True
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
        Run full model with and without using DeepSpeed Evoformer attention kernel
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
        with torch.no_grad():
            with torch.cuda.amp.autocast(dtype=torch.bfloat16):
                model = compare_utils.get_global_pretrained_openfold()
                model.globals.use_deepspeed_evo_attention = False
                out_repro = model(batch)

                # Enable kernel
                model.globals.use_deepspeed_evo_attention = True
                out_repro_ds = model(batch)

                out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)
                out_repro_ds = tensor_tree_map(lambda t: t.cpu(), out_repro_ds)

                out_repro = out_repro["sm"]["positions"][-1].squeeze(0)
                out_repro_ds = out_repro_ds["sm"]["positions"][-1].squeeze(0)

                compare_utils.assert_mean_abs_diff_small(out_repro, out_repro_ds, eps)


if __name__ == "__main__":
    unittest.main()
