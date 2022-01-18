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

import torch
import numpy as np
import unittest
from openfold.model.evoformer import (
    MSATransition,
    EvoformerStack,
    ExtraMSAStack,
)
from openfold.utils.tensor_utils import tree_map
import tests.compare_utils as compare_utils
from tests.config import consts

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestEvoformerStack(unittest.TestCase):
    def test_shape(self):
        batch_size = consts.batch_size
        n_seq = consts.n_seq
        n_res = consts.n_res
        c_m = consts.c_m
        c_z = consts.c_z
        c_hidden_msa_att = 12
        c_hidden_opm = 17
        c_hidden_mul = 19
        c_hidden_pair_att = 14
        c_s = consts.c_s
        no_heads_msa = 3
        no_heads_pair = 7
        no_blocks = 2
        transition_n = 2
        msa_dropout = 0.15
        pair_stack_dropout = 0.25
        inf = 1e9
        eps = 1e-10

        es = EvoformerStack(
            c_m,
            c_z,
            c_hidden_msa_att,
            c_hidden_opm,
            c_hidden_mul,
            c_hidden_pair_att,
            c_s,
            no_heads_msa,
            no_heads_pair,
            no_blocks,
            transition_n,
            msa_dropout,
            pair_stack_dropout,
            blocks_per_ckpt=None,
            inf=inf,
            eps=eps,
        ).eval()

        m = torch.rand((batch_size, n_seq, n_res, c_m))
        z = torch.rand((batch_size, n_res, n_res, c_z))
        msa_mask = torch.randint(0, 2, size=(batch_size, n_seq, n_res))
        pair_mask = torch.randint(0, 2, size=(batch_size, n_res, n_res))

        shape_m_before = m.shape
        shape_z_before = z.shape

        m, z, s = es(
            m, z, chunk_size=4, msa_mask=msa_mask, pair_mask=pair_mask
        )

        self.assertTrue(m.shape == shape_m_before)
        self.assertTrue(z.shape == shape_z_before)
        self.assertTrue(s.shape == (batch_size, n_res, c_s))

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def run_ei(activations, masks):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            ei = alphafold.model.modules.EvoformerIteration(
                c_e, config.model.global_config, is_extra_msa=False
            )
            return ei(activations, masks, is_training=False)

        f = hk.transform(run_ei)

        n_res = consts.n_res
        n_seq = consts.n_seq

        activations = {
            "msa": np.random.rand(n_seq, n_res, consts.c_m).astype(np.float32),
            "pair": np.random.rand(n_res, n_res, consts.c_z).astype(np.float32),
        }

        masks = {
            "msa": np.random.randint(0, 2, (n_seq, n_res)).astype(np.float32),
            "pair": np.random.randint(0, 2, (n_res, n_res)).astype(np.float32),
        }

        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/evoformer_iteration"
        )
        params = tree_map(lambda n: n[0], params, jax.numpy.DeviceArray)

        key = jax.random.PRNGKey(42)
        out_gt = f.apply(params, key, activations, masks)
        jax.tree_map(lambda x: x.block_until_ready(), out_gt)
        out_gt_msa = torch.as_tensor(np.array(out_gt["msa"]))
        out_gt_pair = torch.as_tensor(np.array(out_gt["pair"]))

        model = compare_utils.get_global_pretrained_openfold()
        out_repro_msa, out_repro_pair = model.evoformer.blocks[0](
            torch.as_tensor(activations["msa"]).cuda(),
            torch.as_tensor(activations["pair"]).cuda(),
            torch.as_tensor(masks["msa"]).cuda(),
            torch.as_tensor(masks["pair"]).cuda(),
            chunk_size=4,
            _mask_trans=False,
        )

        out_repro_msa = out_repro_msa.cpu()
        out_repro_pair = out_repro_pair.cpu()

        assert(torch.max(torch.abs(out_repro_msa - out_gt_msa)) < consts.eps)
        assert(torch.max(torch.abs(out_repro_pair - out_gt_pair)) < consts.eps)


class TestExtraMSAStack(unittest.TestCase):
    def test_shape(self):
        batch_size = 2
        s_t = 23
        n_res = 5
        c_m = 7
        c_z = 11
        c_hidden_msa_att = 12
        c_hidden_opm = 17
        c_hidden_mul = 19
        c_hidden_tri_att = 16
        no_heads_msa = 3
        no_heads_pair = 8
        no_blocks = 2
        transition_n = 5
        msa_dropout = 0.15
        pair_stack_dropout = 0.25
        inf = 1e9
        eps = 1e-10

        es = ExtraMSAStack(
            c_m,
            c_z,
            c_hidden_msa_att,
            c_hidden_opm,
            c_hidden_mul,
            c_hidden_tri_att,
            no_heads_msa,
            no_heads_pair,
            no_blocks,
            transition_n,
            msa_dropout,
            pair_stack_dropout,
            ckpt=False,
            inf=inf,
            eps=eps,
        ).eval()

        m = torch.rand((batch_size, s_t, n_res, c_m))
        z = torch.rand((batch_size, n_res, n_res, c_z))
        msa_mask = torch.randint(
            0,
            2,
            size=(
                batch_size,
                s_t,
                n_res,
            ),
        )
        pair_mask = torch.randint(
            0,
            2,
            size=(
                batch_size,
                n_res,
                n_res,
            ),
        )

        shape_z_before = z.shape

        z = es(m, z, chunk_size=4, msa_mask=msa_mask, pair_mask=pair_mask)

        self.assertTrue(z.shape == shape_z_before)


class TestMSATransition(unittest.TestCase):
    def test_shape(self):
        batch_size = 2
        s_t = 3
        n_r = 5
        c_m = 7
        n = 11

        mt = MSATransition(c_m, n)

        m = torch.rand((batch_size, s_t, n_r, c_m))

        shape_before = m.shape
        m = mt(m, chunk_size=4)
        shape_after = m.shape

        self.assertTrue(shape_before == shape_after)

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def run_msa_transition(msa_act, msa_mask):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            msa_trans = alphafold.model.modules.Transition(
                c_e.msa_transition,
                config.model.global_config,
                name="msa_transition",
            )
            act = msa_trans(act=msa_act, mask=msa_mask)
            return act

        f = hk.transform(run_msa_transition)

        n_res = consts.n_res
        n_seq = consts.n_seq

        msa_act = np.random.rand(n_seq, n_res, consts.c_m).astype(np.float32)
        msa_mask = np.ones((n_seq, n_res)).astype(
            np.float32
        )  # no mask here either

        # Fetch pretrained parameters (but only from one block)]
        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/evoformer_iteration/"
            + "msa_transition"
        )
        params = tree_map(lambda n: n[0], params, jax.numpy.DeviceArray)

        out_gt = f.apply(params, None, msa_act, msa_mask).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt))

        model = compare_utils.get_global_pretrained_openfold()
        
        out_repro = (
            model.evoformer.blocks[0].core.msa_transition(
                torch.as_tensor(msa_act, dtype=torch.float32).cuda(),
                mask=torch.as_tensor(msa_mask, dtype=torch.float32).cuda(),
            )
            .cpu()
        )

        print(out_gt)
        print(out_repro)
        
        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)


if __name__ == "__main__":
    unittest.main()
