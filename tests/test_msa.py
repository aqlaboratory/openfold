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
from openfold.model.msa import (
    MSARowAttentionWithPairBias,
    MSAColumnAttention,
    MSAColumnGlobalAttention,
)
from openfold.utils.tensor_utils import tree_map
import tests.compare_utils as compare_utils
from tests.config import consts

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestMSARowAttentionWithPairBias(unittest.TestCase):
    def test_shape(self):
        batch_size = consts.batch_size
        n_seq = consts.n_seq
        n_res = consts.n_res
        c_m = consts.c_m
        c_z = consts.c_z
        c = 52
        no_heads = 4
        chunk_size = None

        mrapb = MSARowAttentionWithPairBias(c_m, c_z, c, no_heads)

        m = torch.rand((batch_size, n_seq, n_res, c_m))
        z = torch.rand((batch_size, n_res, n_res, c_z))

        shape_before = m.shape
        m = mrapb(m, z=z, chunk_size=chunk_size)
        shape_after = m.shape

        self.assertTrue(shape_before == shape_after)

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def run_msa_row_att(msa_act, msa_mask, pair_act):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            msa_row = alphafold.model.modules.MSARowAttentionWithPairBias(
                c_e.msa_row_attention_with_pair_bias, config.model.global_config
            )
            act = msa_row(msa_act=msa_act, msa_mask=msa_mask, pair_act=pair_act)
            return act

        f = hk.transform(run_msa_row_att)

        n_res = consts.n_res
        n_seq = consts.n_seq

        msa_act = np.random.rand(n_seq, n_res, consts.c_m).astype(np.float32)
        msa_mask = np.random.randint(low=0, high=2, size=(n_seq, n_res)).astype(
            np.float32
        )
        pair_act = np.random.rand(n_res, n_res, consts.c_z).astype(np.float32)

        # Fetch pretrained parameters (but only from one block)]
        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/evoformer_iteration/"
            + "msa_row_attention"
        )
        params = tree_map(lambda n: n[0], params, jax.numpy.DeviceArray)

        out_gt = f.apply(
            params, None, msa_act, msa_mask, pair_act
        ).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt))

        model = compare_utils.get_global_pretrained_openfold()
        out_repro = (
            model.evoformer.blocks[0].msa_att_row(
                torch.as_tensor(msa_act).cuda(),
                z=torch.as_tensor(pair_act).cuda(),
                chunk_size=4,
                mask=torch.as_tensor(msa_mask).cuda(),
            )
        ).cpu()

        self.assertTrue(torch.all(torch.abs(out_gt - out_repro) < consts.eps))


class TestMSAColumnAttention(unittest.TestCase):
    def test_shape(self):
        batch_size = consts.batch_size
        n_seq = consts.n_seq
        n_res = consts.n_res
        c_m = consts.c_m
        c = 44
        no_heads = 4

        msaca = MSAColumnAttention(c_m, c, no_heads)

        x = torch.rand((batch_size, n_seq, n_res, c_m))

        shape_before = x.shape
        x = msaca(x, chunk_size=None)
        shape_after = x.shape

        self.assertTrue(shape_before == shape_after)

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def run_msa_col_att(msa_act, msa_mask):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            msa_col = alphafold.model.modules.MSAColumnAttention(
                c_e.msa_column_attention, config.model.global_config
            )
            act = msa_col(msa_act=msa_act, msa_mask=msa_mask)
            return act

        f = hk.transform(run_msa_col_att)

        n_res = consts.n_res
        n_seq = consts.n_seq

        msa_act = np.random.rand(n_seq, n_res, consts.c_m).astype(np.float32)
        msa_mask = np.random.randint(low=0, high=2, size=(n_seq, n_res)).astype(
            np.float32
        )

        # Fetch pretrained parameters (but only from one block)]
        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/evoformer_iteration/"
            + "msa_column_attention"
        )
        params = tree_map(lambda n: n[0], params, jax.numpy.DeviceArray)

        out_gt = f.apply(params, None, msa_act, msa_mask).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt))

        model = compare_utils.get_global_pretrained_openfold()
        out_repro = (
            model.evoformer.blocks[0].msa_att_col(
                torch.as_tensor(msa_act).cuda(),
                chunk_size=4,
                mask=torch.as_tensor(msa_mask).cuda(),
            )
        ).cpu()

        print(torch.mean(torch.abs(out_gt - out_repro)))

        self.assertTrue(torch.all(torch.abs(out_gt - out_repro) < consts.eps))


class TestMSAColumnGlobalAttention(unittest.TestCase):
    def test_shape(self):
        batch_size = consts.batch_size
        n_seq = consts.n_seq
        n_res = consts.n_res
        c_m = consts.c_m
        c = 44
        no_heads = 4

        msagca = MSAColumnGlobalAttention(c_m, c, no_heads)

        x = torch.rand((batch_size, n_seq, n_res, c_m))

        shape_before = x.shape
        x = msagca(x, chunk_size=None)
        shape_after = x.shape

        self.assertTrue(shape_before == shape_after)

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def run_msa_col_global_att(msa_act, msa_mask):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            msa_col = alphafold.model.modules.MSAColumnGlobalAttention(
                c_e.msa_column_attention,
                config.model.global_config,
                name="msa_column_global_attention",
            )
            act = msa_col(msa_act=msa_act, msa_mask=msa_mask)
            return act

        f = hk.transform(run_msa_col_global_att)

        n_res = consts.n_res
        n_seq = consts.n_seq
        c_e = consts.c_e

        msa_act = np.random.rand(n_seq, n_res, c_e)
        msa_mask = np.random.randint(low=0, high=2, size=(n_seq, n_res))

        # Fetch pretrained parameters (but only from one block)]
        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/extra_msa_stack/"
            + "msa_column_global_attention"
        )
        params = tree_map(lambda n: n[0], params, jax.numpy.DeviceArray)

        out_gt = f.apply(params, None, msa_act, msa_mask).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt.block_until_ready()))

        model = compare_utils.get_global_pretrained_openfold()
        out_repro = (
            model.extra_msa_stack.blocks[0].msa_att_col(
                torch.as_tensor(msa_act, dtype=torch.float32).cuda(),
                chunk_size=4,
                mask=torch.as_tensor(msa_mask, dtype=torch.float32).cuda(),
            )
            .cpu()
        )

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro) < consts.eps))


if __name__ == "__main__":
    unittest.main()
