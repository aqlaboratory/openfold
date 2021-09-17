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
from alphafold.model.evoformer import *


class TestEvoformerStack(unittest.TestCase):
    def test_shape(self):  
        batch_size = 5
        s_t = 27
        n_res = 29
        c_m = 7
        c_z = 11
        c_hidden_msa_att = 12
        c_hidden_opm = 17
        c_hidden_mul = 19
        c_hidden_pair_att = 14
        c_s = 23
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
            chunk_size=4,
            inf=inf,
            eps=eps,
        ).eval()

        m = torch.rand((batch_size, s_t, n_res, c_m))
        z = torch.rand((batch_size, n_res, n_res, c_z))
        msa_mask = torch.randint(0, 2, size=(batch_size, s_t, n_res))
        pair_mask = torch.randint(0, 2, size=(batch_size, n_res, n_res))

        shape_m_before = m.shape
        shape_z_before = z.shape

        m, z, s = es(m, z, msa_mask, pair_mask)

        self.assertTrue(m.shape == shape_m_before)
        self.assertTrue(z.shape == shape_z_before)
        self.assertTrue(s.shape == (batch_size, n_res, c_s))


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
            blocks_per_ckpt=None,
            chunk_size=4,
            inf=inf,
            eps=eps,
        ).eval()

        m = torch.rand((batch_size, s_t, n_res, c_m))
        z = torch.rand((batch_size, n_res, n_res, c_z))
        msa_mask = torch.randint(0, 2, size=(batch_size, s_t, n_res,))
        pair_mask = torch.randint(0, 2, size=(batch_size, n_res, n_res,))

        shape_z_before = z.shape

        z = es(m, z, msa_mask, pair_mask)

        self.assertTrue(z.shape == shape_z_before)


class TestMSATransition(unittest.TestCase):
    def test_shape(self): 
        batch_size = 2
        s_t = 3
        n_r = 5
        c_m = 7
        n = 11

        mt = MSATransition(c_m, n, chunk_size=4)

        m = torch.rand((batch_size, s_t, n_r, c_m))

        shape_before = m.shape
        m = mt(m)
        shape_after = m.shape

        self.assertTrue(shape_before == shape_after)


if __name__ == "__main__":
    unittest.main()
