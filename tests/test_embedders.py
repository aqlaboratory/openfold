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

import random
import torch
import unittest
from tests.config import consts
from tests.data_utils import random_asym_ids
from openfold.model.embedders import (
    InputEmbedder,
    InputEmbedderMultimer,
    PreembeddingEmbedder,
    RecyclingEmbedder,
    TemplateSingleEmbedder,
    TemplatePairEmbedder
)


class TestInputEmbedder(unittest.TestCase):
    def test_shape(self):
        tf_dim = 2
        msa_dim = 3
        c_z = 5
        c_m = 7
        relpos_k = 11

        b = 13
        n_res = 17
        n_clust = 19

        max_relative_chain = 2
        max_relative_idx = 32
        use_chain_relative = True

        tf = torch.rand((b, n_res, tf_dim))
        ri = torch.rand((b, n_res))
        msa = torch.rand((b, n_clust, n_res, msa_dim))
        asym_ids_flat = torch.Tensor(random_asym_ids(n_res))
        asym_id = torch.tile(asym_ids_flat.unsqueeze(0), (b, 1))
        entity_id = asym_id
        sym_id = torch.zeros_like(entity_id)

        if consts.is_multimer:
            ie = InputEmbedderMultimer(tf_dim, msa_dim, c_z, c_m,
                                       max_relative_idx=max_relative_idx,
                                       use_chain_relative=use_chain_relative,
                                       max_relative_chain=max_relative_chain)
            batch = {"target_feat": tf, "residue_index": ri, "msa_feat": msa,
                     "asym_id": asym_id, "entity_id": entity_id, "sym_id": sym_id}
            msa_emb, pair_emb = ie(batch)
        else:
            ie = InputEmbedder(tf_dim, msa_dim, c_z, c_m, relpos_k)
            msa_emb, pair_emb = ie(tf=tf, ri=ri, msa=msa, inplace_safe=False)

        self.assertTrue(msa_emb.shape == (b, n_clust, n_res, c_m))
        self.assertTrue(pair_emb.shape == (b, n_res, n_res, c_z))


class TestPreembeddingEmbedder(unittest.TestCase):
    def test_shape(self):
        tf_dim = 22
        preembedding_dim = 1280
        c_z = 4
        c_m = 6
        relpos_k = 10

        batch_size = 4
        num_res = 20

        tf = torch.rand((batch_size, num_res, tf_dim))
        ri = torch.rand((batch_size, num_res))
        preemb = torch.rand((batch_size, num_res, preembedding_dim))

        pe = PreembeddingEmbedder(tf_dim, preembedding_dim, c_z, c_m, relpos_k)

        seq_emb, pair_emb = pe(tf, ri, preemb)
        self.assertTrue(seq_emb.shape == (batch_size, 1, num_res, c_m))
        self.assertTrue(pair_emb.shape == (batch_size, num_res, num_res, c_z))


class TestRecyclingEmbedder(unittest.TestCase):
    def test_shape(self):
        batch_size = 2
        n = 3
        c_z = 5
        c_m = 7
        min_bin = 0
        max_bin = 10
        no_bins = 9

        re = RecyclingEmbedder(c_m, c_z, min_bin, max_bin, no_bins)

        m_1 = torch.rand((batch_size, n, c_m))
        z = torch.rand((batch_size, n, n, c_z))
        x = torch.rand((batch_size, n, 3))

        m_1, z = re(m_1, z, x)

        self.assertTrue(z.shape == (batch_size, n, n, c_z))
        self.assertTrue(m_1.shape == (batch_size, n, c_m))


class TestTemplateAngleEmbedder(unittest.TestCase):
    def test_shape(self):
        template_angle_dim = 51
        c_m = 256
        batch_size = 4
        n_templ = 4
        n_res = 256

        tae = TemplateSingleEmbedder(
            template_angle_dim,
            c_m,
        )

        x = torch.rand((batch_size, n_templ, n_res, template_angle_dim))
        x = tae(x)

        self.assertTrue(x.shape == (batch_size, n_templ, n_res, c_m))


class TestTemplatePairEmbedder(unittest.TestCase):
    def test_shape(self):
        batch_size = 2
        n_templ = 3
        n_res = 5
        template_pair_dim = 7
        c_t = 11

        tpe = TemplatePairEmbedder(
            template_pair_dim,
            c_t,
        )

        x = torch.rand((batch_size, n_templ, n_res, n_res, template_pair_dim))
        x = tpe(x)

        self.assertTrue(x.shape == (batch_size, n_templ, n_res, n_res, c_t))


if __name__ == "__main__":
    unittest.main()
