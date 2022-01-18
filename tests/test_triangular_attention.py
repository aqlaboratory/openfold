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
from openfold.model.triangular_attention import TriangleAttention
from openfold.utils.tensor_utils import tree_map

import tests.compare_utils as compare_utils
from tests.config import consts

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestTriangularAttention(unittest.TestCase):
    def test_shape(self):
        c_z = consts.c_z
        c = 12
        no_heads = 4
        starting = True

        tan = TriangleAttention(c_z, c, no_heads, starting)

        batch_size = consts.batch_size
        n_res = consts.n_res

        x = torch.rand((batch_size, n_res, n_res, c_z))
        shape_before = x.shape
        x = tan(x, chunk_size=None)
        shape_after = x.shape

        self.assertTrue(shape_before == shape_after)

    def _tri_att_compare(self, starting=False):
        name = (
            "triangle_attention_"
            + ("starting" if starting else "ending")
            + "_node"
        )

        def run_tri_att(pair_act, pair_mask):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            tri_att = alphafold.model.modules.TriangleAttention(
                c_e.triangle_attention_starting_node
                if starting
                else c_e.triangle_attention_ending_node,
                config.model.global_config,
                name=name,
            )
            act = tri_att(pair_act=pair_act, pair_mask=pair_mask)
            return act

        f = hk.transform(run_tri_att)

        n_res = consts.n_res

        pair_act = np.random.rand(n_res, n_res, consts.c_z) * 100
        pair_mask = np.random.randint(low=0, high=2, size=(n_res, n_res))

        # Fetch pretrained parameters (but only from one block)]
        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/evoformer_iteration/"
            + name
        )
        params = tree_map(lambda n: n[0], params, jax.numpy.DeviceArray)

        out_gt = f.apply(params, None, pair_act, pair_mask).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt))

        model = compare_utils.get_global_pretrained_openfold()
        module = (
            model.evoformer.blocks[0].core.tri_att_start
            if starting
            else model.evoformer.blocks[0].core.tri_att_end
        )
        out_repro = module(
            torch.as_tensor(pair_act, dtype=torch.float32).cuda(),
            mask=torch.as_tensor(pair_mask, dtype=torch.float32).cuda(),
            chunk_size=None,
        ).cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_tri_att_end_compare(self):
        self._tri_att_compare()

    @compare_utils.skip_unless_alphafold_installed()
    def test_tri_att_start_compare(self):
        self._tri_att_compare(starting=True)


if __name__ == "__main__":
    unittest.main()
