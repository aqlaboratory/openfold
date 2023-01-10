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
from openfold.model.triangular_multiplicative_update import *
from openfold.utils.tensor_utils import tree_map
import tests.compare_utils as compare_utils
from tests.config import consts

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestTriangularMultiplicativeUpdate(unittest.TestCase):
    def test_shape(self):
        c_z = consts.c_z
        c = 11

        tm = TriangleMultiplicationOutgoing(
            c_z,
            c,
        )

        n_res = consts.c_z
        batch_size = consts.batch_size

        x = torch.rand((batch_size, n_res, n_res, c_z))
        mask = torch.randint(0, 2, size=(batch_size, n_res, n_res))
        shape_before = x.shape
        x = tm(x, mask)
        shape_after = x.shape

        self.assertTrue(shape_before == shape_after)

    def _tri_mul_compare(self, incoming=False):
        name = "triangle_multiplication_" + (
            "incoming" if incoming else "outgoing"
        )

        def run_tri_mul(pair_act, pair_mask):
            config = compare_utils.get_alphafold_config()
            c_e = config.model.embeddings_and_evoformer.evoformer
            tri_mul = alphafold.model.modules.TriangleMultiplication(
                c_e.triangle_multiplication_incoming
                if incoming
                else c_e.triangle_multiplication_outgoing,
                config.model.global_config,
                name=name,
            )
            act = tri_mul(act=pair_act, mask=pair_mask)
            return act

        f = hk.transform(run_tri_mul)

        n_res = consts.n_res

        pair_act = np.random.rand(n_res, n_res, consts.c_z).astype(np.float32)
        pair_mask = np.random.randint(low=0, high=2, size=(n_res, n_res))
        pair_mask = pair_mask.astype(np.float32)

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
            model.evoformer.blocks[0].core.tri_mul_in
            if incoming
            else model.evoformer.blocks[0].core.tri_mul_out
        )
        out_repro = module(
            torch.as_tensor(pair_act, dtype=torch.float32).cuda(),
            mask=torch.as_tensor(pair_mask, dtype=torch.float32).cuda(),
            _inplace=True, _inplace_chunk_size=4,
        ).cpu()

        self.assertTrue(torch.mean(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_tri_mul_out_compare(self):
        self._tri_mul_compare()

    @compare_utils.skip_unless_alphafold_installed()
    def test_tri_mul_in_compare(self):
        self._tri_mul_compare(incoming=True)

    def _tri_mul_inplace(self, incoming=False):
        n_res = consts.n_res
        
        pair_act = np.random.rand(n_res, n_res, consts.c_z).astype(np.float32)
        pair_mask = np.random.randint(low=0, high=2, size=(n_res, n_res))
        pair_mask = pair_mask.astype(np.float32)


        model = compare_utils.get_global_pretrained_openfold()
        module = (
            model.evoformer.blocks[0].core.tri_mul_in
            if incoming
            else model.evoformer.blocks[0].core.tri_mul_out
        )
        out_stock = module(
            torch.as_tensor(pair_act, dtype=torch.float32).cuda(),
            mask=torch.as_tensor(pair_mask, dtype=torch.float32).cuda(),
            _inplace=False,
        ).cpu()
        
        # This has to come second because inference mode is in-place
        out_inplace = module(
            torch.as_tensor(pair_act, dtype=torch.float32).cuda(),
            mask=torch.as_tensor(pair_mask, dtype=torch.float32).cuda(),
            _inplace=True, _inplace_chunk_size=2,
        ).cpu()

        self.assertTrue(torch.mean(torch.abs(out_stock - out_inplace)) < consts.eps)

    def test_tri_mul_out_inference(self):
        self._tri_mul_inplace()

    def test_tri_mul_in_inference(self):
        self._tri_mul_inplace(incoming=True)

if __name__ == "__main__":
    unittest.main()
