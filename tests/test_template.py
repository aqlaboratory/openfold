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
from alphafold.model.template import *


class TestTemplatePointwiseAttention(unittest.TestCase):
    def test_shape(self): 
        batch_size = 2
        s_t = 3
        c_t = 5
        c_z = 7
        c = 26
        no_heads = 13
        n = 17
        
        tpa = TemplatePointwiseAttention(c_t, c_z, c, no_heads, chunk_size=4)

        t = torch.rand((batch_size, s_t, n, n, c_t))
        z = torch.rand((batch_size, n, n, c_z))

        z_update = tpa(t, z)

        self.assertTrue(z_update.shape == z.shape)


class TestTemplatePairStack(unittest.TestCase):
    def test_shape(self):
        batch_size = 2
        c_t = 5
        c_hidden_tri_att = 7
        c_hidden_tri_mul = 7
        no_blocks = 2
        no_heads = 4
        pt_inner_dim = 15
        dropout = 0.25
        n_templ = 3
        n_res = 5
        chunk_size = 4

        tpe = TemplatePairStack(
            c_t,
            c_hidden_tri_att=c_hidden_tri_att,
            c_hidden_tri_mul=c_hidden_tri_mul,
            no_blocks=no_blocks,
            no_heads=no_heads,
            pair_transition_n=pt_inner_dim,
            dropout_rate=dropout,
            chunk_size=chunk_size,
        )

        t = torch.rand((batch_size, n_templ, n_res, n_res, c_t))
        mask = torch.randint(0, 2, (batch_size, n_templ, n_res, n_res))
        shape_before = t.shape
        t = tpe(t, mask)
        shape_after = t.shape

        self.assertTrue(shape_before == shape_after)




if __name__ == "__main__":
    unittest.main() 
