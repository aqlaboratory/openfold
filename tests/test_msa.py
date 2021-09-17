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
from alphafold.model.msa import *


class TestMSARowAttentionWithPairBias(unittest.TestCase):
    def test_shape(self): 
        batch_size = 2
        s_t = 3
        n = 5
        c_m = 7
        c_z = 11
        c = 52
        no_heads = 4

        mrapb = MSARowAttentionWithPairBias(c_m, c_z, c, no_heads)

        m = torch.rand((batch_size, s_t, n, c_m))
        z = torch.rand((batch_size, n, n, c_z))

        shape_before = m.shape
        m = mrapb(m, z)
        shape_after = m.shape

        self.assertTrue(shape_before == shape_after)


class TestMSAColumnAttention(unittest.TestCase):
    def test_shape(self): 
        batch_size = 2
        s_t = 3
        n = 5
        c_m = 7
        c = 44
        no_heads = 4

        msaca = MSAColumnAttention(c_m, c, no_heads)

        x = torch.rand((batch_size, s_t, n, c_m))

        shape_before = x.shape
        x = msaca(x)
        shape_after = x.shape

        self.assertTrue(shape_before == shape_after)


class TestMSAColumnGlobalAttention(unittest.TestCase):
    def test_shape(self): 
        batch_size = 2
        s_t = 3
        n = 5
        c_m = 7
        c = 44
        no_heads = 4

        msagca = MSAColumnGlobalAttention(c_m, c, no_heads)

        x = torch.rand((batch_size, s_t, n, c_m))

        shape_before = x.shape
        x = msagca(x)
        shape_after = x.shape

        self.assertTrue(shape_before == shape_after)


if __name__ == "__main__":
    unittest.main()
