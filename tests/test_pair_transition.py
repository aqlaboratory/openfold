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
from alphafold.model.pair_transition import *


class TestPairTransition(unittest.TestCase):
    def test_shape(self): 
        c_z = 5
        n = 4

        pt = PairTransition(c_z, n)

        batch_size = 4
        n_res = 256

        z = torch.rand((batch_size, n_res, n_res, c_z))
        mask = torch.randint(0, 2, size=(batch_size, n_res, n_res))
        shape_before = z.shape
        z = pt(z, mask)
        shape_after = z.shape

        self.assertTrue(shape_before == shape_after)


if __name__ == "__main__":
    unittest.main()

