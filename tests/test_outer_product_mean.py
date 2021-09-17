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
from alphafold.model.outer_product_mean import *


class TestOuterProductMean(unittest.TestCase):
    def test_shape(self):
        batch_size = 2
        s = 5
        n_res = 7
        c_m = 11
        c = 13
        c_z = 17

        opm = OuterProductMean(c_m, c_z, c)

        m = torch.rand((batch_size, s, n_res, c_m))
        mask = torch.randint(0, 2, size=(batch_size, s, n_res))
        m = opm(m, mask)

        self.assertTrue(m.shape == (batch_size, n_res, n_res, c_z))


if __name__ == "__main__":
    unittest.main()

