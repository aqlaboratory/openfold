# Copyright 2023 RIKEN & Fujitsu Limited
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

import unittest
import torch
from openfold.model.outer_product_mean import einsum_bac_dae_bdce

#def einsum_bac_dae_bdce(a, b):
#    dim_notation_n = a.shape[:-3]
#    dim_notation_a = a.shape[-2]
#    dim_notation_b = a.shape[-3]
#    dim_notation_c = a.shape[-1]
#    dim_notation_d = b.shape[-3]
#    dim_notation_e = b.shape[-1]
#
#    # ...bca <- ...bac
#    a = a.transpose(-2, -1)
#    # ...ade <- ...dae
#    b = b.transpose(-3, -2)
#
#    # ...xa <- ...bca
#    a = a.reshape(-1, dim_notation_b * dim_notation_c, dim_notation_a)
#    # ...ay <- ...ade
#    b = b.reshape(-1, dim_notation_a, dim_notation_d * dim_notation_e)
#
#    # ...xy <- ...xa * ...ay
#    if a.shape[0] == 1:
#        # mm
#        r = torch.matmul(a[0,...], b[0,...])
#    else:
#        r = torch.empty([a.shape[0]] + [dim_notation_b * dim_notation_c, dim_notation_d * dim_notation_e])
#        for i in range(a.shape[0]):
#            r[i,...] = torch.matmul(a[i,...], b[i,...])
#
#    # ...bcde <- ...xy
#    r = r.reshape(list(dim_notation_n) + [dim_notation_b, dim_notation_c, dim_notation_d, dim_notation_e])
#    # ...bdce <- ...bcde
#    r = r.transpose(-3, -2)
#
#    return r

patterns = [
    [128, 17, 8],
    [1, 128, 17, 8],
    [2, 128, 17, 8],
    [2, 3, 128, 17, 8],
    [128, 17, 1],
    [128, 1, 8],
    [1, 17, 8],
]

class TestEinsum(unittest.TestCase):
    def test_einsum(self):
        for shape in patterns:
            t_a = torch.randn(shape)
            t_b = torch.randn(shape)
            t_a_2 = t_a.clone()
            t_b_2 = t_b.clone()
            t_r = torch.einsum("...bac,...dae->...bdce", t_a, t_b)
            t_r_2 = einsum_bac_dae_bdce(t_a_2, t_b_2)
            self.assertTrue(t_r.shape == t_r_2.shape)
            self.assertTrue(torch.allclose(t_r, t_r_2))
    
if __name__ == '__main__':
    unittest.main()
