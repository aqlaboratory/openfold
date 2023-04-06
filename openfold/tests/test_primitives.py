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

from openfold.model.primitives import (
    Attention,
)
from tests.config import consts


class TestLMA(unittest.TestCase):
    def test_lma_vs_attention(self):
        batch_size = consts.batch_size
        c_hidden = 32 
        n = 2**12
        no_heads = 4

        q = torch.rand(batch_size, n, c_hidden).cuda()
        kv = torch.rand(batch_size, n, c_hidden).cuda()

        bias = [torch.rand(no_heads, 1, n)]
        bias = [b.cuda() for b in bias]
        
        gating_fill = torch.rand(c_hidden * no_heads, c_hidden)
        o_fill = torch.rand(c_hidden, c_hidden * no_heads)
        
        a = Attention(
            c_hidden, c_hidden, c_hidden, c_hidden, no_heads
        ).cuda()
        
        with torch.no_grad():
            l = a(q, kv, biases=bias, use_lma=True)
            real = a(q, kv, biases=bias)
        
        self.assertTrue(torch.max(torch.abs(l - real)) < consts.eps)

 
if __name__ == "__main__":
    unittest.main()
