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
    _lma,
    _attention,
    DEFAULT_LMA_Q_CHUNK_SIZE,
    DEFAULT_LMA_KV_CHUNK_SIZE
)
from tests.config import consts


class TestLMA(unittest.TestCase):
    def test_lma_vs_attention(self):
        batch_size = consts.batch_size
        c_hidden = 32 
        n = 2**12
        n_seq = 12
        no_heads = 4

        q = torch.rand(batch_size, n_seq, no_heads, n, c_hidden).cuda()
        k = torch.rand(batch_size, n_seq, no_heads, n, c_hidden).cuda()
        v = torch.rand(batch_size, n_seq, no_heads, n, c_hidden).cuda()

        bias = [torch.rand(batch_size, n_seq, 1, 1, n), torch.rand(batch_size, 1, no_heads, n, n)]
        biases = [b.cuda() for b in bias]
        
        with torch.no_grad():
            lma_biases = [
                b.expand(b.shape[:-2] + (q.shape[-2],) + (k.shape[-2],))
                for b in biases
            ]
            l = _lma(q, k, v, lma_biases, DEFAULT_LMA_Q_CHUNK_SIZE, DEFAULT_LMA_KV_CHUNK_SIZE).cpu()
            real = _attention(q, k, v, biases).cpu()
       
        err = torch.max(torch.abs(l - real))
        self.assertTrue(err < consts.eps, f'Error: {err}')

 
if __name__ == "__main__":
    unittest.main()
