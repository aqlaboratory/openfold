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
import unittest

from openfold.model.primitives import (
    lecun_normal_init_,
    Attention,
)
from tests.config import consts
from tests.data_utils import random_attention_inputs


class TestLMA(unittest.TestCase):
    def test_lma_vs_attention(self):
        c_hidden = 32
        no_heads = 4

        q, kv, _, biases = random_attention_inputs(batch_size=consts.batch_size,
                                                   n_seq=consts.n_seq,
                                                   n=2 ** 12,
                                                   no_heads=no_heads,
                                                   c_hidden=c_hidden)

        a = Attention(
            c_hidden, c_hidden, c_hidden, c_hidden, no_heads
        ).cuda()

        with torch.no_grad():
            lecun_normal_init_(a.linear_g.weight)
            lecun_normal_init_(a.linear_o.weight)

            l = a(q, kv, biases=biases, use_lma=True).cpu()
            real = a(q, kv, biases=biases).cpu()

        err = torch.max(torch.abs(l - real))
        self.assertTrue(err < consts.eps, f'Error: {err}')


if __name__ == "__main__":
    unittest.main()