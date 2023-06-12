#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
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

import torch
import unittest

from openfold.model.primitives import _attention
from openfold.utils.kernel.attention_core import attention_core
from tests.config import consts

dtype = torch.float32
device = torch.device('cpu')

n_res = consts.n_res
h = 8
c = consts.c_e

patterns = [
    {'qkv_shape':      [516, h, n_res,     c],
     'bias1_shape':    [516, 1,     1, n_res],
     'bias2_shape':    [  1, h, n_res, n_res],},
    {'qkv_shape':      [516, h, n_res,     c],
     'bias1_shape':    [516, 1,     1, n_res],
     'bias2_shape':    [516, 1,     1, n_res],},
    {'qkv_shape':   [1, 516, h, n_res,     c],
     'bias1_shape': [1, 516, 1,     1, n_res],
     'bias2_shape': [1,   1, h, n_res, n_res],},
    {'qkv_shape':   [1, 516, h, n_res,     c],
     'bias1_shape': [1, 516, 1,     1, n_res],
     'bias2_shape': [1, 516, 1,     1, n_res],},
    ]


class TestAttentionCore(unittest.TestCase):
    def test_attention_core_forward(self):
        for pattern in patterns:
            with self.subTest('test case {}'.format(pattern)):
                q = torch.rand(
                    pattern['qkv_shape'], dtype=dtype, device=device)
                k = torch.rand(
                    pattern['qkv_shape'], dtype=dtype, device=device)
                v = torch.rand(
                    pattern['qkv_shape'], dtype=dtype, device=device)
                mask1 = torch.randint(0, 2, pattern['bias1_shape'], device=device)
                mask2 = torch.randint(0, 2, pattern['bias2_shape'], device=device)
                bias1 = (1e9 * mask1 - 1)
                bias2 = (1e9 * mask2 - 1)
            
                out_repro = attention_core(q, k, v, bias1, bias2)
                out_gt = _attention(q, k, v, [bias1, bias2])

                d = torch.max(torch.abs(out_repro - out_gt))
                print(d, "close={}, d_sum={}".format(torch.allclose(out_repro, out_gt),
                                                     torch.sum(torch.abs(out_repro - out_gt), dtype=torch.float64)))
                self.assertTrue(torch.max(torch.abs(out_repro - out_gt)) < consts.eps)
                self.assertTrue(out_repro.shape == out_gt.shape)

    def test_attention_core_backward(self):
        for pattern in patterns:
            with self.subTest('test case {}'.format(pattern)):
                q = torch.rand(
                    pattern['qkv_shape'], dtype=dtype, requires_grad=True, device=device)
                k = torch.rand(
                    pattern['qkv_shape'], dtype=dtype, requires_grad=True, device=device)
                v = torch.rand(
                    pattern['qkv_shape'], dtype=dtype, requires_grad=True, device=device)
                mask1 = torch.randint(0, 2, pattern['bias1_shape'], device=device)
                mask2 = torch.randint(0, 2, pattern['bias2_shape'], device=device)
                bias1 = (1e9 * mask1 - 1)
                bias2 = (1e9 * mask2 - 1)
                bias1.requires_grad = True
                bias2.requires_grad = True
        
                def clone(t):
                    t = t.clone()
                    if(t.requires_grad):
                        t.retain_grad()
                    return t

                q_repro = clone(q)
                k_repro = clone(k)
                v_repro = clone(v)
                bias1_repro = clone(bias1)
                bias2_repro = clone(bias2)
                out_repro = attention_core(
                    q_repro, k_repro, v_repro, bias1_repro, bias2_repro
                )

                loss_repro = torch.mean(out_repro)
                loss_repro.backward()
        
                q_gt = clone(q)
                k_gt = clone(k)
                v_gt = clone(v)
                bias1_gt = clone(bias1)
                bias2_gt = clone(bias2)
                out_gt = _attention(
                    q_gt, k_gt, v_gt, [bias1_gt, bias2_gt]
                )

                loss_gt = torch.mean(out_gt)
                loss_gt.backward()

                for t_repro, t_gt in zip([q_repro, k_repro, v_repro, bias1_repro, bias2_repro],
                                         [q_gt, k_gt, v_gt, bias1_gt, bias2_gt]):
                    d_grad = torch.max(torch.abs((t_repro.grad - t_gt.grad)))
                    d = torch.max(torch.abs(t_repro - t_gt))
                    print(d, d_grad, "close={}, d_grad_sum={}".format(torch.allclose(t_repro, t_gt), torch.sum(torch.abs(t_repro.grad - t_gt.grad), dtype=torch.float64)))

                    self.assertTrue(
                        torch.max(torch.abs(t_repro.grad - t_gt.grad)) < consts.eps
                    ) 
                    self.assertTrue(
                        torch.max(torch.abs(t_repro - t_gt)) < consts.eps
                    )
                    self.assertTrue(t_repro.shape == t_gt.shape)

if __name__ == '__main__':
    unittest.main()

