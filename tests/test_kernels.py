#!/usr/bin/env python
# -*- coding: utf-8 -*-
import torch
import unittest

from openfold.model.primitives import _attention
from openfold.utils.kernel.attention_core import attention_core
from tests.config import consts


class TestAttentionCore(unittest.TestCase):
    def test_attention_core_forward(self):
        n_res = consts.n_res
        h = consts.n_heads_extra_msa
        n_seq = consts.n_extra
        c = consts.c_e
        dtype = torch.float32
        
        q = torch.rand([n_seq, h, n_res, c], dtype=dtype).cuda()
        k = torch.rand([n_seq, h, n_res, c], dtype=dtype).cuda()
        v = torch.rand([n_seq, h, n_res, c], dtype=dtype).cuda()
        mask = torch.randint(0, 2, [n_seq, n_res]).cuda()
        mask_bias = (1e9 * mask - 1)[..., None, None, :].to(dtype)
        
        out_repro = attention_core(q, k, v, mask_bias, None)
        out_gt = _attention(q, k, v, [mask_bias])
        
        self.assertTrue(torch.max(torch.abs(out_repro - out_gt)) < consts.eps)

    def test_attention_core_backward(self):
        n_res = consts.n_res
        h = consts.n_heads_extra_msa
        n_seq = consts.n_extra
        c = consts.c_e
        dtype = torch.float32
        
        q = torch.rand(
            [n_seq, h, n_res, c], dtype=dtype, requires_grad=True
        ).cuda()
        k = torch.rand(
            [n_seq, h, n_res, c], dtype=dtype, requires_grad=True
        ).cuda()
        v = torch.rand(
            [n_seq, h, n_res, c], dtype=dtype, requires_grad=True
        ).cuda()
        mask = torch.randint(0, 2, [n_seq, n_res]).cuda()
        mask_bias = (1e9 * mask - 1)[..., None, None, :].to(dtype)
        
        def clone(t):
            t = t.clone()
            if(t.requires_grad):
                t.retain_grad()
            return t

        q_repro = clone(q)
        k_repro = clone(k)
        v_repro = clone(v)
        out_repro = attention_core(
            q_repro, k_repro, v_repro, mask_bias, None
        )

        loss_repro = torch.mean(out_repro)
        loss_repro.backward()
        
        q_gt = clone(q)
        k_gt = clone(k)
        v_gt = clone(v)
        out_gt = _attention(
            q_gt, k_gt, v_gt, [mask_bias]
        )

        loss_gt = torch.mean(out_gt)
        loss_gt.backward()

        pairs = zip([q_repro, k_repro, v_repro], [q_gt, k_gt, v_gt])
        for t_repro, t_gt in pairs:
            self.assertTrue(
                torch.max(torch.abs(t_repro.grad - t_gt.grad)) < consts.eps
            ) 


if __name__ == '__main__':
    unittest.main()

