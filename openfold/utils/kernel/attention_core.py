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
import importlib
from functools import reduce
from operator import mul

import torch

attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")


SUPPORTED_DTYPES = [torch.float32, torch.bfloat16]


class AttentionCoreFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, q, k, v, bias_1=None, bias_2=None):
        if(bias_1 is None and bias_2 is not None):
            raise ValueError("bias_1 must be specified before bias_2")
        if(q.dtype not in SUPPORTED_DTYPES):
            raise ValueError("Unsupported datatype")

        q = q.contiguous()
        k = k.contiguous()

        # [*, H, Q, K] 
        attention_logits = torch.matmul(
            q, k.transpose(-1, -2), 
        )

        if(bias_1 is not None):
            attention_logits += bias_1
        if(bias_2 is not None):
            attention_logits += bias_2

        attn_core_inplace_cuda.forward_(
            attention_logits, 
            reduce(mul, attention_logits.shape[:-1]),
            attention_logits.shape[-1],
        )

        o = torch.matmul(attention_logits, v) 

        ctx.bias_1_shape = bias_1.shape if bias_1 is not None else None
        ctx.bias_2_shape = bias_2.shape if bias_2 is not None else None
        ctx.save_for_backward(q, k, v, attention_logits)

        return o

    @staticmethod
    def backward(ctx, grad_output):
        q, k, v, attention_logits = ctx.saved_tensors
        grad_q = grad_k = grad_v = grad_bias_1 = grad_bias_2 = None
       
        grad_v = torch.matmul(
            attention_logits.transpose(-1, -2), 
            grad_output
        )

        attn_core_inplace_cuda.backward_(
            attention_logits,
            grad_output.contiguous(),
            v.contiguous(), # v is implicitly transposed in the kernel
            reduce(mul, attention_logits.shape[:-1]),
            attention_logits.shape[-1],
            grad_output.shape[-1],
        )

        if(ctx.bias_1_shape is not None):
            grad_bias_1 = torch.sum(
                attention_logits,
                dim=tuple(i for i,d in enumerate(ctx.bias_1_shape) if d == 1),
                keepdim=True,
            )

        if(ctx.bias_2_shape is not None):
            grad_bias_2 = torch.sum(
                attention_logits,
                dim=tuple(i for i,d in enumerate(ctx.bias_2_shape) if d == 1),
                keepdim=True,
            )

        grad_q = torch.matmul(
            attention_logits, k
        )
        grad_k = torch.matmul(
            q.transpose(-1, -2), attention_logits,
        ).transpose(-1, -2)

        return grad_q, grad_k, grad_v, grad_bias_1, grad_bias_2

attention_core = AttentionCoreFunction.apply
