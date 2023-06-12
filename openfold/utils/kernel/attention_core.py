# Copyright 2023 RIKEN & Fujitsu Limited
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

if torch.cuda.is_available():
    attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")
else:
    import openfold.utils.kernel.attention_core_cpu as attn_core_inplace_cuda


SUPPORTED_DTYPES = [torch.float32, torch.bfloat16]

_recompute = False
chunk_size = 256
seq_dim = -4

def split_chunks(tensor, dim_size, chunk_size, dim):
    n_chunks = (dim_size - 1) // chunk_size + 1
    if tensor is None:
        return [None] * n_chunks
    elif tensor.shape[dim] == 1:
        return [tensor] * n_chunks
    else:
        return torch.split(tensor, chunk_size, dim=dim)

class AttentionCoreFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, q, k, v, bias_1=None, bias_2=None):
        if(bias_1 is None and bias_2 is not None):
            raise ValueError("bias_1 must be specified before bias_2")
        if(q.dtype not in SUPPORTED_DTYPES):
            raise ValueError("Unsupported datatype")

        # torch.Size([1, 5120, 8, 384, 8]) torch.Size([1, 5120, 8, 384, 8]) torch.Size([1, 5120, 8, 384, 8]) torch.Size([1, 5120, 1, 1, 384]) torch.Size([1, 1, 8, 384, 384])
        q = q.contiguous()
        k = k.contiguous()

        if _recompute:
            attention_logits = None
        else:
            attention_logits = torch.empty(list(q.shape[:-1] + (q.shape[-2],)), dtype=q.dtype)

        o = torch.empty_like(q)

        for q_c, k_c, v_c, b1_c, b2_c, o_c, al_c in zip(torch.split(q, chunk_size, dim=seq_dim),
                                                        torch.split(k, chunk_size, dim=seq_dim),
                                                        torch.split(v, chunk_size, dim=seq_dim),
                                                        split_chunks(bias_1, q.shape[seq_dim], chunk_size, dim=seq_dim),
                                                        split_chunks(bias_2, q.shape[seq_dim], chunk_size, dim=seq_dim),
                                                        torch.split(o, chunk_size, dim=seq_dim),
                                                        split_chunks(attention_logits, q.shape[seq_dim], chunk_size, dim=seq_dim)):
            if al_c is None:
                al_c = torch.empty(list(q_c.shape[:-1] + (q_c.shape[-2],)), dtype=q_c.dtype)

            # [*, H, Q, K]
            torch.matmul(
                q_c, k_c.transpose(-1, -2), out=al_c
            )

            if(b1_c is not None):
                al_c += b1_c
            if(b2_c is not None):
                al_c += b2_c

            attn_core_inplace_cuda.forward_(
                al_c,
                reduce(mul, al_c.shape[:-1]),
                al_c.shape[-1],
            )

            # error by specifying non-contiguous o_c as out argument
            tmp = torch.matmul(al_c, v_c)
            o_c[...] = tmp[...]

        ctx.bias_1_shape = bias_1.shape if bias_1 is not None else None
        ctx.bias_2_shape = bias_2.shape if bias_2 is not None else None

        if _recompute:
            ctx.save_for_backward(q, k, v, bias_1, bias_2)
        else:
            ctx.save_for_backward(q, k, v, attention_logits)

        return o

    @staticmethod
    def backward(ctx, grad_output):
        if _recompute:
            q, k, v, bias_1, bias_2 = ctx.saved_tensors
            attention_logits = None
        else:
            q, k, v, attention_logits = ctx.saved_tensors
            bias_1 = bias_2 = None
        grad_q = grad_k = grad_v = grad_bias_1 = grad_bias_2 = None

        grad_q = torch.empty_like(q)
        grad_k = torch.empty_like(k)
        grad_v = torch.empty_like(v)
        if(ctx.bias_1_shape is not None):
            grad_bias_1 = torch.zeros(ctx.bias_1_shape)
        if(ctx.bias_2_shape is not None):
            grad_bias_2 = torch.zeros(ctx.bias_2_shape)

        for q_c, k_c, v_c, b1_c, b2_c, al_c, g_out_c, g_q_c, g_k_c, g_v_c, g_b1_c, g_b2_c in \
            zip(torch.split(q, chunk_size, dim=seq_dim),
                torch.split(k, chunk_size, dim=seq_dim),
                torch.split(v, chunk_size, dim=seq_dim),
                split_chunks(bias_1, q.shape[seq_dim], chunk_size, dim=seq_dim),
                split_chunks(bias_2, q.shape[seq_dim], chunk_size, dim=seq_dim),
                split_chunks(attention_logits, q.shape[seq_dim], chunk_size, dim=seq_dim),
                torch.split(grad_output, chunk_size, dim=seq_dim),
                torch.split(grad_q, chunk_size, dim=seq_dim),
                torch.split(grad_k, chunk_size, dim=seq_dim),
                torch.split(grad_v, chunk_size, dim=seq_dim),
                split_chunks(grad_bias_1, q.shape[seq_dim], chunk_size, dim=seq_dim),
                split_chunks(grad_bias_2, q.shape[seq_dim], chunk_size, dim=seq_dim),
            ):
            if _recompute:
                al_c = torch.empty(list(q_c.shape[:-1] + (q_c.shape[-2],)), dtype=q_c.dtype)

                # [*, H, Q, K]
                torch.matmul(
                    q_c, k_c.transpose(-1, -2),
                    out=al_c
                )

                if(b1_c is not None):
                    al_c += b1_c
                if(b2_c is not None):
                    al_c += b2_c

                attn_core_inplace_cuda.forward_(
                    al_c,
                    reduce(mul, al_c.shape[:-1]),
                    al_c.shape[-1],
                )

            tmp_g_v_c = torch.matmul(
                al_c.transpose(-1, -2),
                g_out_c,
            )
            g_v_c.copy_(tmp_g_v_c)

            attn_core_inplace_cuda.backward_(
                al_c,
                g_out_c.contiguous(),
                v_c.contiguous(), # v is implicitly transposed in the kernel
                reduce(mul, al_c.shape[:-1]),
                al_c.shape[-1],
                g_out_c.shape[-1],
            )

            if(ctx.bias_1_shape is not None):
                tmp = torch.sum(
                    al_c,
                    dim=tuple(i for i,d in enumerate(ctx.bias_1_shape) if d == 1),
                    keepdim=True,
                )
                g_b1_c.add_(tmp)

            if(ctx.bias_2_shape is not None):
                tmp = torch.sum(
                    al_c,
                    dim=tuple(i for i,d in enumerate(ctx.bias_2_shape) if d == 1),
                    keepdim=True,
                )
                g_b2_c.add_(tmp)

            torch.matmul(
                al_c, k_c,
                out=g_q_c
            )
            tmp_g_k_c = torch.matmul(
                q_c.transpose(-1, -2), al_c,
            ).transpose(-1, -2)
            g_k_c.copy_(tmp_g_k_c)

        return grad_q, grad_k, grad_v, grad_bias_1, grad_bias_2

attention_core = AttentionCoreFunction.apply
