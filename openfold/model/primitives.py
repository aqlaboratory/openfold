# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
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

from functools import partial
import math
from typing import Optional, Callable, List, Tuple, Sequence
import numpy as np

import torch
import torch.nn as nn
from scipy.stats import truncnorm

from openfold.utils.checkpointing import get_checkpoint_fn
from openfold.utils.tensor_utils import (
    permute_final_dims,
    flatten_final_dims,
    _chunk_slice,
)


def _prod(nums):
    out = 1
    for n in nums:
        out = out * n
    return out


def _calculate_fan(linear_weight_shape, fan="fan_in"):
    fan_out, fan_in = linear_weight_shape

    if fan == "fan_in":
        f = fan_in
    elif fan == "fan_out":
        f = fan_out
    elif fan == "fan_avg":
        f = (fan_in + fan_out) / 2
    else:
        raise ValueError("Invalid fan option")

    return f


def trunc_normal_init_(weights, scale=1.0, fan="fan_in"):
    shape = weights.shape
    f = _calculate_fan(shape, fan)
    scale = scale / max(1, f)
    a = -2
    b = 2
    std = math.sqrt(scale) / truncnorm.std(a=a, b=b, loc=0, scale=1)
    size = _prod(shape)
    samples = truncnorm.rvs(a=a, b=b, loc=0, scale=std, size=size)
    samples = np.reshape(samples, shape)
    with torch.no_grad():
        weights.copy_(torch.tensor(samples, device=weights.device))


def lecun_normal_init_(weights):
    trunc_normal_init_(weights, scale=1.0)


def he_normal_init_(weights):
    trunc_normal_init_(weights, scale=2.0)


def glorot_uniform_init_(weights):
    nn.init.xavier_uniform_(weights, gain=1)


def final_init_(weights):
    with torch.no_grad():
        weights.fill_(0.0)


def gating_init_(weights):
    with torch.no_grad():
        weights.fill_(0.0)


def normal_init_(weights):
    torch.nn.init.kaiming_normal_(weights, nonlinearity="linear")


def ipa_point_weights_init_(weights):
    with torch.no_grad():
        softplus_inverse_1 = 0.541324854612918
        weights.fill_(softplus_inverse_1)


class Linear(nn.Linear):
    """
    A Linear layer with built-in nonstandard initializations. Called just
    like torch.nn.Linear.

    Implements the initializers in 1.11.4, plus some additional ones found
    in the code.
    """

    def __init__(
        self,
        in_dim: int,
        out_dim: int,
        bias: bool = True,
        init: str = "default",
        init_fn: Optional[Callable[[torch.Tensor, torch.Tensor], None]] = None,
    ):
        """
        Args:
            in_dim:
                The final dimension of inputs to the layer
            out_dim:
                The final dimension of layer outputs
            bias:
                Whether to learn an additive bias. True by default
            init:
                The initializer to use. Choose from:

                "default": LeCun fan-in truncated normal initialization
                "relu": He initialization w/ truncated normal distribution
                "glorot": Fan-average Glorot uniform initialization
                "gating": Weights=0, Bias=1
                "normal": Normal initialization with std=1/sqrt(fan_in)
                "final": Weights=0, Bias=0

                Overridden by init_fn if the latter is not None.
            init_fn:
                A custom initializer taking weight and bias as inputs.
                Overrides init if not None.
        """
        super(Linear, self).__init__(in_dim, out_dim, bias=bias)

        if bias:
            with torch.no_grad():
                self.bias.fill_(0)

        if init_fn is not None:
            init_fn(self.weight, self.bias)
        else:
            if init == "default":
                lecun_normal_init_(self.weight)
            elif init == "relu":
                he_normal_init_(self.weight)
            elif init == "glorot":
                glorot_uniform_init_(self.weight)
            elif init == "gating":
                gating_init_(self.weight)
                if bias:
                    with torch.no_grad():
                        self.bias.fill_(1.0)
            elif init == "normal":
                normal_init_(self.weight)
            elif init == "final":
                final_init_(self.weight)
            else:
                raise ValueError("Invalid init string.")


def _attention(query, key, value, biases):
    a = torch.matmul(query, key)

    for b in biases:
        a += b

    a = torch.nn.functional.softmax(a, dim=-1)

    # [*, H, Q, C_hidden]
    o = torch.matmul(a, value)

    # [*, Q, H, C_hidden]
    o = o.transpose(-2, -3)

    return o


@torch.jit.ignore
def _attention_chunk_and_checkpoint(query, key, value, biases, chunk_size):
    if(len(biases) > 2):
        raise ValueError(
            "_chunk_and_checkpoint only permits two bias terms"
        )

    biases = biases + [None, None]
    bias_1, bias_2 = biases[:2]

    def _checkpointable_attention(q, k, v, b1, b2):
        bs = [b1, b2]
        return _attention(q, k, v, bs)

    batch_dims = query.shape[:-3]
    no_batch_dims = len(query.shape[:-3])

    # q, k, and v are assumed to have no singleton dimensions
    flat_q = query.reshape(-1, *query.shape[-3:])
    flat_k = key.reshape(-1, *key.shape[-3:])
    flat_v = value.reshape(-1, *value.shape[-3:])

    o_chunks = []
    checkpoint_fn = get_checkpoint_fn()
    count = flat_q.shape[0]
    for start in range(0, count, chunk_size):
        end = start + chunk_size
        q_chunk = flat_q[start: end, ...]
        k_chunk = flat_k[start: end, ...]
        v_chunk = flat_v[start: end, ...]
        bias_1_chunk = _chunk_slice(bias_1, start, end, no_batch_dims)
        bias_2_chunk = _chunk_slice(bias_2, start, end, no_batch_dims)

        o_chunk = checkpoint_fn(_checkpointable_attention,
            q_chunk, k_chunk, v_chunk, bias_1_chunk, bias_2_chunk
        )

        o_chunks.append(o_chunk)

    o_flat = torch.cat(o_chunks, dim=0)

    return o_flat.reshape(batch_dims + o_flat.shape[1:])


class Attention(nn.Module):
    """
    Standard multi-head attention using AlphaFold's default layer
    initialization. Allows multiple bias vectors.
    """
    def __init__(
        self,
        c_q: int,
        c_k: int,
        c_v: int,
        c_hidden: int,
        no_heads: int,
        gating: bool = True,
    ):
        """
        Args:
            c_q:
                Input dimension of query data
            c_k:
                Input dimension of key data
            c_v:
                Input dimension of value data
            c_hidden:
                Per-head hidden dimension
            no_heads:
                Number of attention heads
            gating:
                Whether the output should be gated using query data
        """
        super(Attention, self).__init__()

        self.c_q = c_q
        self.c_k = c_k
        self.c_v = c_v
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.gating = gating

        # DISCREPANCY: c_hidden is not the per-head channel dimension, as
        # stated in the supplement, but the overall channel dimension.

        self.linear_q = Linear(
            self.c_q, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )
        self.linear_k = Linear(
            self.c_k, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )
        self.linear_v = Linear(
            self.c_v, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )
        self.linear_o = Linear(
            self.c_hidden * self.no_heads, self.c_q, init="final"
        )

        self.linear_g = None
        if self.gating:
            self.linear_g = Linear(
                self.c_q, self.c_hidden * self.no_heads, init="gating"
            )

        self.sigmoid = nn.Sigmoid()

    def forward(
        self,
        q_x: torch.Tensor,
        k_x: torch.Tensor,
        v_x: torch.Tensor,
        biases: Optional[List[torch.Tensor]] = None,
        use_lma: bool = False,
        q_chunk_size: Optional[int] = None,
        kv_chunk_size: Optional[int] = None,
        _chunk_and_checkpoint: Optional[int] = None
    ) -> torch.Tensor:
        """
        Args:
            q_x:
                [*, Q, C_q] query data
            k_x:
                [*, K, C_k] key data
            v_x:
                [*, V, C_v] value data
        Returns
            [*, Q, C_q] attention update
        """
        if(biases is None):
            biases = []
        if(use_lma and (q_chunk_size is None or kv_chunk_size is None)):
            raise ValueError(
                "If use_lma is specified, q_chunk_size and kv_chunk_size must "
                "be provided"
            )
        if(use_lma and _chunk_and_checkpoint is not None):
            raise ValueError(
                "use_lma and _chunk_and_checkpoint are mutually exclusive"
            )

        # [*, Q/K/V, H * C_hidden]
        q = self.linear_q(q_x)
        k = self.linear_k(k_x)
        v = self.linear_v(v_x)

        # [*, Q/K, H, C_hidden]
        q = q.view(q.shape[:-1] + (self.no_heads, -1))
        k = k.view(k.shape[:-1] + (self.no_heads, -1))
        v = v.view(v.shape[:-1] + (self.no_heads, -1))

        q = q / math.sqrt(self.c_hidden) 

        if(use_lma):
            biases = [
                b.expand(b.shape[:-2] + (q_x.shape[-2],) + (k_x.shape[-2],)) 
                for b in biases
            ]
            o = _lma(q, k, v, biases, q_chunk_size, kv_chunk_size)
        else:
            # [*, H, Q, C_hidden]
            q = permute_final_dims(q, (1, 0, 2))
            
            # [*, H, C_hidden, K]
            k = permute_final_dims(k, (1, 2, 0))
    
            # [*, H, V, C_hidden]
            v = permute_final_dims(v, (1, 0, 2))

            if(_chunk_and_checkpoint):
                # REMEMBER THAT THE K, Q, V COMPUTATION AND GATING ARE *NOT* 
                # CHECKPOINTED HERE
                o = _attention_chunk_and_checkpoint(
                    q, k, v, biases, _chunk_and_checkpoint
                )
            else:
                o = _attention(q, k, v, biases)

        if(self.linear_g is not None):
            g = self.sigmoid(self.linear_g(q_x))
            # [*, Q, H, C_hidden]
            g = g.view(g.shape[:-1] + (self.no_heads, -1))
            o = o * g

        # [*, Q, H * C_hidden]
        o = flatten_final_dims(o, 2)

        # [*, Q, C_q]
        o = self.linear_o(o)

        return o


class GlobalAttention(nn.Module):
    def __init__(self, c_in, c_hidden, no_heads, inf, eps):
        super(GlobalAttention, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.inf = inf
        self.eps = eps

        self.linear_q = Linear(
            c_in, c_hidden * no_heads, bias=False, init="glorot"
        )

        self.linear_k = Linear(
            c_in, c_hidden, bias=False, init="glorot",
        )
        self.linear_v = Linear(
            c_in, c_hidden, bias=False, init="glorot",
        )
        self.linear_g = Linear(c_in, c_hidden * no_heads, init="gating")
        self.linear_o = Linear(c_hidden * no_heads, c_in, init="final")

        self.sigmoid = nn.Sigmoid()
        self.softmax = nn.Softmax(dim=-1)

    def forward(self, m: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        # [*, N_res, C_in]
        q = torch.sum(m * mask.unsqueeze(-1), dim=-2) / (
            torch.sum(mask, dim=-1)[..., None] + self.eps
        )

        # [*, N_res, H * C_hidden]
        q = self.linear_q(q)
        q *= (self.c_hidden ** (-0.5))

        # [*, N_res, H, C_hidden]
        q = q.view(q.shape[:-1] + (self.no_heads, -1))

        # [*, N_res, N_seq, C_hidden]
        k = self.linear_k(m)
        v = self.linear_v(m)

        # [*, N_res, H, N_seq]
        a = torch.matmul(
            q,
            k.transpose(-1, -2),  # [*, N_res, C_hidden, N_seq]
        )
        bias = (self.inf * (mask - 1))[..., :, None, :]
        a += bias
        a = self.softmax(a)

        # [*, N_res, H, C_hidden]
        o = torch.matmul(
            a,
            v,
        )

        # [*, N_res, N_seq, C_hidden]
        g = self.sigmoid(self.linear_g(m))

        # [*, N_res, N_seq, H, C_hidden]
        g = g.view(g.shape[:-1] + (self.no_heads, -1))

        # [*, N_res, N_seq, H, C_hidden]
        o = o.unsqueeze(-3) * g

        # [*, N_res, N_seq, H * C_hidden]
        o = o.reshape(o.shape[:-2] + (-1,))

        # [*, N_res, N_seq, C_in]
        m = self.linear_o(o)

        return m


def _lma(
    q: torch.Tensor, 
    k: torch.Tensor, 
    v: torch.Tensor, 
    biases: List[torch.Tensor], 
    q_chunk_size: int, 
    kv_chunk_size: int,
):
    no_q, no_kv = q.shape[-3], k.shape[-3]

    # [*, Q, H, C_hidden]
    o = q.new_zeros(q.shape)
    for q_s in range(0, no_q, q_chunk_size):
        q_chunk = q[..., q_s: q_s + q_chunk_size, :, :]
        large_bias_chunks = [
            b[..., q_s: q_s + q_chunk_size, :] for b in biases
        ]

        maxes = []
        weights = []
        values = []
        for kv_s in range(0, no_kv, kv_chunk_size):
            k_chunk = k[..., kv_s: kv_s + kv_chunk_size, :, :]
            v_chunk = v[..., kv_s: kv_s + kv_chunk_size, :, :]
            small_bias_chunks = [
                b[..., kv_s: kv_s + kv_chunk_size] for b in large_bias_chunks
            ]

            a = torch.einsum(
                "...qhd,...khd->...hqk", query, key
            )
        
            for b in small_bias_chunks:
                a += b
        
            a = a.transpose(-2, -3)
        
            max_a = torch.max(a, dim=-1, keepdim=True)[0]
            exp_a = torch.exp(a - max_a)
            exp_v = torch.einsum("...vhf,...qhv->...qhf", value, exp_a)
 
            maxes.append(max_a.detach().squeeze(-1))
            weights.append(torch.sum(exp_a, dim=-1))
            values.append(exp_v)

        chunk_max = torch.stack(maxes, dim=-3)
        chunk_weights = torch.stack(weights, dim=-3)
        chunk_values = torch.stack(values, dim=-4)

        global_max = torch.max(chunk_max, dim=-3, keepdim=True)[0]
        max_diffs = torch.exp(chunk_max - global_max)
        chunk_values *= max_diffs.unsqueeze(-1)
        chunk_weights *= max_diffs

        all_values = torch.sum(chunk_values, dim=-4)
        all_weights = torch.sum(chunk_weights.unsqueeze(-1), dim=-4)

        q_chunk_out = all_values / all_weights

        o[..., q_s: q_s + q_chunk_size, :, :] = q_chunk_out

    return o
