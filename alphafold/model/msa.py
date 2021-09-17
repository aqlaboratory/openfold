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

import math
import torch
import torch.nn as nn

from alphafold.model.primitives import Linear, Attention
from alphafold.utils.tensor_utils import (
    chunk_layer,
    permute_final_dims, 
    flatten_final_dims,
)


class MSAAttention(nn.Module):
    def __init__(self, 
        c_in, 
        c_hidden, 
        no_heads, 
        pair_bias=False, 
        c_z=None, 
        chunk_size=4,
        inf=1e9,
    ):
        """
            Args:
                c_in:
                    Input channel dimension
                c_hidden:
                    Per-head hidden channel dimension
                no_heads:
                    Number of attention heads
                pair_bias:
                    Whether to use pair embedding bias
                c_z:
                    Pair embedding channel dimension. Ignored unless pair_bias
                    is true
                inf:
                    A large number to be used in computing the attention mask
        """
        super(MSAAttention, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.pair_bias = pair_bias
        self.c_z = c_z
        self.chunk_size = chunk_size
        self.inf = inf

        self.layer_norm_m = nn.LayerNorm(self.c_in)

        if(self.pair_bias):
            self.layer_norm_z = nn.LayerNorm(self.c_z)
            self.linear_z = Linear(
                self.c_z, self.no_heads, bias=False, init="normal"
            )

        self.mha = Attention(
            self.c_in, self.c_in, self.c_in, 
            self.c_hidden, 
            self.no_heads
        )

    def forward(self, m, z=None, mask=None):
        """
            Args:
                m:
                    [*, N_seq, N_res, C_m] MSA embedding
                z:
                    [*, N_res, N_res, C_z] pair embedding. Required only if 
                    pair_bias is True
                mask:
                    [*, N_seq, N_res] MSA mask
        """
        # [*, N_seq, N_res, C_m]
        m = self.layer_norm_m(m)

        n_seq, n_res = m.shape[-3:-1]
        if(mask is None):
            # [*, N_seq, N_res]
            mask = torch.ones(
                (*m.shape[:-3], n_seq, n_res), 
                device=m.device, 
                requires_grad=False
            )

        # [*, N_seq, 1, 1, N_res]
        bias = (self.inf * (mask - 1))[..., :, None, None, :]
        
        # [*, N_seq, no_heads, N_res, N_res]
        bias = bias.expand(
            (*((-1,) * len(bias.shape[:-4])), -1, self.no_heads, n_res, -1)
        )

        if(self.pair_bias):
            # [*, N_res, N_res, C_z]
            z = self.layer_norm_z(z)

            # [*, N_res, N_res, no_heads]
            z = self.linear_z(z)

            # [*, 1, no_heads, N_res, N_res]
            z = permute_final_dims(z, 2, 0, 1).unsqueeze(-4)

            # [*, N_seq, no_heads, N_res, N_res]
            bias = bias + z

        mha_inputs = {
            "q_x": m, 
            "k_x": m, 
            "v_x": m, 
            "biases": [bias]
        }
        if(not self.training and self.chunk_size is not None):
            m = chunk_layer(
                self.mha,
                mha_inputs,
                chunk_size=self.chunk_size,
                no_batch_dims=len(m.shape[:-2])
            )
        else:
            m = self.mha(**mha_inputs)

        return m


class MSARowAttentionWithPairBias(MSAAttention):
    """
        Implements Algorithm 7.
    """
    def __init__(self, c_m, c_z, c_hidden, no_heads, inf=1e9):
        """
            Args:
                c_m:
                    Input channel dimension
                c_z:
                    Pair embedding channel dimension
                c_hidden:
                    Per-head hidden channel dimension
                no_heads:
                    Number of attention heads
                inf:
                    Large number used to construct attention masks
        """
        super(MSARowAttentionWithPairBias, self).__init__(
            c_m, 
            c_hidden, 
            no_heads, 
            pair_bias=True, 
            c_z=c_z, 
            inf=inf,
        )


class MSAColumnAttention(MSAAttention):
    """
        Implements Algorithm 8.
    """
    def __init__(self, c_m, c_hidden, no_heads, chunk_size=4, inf=1e9):
        """
            Args:
                c_m:
                    MSA channel dimension
                c_hidden:
                    Per-head hidden channel dimension
                no_heads:
                    Number of attention heads
                inf:
                    Large number used to construct attention masks
        """
        super(MSAColumnAttention, self).__init__(
            c_in=c_m,
            c_hidden=c_hidden,
            no_heads=no_heads,
            pair_bias=False,
            c_z=None,
            chunk_size=chunk_size,
            inf=inf,
        )

        
    def forward(self, m, mask=None):
        """
            Args:
                m:
                    [*, N_seq, N_res, C_m] MSA embedding
                mask:
                    [*, N_seq, N_res] MSA mask
        """
        # [*, N_res, N_seq, C_in]
        m = m.transpose(-2, -3)
        if(mask is not None):
            mask = mask.transpose(-1, -2)

        m = super().forward(m, mask=mask)

        # [*, N_seq, N_res, C_in]
        m = m.transpose(-2, -3)
        if(mask is not None):
            mask = mask.transpose(-1, -2)
        return m


class MSAColumnGlobalAttention(nn.Module):
    def __init__(self, 
        c_in, 
        c_hidden, 
        no_heads, 
        chunk_size=4, 
        inf=1e9, 
        eps=1e-10
    ):
        super(MSAColumnGlobalAttention, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.chunk_size = chunk_size
        self.inf = inf
        self.eps = eps

        self.layer_norm_m = nn.LayerNorm(self.c_in)

        self.linear_q = Linear(
            self.c_in, self.c_hidden * self.no_heads, bias=False, init="glorot"
        )

        C_hidden = self.c_hidden
        self.linear_k = Linear(
            self.c_in, C_hidden, bias=False, init="glorot",
        )
        self.linear_v = Linear(
            self.c_in, C_hidden, bias=False, init="glorot",
        )
        self.linear_g = Linear(self.c_in, self.c_hidden * self.no_heads, init="gating")
        self.linear_o = Linear(self.c_hidden * self.no_heads, self.c_in, init="final")

        self.sigmoid = nn.Sigmoid()
        self.softmax = nn.Softmax(dim=-1)

    def global_attention(self, m, mask):
        # [*, N_res, C_in]
        q = (torch.sum(m * mask.unsqueeze(-1), dim=-2) / 
               (torch.sum(mask, dim=-1)[..., None] + self.eps))

        # [*, N_res, H * C_hidden]
        q = self.linear_q(q)
        q *= self.c_hidden ** (-0.5)

        # [*, N_res, H, C_hidden]
        q = q.view(*q.shape[:-1], self.no_heads, -1)

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
        g = g.view(*g.shape[:-1], self.no_heads, -1)

        # [*, N_res, N_seq, H, C_hidden]
        o = o.unsqueeze(-3) * g
        
        # [*, N_res, N_seq, H * C_hidden]
        o = o.reshape(*o.shape[:-2], -1)

        # [*, N_res, N_seq, C_in]
        m = self.linear_o(o)

        return m

    def forward(self, m, mask=None):
        n_seq, n_res, c_in = m.shape[-3:]

        if(mask is None):
            # [*, N_seq, N_res]
            mask = m.new_ones(m.shape[:-1], requires_grad=False)

        # [*, N_res, N_seq, C_in]
        m = m.transpose(-2, -3)
        mask = mask.transpose(-1, -2)

        # [*, N_res, N_seq, C_in]
        m = self.layer_norm_m(m)

        mha_input = {
            "m": m,
            "mask": mask,
        }
        if(not self.training and self.chunk_size is not None):
            m = chunk_layer(
                self.global_attention,
                mha_input,
                chunk_size=self.chunk_size,
                no_batch_dims=len(m.shape[:-2])
            )
        else:
            m = self.global_attention(**mha_input)
 
        # [*, N_seq, N_res, C_in]
        m = m.transpose(-2, -3)

        return m


if __name__ == "__main__":
    batch_size = 2
    s_t = 3
    n = 100
    c_in = 128
    c = 32
    no_heads = 4

    msaca = MSAColumnAttention(c_in, c, no_heads)

    x = torch.rand((batch_size, s_t, n, c_in))

    shape_before = x.shape
    x = msaca(x)
    shape_after = x.shape

    assert(shape_before == shape_after)
