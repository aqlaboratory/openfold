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
from typing import Optional

from openfold.model.primitives import Linear, Attention, GlobalAttention 
from openfold.utils.tensor_utils import (
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
            mask = m.new_ones(
                m.shape[:-3] + (n_seq, n_res), 
            )

        # [*, N_seq, 1, 1, N_res]
        bias = (self.inf * (mask - 1))[..., :, None, None, :]
        
        # [*, N_seq, no_heads, N_res, N_res]
        bias = bias.expand(
            ((-1,) * len(bias.shape[:-4])) + (-1, self.no_heads, n_res, -1)
        )
        biases = [bias]
        if(self.pair_bias):
            # [*, N_res, N_res, C_z]
            z = self.layer_norm_z(z)

            # [*, N_res, N_res, no_heads]
            z = self.linear_z(z)

            # [*, 1, no_heads, N_res, N_res]
            z = permute_final_dims(z, (2, 0, 1)).unsqueeze(-4)

            biases.append(z)

        mha_inputs = {
            "q_x": m, 
            "k_x": m, 
            "v_x": m, 
            "biases": biases
        }
        if(self.chunk_size is not None):
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
    def __init__(self, c_m, c_z, c_hidden, no_heads, chunk_size, inf=1e9):
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
            chunk_size=chunk_size,
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

        self.layer_norm_m = nn.LayerNorm(c_in)

        self.global_attention = GlobalAttention(
            c_in=c_in,
            c_hidden=c_hidden,
            no_heads=no_heads,
            inf=inf,
            eps=eps,
        )

    def forward(self, 
        m: torch.Tensor, 
        mask: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        n_seq, n_res, c_in = m.shape[-3:]

        if(mask is None):
            # [*, N_seq, N_res]
            mask = torch.ones(
                m.shape[:-1],
                dtype=m.dtype,
                device=m.device,
            ).detach()

        # [*, N_res, N_seq, C_in]
        m = m.transpose(-2, -3)
        mask = mask.transpose(-1, -2)

        # [*, N_res, N_seq, C_in]
        m = self.layer_norm_m(m)

        mha_input = {
            "m": m,
            "mask": mask,
        }
        if(self.chunk_size is not None):
            m = chunk_layer(
                self.global_attention,
                mha_input,
                chunk_size=self.chunk_size,
                no_batch_dims=len(m.shape[:-2])
            )
        else:
            m = self.global_attention(m=mha_input["m"], mask=mha_input["mask"])
 
        # [*, N_seq, N_res, C_in]
        m = m.transpose(-2, -3)

        return m
