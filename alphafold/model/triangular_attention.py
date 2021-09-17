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

from functools import partialmethod
import math
import torch
import torch.nn as nn

from alphafold.model.primitives import Linear, Attention
from alphafold.utils.tensor_utils import (
    chunk_layer, 
    permute_final_dims, 
    flatten_final_dims,
)


class TriangleAttention(nn.Module):
    def __init__(self, 
        c_in, 
        c_hidden, 
        no_heads, 
        starting, 
        chunk_size=4, 
        inf=1e9
    ):
        """
            Args:
                c_in:
                    Input channel dimension
                c_hidden:
                    Overall hidden channel dimension (not per-head)
                no_heads:
                    Number of attention heads
        """
        super(TriangleAttention, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.starting = starting
        self.chunk_size = chunk_size
        self.inf = inf

        self.layer_norm = nn.LayerNorm(self.c_in)
        
        self.linear = Linear(c_in, self.no_heads, bias=False, init="normal")

        self.mha = Attention(
            self.c_in, self.c_in, self.c_in,
            self.c_hidden, 
            self.no_heads
        )

    def forward(self, x, mask=None):
        """
            Args:
                x:
                    [*, I, J, C_in] input tensor (e.g. the pair representation)
            Returns:
                [*, I, J, C_in] output tensor
        """
        if(mask is None):
            # [*, I, J]
            mask = torch.ones(
                x.shape[:-1], 
                device=x.device,
                requires_grad=False,
            )

        # Shape annotations assume self.starting. Else, I and J are flipped
        if(not self.starting):
            x = x.transpose(-2, -3)
            mask = mask.transpose(-1, -2)

        # [*, I, J, C_in]
        x = self.layer_norm(x)
    
        # [*, I, 1, 1, J]
        mask_bias = (self.inf * (mask - 1))[..., :, None, None, :]
        
        # [*, H, I, J]
        triangle_bias = permute_final_dims(self.linear(x), 2, 0, 1)

        # [*, 1, H, I, J]
        triangle_bias = triangle_bias.unsqueeze(-4)

        # Broadcasting and chunking doesn't really work yet (TODO)
        # [*, I, H, I, J]
        i = x.shape[-3]
        triangle_bias = triangle_bias.expand(
            (*((-1,) * len(triangle_bias.shape[:-4])), i, -1, -1, -1)
        )

        #print(x.shape)
        #print(mask_bias.shape)
        #print(triangle_bias.shape)

        mha_inputs = {
            "q_x": x,
            "k_x": x,
            "v_x": x,
            "biases": [mask_bias, triangle_bias],
        }
        if(not self.training and self.chunk_size is not None):
            x = chunk_layer(
                self.mha,
                mha_inputs, 
                chunk_size=self.chunk_size,
                no_batch_dims=len(x.shape[:-2])
            )
        else:
            x = self.mha(**mha_inputs)

        if(not self.starting):
            x = x.transpose(-2, -3)

        return x


class TriangleAttentionStartingNode(TriangleAttention):
    """
        Implements Algorithm 13.
    """
    __init__ = partialmethod(TriangleAttention.__init__, starting=True)


class TriangleAttentionEndingNode(TriangleAttention):
    """
        Implements Algorithm 14.
    """
    __init__ = partialmethod(TriangleAttention.__init__, starting=False)


if __name__ == "__main__":
    c_in = 256
    c = 32
    no_heads = 4
    starting = True

    tan = TriangleAttention(
        c_in,
        c,
        no_heads,
        starting
    )

    batch_size = 16
    n_res = 256

    x = torch.rand((batch_size, n_res, n_res, c_in))
    shape_before = x.shape
    x = tan(x)
    shape_after = x.shape

    assert(shape_before == shape_after)
