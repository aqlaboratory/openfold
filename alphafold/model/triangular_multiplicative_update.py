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
import torch
import torch.nn as nn

from alphafold.model.primitives import Linear
from alphafold.utils.tensor_utils import permute_final_dims


class TriangleMultiplicativeUpdate(nn.Module):
    """
        Implements Algorithms 11 and 12.
    """
    def __init__(self, c_z, c_hidden, _outgoing=True):
        """
            Args:
                c_z:
                    Input channel dimension
                c:
                    Hidden channel dimension
        """ 
        super(TriangleMultiplicativeUpdate, self).__init__()
        self.c_z = c_z
        self.c_hidden = c_hidden
        self._outgoing = _outgoing

        self.linear_a_p = Linear(self.c_z, self.c_hidden)
        self.linear_a_g = Linear(self.c_z, self.c_hidden, init="gating")
        self.linear_b_p = Linear(self.c_z, self.c_hidden)
        self.linear_b_g = Linear(self.c_z, self.c_hidden, init="gating")
        self.linear_g = Linear(self.c_z, self.c_z, init="gating")
        self.linear_z = Linear(self.c_hidden, self.c_z, init="final")

        self.layer_norm_in = nn.LayerNorm(self.c_z)
        self.layer_norm_out = nn.LayerNorm(self.c_hidden)

        self.sigmoid = nn.Sigmoid()

        cp = self._outgoing_matmul if self._outgoing else self._incoming_matmul
        self.combine_projections = cp

    def _outgoing_matmul(self, 
        a: torch.Tensor,    # [*, N_i, N_k, C] 
        b: torch.Tensor,    # [*, N_j, N_k, C]
    ):
        # [*, C, N_i, N_j]
        p = torch.matmul(
            permute_final_dims(a, 2, 0, 1),
            permute_final_dims(b, 2, 1, 0),
        )
        
        # [*, N_i, N_j, C]
        return permute_final_dims(p, 1, 2, 0)

    def _incoming_matmul(self, 
        a: torch.Tensor,    # [*, N_k, N_i, C] 
        b: torch.Tensor,    # [*, N_k, N_j, C]
    ):

        # [*, C, N_i, N_j]
        p = torch.matmul(
            permute_final_dims(a, 2, 1, 0),
            permute_final_dims(b, 2, 0, 1),
        )
       
        # [*, N_i, N_j, C]
        return permute_final_dims(p, 1, 2, 0)
    
    def forward(self, z, mask=None):
        """
            Args:
                x:
                    [*, N_res, N_res, C_z] input tensor
                mask:
                    [*, N_res, N_res] input mask
            Returns:
                [*, N_res, N_res, C_z] output tensor
        """
        if(mask is None):
            mask = z.new_ones(z.shape[:-1], requires_grad=False)

        mask = mask.unsqueeze(-1)

        z = self.layer_norm_in(z)
        a = self.linear_a_p(z) * self.sigmoid(self.linear_a_g(z))
        a = a * mask
        b = self.linear_b_p(z) * self.sigmoid(self.linear_b_g(z))
        b = b * mask
        x = self.combine_projections(a, b)
        x = self.layer_norm_out(x)
        x = self.linear_z(x)
        g = self.sigmoid(self.linear_g(z))
        z = x * g

        return z


class TriangleMultiplicationOutgoing(TriangleMultiplicativeUpdate):
    """
        Implements Algorithm 11.
    """
    __init__ = partialmethod(
        TriangleMultiplicativeUpdate.__init__, _outgoing=True,
    )
   

class TriangleMultiplicationIncoming(TriangleMultiplicativeUpdate):
    """
        Implements Algorithm 12.
    """
    __init__ = partialmethod(
        TriangleMultiplicativeUpdate.__init__, _outgoing=False,
    )


if __name__ == "__main__":
    c_in = 256 # doubled to make shape changes more apparent
    c = 128
    outgoing = True

    tm = TriangleMultiplication(
        c_in,
        c,
        outgoing,
    )

    n_res = 300
    batch_size = 16
    x = torch.rand((batch_size, n_res, n_res, c_in))
    shape_before = x.shape
    x = tm(x)
    shape_after = x.shape

    assert(shape_before == shape_after)
