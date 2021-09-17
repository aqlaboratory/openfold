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
import torch
import torch.nn as nn

from alphafold.model.primitives import Linear
from alphafold.utils.tensor_utils import chunk_layer


class OuterProductMean(nn.Module):
    """
        Implements Algorithm 10.
    """
    def __init__(self, c_m, c_z, c_hidden, chunk_size=4, eps=1e-3):
        """
            Args:
                c_m:
                    MSA embedding channel dimension
                c_z:
                    Pair embedding channel dimension
                c_hidden:
                    Hidden channel dimension
        """
        super(OuterProductMean, self).__init__()

        self.c_z = c_z
        self.c_hidden = c_hidden
        self.chunk_size = chunk_size
        self.eps = eps

        self.layer_norm = nn.LayerNorm(c_m)
        self.linear_1 = Linear(c_m, c_hidden)
        self.linear_2 = Linear(c_m, c_hidden)
        self.linear_out = Linear(c_hidden**2, c_z, init="final")

    def _opm(self, a, b):
        # [*, N_res, N_res, C, C]
        outer = torch.einsum("...bac,...dae->...bdce", a, b)
        
        # [*, N_res, N_res, C * C]
        outer = outer.reshape(*outer.shape[:-2], -1)

        # [*, N_res, N_res, C_z]
        outer = self.linear_out(outer)

        return outer

    def forward(self, m, mask=None):
        """
            Args:
                m:
                    [*, N_seq, N_res, C_m] MSA embedding
                mask:
                    [*, N_seq, N_res] MSA mask
            Returns:
                [*, N_res, N_res, C_z] pair embedding update
        """
        if(mask is None):
            mask = m.new_ones(m.shape[:-1], requires_grad=False)

        # [*, N_seq, N_res, C_m]
        m = self.layer_norm(m)
        
        # [*, N_seq, N_res, C]
        mask = mask.unsqueeze(-1)
        a = self.linear_1(m) * mask
        b = self.linear_2(m) * mask

        a = a.transpose(-2, -3)
        b = b.transpose(-2, -3)

        if(not self.training and self.chunk_size is not None):
            # Since the "batch dim" in this case is not a true batch dimension
            # (in that the shape of the output depends on it), we need to
            # iterate over it ourselves
            a_reshape = a.reshape(-1, *a.shape[-3:])
            b_reshape = b.reshape(-1, *b.shape[-3:])
            out = []
            for a_prime, b_prime in zip(a_reshape, b_reshape):
                outer = chunk_layer(
                    partial(self._opm, b=b_prime),
                    {"a": a_prime},
                    chunk_size=self.chunk_size,
                    no_batch_dims=1,
                )
                out.append(outer)
            outer = torch.stack(out, dim=0)
            outer = outer.reshape(*a.shape[:-3], *outer.shape[1:])
        else:
            outer = self._opm(a, b)

        # [*, N_res, N_res, 1]
        norm = torch.einsum("...abc,...adc->...bdc", mask, mask)

        # [*, N_res, N_res, C_z]
        outer /= self.eps + norm

        return outer


if __name__ == "__main__":
    batch_size = 2
    s = 5
    n_res = 100
    c_m = 256
    c = 32
    c_z = 128

    opm = OuterProductMean(c_m, c_z, c)

    m = torch.rand((batch_size, s, n_res, c_m))
    m = opm(m)

    assert(m.shape == (batch_size, n_res, n_res, c_z))
