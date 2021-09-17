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

import torch
import torch.nn as nn

from alphafold.model.primitives import Linear
from alphafold.utils.tensor_utils import chunk_layer


class PairTransition(nn.Module):
    """
        Implements Algorithm 15.
    """
    def __init__(self, c_z, n, chunk_size=4):
        """
            Args:
                c_z:
                    Pair transition channel dimension
                n:
                    Factor by which c_z is multiplied to obtain hidden channel 
                    dimension
        """
        super(PairTransition, self).__init__()

        self.c_z = c_z
        self.n = n
        self.chunk_size = chunk_size

        self.layer_norm = nn.LayerNorm(self.c_z)
        self.linear_1 = Linear(self.c_z, self.n * self.c_z, init="relu")
        self.relu = nn.ReLU()
        self.linear_2 = Linear(self.n * self.c_z, c_z, init="final")

    def _transition(self, z, mask):
        # [*, N_res, N_res, C_hidden]
        z = self.linear_1(z)
        z = self.relu(z)

        # [*, N_res, N_res, C_z]
        z = self.linear_2(z) * mask

        return z

    def forward(self, z, mask=None):
        """
            Args:
                z:
                    [*, N_res, N_res, C_z] pair embedding
            Returns:
                [*, N_res, N_res, C_z] pair embedding update
        """
        # DISCREPANCY: DeepMind forgets to apply the mask in this module.
        if(mask is None):
            mask = z.new_ones(z.shape[:-1], requires_grad=False)

        # [*, N_res, N_res, 1]
        mask = mask.unsqueeze(-1)

        # [*, N_res, N_res, C_z]
        z = self.layer_norm(z)

        inp = {"z": z, "mask": mask}
        if(not self.training and self.chunk_size is not None):
            z = chunk_layer(
                self._transition,
                inp,
                chunk_size=self.chunk_size,
                no_batch_dims=len(z.shape[:-2]), 
            )
        else:
            z = self._transition(**inp)

        return z


if __name__ == "__main__":
    n = 4
    c_in = 128

    pt = PairTransition(n, c_in)

    batch_size = 4
    n_res = 256

    z = torch.rand((batch_size, n_res, n_res, c_in))
    shape_before = z.shape
    z = pt(z)
    shape_after = z.shape

    assert(shape_before == shape_after)
