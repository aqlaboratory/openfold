# Copyright 2023 RIKEN & Fujitsu Limited
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
from openfold.utils.tensor_utils import permute_final_dims

def flat(t, keep_dims):
    return t.reshape([-1] + list(t.shape[-keep_dims:]))

def forward_(attention_logits, rows, cols):
    chunk_size = 512 * 1024

    for al_c in torch.split(flat(attention_logits, keep_dims=1), chunk_size, dim=0):
        tmp = torch.nn.functional.softmax(al_c, dim=-1)
        al_c.copy_(tmp)

def backward_(attention_logits, #output
              grad_output, #d_ov
              v, #values
              rows,
              cols_output,
              cols_values):
    chunk_size = 2048

    for al_c, go_c, v_c in zip(torch.split(flat(attention_logits, keep_dims=2), chunk_size, dim=0),
                               torch.split(flat(grad_output, keep_dims=2), chunk_size, dim=0),
                               torch.split(flat(v, keep_dims=2), chunk_size, dim=0)):
        v_t = permute_final_dims(v_c, (1, 0))
        dy_buf = torch.matmul(go_c, v_t)
        tmp_sum = torch.sum(al_c * dy_buf, dim=-1)
        al_c.mul_(dy_buf.sub_(tmp_sum.unsqueeze(-1)))
