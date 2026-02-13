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

from random import randint
import torch
import numpy as np
from scipy.spatial.transform import Rotation

from tests.config import consts


def random_asym_ids(n_res, split_chains=True, min_chain_len=4):
    n_chain = randint(1, n_res // min_chain_len) if consts.is_multimer else 1

    if not split_chains:
        return [0] * n_res

    assert n_res >= n_chain

    pieces = []
    asym_ids = []
    final_idx = n_chain - 1
    for idx in range(n_chain - 1):
        n_stop = (n_res - sum(pieces) - n_chain + idx - min_chain_len)
        if n_stop <= min_chain_len:
            final_idx = idx
            break
        piece = randint(min_chain_len, n_stop)
        pieces.append(piece)
        asym_ids.extend(piece * [idx])
    asym_ids.extend((n_res - sum(pieces)) * [final_idx])

    return np.array(asym_ids).astype(np.float32) + 1


def random_template_feats(n_templ, n, batch_size=None):
    b = []
    if batch_size is not None:
        b.append(batch_size)
    batch = {
        "template_mask": np.random.randint(0, 2, (*b, n_templ)),
        "template_pseudo_beta_mask": np.random.randint(0, 2, (*b, n_templ, n)),
        "template_pseudo_beta": np.random.rand(*b, n_templ, n, 3),
        "template_aatype": np.random.randint(0, 22, (*b, n_templ, n)),
        "template_all_atom_mask": np.random.randint(
            0, 2, (*b, n_templ, n, 37)
        ),
        "template_all_atom_positions": 
            np.random.rand(*b, n_templ, n, 37, 3) * 10,
        "template_torsion_angles_sin_cos": 
            np.random.rand(*b, n_templ, n, 7, 2),
        "template_alt_torsion_angles_sin_cos": 
            np.random.rand(*b, n_templ, n, 7, 2),
        "template_torsion_angles_mask": 
            np.random.rand(*b, n_templ, n, 7),
    }
    batch = {k: v.astype(np.float32) for k, v in batch.items()}
    batch["template_aatype"] = batch["template_aatype"].astype(np.int64)

    if consts.is_multimer:
        asym_ids = np.array(random_asym_ids(n))
        batch["asym_id"] = np.tile(asym_ids[np.newaxis, :], (*b, n_templ, 1))

    return batch


def random_extra_msa_feats(n_extra, n, batch_size=None):
    b = []
    if batch_size is not None:
        b.append(batch_size)
    batch = {
        "extra_msa": np.random.randint(0, 22, (*b, n_extra, n)).astype(
            np.int64
        ),
        "extra_has_deletion": np.random.randint(0, 2, (*b, n_extra, n)).astype(
            np.float32
        ),
        "extra_deletion_value": np.random.rand(*b, n_extra, n).astype(
            np.float32
        ),
        "extra_msa_mask": np.random.randint(0, 2, (*b, n_extra, n)).astype(
            np.float32
        ),
    }
    return batch


def random_affines_vector(dim):
    prod_dim = 1
    for d in dim:
        prod_dim *= d

    affines = np.zeros((prod_dim, 7)).astype(np.float32)

    for i in range(prod_dim):
        affines[i, :4] = Rotation.random(random_state=42).as_quat()
        affines[i, 4:] = np.random.rand(
            3,
        ).astype(np.float32)

    return affines.reshape(*dim, 7)


def random_affines_4x4(dim):
    prod_dim = 1
    for d in dim:
        prod_dim *= d

    affines = np.zeros((prod_dim, 4, 4)).astype(np.float32)

    for i in range(prod_dim):
        affines[i, :3, :3] = Rotation.random(random_state=42).as_matrix()
        affines[i, :3, 3] = np.random.rand(
            3,
        ).astype(np.float32)

    affines[:, 3, 3] = 1

    return affines.reshape(*dim, 4, 4)


def random_attention_inputs(batch_size, n_seq, n, no_heads, c_hidden, inf=1e9,
                            dtype=torch.float32, requires_grad=False):
    q = torch.rand(batch_size, n_seq, n, c_hidden, dtype=dtype, requires_grad=requires_grad).cuda()
    kv = torch.rand(batch_size, n_seq, n, c_hidden, dtype=dtype, requires_grad=requires_grad).cuda()

    mask = torch.randint(0, 2, (batch_size, n_seq, 1, 1, n), dtype=dtype, requires_grad=False).cuda()
    z_bias = torch.rand(batch_size, 1, no_heads, n, n, dtype=dtype, requires_grad=requires_grad).cuda()
    mask_bias = inf * (mask - 1)

    biases = [mask_bias, z_bias]

    return q, kv, mask, biases
