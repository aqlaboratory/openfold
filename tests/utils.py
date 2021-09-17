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

import numpy as np


def random_template_feats(n_templ, n, batch_size=None):
    b = []
    if(batch_size is not None):
        b.append(batch_size)
    batch = {
        "template_mask": np.random.randint(0, 2, (*b, n_templ)),
        "template_pseudo_beta_mask": np.random.randint(0, 2, (*b, n_templ, n)),
        "template_pseudo_beta": np.random.rand(*b, n_templ, n, 3),
        "template_aatype": np.random.randint(0, 22, (*b, n_templ, n)),
        "template_all_atom_masks": np.random.randint(
            0, 2, (*b, n_templ, n, 37)
        ),
        "template_all_atom_positions": np.random.rand(
            *b, n_templ, n, 37, 3
        ) * 10,
    }
    batch = {k:v.astype(np.float32) for k,v in batch.items()}
    batch["template_aatype"] = batch["template_aatype"].astype(np.int64)
    return batch

def random_extra_msa_feats(n_extra, n, batch_size=None):
    b = []
    if(batch_size is not None):
        b.append(batch_size)
    batch = {
        "extra_msa": 
            np.random.randint(0, 22, (*b, n_extra, n)).astype(np.int64),
        "extra_has_deletion": 
            np.random.randint(0, 2, (*b, n_extra, n)).astype(np.float32),
        "extra_deletion_value": 
            np.random.rand(*b, n_extra, n).astype(np.float32),
        "extra_msa_mask":
            np.random.randint(0, 2, (*b, n_extra, n)).astype(np.float32),
    }
    return batch
