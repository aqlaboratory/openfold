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
import numpy as np
from scipy.spatial.transform import Rotation
import pickle,os
from tests.config import consts
import gzip

def process_label(all_atom_positions: np.ndarray, operation) -> np.ndarray:
    """
    Adapted from unifold
    https://github.com/dptech-corp/Uni-Fold/blob/b1c89a2cebd4e4ee4c47b4e443f92beeb9138fbb/unifold/dataset.py#L55-L61
    """
    if operation == "I":
        return all_atom_positions
    rot, trans = operation
    rot = np.array(rot).reshape(3, 3)
    trans = np.array(trans).reshape(3)
    return all_atom_positions @ rot.T + trans

def load_single_label(
    label_id,label_dir,
    symmetry_operation=None,
):
    """
    Adapted from unifold
    https://github.com/dptech-corp/Uni-Fold/blob/b1c89a2cebd4e4ee4c47b4e443f92beeb9138fbb/unifold/dataset.py#L101-L116
    
    args:
    label: is the dictionary of numpy arrays. created by loading the label pickle file
    """
    label_path = os.path.join(label_dir,f"{label_id}.label.pkl.gz")
    label = pickle.load(gzip.open(label_path,"rb"))
    if symmetry_operation is not None:
        label["all_atom_positions"] = process_label(
            label["all_atom_positions"], symmetry_operation
        )
    label = {
        k: v
        for k, v in label.items()
        if k in ["aatype", "all_atom_positions", "all_atom_mask", "resolution"]
    }
    return label

def load_labels(label_dir,label_ids:list):
    symmetry_operations = ["I" for _ in label_ids] # for now suppose there are NO symmetry operations
    all_chain_labels = [
        load_single_label(l, label_dir, o)
        for l, o in zip(label_ids, symmetry_operations)
    ]



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

    return np.array(asym_ids).astype(np.int64)


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
