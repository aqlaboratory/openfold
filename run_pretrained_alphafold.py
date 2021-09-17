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

import os
import sys
sys.path.append("lib/conda/lib/python3.9/site-packages")

import math
import pickle
import time
import torch
import torch.nn as nn
import numpy as np

from config import model_config
from alphafold.model.model import AlphaFold
import alphafold.np.protein as protein
import alphafold.np.relax.relax as relax
from alphafold.np import residue_constants
from alphafold.utils.import_weights import (
    import_jax_weights_,
)
from alphafold.utils.tensor_utils import (
    tree_map,
    tensor_tree_map,
)


MODEL_NAME = "model_1"
MODEL_DEVICE = "cuda:1"
PARAM_PATH = "alphafold/resources/params/params_model_1.npz"
FEAT_PATH = "tests/test_data/sample_feats.pickle"


config = model_config(MODEL_NAME)
model = AlphaFold(config.model)
model = model.eval()
import_jax_weights_(model, PARAM_PATH)
model_device = 'cuda:1'
model = model.to(model_device)

with open(FEAT_PATH, "rb") as f:
    batch = pickle.load(f)

batch = {k:torch.as_tensor(v, device=model_device) for k,v in batch.items()}

longs = [
    "aatype", 
    "template_aatype", 
    "extra_msa", 
    "residx_atom37_to_atom14",
    "residx_atom14_to_atom37",
]
for l in longs:
    batch[l] = batch[l].long()


# Move the recycling dimension to the end
move_dim = lambda t: t.permute(*range(len(t.shape))[1:], 0).contiguous()
batch = tensor_tree_map(move_dim, batch)

with torch.no_grad():
    t = time.time()
    out = model(batch)
    print(f"Inference time: {time.time() - t}")

# Toss out the recycling dimensions --- we don't need them anymore
batch = tensor_tree_map(lambda x: np.array(x[..., -1].cpu()), batch)
out = tensor_tree_map(lambda x: np.array(x.cpu()), out)

plddt = out["plddt"]
mean_plddt = np.mean(plddt)

plddt_b_factors = np.repeat(
    plddt[..., None], residue_constants.atom_type_num, axis=-1
)

unrelaxed_protein = protein.from_prediction(
    features=batch,
    result=out,
    b_factors=plddt_b_factors
)

amber_relaxer = relax.AmberRelaxation(
    **config.relax
)

# Relax the prediction.
t = time.time()
relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
print(f"Relaxation time: {time.time() - t}")

# Save the relaxed PDB.
output_dir = '.'
relaxed_output_path = os.path.join(output_dir, f'relaxed_{MODEL_NAME}.pdb')
with open(relaxed_output_path, 'w') as f:
    f.write(relaxed_pdb_str)
