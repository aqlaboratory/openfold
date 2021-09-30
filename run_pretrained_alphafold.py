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

import argparse
import math
import pickle
import os

# A hack to get OpenMM and PyTorch to peacefully coexist
os.environ["OPENMM_DEFAULT_PLATFORM"] = "OpenCL"

import time

import numpy as np
import torch
import torch.nn as nn

from config import model_config
from openfold.model.model import AlphaFold
from openfold.np import residue_constants, protein
import openfold.np.relax.relax as relax
from openfold.utils.import_weights import (
    import_jax_weights_,
)
from openfold.utils.tensor_utils import (
    tree_map,
    tensor_tree_map,
)

FEAT_PATH = "tests/test_data/sample_feats.pickle"

def main(args):
    config = model_config(args.model_name)
    model = AlphaFold(config.model)
    model = model.eval()
    import_jax_weights_(model, args.param_path)
    model = model.to(args.device)
    
    with open(FEAT_PATH, "rb") as f:
        batch = pickle.load(f)
    
    with torch.no_grad():
        batch = {
            k:torch.as_tensor(v, device=args.device) 
            for k,v in batch.items()
        }
        
        longs = [
            "aatype", 
            "template_aatype", 
            "extra_msa", 
            "residx_atom37_to_atom14",
            "residx_atom14_to_atom37",
            "true_msa",
            "residue_index",
        ]
        for l in longs:
            batch[l] = batch[l].long()
        
        # Move the recycling dimension to the end
        move_dim = lambda t: t.permute(*range(len(t.shape))[1:], 0)
        batch = tensor_tree_map(move_dim, batch)
        make_contig = lambda t: t.contiguous()
        batch = tensor_tree_map(make_contig, batch)
    
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
    
    os.environ["CUDA_VISIBLE_DEVICES"] = "7"
    
    amber_relaxer = relax.AmberRelaxation(
        **config.relax
    )
    
    # Relax the prediction.
    t = time.time()
    relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
    print(f"Relaxation time: {time.time() - t}")
    
    # Save the relaxed PDB.
    relaxed_output_path = os.path.join(
        args.output_dir, f'relaxed_{args.model_name}.pdb'
    )
    with open(relaxed_output_path, 'w') as f:
        f.write(relaxed_pdb_str)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--device", type=str, default="cpu",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--model_name", type=str, default="model_1",
        help="""Name of a model config. Choose one of model_{1-5} or 
             model_{1-5}_ptm, as defined on the AlphaFold GitHub."""
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction"""
    )
    parser.add_argument(
        "--param_path", type=str, default=None,
        help="""Path to model parameters. If None, parameters are selected
             automatically according to the model name from 
             openfold/resources/params"""
    )

    args = parser.parse_args()

    if(args.param_path is None):
        args.param_path = os.path.join(
            "openfold", "resources", "params", 
            "params_" + args.model_name + ".npz"
        )

    main(args)
