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
import pickle
import os

# A hack to get OpenMM and PyTorch to peacefully coexist
import random
import sys

from openfold.features import templates, feature_pipeline
from openfold.features.np import data_pipeline

os.environ["OPENMM_DEFAULT_PLATFORM"] = "OpenCL"

import time

import numpy as np
import torch

from config import model_config
from openfold.model.model import AlphaFold
from openfold.np import residue_constants, protein
import openfold.np.relax.relax as relax
from openfold.utils.import_weights import (
    import_jax_weights_,
)
from openfold.utils.tensor_utils import (
    tensor_tree_map,
)

FEAT_PATH = "tests/test_data/sample_feats.pickle"

MAX_TEMPLATE_HITS = 20

def main(args):
    config = model_config(args.model_name)
    model = AlphaFold(config.model)
    model = model.eval()
    import_jax_weights_(model, args.param_path)
    model = model.to(args.device)
    
    # FEATURE COLLECTION AND PROCESSING
    use_small_bfd = args.preset == "reduced_dbs"
    num_ensemble = 1

    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=args.template_mmcif_dir,
        max_template_date=args.max_template_date,
        max_hits=MAX_TEMPLATE_HITS,
        kalign_binary_path=args.kalign_binary_path,
        release_dates_path=None,
        obsolete_pdbs_path=args.obsolete_pdbs_path)

    data_processor = data_pipeline.DataPipeline(
        jackhammer_binary_path=args.jackhmmer_binary_path,
        hhblits_binary_path=args.hhblits_binary_path,
        hhsearch_binary_path=args.hhsearch_binary_path,
        uniref90_database_path=args.uniref90_database_path,
        mgnify_database_path=args.mgnify_database_path,
        bfd_database_path=args.bfd_database_path,
        uniclust30_database_path=args.uniclust30_database_path,
        small_bfd_database_path=args.small_bfd_database_path,
        pdb70_database_path=args.pdb70_database_path,
        template_featurizer=template_featurizer,
        use_small_bfd=use_small_bfd
    )

    output_dir_base = args.output_dir
    random_seed = args.random_seed
    if random_seed is None:
        random_seed = random.randrange(sys.maxsize)
    config.data.eval.num_ensemble = num_ensemble
    feature_processor = feature_pipeline.FeaturePipeline(config)
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)
    msa_output_dir = os.path.join(output_dir_base, "msas")
    if not os.path.exists(msa_output_dir):
        os.makedirs(msa_output_dir)

    print("Collecting data...")
    feature_dict = data_processor.process(
        input_fasta_path=args.fasta_path, msa_output_dir=msa_output_dir)
    # Output the features
    features_output_path = os.path.join(output_dir_base, 'features.pkl')
    with open(features_output_path, 'wb') as f:
        pickle.dump(feature_dict, f, protocol=4)

    print("Generating features...")
    processed_feature_dict = feature_processor.process_features(feature_dict, random_seed)

    with open(os.path.join(output_dir_base, 'processed_feats.pkl'), 'wb') as f:
        pickle.dump(processed_feature_dict, f, protocol=4)

    print("Executing model...")
    batch = processed_feature_dict
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
        "--fasta_path", type=str, default=None, required=True
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction""",
        required=True
    )
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
        "--param_path", type=str, default=None,
        help="""Path to model parameters. If None, parameters are selected
             automatically according to the model name from 
             openfold/resources/params"""
    )
    parser.add_argument(
        '--jackhmmer_binary_path', type=str, default='/usr/bin/jackhmmer'
    )
    parser.add_argument(
        '--hhblits_binary_path', type=str, default='/usr/bin/hhblits'
    )
    parser.add_argument(
        '--hhsearch_binary_path', type=str, default='/usr/bin/hhsearch'
    )
    parser.add_argument(
        '--kalign_binary_path', type=str, default='/usr/bin/kalign'
    )
    parser.add_argument('--uniref90_database_path', type=str, default=None, required=True
    )
    parser.add_argument(
        '--mgnify_database_path', type=str, default=None, required=True
    )
    parser.add_argument(
        '--bfd_database_path', type=str, default=None, required=True
    )
    parser.add_argument(
        '--small_bfd_database_path', type=str, default=None
    )
    parser.add_argument(
        '--uniclust30_database_path', type=str, default=None
    )
    parser.add_argument(
        '--pdb70_database_path', type=str, default=None, required=True
    )
    parser.add_argument(
        '--template_mmcif_dir', type=str, default=None, required=True
    )
    parser.add_argument(
        '--max_template_date', type=str, default=None, required=True
    )
    parser.add_argument(
        '--obsolete_pdbs_path', type=str, default=None
    )
    parser.add_argument(
        '--preset', type=str, default='full_dbs', required=True,
        choices=('reduced_dbs', 'full_dbs')
    )
    parser.add_argument(
        '--random_seed', type=str, default=None
    )

    args = parser.parse_args()

    if(args.param_path is None):
        args.param_path = os.path.join(
            "openfold", "resources", "params", 
            "params_" + args.model_name + ".npz"
        )

    main(args)
