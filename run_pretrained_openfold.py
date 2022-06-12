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
from datetime import date
import gc
import logging
import numpy as np
import os

import pickle
from pytorch_lightning.utilities.deepspeed import (
    convert_zero_checkpoint_to_fp32_state_dict
)
import random
import sys
import time
import torch

from openfold.config import model_config
from openfold.data import templates, feature_pipeline, data_pipeline
from openfold.model.model import AlphaFold
from openfold.model.torchscript import script_preset_
from openfold.np import residue_constants, protein
import openfold.np.relax.relax as relax
from openfold.utils.import_weights import (
    import_jax_weights_,
)
from openfold.utils.tensor_utils import (
    tensor_tree_map,
)

from scripts.utils import add_data_args


def precompute_alignments(tags, seqs, alignment_dir, args):
    for tag, seq in zip(tags, seqs):
        tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)
        if(args.use_precomputed_alignments is None):
            logging.info(f"Generating alignments for {tag}...") 
            if not os.path.exists(local_alignment_dir):
                os.makedirs(local_alignment_dir)
            
            use_small_bfd=(args.bfd_database_path is None)
            alignment_runner = data_pipeline.AlignmentRunner(
                jackhmmer_binary_path=args.jackhmmer_binary_path,
                hhblits_binary_path=args.hhblits_binary_path,
                hhsearch_binary_path=args.hhsearch_binary_path,
                uniref90_database_path=args.uniref90_database_path,
                mgnify_database_path=args.mgnify_database_path,
                bfd_database_path=args.bfd_database_path,
                uniclust30_database_path=args.uniclust30_database_path,
                pdb70_database_path=args.pdb70_database_path,
                use_small_bfd=use_small_bfd,
                no_cpus=args.cpus,
            )
            alignment_runner.run(
                tmp_fasta_path, local_alignment_dir
            )

        # Remove temporary FASTA file
        os.remove(tmp_fasta_path)


def run_model(model, batch, tag, args):
    logging.info("Executing model...")
    with torch.no_grad():
        batch = {
            k:torch.as_tensor(v, device=args.model_device) 
            for k,v in batch.items()
        }
 
        # Disable templates if there aren't any in the batch
        model.config.template.enabled = any([
            "template_" in k for k in batch
        ])

        logging.info(f"Running inference for {tag}...")
        t = time.perf_counter()
        out = model(batch)
        logging.info(f"Inference time: {time.perf_counter() - t}")
    
    return out


def prep_output(out, batch, feature_dict, feature_processor, args):
    plddt = out["plddt"]
    mean_plddt = np.mean(plddt)
    
    plddt_b_factors = np.repeat(
        plddt[..., None], residue_constants.atom_type_num, axis=-1
    )

    # Prep protein metadata
    template_domain_names = []
    template_chain_index = None
    if(feature_processor.config.common.use_templates and "template_domain_names" in feature_dict):
        template_domain_names = [
            t.decode("utf-8") for t in feature_dict["template_domain_names"]
        ]

        # This works because templates are not shuffled during inference
        template_domain_names = template_domain_names[
            :feature_processor.config.predict.max_templates
        ]

        if("template_chain_index" in feature_dict):
            template_chain_index = feature_dict["template_chain_index"]
            template_chain_index = template_chain_index[
                :feature_processor.config.predict.max_templates
            ]

    no_recycling = feature_processor.config.common.max_recycling_iters
    remark = ', '.join([
        f"no_recycling={no_recycling}",
        f"max_templates={feature_processor.config.predict.max_templates}",
        f"config_preset={args.model_name}",
    ])

    # For multi-chain FASTAs
    ri = feature_dict["residue_index"]
    chain_index = (ri - np.arange(ri.shape[0])) / args.multimer_ri_gap
    chain_index = chain_index.astype(np.int64)
    cur_chain = 0
    prev_chain_max = 0
    for i, c in enumerate(chain_index):
        if(c != cur_chain):
            cur_chain = c
            prev_chain_max = i + cur_chain * args.multimer_ri_gap

        batch["residue_index"][i] -= prev_chain_max

    unrelaxed_protein = protein.from_prediction(
        features=batch,
        result=out,
        b_factors=plddt_b_factors,
        chain_index=chain_index,
        remark=remark,
        parents=template_domain_names,
        parents_chain_index=template_chain_index,
    )

    return unrelaxed_protein


def main(args):
    # Create the output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Prep the model
    config = model_config(args.model_name)
    model = AlphaFold(config)
    model = model.eval()

    if(args.jax_param_path):
        import_jax_weights_(
            model, args.jax_param_path, version=args.model_name
        )
    elif(args.openfold_checkpoint_path):
        if(os.path.isdir(args.openfold_checkpoint_path)):
            checkpoint_basename = os.path.splitext(
                os.path.basename(
                    os.path.normpath(args.openfold_checkpoint_path)
                )
            )[0]
            ckpt_path = os.path.join(
                args.output_dir,
                checkpoint_basename + ".pt",
            ) 

            if(not os.path.isfile(ckpt_path)):
                convert_zero_checkpoint_to_fp32_state_dict(
                    args.openfold_checkpoint_path,
                    ckpt_path,
                )
        else:
            ckpt_path = args.openfold_checkpoint_path

        d = torch.load(ckpt_path)
        model.load_state_dict(d["ema"]["params"])
    else:
        raise ValueError(
            "At least one of jax_param_path or openfold_checkpoint_path must "
            "be specified."
        )

    model = model.to(args.model_device)
 
    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=args.template_mmcif_dir,
        max_template_date=args.max_template_date,
        max_hits=config.data.predict.max_templates,
        kalign_binary_path=args.kalign_binary_path,
        release_dates_path=args.release_dates_path,
        obsolete_pdbs_path=args.obsolete_pdbs_path
    )

    data_processor = data_pipeline.DataPipeline(
        template_featurizer=template_featurizer,
    )

    output_dir_base = args.output_dir
    random_seed = args.data_random_seed
    if random_seed is None:
        random_seed = random.randrange(sys.maxsize)
    feature_processor = feature_pipeline.FeaturePipeline(config.data)
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)
    if(args.use_precomputed_alignments is None):
        alignment_dir = os.path.join(output_dir_base, "alignments")
    else:
        alignment_dir = args.use_precomputed_alignments

    prediction_dir = os.path.join(args.output_dir, "predictions")
    os.makedirs(prediction_dir, exist_ok=True)

    for fasta_file in os.listdir(args.fasta_dir):
        # Gather input sequences
        with open(os.path.join(args.fasta_dir, fasta_file), "r") as fp:
            data = fp.read()

        lines = [
            l.replace('\n', '') 
            for prot in data.split('>') for l in prot.strip().split('\n', 1)
        ][1:]
        tags, seqs = lines[::2], lines[1::2]

        tags = [t.split()[0] for t in tags]
        assert len(tags) == len(set(tags)), "All FASTA tags must be unique"
        tag = '-'.join(tags)

        precompute_alignments(tags, seqs, alignment_dir, args)

        tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
        if(len(seqs) == 1):
            seq = seqs[0]
            with open(tmp_fasta_path, "w") as fp:
                fp.write(f">{tag}\n{seq}")

            local_alignment_dir = os.path.join(alignment_dir, tag)
            feature_dict = data_processor.process_fasta(
                fasta_path=tmp_fasta_path, alignment_dir=local_alignment_dir
            )
        else:
            with open(tmp_fasta_path, "w") as fp:
                fp.write(
                    '\n'.join([f">{tag}\n{seq}" for tag, seq in zip(tags, seqs)])
                )
            feature_dict = data_processor.process_multiseq_fasta(
                fasta_path=tmp_fasta_path, super_alignment_dir=alignment_dir, 
            )
 
        # Remove temporary FASTA file
        os.remove(tmp_fasta_path)
    
        processed_feature_dict = feature_processor.process_features(
            feature_dict, mode='predict',
        )

        batch = processed_feature_dict
        out = run_model(model, batch, tag, args)

        # Toss out the recycling dimensions --- we don't need them anymore
        batch = tensor_tree_map(lambda x: np.array(x[..., -1].cpu()), batch)
        out = tensor_tree_map(lambda x: np.array(x.cpu()), out)
        
        unrelaxed_protein = prep_output(
            out, batch, feature_dict, feature_processor, args
        )

        output_name = f'{tag}_{args.model_name}'
        if(args.output_postfix is not None):
            output_name = f'{output_name}_{args.output_postfix}'

        # Save the unrelaxed PDB.
        unrelaxed_output_path = os.path.join(
            prediction_dir, f'{output_name}_unrelaxed.pdb'
        )
        with open(unrelaxed_output_path, 'w') as fp:
            fp.write(protein.to_pdb(unrelaxed_protein))

        if(not args.skip_relaxation):
            amber_relaxer = relax.AmberRelaxation(
                use_gpu=(args.model_device != "cpu"),
                **config.relax,
            )
            
            # Relax the prediction.
            t = time.perf_counter()
            visible_devices = os.getenv("CUDA_VISIBLE_DEVICES", default="")
            if("cuda" in args.model_device):
                device_no = args.model_device.split(":")[-1]
                os.environ["CUDA_VISIBLE_DEVICES"] = device_no
            relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
            os.environ["CUDA_VISIBLE_DEVICES"] = visible_devices
            logging.info(f"Relaxation time: {time.perf_counter() - t}")
            
            # Save the relaxed PDB.
            relaxed_output_path = os.path.join(
                prediction_dir, f'{output_name}_relaxed.pdb'
            )
            with open(relaxed_output_path, 'w') as fp:
                fp.write(relaxed_pdb_str)

        if(args.save_outputs):
            output_dict_path = os.path.join(
                args.output_dir, f'{output_name}_output_dict.pkl'
            )
            with open(output_dict_path, "wb") as fp:
                pickle.dump(out, fp, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_dir", type=str,
        help="Path to directory containing FASTA files, one sequence per file"
    )
    parser.add_argument(
        "template_mmcif_dir", type=str,
    )
    parser.add_argument(
        "--use_precomputed_alignments", type=str, default=None,
        help="""Path to alignment directory. If provided, alignment computation 
                is skipped and database path arguments are ignored."""
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction""",
    )
    parser.add_argument(
        "--model_device", type=str, default="cpu",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--model_name", type=str, default="model_1",
        help="""Name of a model config. Choose one of model_{1-5} or 
             model_{1-5}_ptm, as defined on the AlphaFold GitHub."""
    )
    parser.add_argument(
        "--jax_param_path", type=str, default=None,
        help="""Path to JAX model parameters. If None, and openfold_checkpoint_path
             is also None, parameters are selected automatically according to 
             the model name from openfold/resources/params"""
    )
    parser.add_argument(
        "--openfold_checkpoint_path", type=str, default=None,
        help="""Path to OpenFold checkpoint. Can be either a DeepSpeed 
             checkpoint directory or a .pt file"""
    )
    parser.add_argument(
        "--save_outputs", action="store_true", default=False,
        help="Whether to save all model outputs, including embeddings, etc."
    )
    parser.add_argument(
        "--cpus", type=int, default=4,
        help="""Number of CPUs with which to run alignment tools"""
    )
    parser.add_argument(
        "--preset", type=str, default='full_dbs',
        choices=('reduced_dbs', 'full_dbs')
    )
    parser.add_argument(
        "--output_postfix", type=str, default=None,
        help="""Postfix for output prediction filenames"""
    )
    parser.add_argument(
        "--data_random_seed", type=str, default=None
    )
    parser.add_argument(
        "--skip_relaxation", action="store_true", default=False,
    )
    parser.add_argument(
        "--multimer_ri_gap", type=int, default=200,
        help="""Residue index offset between multiple sequences, if provided"""
    )
    add_data_args(parser)
    args = parser.parse_args()

    if(args.jax_param_path is None and args.openfold_checkpoint_path is None):
        args.jax_param_path = os.path.join(
            "openfold", "resources", "params", 
            "params_" + args.model_name + ".npz"
        )

    if(args.model_device == "cpu" and torch.cuda.is_available()):
        logging.warning(
            """The model is being run on CPU. Consider specifying 
            --model_device for better performance"""
        )

    main(args)
