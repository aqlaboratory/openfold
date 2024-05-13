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
import logging
import math
import numpy as np
import os
import pickle
import random
import time
import json

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

import torch
torch_versions = torch.__version__.split(".")
torch_major_version = int(torch_versions[0])
torch_minor_version = int(torch_versions[1])
if (
    torch_major_version > 1 or
    (torch_major_version == 1 and torch_minor_version >= 12)
):
    # Gives a large speedup on Ampere-class GPUs
    torch.set_float32_matmul_precision("high")

torch.set_grad_enabled(False)

from openfold.config import model_config
from openfold.data import templates, feature_pipeline, data_pipeline
from openfold.data.tools import hhsearch, hmmsearch
from openfold.np import protein
from openfold.utils.script_utils import (load_models_from_command_line, parse_fasta, run_model,
                                         prep_output, relax_protein)
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.trace_utils import (
    pad_feature_dict_seq,
    trace_model_,
)

from scripts.precompute_embeddings import EmbeddingGenerator
from scripts.utils import add_data_args


TRACING_INTERVAL = 50


def precompute_alignments(tags, seqs, alignment_dir, args):
    for tag, seq in zip(tags, seqs):
        tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)

        if args.use_precomputed_alignments is None:
            logger.info(f"Generating alignments for {tag}...")

            os.makedirs(local_alignment_dir, exist_ok=True)

            if "multimer" in args.config_preset:
                template_searcher = hmmsearch.Hmmsearch(
                    binary_path=args.hmmsearch_binary_path,
                    hmmbuild_binary_path=args.hmmbuild_binary_path,
                    database_path=args.pdb_seqres_database_path,
                )
            else:
                template_searcher = hhsearch.HHSearch(
                    binary_path=args.hhsearch_binary_path,
                    databases=[args.pdb70_database_path],
                )

            # In seqemb mode, use AlignmentRunner only to generate templates
            if args.use_single_seq_mode:
                alignment_runner = data_pipeline.AlignmentRunner(
                    jackhmmer_binary_path=args.jackhmmer_binary_path,
                    uniref90_database_path=args.uniref90_database_path,
                    template_searcher=template_searcher,
                    no_cpus=args.cpus,
                )
                embedding_generator = EmbeddingGenerator()
                embedding_generator.run(tmp_fasta_path, alignment_dir)
            else:
                alignment_runner = data_pipeline.AlignmentRunner(
                    jackhmmer_binary_path=args.jackhmmer_binary_path,
                    hhblits_binary_path=args.hhblits_binary_path,
                    uniref90_database_path=args.uniref90_database_path,
                    mgnify_database_path=args.mgnify_database_path,
                    bfd_database_path=args.bfd_database_path,
                    uniref30_database_path=args.uniref30_database_path,
                    uniclust30_database_path=args.uniclust30_database_path,
                    uniprot_database_path=args.uniprot_database_path,
                    template_searcher=template_searcher,
                    use_small_bfd=args.bfd_database_path is None,
                    no_cpus=args.cpus
                )

            alignment_runner.run(
                tmp_fasta_path, local_alignment_dir
            )
        else:
            logger.info(
                f"Using precomputed alignments for {tag} at {alignment_dir}..."
            )

        # Remove temporary FASTA file
        os.remove(tmp_fasta_path)


def round_up_seqlen(seqlen):
    return int(math.ceil(seqlen / TRACING_INTERVAL)) * TRACING_INTERVAL


def generate_feature_dict(
    tags,
    seqs,
    alignment_dir,
    data_processor,
    args,
):
    tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")

    if "multimer" in args.config_preset:
        with open(tmp_fasta_path, "w") as fp:
            fp.write(
                '\n'.join([f">{tag}\n{seq}" for tag, seq in zip(tags, seqs)])
            )
        feature_dict = data_processor.process_fasta(
            fasta_path=tmp_fasta_path, alignment_dir=alignment_dir,
        )
    elif len(seqs) == 1:
        tag = tags[0]
        seq = seqs[0]
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)
        feature_dict = data_processor.process_fasta(
            fasta_path=tmp_fasta_path,
            alignment_dir=local_alignment_dir,
            seqemb_mode=args.use_single_seq_mode,
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

    return feature_dict


def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]


def main(args):
    # Create the output directory
    os.makedirs(args.output_dir, exist_ok=True)

    if args.config_preset.startswith("seq"):
        args.use_single_seq_mode = True

    config = model_config(
        args.config_preset, 
        long_sequence_inference=args.long_sequence_inference,
        use_deepspeed_evoformer_attention=args.use_deepspeed_evoformer_attention,
        )

    if args.experiment_config_json: 
        with open(args.experiment_config_json, 'r') as f:
            custom_config_dict = json.load(f)
        config.update_from_flattened_dict(custom_config_dict)

    if args.experiment_config_json: 
        with open(args.experiment_config_json, 'r') as f:
            custom_config_dict = json.load(f)
        config.update_from_flattened_dict(custom_config_dict)

    if args.trace_model:
        if not config.data.predict.fixed_size:
            raise ValueError(
                "Tracing requires that fixed_size mode be enabled in the config"
            )

    is_multimer = "multimer" in args.config_preset

    if is_multimer:
        template_featurizer = templates.HmmsearchHitFeaturizer(
            mmcif_dir=args.template_mmcif_dir,
            max_template_date=args.max_template_date,
            max_hits=config.data.predict.max_templates,
            kalign_binary_path=args.kalign_binary_path,
            release_dates_path=args.release_dates_path,
            obsolete_pdbs_path=args.obsolete_pdbs_path
        )
    else:
        template_featurizer = templates.HhsearchHitFeaturizer(
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

    if is_multimer:
        data_processor = data_pipeline.DataPipelineMultimer(
            monomer_data_pipeline=data_processor,
        )

    output_dir_base = args.output_dir
    random_seed = args.data_random_seed
    if random_seed is None:
        random_seed = random.randrange(2 ** 32)

    np.random.seed(random_seed)
    torch.manual_seed(random_seed + 1)

    feature_processor = feature_pipeline.FeaturePipeline(config.data)
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)
    if args.use_precomputed_alignments is None:
        alignment_dir = os.path.join(output_dir_base, "alignments")
    else:
        alignment_dir = args.use_precomputed_alignments

    tag_list = []
    seq_list = []
    for fasta_file in list_files_with_extensions(args.fasta_dir, (".fasta", ".fa")):
        # Gather input sequences
        fasta_path = os.path.join(args.fasta_dir, fasta_file)
        with open(fasta_path, "r") as fp:
            data = fp.read()

        tags, seqs = parse_fasta(data)

        if not is_multimer and len(tags) != 1:
            print(
                f"{fasta_path} contains more than one sequence but "
                f"multimer mode is not enabled. Skipping..."
            )
            continue

        # assert len(tags) == len(set(tags)), "All FASTA tags must be unique"
        tag = '-'.join(tags)

        tag_list.append((tag, tags))
        seq_list.append(seqs)

    seq_sort_fn = lambda target: sum([len(s) for s in target[1]])
    sorted_targets = sorted(zip(tag_list, seq_list), key=seq_sort_fn)
    feature_dicts = {}

    if is_multimer and args.openfold_checkpoint_path:
        raise ValueError(
            '`openfold_checkpoint_path` was specified, but no OpenFold checkpoints are available for multimer mode')

    model_generator = load_models_from_command_line(
        config,
        args.model_device,
        args.openfold_checkpoint_path,
        args.jax_param_path,
        args.output_dir)

    for model, output_directory in model_generator:
        cur_tracing_interval = 0
        for (tag, tags), seqs in sorted_targets:
            output_name = f'{tag}_{args.config_preset}'
            if args.output_postfix is not None:
                output_name = f'{output_name}_{args.output_postfix}'

            # Does nothing if the alignments have already been computed
            precompute_alignments(tags, seqs, alignment_dir, args)

            feature_dict = feature_dicts.get(tag, None)
            if feature_dict is None:
                feature_dict = generate_feature_dict(
                    tags,
                    seqs,
                    alignment_dir,
                    data_processor,
                    args,
                )

                if args.trace_model:
                    n = feature_dict["aatype"].shape[-2]
                    rounded_seqlen = round_up_seqlen(n)
                    feature_dict = pad_feature_dict_seq(
                        feature_dict, rounded_seqlen,
                    )

                feature_dicts[tag] = feature_dict

            processed_feature_dict = feature_processor.process_features(
                feature_dict, mode='predict', is_multimer=is_multimer
            )

            processed_feature_dict = {
                k: torch.as_tensor(v, device=args.model_device)
                for k, v in processed_feature_dict.items()
            }

            if args.trace_model:
                if rounded_seqlen > cur_tracing_interval:
                    logger.info(
                        f"Tracing model at {rounded_seqlen} residues..."
                    )
                    t = time.perf_counter()
                    trace_model_(model, processed_feature_dict)
                    tracing_time = time.perf_counter() - t
                    logger.info(
                        f"Tracing time: {tracing_time}"
                    )
                    cur_tracing_interval = rounded_seqlen

            out = run_model(model, processed_feature_dict, tag, args.output_dir)

            # Toss out the recycling dimensions --- we don't need them anymore
            processed_feature_dict = tensor_tree_map(
                lambda x: np.array(x[..., -1].cpu()),
                processed_feature_dict
            )
            out = tensor_tree_map(lambda x: np.array(x.cpu()), out)

            unrelaxed_protein = prep_output(
                out,
                processed_feature_dict,
                feature_dict,
                feature_processor,
                args.config_preset,
                args.multimer_ri_gap,
                args.subtract_plddt
            )

            unrelaxed_file_suffix = "_unrelaxed.pdb"
            if args.cif_output:
                unrelaxed_file_suffix = "_unrelaxed.cif"
            unrelaxed_output_path = os.path.join(
                output_directory, f'{output_name}{unrelaxed_file_suffix}'
            )

            with open(unrelaxed_output_path, 'w') as fp:
                if args.cif_output:
                    fp.write(protein.to_modelcif(unrelaxed_protein))
                else:
                    fp.write(protein.to_pdb(unrelaxed_protein))

            logger.info(f"Output written to {unrelaxed_output_path}...")

            if not args.skip_relaxation:
                # Relax the prediction.
                logger.info(f"Running relaxation on {unrelaxed_output_path}...")
                relax_protein(config, args.model_device, unrelaxed_protein, output_directory, output_name,
                              args.cif_output)

            if args.save_outputs:
                output_dict_path = os.path.join(
                    output_directory, f'{output_name}_output_dict.pkl'
                )
                with open(output_dict_path, "wb") as fp:
                    pickle.dump(out, fp, protocol=pickle.HIGHEST_PROTOCOL)

                logger.info(f"Model output written to {output_dict_path}...")


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
        "--use_single_seq_mode", action="store_true", default=False,
        help="""Use single sequence embeddings instead of MSAs."""
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
        "--config_preset", type=str, default="model_1",
        help="""Name of a model config preset defined in openfold/config.py"""
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
        "--data_random_seed", type=int, default=None
    )
    parser.add_argument(
        "--skip_relaxation", action="store_true", default=False,
    )
    parser.add_argument(
        "--multimer_ri_gap", type=int, default=200,
        help="""Residue index offset between multiple sequences, if provided"""
    )
    parser.add_argument(
        "--trace_model", action="store_true", default=False,
        help="""Whether to convert parts of each model to TorchScript.
                Significantly improves runtime at the cost of lengthy
                'compilation.' Useful for large batch jobs."""
    )
    parser.add_argument(
        "--subtract_plddt", action="store_true", default=False,
        help=""""Whether to output (100 - pLDDT) in the B-factor column instead
                 of the pLDDT itself"""
    )
    parser.add_argument(
        "--long_sequence_inference", action="store_true", default=False,
        help="""enable options to reduce memory usage at the cost of speed, helps longer sequences fit into GPU memory, see the README for details"""
    )
    parser.add_argument(
        "--cif_output", action="store_true", default=False,
        help="Output predicted models in ModelCIF format instead of PDB format (default)"
    )
    parser.add_argument(
        "--experiment_config_json", default="", help="Path to a json file with custom config values to overwrite config setting",
    )
    parser.add_argument(
        "--use_deepspeed_evoformer_attention", action="store_true", default=False, 
        help="Whether to use the DeepSpeed evoformer attention layer. Must have deepspeed installed in the environment.",
    )
    add_data_args(parser)
    args = parser.parse_args()

    if args.jax_param_path is None and args.openfold_checkpoint_path is None:
        args.jax_param_path = os.path.join(
            "openfold", "resources", "params",
            "params_" + args.config_preset + ".npz"
        )

    if args.model_device == "cpu" and torch.cuda.is_available():
        logging.warning(
            """The model is being run on CPU. Consider specifying 
            --model_device for better performance"""
        )

    main(args)
