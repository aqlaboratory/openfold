import argparse
from functools import partial
import json
import logging
import os
import threading
from multiprocessing import cpu_count
from shutil import copyfile
import tempfile
import traceback

import openfold.data.mmcif_parsing as mmcif_parsing
from openfold.data.data_pipeline import AlignmentRunner
from openfold.data.parsers import parse_fasta
from openfold.np import protein, residue_constants

from scripts.utils import add_data_args
import glob

from scripts.openfold_runner import OpenFoldInference


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s :%(message)s")


def pdb_path(pred_dir, name, model, relaxed):
    if relaxed:
        return os.path.join(pred_dir, '{}_{}_relaxed.pdb'.format(name, model))
    else:
        return os.path.join(pred_dir, '{}_{}_unrelaxed.pdb'.format(name, model))

def is_inferred(name, pred_dir, args):
    # e.g. 4X96_D_model_1_relaxed.pdb
    unrelaxed_pdb_path = pdb_path(pred_dir, name, args.config_preset, False)
    relaxed_pdb_path = pdb_path(pred_dir, name, args.config_preset, True)

    return os.path.isfile(unrelaxed_pdb_path) and (args.skip_relaxation or os.path.isfile(relaxed_pdb_path))

def run_inference(runner, name, seq, args):
    with tempfile.TemporaryDirectory() as fasta_dir:
        logging.info(f"Temporal fasta dir {fasta_dir}")
        fd, fasta_path = tempfile.mkstemp(dir=fasta_dir, suffix=".fasta")
        with os.fdopen(fd, 'w') as fp:
            fp.write(f'>{name}\n{seq}')

        logging.info(f"Processing for {name} on {fasta_path}")
        runner.run(fasta_dir, args.template_mmcif_dir, args, timeout=args.timeout)
        logging.info(f"Processing for {name} done!")

    generated_files = [pdb_path(pred_dir, name, args.config_preset, False)]
    if not args.skip_relaxation:
        generated_files.append(pdb_path(pred_dir, name, args.config_preset, True))

    for gen_file in generated_files:
        if os.path.isfile(gen_file):
            logging.info(f"Generated pdb file: {gen_file}")
        else:
            raise Exception(f'{gen_File} is not exist')
    return generated_files

def run_seq_group_inference(seq_groups, args):
    dirs = set(os.listdir(args.output_dir))
    pred_dir = os.path.join(args.output_dir, 'predictions')
    runner = OpenFoldInference('run_pretrained_openfold.py', os.environ.get('OPENFOLDDIR'))

    for seq, names in seq_groups:
        print("seq, names", seq, names)
        first_name = names[0]

        if not is_inferred(first_name, pred_dir, args):
            try:
                generated_files = run_inference(runner, first_name, seq, args)
            except:
                traceback.print_exc()
                logging.warning(f"Failed to run inference for {first_name}. Skipping...")
                continue

        for name in names[1:]:
            if not is_inferred(name, pred_Dir, args):
                for f in generated_files:
                        copy_file = os.path.join(pred_dir, '{}{}'.format(name, os.path.basename(f)[len(first_name):]))
                        logging.info(f"Copying result from {f} to {copy_file}")
                        copyfile(f, copy_file)


def make_uniq_seq_groups(input_seqs, input_chains):
    assert len(input_seqs) == len(input_chains)
    assert len(input_chains) == len(set(input_chains)) # Chain IDs must be unique

    dic = {}
    for seq, chain in zip(input_seqs, input_chains):
        if seq not in dic.keys():
            dic[seq] = []

        dic[seq].append(chain)

    for k in dic.keys():
        dic[k] = sorted(dic[k])

    items = list(dic.items())

    # items must be inter-process consistent as it is divided by processes
    items = sorted(items, key=lambda x: x[1])

    return [x[0] for x in items], [x[1] for x in items]

def main(args):
    input_file = args.input_file
    with open(input_file, 'r') as fp:
        fasta_str = fp.read()
    input_seqs, input_chains = parse_fasta(fasta_str)
    orig_total_count = len(input_seqs)

    if not args.ignore_unique:
        input_seqs, input_chains = make_uniq_seq_groups(input_seqs, input_chains)
    else:
        logging.warning(f"--ignore_unique is enabled. The process might be redundant")
        input_chains = [[x] for x in input_chains]

    # input_seqs   = [AAA, BBB, ...]
    # input_chains = [[A_1, A_2], [B_1], ...]

    if "OMPI_COMM_WORLD_RANK" in os.environ:
        # ABCI (OpenMPI)
        mpi_rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
        mpi_size = int(os.environ["OMPI_COMM_WORLD_SIZE"])

    elif "PMIX_RANK" in os.environ:
        # Fugaku (Fujitsu MPI)
        mpi_rank = int(os.environ["PMIX_RANK"])
        mpi_size = int(os.environ["OMPI_UNIVERSE_SIZE"])

    else:
        logging.warning("MPI rank/size environment variables not found")
        mpi_rank = 0
        mpi_size = 1

    if args.weak_scale:
        logging.warning(f"--weak_scale is enabled. The process might be redundant")
        assert len(input_seqs) == 1
        assert len(input_chains) == 1
        assert len(input_chains[0]) == 1
        input_seqs = [input_seqs[0]]*mpi_size
        input_chains = [[f"{input_chains[0][0]}_{i}"] for i in range(mpi_size)]

    assert mpi_size > 0
    assert mpi_rank >= 0 and mpi_rank < mpi_size
    total_count  = len(input_seqs)
    input_seqs   =   input_seqs[mpi_rank::mpi_size]
    input_chains = input_chains[mpi_rank::mpi_size]

    logging.info(f"mpi_rank={mpi_rank}, mpi_size={mpi_size}, orig_total_count={orig_total_count}, total_count={total_count}, my_count={len(input_seqs)}")
    logging.info(f"my chains: {input_chains}")

    run_seq_group_inference(
        zip(input_seqs, input_chains),
        args)

    logging.info("DONE!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file", type=str,
        help="""The input FASTA file"""
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Directory in which to output pdb"
    )
    parser.add_argument(
        "template_mmcif_dir", type=str,
    )

    add_data_args(parser)
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )
    parser.add_argument(
        "--ignore_unique",
        dest="ignore_unique", action="store_const",
        const=True, default=False,
        help="Do not merge the same sequences (use only for testing)"
    )
    parser.add_argument(
        "--weak_scale",
        dest="weak_scale", action="store_const",
        const=True, default=False,
        help="Duplicate the single input sequence for every process (use only for testing)"
    )
    parser.add_argument(
        '--timeout', type=float, default=None
    )

    
    # for inference
    parser.add_argument(
        "--use_precomputed_alignments", type=str, default=None,
        help="""Path to alignment directory. If provided, alignment computation 
                is skipped and database path arguments are ignored."""
    )
    parser.add_argument(
        "--model_device", type=str, default="cpu",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--config_preset", type=str, default="model_1",
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

    args = parser.parse_args()

    main(args)
