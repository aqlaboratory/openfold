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

import argparse
import logging
import os
from shutil import copyfile
import tempfile
import traceback
import time
import subprocess

import glob
from mpi4py import MPI

from openfold.data.parsers import parse_fasta
from scripts.utils import add_data_args
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

def has_alignment(name, alignment_dir):
    alignment_files = [
        'mgnify_hits.a3m',
        'pdb70_hits.hhr',
        'small_bfd_hits.sto',
        'uniref90_hits.a3m'
    ]
    return all([os.path.isfile(os.path.join(alignment_dir, name, f)) for f in alignment_files])

def run_inference(runner, name, seq, subdir, args):
    with tempfile.TemporaryDirectory() as fasta_dir:
        logging.info(f"Temporal fasta dir {fasta_dir}")
        fd, fasta_path = tempfile.mkstemp(dir=fasta_dir, suffix=".fasta")
        with os.fdopen(fd, 'w') as fp:
            fp.write(f'>{name}\n{seq}')

        logging.info(f"Processing for {name} on {fasta_path}")
        ret = runner.run(fasta_dir, args.template_mmcif_dir, subdir, args, timeout=args.timeout)
        logging.info(f"Processing for {name} done!")
    return ret

def run_seq_group_inference(seq_groups, subdir_map, args):
    dirs = set(os.listdir(args.output_dir))
    pred_dir_base = os.path.join(args.output_dir, 'predictions')
    runner = OpenFoldInference(os.path.join(os.environ.get('OPENFOLDDIR'), 'run_pretrained_openfold.py'))

    comm = MPI.COMM_WORLD
    jobid = int(os.environ.get('PJM_JOBID', '0'))
    result_dir = os.path.join(args.output_dir, 'result')
    os.makedirs(result_dir, exist_ok=True)
    result_file_path = os.path.join(result_dir, f'job_{jobid}.csv')
    result_file = MPI.File.Open(comm, result_file_path, MPI.MODE_CREATE | MPI.MODE_WRONLY | MPI.MODE_APPEND)
    result_file.Set_atomicity(True)

    for seq, names in seq_groups:
        print("seq, names", seq, names)
        first_name = names[0]

        if subdir_map is None:
            pred_dir = pred_dir_base
            subdir = None
            alignment_dir = args.use_precomputed_alignments
        else:
            subdir = subdir_map[first_name]
            logging.info(f"Sub-directory of {first_name}: {subdir}")
            pred_dir = os.path.join(pred_dir_base, subdir)
            alignment_dir = os.path.join(args.use_precomputed_alignments, subdir) if args.use_precomputed_alignments else None

        if is_inferred(first_name, pred_dir, args):
            state = 'OK'
            duration = time_inference = time_relaxation = 0
        elif alignment_dir and (not has_alignment(first_name, alignment_dir)):
            state = 'NG_noalignment'
            duration = time_inference = time_relaxation = 0
        else:
            begin_time = time.time()
            try:
                ret = run_inference(runner, first_name, seq, subdir, args)
            except Exception as e:
                duration = time.time() - begin_time
                traceback.print_exc()
                logging.warning(f"Failed to run inference for {first_name}. Skipping...")
                if isinstance(e, subprocess.TimeoutExpired):
                    state = 'NG_timeout'
                else:
                    state = 'NG_unknown'
                time_inference = 0
                time_relaxation = 0
            else:
                duration = time.time() - begin_time
                state = 'OK'
                time_inference = ret['inference_time']
                time_relaxation = ret['relaxation_time']

        logging.info(f"inference_stat {first_name} {len(seq)} {state} {duration:.1f} {time_inference:.1f} {time_relaxation:.1f}")

        write_line = f'{first_name},{len(seq)},{state},{duration:.1f},{time_inference:.1f},{time_relaxation:.1f}\n'
        result_file.Write_shared(write_line.encode('utf-8'))
        result_file.Sync()

        if state != 'OK':
            continue

        generated_pdbs = [pdb_path(pred_dir, first_name, args.config_preset, False)]
        if not args.skip_relaxation:
            generated_pdbs.append(pdb_path(pred_dir, first_name, args.config_preset, True))

        for gen_file in generated_pdbs:
            if not os.path.isfile(gen_file):
                raise Exception(f'{gen_file} is not exist')

        for name in names[1:]:
            another_pred_dir = os.path.join(pred_dir_base, subdir_map[name])
            if not is_inferred(name, another_pred_dir, args):
                for f in generated_pdbs:
                        copy_file = os.path.join(another_pred_dir, '{}{}'.format(name, os.path.basename(f)[len(first_name):]))
                        logging.info(f"Copying result from {f} to {copy_file}")
                        os.makedirs(another_pred_dir, exist_ok=True)
                        copyfile(f, copy_file)

    result_file.Close()

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

def remove_non_target_seqs(input_seqs, input_chains, args):

    non_targets = set()

    if args.ignore_file:
        with open(args.ignore_file, 'r') as ignore_file:
            lines = ignore_file.readlines()
            non_targets |= set( [x.strip() for x in lines] )

    result_file_paths = glob.glob(os.path.join(args.output_dir, 'result', f'*.csv'))
    for result_file_path in result_file_paths:
        with open(result_file_path, 'r') as result_file:
            lines = result_file.readlines()
            non_targets |= set( [x.strip().split(',')[0] for x in lines] )

    logging.info(f'non_targets: {non_targets}')
    items = [(seq, chain) for seq, chain in zip(input_seqs, input_chains) if chain not in non_targets]

    return [x[0] for x in items], [x[1] for x in items]

def main(args):
    input_file = args.input_file
    with open(input_file, 'r') as fp:
        fasta_str = fp.read()
    input_seqs, input_chains = parse_fasta(fasta_str)
    orig_total_count = len(input_seqs)

    def to_first_lower(s):
        x = s.split(sep='_')
        x[0] = x[0].lower()
        return '_'.join(x)

    if args.first_lower:
        input_chains = [list(map(to_first_lower, g)) for g in input_chains]

    # Compute sub-directory mapping
    if args.sub_directory_size > 0:
        logging.info(f"Sub directory input/output is enabled")
        subdir_map = {}
        for i, name in enumerate(input_chains):
            subdir_map[name] = str(i//args.sub_directory_size)
    else:
        subdir_map = None

    input_seqs, input_chains = remove_non_target_seqs(input_seqs, input_chains, args)

    if not args.ignore_unique:
        input_seqs, input_chains = make_uniq_seq_groups(input_seqs, input_chains)
    else:
        logging.warning(f"--ignore_unique is enabled. The process might be redundant")
        input_chains = [[x] for x in input_chains]

    # sort by sequence length
    zip_seqs_chains = zip(input_seqs, input_chains)
    zip_seqs_chains_sorted = sorted(zip_seqs_chains, key=lambda x: len(x[0]))
    input_seqs, input_chains = zip(*zip_seqs_chains_sorted)

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
        subdir_map,
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
    parser.add_argument(
        "--first_lower",
        dest="first_lower", action="store_const",
        const=True, default=False,
        help="Convert first part of chain name to lower (e.g. 4X96_D -> 4x96_D)"
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
    parser.add_argument(
        "--max_memory", type=int, default=None,
        help="""Limit memory consumption"""
    )
    parser.add_argument(
        '--sub_directory_size', type=int, default=0,
        help="If this is set, create subdirectories for each number of sequences specified by this (default: 0)",
    )
    parser.add_argument(
        "--ignore_file", type=str, default=None,
        help="""The file of chain name list to ignore"""
    )

    args = parser.parse_args()

    main(args)
