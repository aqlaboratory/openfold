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
import re

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
    other_basenames = [
        'mgnify_hits',
        'pdb70_hits',
        'uniref90_hits'
    ]
    bfd_basenames = [
        'small_bfd_hits',
        'bfd_uniclust_hits'
    ]

    basenames = [os.path.splitext(x)[0] for x in os.listdir(os.path.join(alignment_dir, name)) \
                 if os.path.isfile(os.path.join(alignment_dir, name, x))]

    return all([x in basenames for x in other_basenames]) and any([x in basenames for x in bfd_basenames])

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
    jobid = os.environ.get('PJM_JOBID', '0')
    result_dir = os.path.join(args.output_dir, 'result', jobid)
    os.makedirs(result_dir, exist_ok=True)
    result_file_path = os.path.join(result_dir, f'processed.csv')
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
                write_line = f'{name},{len(seq)},OK,0,0,0\n'
                result_file.Write_shared(write_line.encode('utf-8'))
                result_file.Sync()


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


def intersection_of_sets(set_map):
    if not set_map:
        return set()

    set_list = list(set_map.values())
    result_set = set_list[0]
    for s in set_list[1:]:
        result_set = result_set.intersection(s)
    return result_set

def get_chains(csv_file):
    with open(csv_file) as f:
        chains = f.read().strip().split("\n")

    return set(filter(None, chains))

def get_success_chains(root_dir, search_task):
    log_dirs = [x for x in os.listdir(root_dir) \
                if os.path.isdir(os.path.join(root_dir, x)) and x.endswith(f'.{search_task}')]
    log_dirs.sort(reverse=True)

    for log_dir in log_dirs:
        matcher = re.compile(r"chains_(.+)_(\d+)_(.+)\.csv")
        chain_files = [(x, matcher.match(x)) for x in os.listdir(os.path.join(root_dir, log_dir))]
        chain_files = [(x[0], x[1].groups()) for x in chain_files if x[1] is not None]

        # 何もなければこのディレクトリはスキップ
        if not chain_files:
            continue

        steps = {}
        for chain_file, (_, step, mode) in chain_files:
            info = steps.get(step, dict())
            info[mode] = chain_file
            steps[step] = info

        max_step = max(steps.keys())

        # 最新ステップの after_complete があれば、それを利用
        if 'after_complete' in steps[max_step]:
            return get_chains(os.path.join(root_dir, log_dir, steps[max_step]['after_complete']))

        # 最新ステップの after_complete が無ければ、before_complete + (あれば)success
        if 'before_complete' in steps[max_step]:
            chains = get_chains(os.path.join(root_dir, log_dir, steps[max_step]['before_complete']))
            if 'success' in steps[max_step]:
                chains |= get_chains(os.path.join(root_dir, log_dir, steps[max_step]['success']))
            return chains

    return set()

def get_subdir_map(root_dir):
    if not root_dir:
        return None

    log_dirs = [x for x in os.listdir(root_dir) \
                if os.path.isdir(os.path.join(root_dir, x))]
    log_dirs.sort(reverse=True)

    for log_dir in log_dirs:
        subdir_map_file = os.path.join(root_dir, log_dir, 'subdir_map.csv')
        if os.path.isfile(subdir_map_file):
            with open(subdir_map_file, 'r') as f:
                lines = f.read().strip().split("\n")
            lines = list(filter(None, lines))
            subdir_map = {k: v for k, v in [l.split(',') for l in lines]}

            logging.info(f'subdir_map.csv is found in alignment log directory. ({subdir_map_file})')
            return subdir_map
    return None

def get_alignment_completed_chains(root_dir):
    search_tasks = ['uniref90', 'small_bfd', 'pdb70', 'mgnify']
    completed_chains_map = {}

    for search_task in search_tasks:
        completed_chains_map[search_task] = get_success_chains(root_dir, search_task)

    all_completed_chains = intersection_of_sets(completed_chains_map)
    return all_completed_chains

def write_chains(args, timing, kind, chains):
    jobid = os.environ.get('PJM_JOBID', '0')
    result_dir = os.path.join(args.output_dir, 'result', jobid)
    os.makedirs(result_dir, exist_ok=True)
    result_file_path = os.path.join(result_dir, f'{timing}_{kind}.csv')

    with open(result_file_path, 'w') as f:
        f.writelines([f'{chain}\n' for chain in chains])

def remove_non_target_seqs(input_seqs, input_chains, args, rank):
    assert rank == 0

    ignore_chains = set()
    if args.ignore_file:
        with open(args.ignore_file, 'r') as ignore_file:
            lines = ignore_file.readlines()
            ignore_chains |= set( [x.strip() for x in lines] )

    result_dir = os.path.join(args.output_dir, 'result')
    os.makedirs(result_dir, exist_ok=True)

    def is_jobid(s):
        return True if re.fullmatch('[0-9]+', s, re.ASCII) else False
    job_dirs = [x for x in os.listdir(result_dir) \
                if os.path.isdir(os.path.join(result_dir, x)) and is_jobid(x)]

    jobids_with_result = [int(job_dir) for job_dir in job_dirs \
                          if os.path.isfile(os.path.join(result_dir, job_dir, 'processed.csv'))]

    completed_chains = set()
    skip_chains = set()

    for job_id in jobids_with_result:
        result_file_path = os.path.join(result_dir, str(job_id), 'processed.csv')
        proc_chains = {'OK': set(),
                       'NG_timeout': set(),
                       'NG_unknown': set(),
                       'NG_noalignment': set(),}

        with open(result_file_path, 'r') as result_file:
            lines = result_file.readlines()
        for l in lines:
            cols = l.strip().split(',')
            proc_chains[cols[2]].add(cols[0])

        completed_chains |= proc_chains['OK']
        if job_id > args.ignore_timeout_chain_history:
            skip_chains |= proc_chains['NG_timeout']
        if job_id > args.ignore_failed_chain_history:
            skip_chains |= proc_chains['NG_unknown']

    # Get alignemnt status
    if args.alignment_log_dir:
        alignment_completed_chains = get_alignment_completed_chains(args.alignment_log_dir)
    else:
        alignment_completed_chains = set()

    target_chains = set([c for c in input_chains if c not in ignore_chains])
    incompleted_chains = target_chains - completed_chains
    noalignment_chains = incompleted_chains - alignment_completed_chains
    write_chains(args, 'before', 'complete', completed_chains)
    write_chains(args, 'before', 'incomplete', incompleted_chains)
    write_chains(args, 'before', 'noalign', noalignment_chains)
    write_chains(args, 'before', 'skip', skip_chains)

    non_targets = set()
    non_targets |= ignore_chains
    non_targets |= completed_chains
    non_targets |= skip_chains
    non_targets |= noalignment_chains

    items = [(seq, chain) for seq, chain in zip(input_seqs, input_chains) \
             if chain not in non_targets]

    return [x[0] for x in items], [x[1] for x in items]

def main(args):
    mpi_rank = MPI.COMM_WORLD.Get_rank()
    mpi_size = MPI.COMM_WORLD.Get_size()

    if mpi_rank == 0:
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

        # Get or compute sub-directory mapping
        subdir_map = get_subdir_map(args.alignment_log_dir)

        if args.sub_directory_size > 0:
            if subdir_map:
                logging.info('subdir_map.csv in alignment log directory is used instead of --sub_directory_size')
            else:
                logging.info(f"Sub directory input/output is enabled (size: {args.sub_directory_size})")
                subdir_map = {name: str(i//args.sub_directory_size) \
                              for i, name in enumerate(input_chains)}

        if subdir_map is None:
            logging.info("no subdir_map")
        input_seqs, input_chains = remove_non_target_seqs(input_seqs, input_chains, args, mpi_rank)

        n_input_seqs = len(input_seqs)
        is_valid_subdir_map = set(input_chains).issubset(set(subdir_map.keys())) if subdir_map is not None else True
    else:
        orig_total_count = None
        n_input_seqs = None
        is_valid_subdir_map = None
        subdir_map = None

    orig_total_count = MPI.COMM_WORLD.bcast(orig_total_count, root=0)
    n_input_seqs = MPI.COMM_WORLD.bcast(n_input_seqs, root=0)
    subdir_map = MPI.COMM_WORLD.bcast(subdir_map, root=0)

    # there is no target seqs, exit
    if n_input_seqs == 0:
        logging.info("There is no sequences to be processed. DONE!")
        return

    is_valid_subdir_map = MPI.COMM_WORLD.bcast(is_valid_subdir_map, root=0)

    # check whehter subdir_map contains all input_chains
    if not is_valid_subdir_map:
        raise Exception(
            "The sub directory map must contain all input_chains sub directory."
        )

    if mpi_rank == 0:
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

        if args.weak_scale:
            logging.warning(f"--weak_scale is enabled. The process might be redundant")
            assert len(input_seqs) == 1
            assert len(input_chains) == 1
            assert len(input_chains[0]) == 1
            input_seqs = [input_seqs[0]]*mpi_size
            input_chains = [[f"{input_chains[0][0]}_{i}"] for i in range(mpi_size)]
    else:
        input_seqs = None
        input_chains = None

    input_seqs = MPI.COMM_WORLD.bcast(input_seqs, root=0)
    input_chains = MPI.COMM_WORLD.bcast(input_chains, root=0)

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
        help="If this is set, create subdirectories for each number of sequences specified by this (default: 0).",
    )
    parser.add_argument(
        '--alignment_log_dir', type=str, default=None,
        help="The log directory of alignment",
    )
    parser.add_argument(
        "--ignore_file", type=str, default=None,
        help="""The file of chain name list to ignore"""
    )
    parser.add_argument(
        '--ignore_timeout_chain_history', type=int, default=0,
        help="Ignores the history of timed out chains in jobs before the specified job ID. (default: 0).",
    )
    parser.add_argument(
        '--ignore_failed_chain_history', type=int, default=0,
        help="Ignores the history of failed chains in jobs before the specified job ID. (default: 0).",
    )

    args = parser.parse_args()

    main(args)
