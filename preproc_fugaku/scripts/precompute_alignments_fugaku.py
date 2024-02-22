# Copyright 2023 RIKEN & Fujitsu Limited
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
from functools import partial
import json
import logging
import os
import threading
from multiprocessing import cpu_count
from shutil import copyfile
import tempfile
import traceback
from mpi4py import MPI
import numpy as np

import os
os.environ["OPENFOLD_IGNORE_IMPORT"] = "1"

import openfold.data.mmcif_parsing as mmcif_parsing
from openfold.data.data_pipeline import AlignmentRunner
from openfold.data.parsers import parse_fasta
from openfold.np import protein, residue_constants

from utils import add_data_args


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s :%(message)s")


UNCOMPLETED_FLAG_DTYPE = bool
UNCOMPLETED_FLAG_AR_OP = MPI.LOR


def run_seq_group_alignments(seqs, alignment_runner, args, subdir_map=None):
    completed_count = 0
    total_count = 0
    for seq, names in seqs:
        if isinstance(names, str):
            names = [names]

        first_generated = None
        for i_name, name in enumerate(names):
            total_count += 1

            if subdir_map is None:
                alignment_dir = os.path.join(args.output_dir, name)
            else:
                alignment_dir = os.path.join(args.output_dir, subdir_map[name], name)
                logging.info(f"Sub-directory of {name}: {subdir_map[name]}")

            if not args.create_dir_on_demand:
                if not os.path.exists(alignment_dir):
                    try:
                        os.makedirs(alignment_dir, exist_ok=True)
                    except Exception as e:
                        if not os.path.exists(alignment_dir):
                            logging.warning(f"Failed to create directory for {name} with exception {e}...")
                            continue

            if i_name == 0:
                fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
                with os.fdopen(fd, 'w') as fp:
                    fp.write(f'>query\n{seq}')

                logging.info(f"Processing for {name} on {fasta_path}")
                try:
                    generated = alignment_runner.run(
                        fasta_path,
                        alignment_dir,
                        input_label=name,
                        ignore_if_exists=True,
                        max_memory=args.max_memory,
                        create_dir_on_demand=args.create_dir_on_demand,
                    )
                    if i_name == 0:
                        first_generated = generated

                    logging.info(f"Processing for {name} done!")
                    completed_count += 1

                except:
                    traceback.print_exc()
                    logging.warning(f"Failed to run alignments for {name}. Skipping...")

                os.remove(fasta_path)

            else:
                if first_generated is None:
                    logging.warning(f"Failed to run alignments for {name}. Skipping...")
                    continue

                first_name = names[0]
                logging.info(f"Linking already generated alignment for {name} from {first_name}")
                for f in first_generated:
                    os.symlink(
                        os.path.join("..", first_name, f),
                        os.path.join(alignment_dir, f))

                logging.info(f"Processing for {name} done!")
                completed_count += 1

    return completed_count, total_count


def add_unique_suffix(input_chains):
    """
    Returns a list of chain names with additional suffix so that every name become unique.

    Args:
        input_chains:
            A list of chain names
    Returns:
        A list of chain names
    """

    ret_chains = []
    known_set = set()
    input_set = set(input_chains)
    for c in input_chains:
       uc  = c
       i = 0
       while True:
           if uc not in known_set and not (i > 0 and uc in input_set):
               ret_chains.append(uc)
               known_set.add(uc)
               break

           else:
               uc = f"{c}_{i}"
               i += 1

    return ret_chains


def get_unique_seqs(input_seq_chains):
    """
    Returns a set of unique sequences.

    Args:
        input_seq_chains:
            A list of (seq., chain_name) tuples
    Returns:
        A list of (seq., [chain_name, ...]) tuples
    """

    input_chains = [x[1] for x in input_seq_chains]
    assert len(input_chains) == len(set(input_chains)) # Chain IDs must be unique

    s2c = {}
    for seq, chain in input_seq_chains:
        if seq not in s2c.keys():
            s2c[seq] = []

        s2c[seq].append(chain)

    for seq in s2c.keys():
        s2c[seq] = list(sorted(s2c[seq]))

    # items must be inter-process consistent as it is divided by processes
    return list(sorted(s2c.items(), key=lambda x: x[0]))


def get_uncompleted_flags(input_seq_chains, output_dir, alignment_runner):
    """
    Returns flags each of which means the search for the corresponding input sequence is already completed.

    Args:
        input_seq_chains:
            A list of (seq., chain_name) tuples. Must be identical among all ranks
        output_dir:
            Path to the root output directory
        alignment_runner:
            The alignment runner
    Returns:
        A flag tensor, whose shape is equivalent to that of input_seq_chains
    """
    flags = np.zeros([len(input_seq_chains)], dtype=UNCOMPLETED_FLAG_DTYPE)
    for i, (seq, chain) in enumerate(input_seq_chains):
        alignment_dir = os.path.join(output_dir, chain)
        dry_run = alignment_runner.dry_run(
            alignment_dir,
            input_label=chain,
            ignore_if_exists=True
        )
        if dry_run:
            flags[i] = 1

    return flags


def get_uncompleted_seqs(input_seq_chains, comm, alignment_runner):
    """
    Check whether search for each input sequence is already completed,
    and returns equally-split uncompleted sequences.

    Args:
        input_seq_chains:
            A list of (seq., chain_name) tuples. Must be identical among all ranks
        comm:
            mpi4py communicator
        alignment_runner:
            The alignment runner
    Returns:
        A list of (seq., chain_name) tuples
    """
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()

    proc_begin = int(len(input_seq_chains)*mpi_rank/mpi_size)
    proc_end   = int(len(input_seq_chains)*(mpi_rank+1)/mpi_size)
    proc_uncompleted_flags = get_uncompleted_flags(
        input_seq_chains[proc_begin:proc_end],
        args.output_dir,
        alignment_runner)
    send_uncomplted_flags = np.zeros([len(input_seq_chains)], dtype=UNCOMPLETED_FLAG_DTYPE)
    recv_uncomplted_flags = np.zeros([len(input_seq_chains)], dtype=UNCOMPLETED_FLAG_DTYPE)
    send_uncomplted_flags[proc_begin:proc_end] = proc_uncompleted_flags
    comm.Allreduce(
        send_uncomplted_flags,
        recv_uncomplted_flags,
        op=UNCOMPLETED_FLAG_AR_OP,
    )
    return [x[0] for x in zip(input_seq_chains, recv_uncomplted_flags) \
            if x[1] > 0]


def main(args):

    comm = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    assert mpi_size > 0
    assert mpi_rank >= 0 and mpi_rank < mpi_size

    assert len(args.temp_dir) > 0
    my_temp_dir = os.path.join(args.temp_dir, f"{mpi_rank}")
    os.makedirs(my_temp_dir, exist_ok=True)

    # Build the alignment tool runner
    alignment_runner = AlignmentRunner(
        jackhmmer_binary_path=args.jackhmmer_binary_path,
        hhblits_binary_path=args.hhblits_binary_path,
        hhsearch_binary_path=args.hhsearch_binary_path,
        uniref90_database_path=args.uniref90_database_path,
        mgnify_database_path=args.mgnify_database_path,
        bfd_database_path=args.bfd_database_path,
        uniclust30_database_path=None,
        pdb70_database_path=args.pdb70_database_path,
        use_small_bfd=True,
        convert_small_bfd_to_a3m=args.convert_small_bfd_to_a3m,
        no_cpus=args.cpus_per_task,
        disable_write_permission=args.disable_write_permission,
        timeout=args.timeout,
        stream_sto_size=args.stream_sto_size,
        uniref_max_hits=args.uniref90_max_hits,
        mgnify_max_hits=args.mgnify_max_hits,
        small_bfd_max_hits=args.small_bfd_max_hits,
        temp_dir=my_temp_dir,
    )

    input_file = args.input_file
    with open(input_file, 'r') as fp:
        fasta_str = fp.read()

    input_seqs, input_chains = parse_fasta(fasta_str)
    input_chains = add_unique_suffix(input_chains)
    input_seq_chains = list(zip(input_seqs, input_chains)) # [(AAAAA, name1), (BBBB, name2), ..]
    orig_total_count = len(input_seq_chains)

    # Compute sub-directory mapping
    if args.sub_directory_size > 0:
        subdir_map = {}
        for i, (_, name) in enumerate(input_seq_chains):
            subdir_map[name] = str(i//args.sub_directory_size)
            
    else:
        subdir_map = None

    # Remove completed chains
    input_seq_chains = get_uncompleted_seqs(input_seq_chains, comm, alignment_runner)

    # Remove duplicated seqs.
    if args.unique:
        if mpi_rank == 0:
            logging.warning(
                "The --unique option might be slow because all ranks compute unique sequences. "
                "Consider removing duplicated sequences from the input file manually.")
        input_seq_chains = get_unique_seqs(input_seq_chains)

    uncompleted_total_count = len(input_seq_chains)

    # Distribute uncompleted chains
    input_seq_chains = input_seq_chains[mpi_rank::mpi_size]

    host = os.environ["HOSTNAME"]
    logging.info(f"host={host}, rank={mpi_rank}/{mpi_size}, "
                 f"total_count={orig_total_count}, "
                 f"total_uncompleted_count={uncompleted_total_count}, "
                 f"my_count={len(input_seq_chains)}, "
                 f"my_temp_dir={my_temp_dir}")

    completed_count, total_count = run_seq_group_alignments(
        input_seq_chains,
        alignment_runner,
        args,
        subdir_map=subdir_map)

    logging.info(f"DONE! "
                 f"host={host}, rank={mpi_rank}/{mpi_size}, "
                 f"my_completed_count={completed_count}, "
                 f"my_count={total_count}")

    completed_count = comm.allreduce(completed_count)
    total_count = comm.allreduce(total_count)

    if mpi_rank == 0:
        if args.report_out_path is not None:
            remaining = total_count-completed_count
            assert remaining >= 0
            with open(args.report_out_path, "w") as f:
                f.write(str(remaining))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file", type=str,
        help="""The input FASTA file"""
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Directory in which to output alignments"
    )
    add_data_args(parser)
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )
    parser.add_argument(
        "--cpus_per_task", type=int, default=cpu_count(),
        help="Number of CPUs to use"
    )
    parser.add_argument(
        "--mmcif_cache", type=str, default=None,
        help="Path to mmCIF cache. Used to filter files to be parsed"
    )
    parser.add_argument(
        "--no_tasks", type=int, default=1,
    )
    parser.add_argument(
        "--filter", type=bool, default=True,
    )

    parser.add_argument(
        '--timeout', type=float, default=None,
        help="Time limit for each search tool in seconds (default: None)",
    )
    parser.add_argument(
        '--max_memory', type=int, default=None,
        help="The RLIMIT_AS memory limit for each search tool in bytes (default: None)",
    )
    parser.add_argument(
        '--stream-sto-size', type=int, default=1024*1024*1024,
        help="Use the stream version of sto-to-a3m conversion if sto file size is larger than this size (default: 1 GiB)",
    )
    parser.add_argument(
        '--sub-directory-size', type=int, default=0,
        help="If this is set, create subdirectories for each number of sequences specified by this (default: 0)",
    )
    parser.add_argument(
        '--report_out_path', type=str, default=None,
        help="Path to output the number of uncompleted seqs. (default: None)",
    )
    parser.add_argument(
        "--unique",
        dest="unique",
        default=False,
        const=True,
        action="store_const",
        help="Find duplicated sequences and create symlinks to existing alignment files "
        "instead of running search tools (default: False)",
    )
    parser.add_argument(
        "--disable-write-permission",
        dest="disable_write_permission",
        default=False,
        const=True,
        action="store_const",
        help="Set permission 440 to output MSA and template files",
    )
    parser.add_argument(
        "--convert-small-bfd-to-a3m",
        dest="convert_small_bfd_to_a3m",
        default=False,
        const=True,
        action="store_const",
        help="Convert small BFD MSAs from STO to A3M",
    )
    parser.add_argument(
        "--create-dir-on-demand",
        dest="create_dir_on_demand",
        default=False,
        const=True,
        action="store_const",
        help="Create output directories only if alignment is succeeded",
    )
    parser.add_argument(
        "--uniref90-max-hits", type=int, default=10000,
        help="The maximum number of MSA hits on UniRef90 (default: 10000)",
    )
    parser.add_argument(
        "--mgnify-max-hits", type=int, default=5000,
        help="The maximum number of MSA hits on MGnify (default: 5000)",
    )
    parser.add_argument(
        "--small-bfd-max-hits", type=int, default=None,
        help="The maximum number of MSA hits on small BFD (default: unlimited)",
    )
    parser.add_argument(
        "--temp-dir", type=str, default="/tmp",
        help="Path to the temporary directory (default: /tmp)",
    )

    args = parser.parse_args()

    main(args)
