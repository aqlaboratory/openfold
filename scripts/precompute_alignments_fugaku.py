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


def run_seq_group_alignments(seq_groups, alignment_runner, args):
    dirs = set(os.listdir(args.output_dir))
    for seq, name in seq_groups:
        alignment_dir = os.path.join(args.output_dir, name)

        if not os.path.exists(alignment_dir):
            try:
                os.makedirs(alignment_dir)
            except Exception as e:
                logging.warning(f"Failed to create directory for {name} with exception {e}...")
                continue

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
            )
            logging.info(f"Processing for {name} done!")
        except:
            traceback.print_exc()
            logging.warning(f"Failed to run alignments for {name}. Skipping...")
            os.remove(fasta_path)
            continue

        os.remove(fasta_path)


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


def make_uncompleted_seq_groups(input_seqs, input_chains, output_dir, alignment_runner):
    assert len(input_seqs) == len(input_chains)
    items = []
    for seq, chain in zip(input_seqs, input_chains):
        alignment_dir = os.path.join(output_dir, chain)
        dry_run = alignment_runner.dry_run(
            alignment_dir,
            input_label=chain,
            ignore_if_exists=True
        )
        if dry_run:
            items.append((seq, chain))

    return [x[0] for x in items], [x[1] for x in items]


def split_chains(input_seqs, input_chains, mpi_rank, mpi_size):
    return input_seqs[mpi_rank::mpi_size], input_chains[mpi_rank::mpi_size]


def allgather_chains(input_seqs, input_chains, comm):
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()

    assert len(input_seqs) == len(input_chains)
    my_count = np.array([len(input_seqs)], dtype=np.int32)
    proc_counts = np.zeros([mpi_size], dtype=np.int32)
    comm.Allgather(my_count, proc_counts)

    send_seqs = np.array(input_seqs, dtype=str)
    send_chains = np.array(input_chains, dtype=str)

    recv_seqs = []
    recv_chains = []
    # TODO: optimize
    for i, count in enumerate(proc_counts):
        recv_seq_buf, recv_chain_buf = \
            comm.bcast([send_seqs, send_chains] if mpi_rank == i else None, root=i)
        recv_seqs.extend(  list(recv_seq_buf  ))
        recv_chains.extend(list(recv_chain_buf))

    return recv_seqs, recv_chains


def main(args):
    # Build the alignment tool runner
    alignment_runner = AlignmentRunner(
        jackhmmer_binary_path=args.jackhmmer_binary_path,
        hhblits_binary_path=args.hhblits_binary_path,
        hhsearch_binary_path=args.hhsearch_binary_path,
        uniref90_database_path=args.uniref90_database_path,
        mgnify_database_path=args.mgnify_database_path,
        bfd_database_path=args.bfd_database_path,
        uniclust30_database_path=args.uniclust30_database_path,
        pdb70_database_path=args.pdb70_database_path,
        use_small_bfd=args.use_small_bfd,
        no_cpus=args.cpus_per_task,
        use_ok_files=True,
        disable_write_permission=True,
        timeout=args.timeout,
        bfd_deflate_a3m_hhm_path=args.bfd_deflate_a3m_hhm_path,
        uniclust30_deflate_a3m_hhm_path=args.uniclust30_deflate_a3m_hhm_path,
    )

    comm = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    assert mpi_size > 0
    assert mpi_rank >= 0 and mpi_rank < mpi_size

    input_file = args.input_file
    with open(input_file, 'r') as fp:
        fasta_str = fp.read()
    input_seqs, input_chains = parse_fasta(fasta_str)
    orig_total_count = len(input_seqs)

    # input_seqs   = [AAA, BBB, ...]
    # input_chains = [A, B, ...]

    # Remove completed chains
    input_seqs, input_chains = split_chains(input_seqs, input_chains, mpi_rank, mpi_size)
    input_seqs, input_chains = make_uncompleted_seq_groups(
        input_seqs,
        input_chains,
        args.output_dir,
        alignment_runner)

    input_seqs, input_chains = allgather_chains(input_seqs, input_chains, comm)

    uncompleted_total_count = len(input_seqs)

    if args.weak_scale:
        logging.warning(f"--weak_scale is enabled. The process might be redundant")
        assert len(input_seqs) == 1
        assert len(input_chains) == 1
        assert len(input_chains[0]) == 1
        input_seqs = [input_seqs[0]]*mpi_size
        input_chains = [[f"{input_chains[0][0]}_{i}"] for i in range(mpi_size)]

    total_count = len(input_seqs)
    input_seqs, input_chains = split_chains(input_seqs, input_chains, mpi_rank, mpi_size)

    logging.info(f"mpi_rank={mpi_rank}, mpi_size={mpi_size}, orig_total_count={orig_total_count}, uncompleted_count={uncompleted_total_count}, total_count={total_count}, my_count={len(input_seqs)}")
    logging.info(f"max memory: {args.max_memory}")
    logging.info(f"my chains: {input_chains}")

    run_seq_group_alignments(
        zip(input_seqs, input_chains),
        alignment_runner,
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
        "--use_small_bfd", type=bool, default=False,
    )

    parser.add_argument(
        "--weak_scale",
        dest="weak_scale", action="store_const",
        const=True, default=False,
        help="Duplicate the single input sequence for every process (use only for testing)"
    )
    parser.add_argument(
        "--bfd_deflate_a3m_hhm_path",
        type=str, default=None,
    )
    parser.add_argument(
        "--uniclust30_deflate_a3m_hhm_path",
        type=str, default=None,
    )
    parser.add_argument(
        '--timeout', type=float, default=None
    )
    parser.add_argument(
        '--max_memory', type=int, default=None
    )

    args = parser.parse_args()

    main(args)
