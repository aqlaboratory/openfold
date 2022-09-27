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

from utils import add_data_args


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s :%(message)s")


def run_seq_group_alignments(seq_groups, alignment_runner, args):
    dirs = set(os.listdir(args.output_dir))
    for seq, names in seq_groups:
        first_name = names[0]
        alignment_dir = os.path.join(args.output_dir, first_name)

        if not os.path.exists(alignment_dir):
            try:
                os.makedirs(alignment_dir)
            except Exception as e:
                logging.warning(f"Failed to create directory for {first_name} with exception {e}...")
                continue

        fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
        with os.fdopen(fd, 'w') as fp:
            fp.write(f'>query\n{seq}')

        try:
            alignment_runner.run(
                fasta_path, alignment_dir
            )
        except:
            traceback.print_exc()
            logging.warning(f"Failed to run alignments for {first_name}. Skipping...")
            os.remove(fasta_path)
            os.rmdir(alignment_dir)
            continue

        os.remove(fasta_path)

        for name in names[1:]:
            if(name in dirs):
                logging.warning(
                    f'{name} has already been processed. Skipping...'
                )
                continue

            cp_dir = os.path.join(args.output_dir, name)
            os.makedirs(cp_dir, exist_ok=True)

            for f in os.listdir(alignment_dir):
                copyfile(os.path.join(alignment_dir, f), os.path.join(cp_dir, f))


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
    items = sorted(items, key=lambda x: x[0])

    return [x[0] for x in items], [x[1] for x in items]


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
    )

    input_file = args.input_file
    with open(input_file, 'r') as fp:
        fasta_str = fp.read()
    input_seqs, input_chains = parse_fasta(fasta_str)
    input_seqs, input_chains = make_uniq_seq_groups(input_seqs, input_chains)

    if "OMPI_COMM_WORLD_RANK" in os.environ:
        # ABCI (OpenMPI)
        mpi_rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
        mpi_size = int(os.environ["OMPI_COMM_WORLD_SIZE"])

    else:
        # Fugaku (Fujitsu MPI)
        mpi_rank = int(os.environ["PMIX_RANK"])
        mpi_size = int(os.environ["OMPI_UNIVERSE_SIZE"])

    assert mpi_size > 0
    assert mpi_rank >= 0 and mpi_rank < mpi_size
    input_seqs   =   input_seqs[mpi_rank::mpi_size]
    input_chains = input_chains[mpi_rank::mpi_size]

    run_seq_group_alignments(
        zip(input_seqs, input_chains),
        alignment_runner,
        args)


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

    args = parser.parse_args()

    main(args)
