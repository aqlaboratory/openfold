"""
This script takes a .fasta file as input and then clusters it on a given
sequence identity threshold using mmseqs2. The mmseqs2 flags are identical to
what PDB officially uses to provide their official sequence clusters
(https://github.com/soedinglab/MMseqs2/issues/452).
"""
import shutil
import subprocess
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path


def reformat_cluster_file(cluster_file: Path, output_file: Path):
    """
    This function takes a mmseqs2 output cluster file and reformats it to a text
    file where each line contains a space-separated list of {PDB_ID}_{CHAIN_ID}
    belonging to the same cluster.
    """
    cluster_to_chains = defaultdict(list)

    # extract all chains belonging to each cluster
    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()
            cluster_name, chain_id = line.split()
            cluster_to_chains[cluster_name].append(chain_id)

    # write all chains belonging to the same cluster on the same line
    with open(output_file, "w") as f:
        for chains in cluster_to_chains.values():
            f.write(f"{' '.join(chains)}\n")


def main(args):
    input_file = args.input_fasta.absolute()
    output_file = args.output_file.absolute()
    output_dir = args.output_file.parent

    # prefix that all output files get
    mmseqs_prefix = "_mmseqs_out"

    # temporary directory that mmseqs2 uses
    tmp_name = f"{mmseqs_prefix}_temp"
    tmp_dir = output_dir / tmp_name

    mmseqs_command = [
        args.mmseqs_binary_path,
        "easy-cluster",
        input_file,
        mmseqs_prefix,
        tmp_name,
        "--min-seq-id",
        str(args.seq_id),
        "-c",
        "0.9",
        "-s",
        "8",
        "--max-seqs",
        "1000",
        "--cluster-mode",
        "1",
    ]

    # run mmseqs with PDB settings
    print("Running mmseqs2...")
    subprocess.run(mmseqs_command, check=True, cwd=output_dir)

    cluster_file = output_dir / "_mmseqs_out_cluster.tsv"

    print("Reformatting output file...")
    reformat_cluster_file(cluster_file, output_file)

    print("Cleaning up mmseqs2 output...")
    mmseqs_outputs = [
        output_dir / f"{mmseqs_prefix}_{suffix}"
        for suffix in ["cluster.tsv", "rep_seq.fasta", "all_seqs.fasta"]
    ]
    for file in mmseqs_outputs:
        file.unlink()
    shutil.rmtree(tmp_dir)

    print("Done!")


if __name__ == "__main__":
    parser = ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        "input_fasta",
        type=Path,
        help="Input .fasta file. Sequence names should be in format >{PDB_ID}_{CHAIN_ID}",
    )
    parser.add_argument(
        "output_file",
        type=Path,
        help="Output file. Each line will contain a space-separated list of {PDB_ID}_{CHAIN_ID} belonging to the same cluster.",
    )
    parser.add_argument("mmseqs_binary_path", type=str, help="Path to mmseqs binary")
    parser.add_argument("--seq-id", type=float, default=0.4, help="Sequence identity threshold for clustering.")

    args = parser.parse_args()

    main(args)
