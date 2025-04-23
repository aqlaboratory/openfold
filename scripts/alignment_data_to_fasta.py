"""
This script generates a FASTA file for all chains in an alignment directory or
alignment DB.
"""

import json
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from tqdm import tqdm


def chain_dir_to_fasta(dir: Path) -> str:
    """
    Generates a FASTA string from a chain directory.
    """
    # take some alignment file
    for alignment_file_type in [
        "mgnify_hits.a3m",
        "uniref90_hits.a3m",
        "bfd_uniclust_hits.a3m",
    ]:
        alignment_file = dir / alignment_file_type
        if alignment_file.exists():
            break

    with open(alignment_file, "r") as f:
        next(f)  # skip the first line
        seq = next(f).strip()

        try:
            next_line = next(f)
        except StopIteration:
            pass
        else:
            assert next_line.startswith(">")  # ensure that sequence ended

    chain_id = dir.name

    return f">{chain_id}\n{seq}\n"


def index_entry_to_fasta(index_entry: dict, db_dir: Path, chain_id: str) -> str:
    """
    Generates a FASTA string from an alignment-db index entry.
    """
    db_file = db_dir / index_entry["db"]

    # look for an alignment file
    for alignment_file_type in [
        "mgnify_hits.a3m",
        "uniref90_hits.a3m",
        "bfd_uniclust_hits.a3m",
    ]:
        for file_info in index_entry["files"]:
            if file_info[0] == alignment_file_type:
                start, size = file_info[1], file_info[2]
                break

    with open(db_file, "rb") as f:
        f.seek(start)
        msa_lines = f.read(size).decode("utf-8").splitlines()
        seq = msa_lines[1]

        try:
            next_line = msa_lines[2]
        except IndexError:
            pass
        else:
            assert next_line.startswith(">")  # ensure that sequence ended

    return f">{chain_id}\n{seq}\n"


def main(
    output_path: Path, alignment_db_index: Optional[Path], alignment_dir: Optional[Path]
) -> None:
    """
    Generate a FASTA file from either an alignment-db index or a chain directory using multi-threading.
    """
    fasta = []

    if alignment_dir and alignment_db_index:
        raise ValueError(
            "Only one of alignment_db_index and alignment_dir can be provided."
        )

    if alignment_dir:
        print("Creating FASTA from alignment directory...")
        chain_dirs = list(alignment_dir.iterdir())

        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(chain_dir_to_fasta, chain_dir)
                for chain_dir in chain_dirs
            ]
            for future in tqdm(as_completed(futures), total=len(chain_dirs)):
                fasta.append(future.result())

    elif alignment_db_index:
        print("Creating FASTA from alignment dbs...")

        with open(alignment_db_index, "r") as f:
            index = json.load(f)

        db_dir = alignment_db_index.parent

        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(index_entry_to_fasta, index_entry, db_dir, chain_id)
                for chain_id, index_entry in index.items()
            ]
            for future in tqdm(as_completed(futures), total=len(index)):
                fasta.append(future.result())
    else:
        raise ValueError("Either alignment_db_index or alignment_dir must be provided.")

    with open(output_path, "w") as f:
        f.write("".join(fasta))
    print(f"FASTA file written to {output_path}.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_path",
        type=Path,
        help="Path to output FASTA file.",
    )
    parser.add_argument(
        "--alignment_db_index",
        type=Path,
        help="Path to alignment-db index file.",
    )
    parser.add_argument(
        "--alignment_dir",
        type=Path,
        help="Path to alignment directory.",
    )

    args = parser.parse_args()
    main(args.output_path, args.alignment_db_index, args.alignment_dir)
