"""
This is a modified version of the create_alignment_db.py script in OpenFold
which supports sharding into multiple files. The created index is already a
super index, meaning that "unify_alignment_db_indices.py" does not need to be
run on the output index. Additionally this script uses threading and
multiprocessing and is much faster than the old version.
"""

import argparse
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from math import ceil
from multiprocessing import cpu_count
from pathlib import Path

from tqdm import tqdm


def split_file_list(file_list: list[Path], n_shards: int):
    """
    Split up the total file list into n_shards sublists.
    """
    split_list = []

    for i in range(n_shards):
        split_list.append(file_list[i::n_shards])

    assert len([f for sublist in split_list for f in sublist]) == len(file_list)

    return split_list


def chunked_iterator(lst: list, chunk_size: int):
    """Iterate over a list in chunks of size chunk_size."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]


def read_chain_dir(chain_dir: Path) -> dict:
    """
    Read all alignment files in a single chain directory and return a dict
    mapping chain name to file names and bytes.
    """
    if not chain_dir.is_dir():
        raise ValueError(f"chain_dir must be a directory, but is {chain_dir}")

    # ensure that PDB IDs are all lowercase
    pdb_id, chain = chain_dir.name.split("_")
    pdb_id = pdb_id.lower()
    chain_name = f"{pdb_id}_{chain}"

    file_data = []

    for file_path in sorted(chain_dir.iterdir()):
        file_name = file_path.name

        with open(file_path, "rb") as file:
            file_bytes = file.read()

        file_data.append((file_name, file_bytes))

    return {chain_name: file_data}


def process_chunk(chain_files: list[Path]) -> dict:
    """
    Returns the file names and bytes for all chains in a chunk of files.
    """
    chunk_data = {}

    with ThreadPoolExecutor() as executor:
        for file_data in executor.map(read_chain_dir, chain_files):
            chunk_data.update(file_data)

    return chunk_data


def create_index_default_dict() -> dict:
    """
    Returns a default dict for the index entries).
    """
    return {"db": None, "files": []}


def create_shard(
    shard_files: list[Path], output_dir: Path, output_name: str, shard_num: int
) -> dict:
    """
    Creates a single shard of the alignment database, and returns the
    corresponding indices for the super index.
    """
    CHUNK_SIZE = 200
    shard_index = defaultdict(
        create_index_default_dict
    )  # e.g. {chain_name: {db: str, files: [(file_name, db_offset, file_length)]}, ...}
    chunk_iter = chunked_iterator(shard_files, CHUNK_SIZE)

    pbar_desc = f"Shard {shard_num}"
    output_path = output_dir / f"{output_name}_{shard_num}.db"

    db_offset = 0
    db_file = open(output_path, "wb")
    for files_chunk in tqdm(
        chunk_iter,
        total=ceil(len(shard_files) / CHUNK_SIZE),
        desc=pbar_desc,
        position=shard_num,
        leave=False,
    ):
        # get processed files for one chunk
        chunk_data = process_chunk(files_chunk)

        # write to db and store info in index
        for chain_name, file_data in chunk_data.items():
            shard_index[chain_name]["db"] = output_path.name

            for file_name, file_bytes in file_data:
                file_length = len(file_bytes)
                shard_index[chain_name]["files"].append(
                    (file_name, db_offset, file_length)
                )
                db_file.write(file_bytes)
                db_offset += file_length
    db_file.close()

    return shard_index


def main(args):
    alignment_dir = args.alignment_dir
    output_dir = args.output_db_path
    output_dir.mkdir(exist_ok=True, parents=True)
    output_db_name = args.output_db_name
    n_shards = args.n_shards

    n_cpus = cpu_count()
    if n_shards > n_cpus:
        print(
            f"Warning: Your number of shards ({n_shards}) is greater than the number of cores on your machine ({n_cpus}). "
            "This may result in slower performance. Consider using a smaller number of shards."
        )

    # get all chain dirs in alignment_dir
    print("Getting chain directories...")
    all_chain_dirs = sorted([f for f in tqdm(alignment_dir.iterdir())])

    # split chain dirs into n_shards sublists
    chain_dir_shards = split_file_list(all_chain_dirs, n_shards)

    # total index for all shards
    super_index = {}

    # create a shard for each sublist
    print(f"Creating {n_shards} alignment-db files...")
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                create_shard, shard_files, output_dir, output_db_name, shard_index
            )
            for shard_index, shard_files in enumerate(chain_dir_shards)
        ]

        for future in as_completed(futures):
            shard_index = future.result()
            super_index.update(shard_index)
    print("\nCreated all shards.")

    if args.duplicate_chains_file:
        print("Extending super index with duplicate chains...")
        duplicates_added = 0
        with open(args.duplicate_chains_file, "r") as fp:
            duplicate_chains = [line.strip().split() for line in fp]

        for chains in duplicate_chains:
            # find representative with alignment
            for chain in chains:
                if chain in super_index:
                    representative_chain = chain
                    break
            else:
                print(f"No representative chain found for {chains}, skipping...")
                continue

            # add duplicates to index
            for chain in chains:
                if chain != representative_chain:
                    super_index[chain] = super_index[representative_chain]
                    duplicates_added += 1

        print(f"Added {duplicates_added} duplicate chains to index.")

    # write super index to file
    print("\nWriting super index...")
    index_path = output_dir / f"{output_db_name}.index"
    with open(index_path, "w") as fp:
        json.dump(super_index, fp, indent=4)

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This script creates an alignment database format from a directory of
        precomputed alignments. For better file system health, the total
        database is split into n_shards files, where each shard contains a
        subset of the total alignments. The output is a directory containing the
        n_shards database files, and a single index file mapping chain names to
        the database file and byte offsets for each alignment file.
        
        Note: For optimal performance, your machine should have at least as many
        cores as shards you want to create.
        """
    )
    parser.add_argument(
        "alignment_dir",
        type=Path,
        help="""Path to precomputed flattened alignment directory, with one
                subdirectory per chain.""",
    )
    parser.add_argument("output_db_path", type=Path)
    parser.add_argument("output_db_name", type=str)
    parser.add_argument(
        "--n_shards",
        type=int,
        help="Number of shards to split the database into",
        default=10,
    )
    parser.add_argument(
        "--duplicate_chains_file",
        type=Path,
        help="""
        Optional path to file containing duplicate chain information, where each
        line contains chains that are 100% sequence identical. If provided,
        duplicate chains will be added to the index and point to the same
        underlying database entry as their representatives in the alignment dir.
        """,
        default=None,
    )

    args = parser.parse_args()

    main(args)
