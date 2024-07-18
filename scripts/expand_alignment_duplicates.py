"""
The OpenProteinSet alignment database is non-redundant, meaning that it only
stores one explicit representative alignment directory for all PDB chains in a
100% sequence identity cluster. In order to add explicit alignments for all PDB
chains, this script will add the missing chain directories and symlink them to
their representative alignment directories. This is required in order to train
OpenFold on the full PDB, not just one representative chain per cluster.
"""

from argparse import ArgumentParser
from pathlib import Path

from tqdm import tqdm


def create_duplicate_dirs(duplicate_chains: list[list[str]], alignment_dir: Path):
    """
    Create duplicate directory symlinks for all chains in the given duplicate lists.

    Args:
        duplicate_lists (list[list[str]]): A list of lists, where each inner list
            contains chains that are 100% sequence identical.
        alignment_dir (Path): Path to flattened alignment directory, with one
            subdirectory per chain.
    """
    print("Creating duplicate directory symlinks...")
    dirs_created = 0
    for chains in tqdm(duplicate_chains):
        # find the chain that has an alignment
        for chain in chains:
            if (alignment_dir / chain).exists():
                representative_chain = chain
                break
        else:
            print(f"No representative chain found for {chains}, skipping...")
            continue

        # create symlinks for all other chains
        for chain in chains:
            if chain != representative_chain:
                target_path = alignment_dir / chain
                if target_path.exists():
                    print(f"Chain {chain} already exists, skipping...")
                else:
                    (target_path).symlink_to(alignment_dir / representative_chain)
                    dirs_created += 1

    print(f"Created directories for {dirs_created} duplicate chains.")


def main(alignment_dir: Path, duplicate_chains_file: Path):
    # read duplicate chains file
    with open(duplicate_chains_file, "r") as fp:
        duplicate_chains = [list(line.strip().split()) for line in fp]

    # convert to absolute path for symlink creation
    alignment_dir = alignment_dir.resolve()

    create_duplicate_dirs(duplicate_chains, alignment_dir)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "alignment_dir",
        type=Path,
        help="""Path to flattened alignment directory, with one subdirectory 
                per chain.""",
    )
    parser.add_argument(
        "duplicate_chains_file",
        type=Path,
        help="""Path to file containing duplicate chains, where each line
                contains a space-separated list of chains that are 100%%
                sequence identical.
                """,
    )
    args = parser.parse_args()
    main(args.alignment_dir, args.duplicate_chains_file)
