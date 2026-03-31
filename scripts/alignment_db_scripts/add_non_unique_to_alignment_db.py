from argparse import ArgumentParser
from pathlib import Path
import json


def main(args):
    # get the super index
    with open(args.alignment_db_super_index_path, "r") as fp:
        super_index = json.load(fp)

    # get all chains and sequences
    chains_to_seqs = {}
    with open(args.all_chains_fasta, "r") as fp:
        lines = fp.readlines()

        # iterate through chain-sequence pairs
        for chain_idx in range(0, len(lines), 2):
            chain = lines[chain_idx][1:].strip()
            seq = lines[chain_idx + 1].strip()
            chains_to_seqs[chain] = seq

    chains_w_alignments = set(super_index.keys())
    chains_wo_alignments = set(chains_to_seqs.keys()) - chains_w_alignments

    seq_to_chain_w_alignment = {
        chains_to_seqs[chain]: chain for chain in chains_w_alignments
    }

    print("Unique sequences with alignments:", len(seq_to_chain_w_alignment))

    # map chain without alignment to alignment entry of another chain with the
    # same sequence
    remaining_unaligned_chains = []
    for chain in chains_wo_alignments:
        seq = chains_to_seqs[chain]

        try:
            corresponding_alignment = super_index[seq_to_chain_w_alignment[seq]]
        # no corresponding chain with alignment found
        except KeyError:
            remaining_unaligned_chains.append(chain)
            continue

        super_index[chain] = corresponding_alignment

    with open(args.output_path, "w") as fp:
        json.dump(super_index, fp)

    print(
        f"No corresponding alignment found for the following {len(remaining_unaligned_chains)} chains:",
        remaining_unaligned_chains,
    )


if __name__ == "__main__":
    parser = ArgumentParser(
        description="""
        If the alignment-db index was created on unique-chain alignments only,
        this will add the missing chain entries to the super-index file based on
        a .fasta file that contains sequences for all chains.
        
        Note that this only modifies the index and not the database itself, as
        the duplicate sequences will just point to the same alignments. 
        """
    )
    parser.add_argument(
        "alignment_db_super_index_path",
        type=Path,
        help="Path to alignment-db super index file.",
    )
    parser.add_argument(
        "output_path", type=Path, help="Write the output super index to this path."
    )
    parser.add_argument(
        "all_chains_fasta",
        type=Path,
        help="Path to the fasta file containing sequences for all chains.",
    )

    args = parser.parse_args()

    main(args)
