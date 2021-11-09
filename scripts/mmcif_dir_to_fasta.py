import argparse
import logging
import os

from openfold.data import mmcif_parsing


def main(args):
    fasta = []
    for fname in os.listdir(args.mmcif_dir):
        basename, ext = os.path.splitext(fname)
        basename = basename.upper()

        if(not ext == ".cif"):
            continue

        fpath = os.path.join(args.mmcif_dir, fname)
        with open(fpath, 'r') as fp:
            mmcif_str = fp.read()
        
        mmcif = mmcif_parsing.parse(
            file_id=basename, mmcif_string=mmcif_str
        )
        if(mmcif.mmcif_object is None):
            logging.warning(f'Failed to parse {fname}...')
            if(args.raise_errors):
                raise list(mmcif.errors.values())[0]
            else:
                continue

        mmcif = mmcif.mmcif_object
        for chain, seq in mmcif.chain_to_seqres.items():
            chain_id = '_'.join([basename, chain])
            fasta.append(f">{chain_id}")
            fasta.append(seq)

    with open(args.output_path, "w") as fp:
        fp.write('\n'.join(fasta))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mmcif_dir", type=str,
        help="Path to a directory containing mmCIF files"
    )
    parser.add_argument(
        "output_path", type=str,
        help="Path to output FASTA file"
    )
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )

    args = parser.parse_args()

    main(args)
