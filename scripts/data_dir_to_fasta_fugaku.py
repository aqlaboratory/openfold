import argparse
import logging
import os
import threading

from openfold.data import mmcif_parsing
from openfold.np import protein, residue_constants


def parse_all(data_dir, fnames, raise_errors, fasta):
    for fname in fnames:
        basename, ext = os.path.splitext(fname)
        basename = basename.upper()
        fpath = os.path.join(data_dir, fname)
        if(ext == ".cif"):
            with open(fpath, 'r') as fp:
                mmcif_str = fp.read()

            mmcif = mmcif_parsing.parse(
                file_id=basename, mmcif_string=mmcif_str
            )
            if(mmcif.mmcif_object is None):
                logging.warning(f'Failed to parse {fname}...')
                if(raise_errors):
                    raise list(mmcif.errors.values())[0]
                else:
                    continue

            mmcif = mmcif.mmcif_object
            for chain, seq in mmcif.chain_to_seqres.items():
                chain_id = '_'.join([basename, chain])
                fasta.append(f">{chain_id}")
                fasta.append(seq)
        elif(ext == ".core"):
            with open(fpath, 'r') as fp:
                core_str = fp.read()

            core_protein = protein.from_proteinnet_string(core_str)
            aatype = core_protein.aatype
            seq = ''.join([
                residue_constants.restypes_with_x[aatype[i]]
                for i in range(len(aatype))
            ])
            fasta.append(f">{basename}")
            fasta.append(seq)

    return fasta


def main(args):
    fnames = os.listdir(args.data_dir)
    # fasta = parse_all(args.data_dir, fnames, args.raise_errors)

    threads = []
    ret_fasta = [None]*args.no_tasks
    for i in range(args.no_tasks):
        ret_fasta[i] = []
        print(f"Started thread {i}...")
        targs = [
            args.data_dir,
            fnames[i::args.no_tasks],
            args.raise_errors,
            ret_fasta[i],
        ]
        t = threading.Thread(target=parse_all, args=targs)
        threads.append(t)
        t.start()

    for t in threads:
        t.join()

    fasta = []
    for ret in ret_fasta:
        fasta.extend(ret)

    with open(args.output_path, "w") as fp:
        fp.write('\n'.join(fasta))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_dir", type=str,
        help="Path to a directory containing mmCIF or .core files"
    )
    parser.add_argument(
        "output_path", type=str,
        help="Path to output FASTA file"
    )
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )
    parser.add_argument(
        "--no_tasks", type=int, default=1,
    )

    args = parser.parse_args()
    if args.no_tasks < 1:
        raise ValueError("Use data_dir_to_fasta.py if you don't use multiple threads.")

    main(args)
