import argparse
import logging
import os
import shutil


def main(args):
    count = 0
    max_count = args.max_count if args.max_count is not None else -1
    msas = sorted(f for f in os.listdir(args.msa_dir))
    mmcifs = sorted(f for f in os.listdir(args.mmcif_dir))
    mmcif_idx = 0
    for f in msas:
        if(count == max_count):
            break

        path = os.path.join(args.msa_dir, f)
        name = os.path.splitext(f)[0]
        spl = name.upper().split('_')
        if(len(spl) != 3):
            continue
         
        pdb_id, _, chain_id = spl
       
        while pdb_id > os.path.splitext(mmcifs[mmcif_idx])[0].upper():
            mmcif_idx += 1

        # Only consider files with matching mmCIF files
        if(pdb_id == os.path.splitext(mmcifs[mmcif_idx])[0].upper()):
            dirname = os.path.join(args.out_dir, '_'.join([pdb_id, chain_id]))
            os.makedirs(dirname, exist_ok=True)
            dest = os.path.join(dirname, f)
            if(args.copy):
                shutil.copyfile(path, dest)
            else:
                os.rename(path, dest)

            count += 1
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
        "Converts raw ProteinNet MSAs into a format recognized by the parser"
    )
    parser.add_argument(
        "msa_dir", type=str, help="Directory containing ProteinNet MSAs"
    )
    parser.add_argument(
        "mmcif_dir", type=str, help="Directory containing PDB mmCIFs"
    )
    parser.add_argument(
        "out_dir", type=str,
        help="Directory to which output should be saved"
    )
    parser.add_argument(
        "--copy", type=bool, default=True,
        help="Whether to copy the MSAs to out_dir rather than moving them"
    )
    parser.add_argument(
        "--max_count", type=int, default=None,
        help="A bound on the number of MSAs to process"
    )

    args = parser.parse_args()

    main(args)
