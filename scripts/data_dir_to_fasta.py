import argparse
import logging
import os
import string
from collections import defaultdict
from openfold.data import mmcif_parsing
from openfold.np import protein, residue_constants


def main(args):
    fasta = []
    for fname in os.listdir(args.data_dir):
        basename, ext = os.path.splitext(fname)
        basename = basename.upper()
        fpath = os.path.join(args.data_dir, fname)
        if(ext == ".cif"):
            with open(fpath, 'r') as fp:
                mmcif_str = fp.read()
            
            mmcif = mmcif_parsing.parse(
                file_id=basename, mmcif_string=mmcif_str
            )
            if(mmcif.mmcif_object is None):
                logging.warning(f'Failed to parse {fname}...')
                if(args.raise_errors):
                    raise Exception(list(mmcif.errors.values())[0])
                else:
                    continue

            mmcif = mmcif.mmcif_object
            for chain, seq in mmcif.chain_to_seqres.items():
                chain_id = '_'.join([basename, chain])
                fasta.append(f">{chain_id}")
                fasta.append(seq)
        elif(ext == ".pdb"):
            with open(fpath, 'r') as fp:
                pdb_str = fp.read()
            
            protein_object = protein.from_pdb_string(pdb_str)
            aatype = protein_object.aatype
            chain_index = protein_object.chain_index

            last_chain_index = chain_index[0]
            chain_dict = defaultdict(list)
            for i in range(aatype.shape[0]):
                chain_dict[chain_index[i]].append(residue_constants.restypes_with_x[aatype[i]])
            
            chain_tags = string.ascii_uppercase
            for chain, seq in chain_dict.items():
                chain_id = '_'.join([basename, chain_tags[chain]])
                fasta.append(f">{chain_id}")
                fasta.append(''.join(seq))

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

    args = parser.parse_args()

    main(args)
