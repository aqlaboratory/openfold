import argparse
import logging
import os
import tempfile

import openfold.features.mmcif_parsing as mmcif_parsing
from openfold.features.data_pipeline import AlignmentRunner
from scripts.utils import add_data_args


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
        small_bfd_database_path=args.small_bfd_database_path,
        pdb70_database_path=args.pdb70_database_path,
        use_small_bfd=args.bfd_database_path is None,
        no_cpus=args.cpus,
    )

    for f in os.listdir(args.input_dir):
        path = os.path.join(args.input_dir, f)
        is_mmcif = f.endswith('.cif')
        is_fasta = f.endswith('.fasta')
        file_id = os.path.splitext(f)[0]
        seqs = {}
        if(is_mmcif):
            with open(path, 'r') as fp:
                mmcif_str = fp.read()
            mmcif = mmcif_parsing.parse(
                file_id=file_id, mmcif_string=mmcif_str
            )
            if(mmcif.mmcif_object is None):
                logging.warning(f'Failed to parse {f}...')
                if(args.raise_errors):
                    raise list(mmcif.errors.values())[0]
                else:
                    continue
            mmcif = mmcif.mmcif_object
            for k,v in mmcif.chain_to_seqres.items():
                chain_id = '_'.join([file_id, k])
                seqs[chain_id] = v
        elif(is_fasta):
            with open(path, 'r') as fp:
                fasta_str = fp.read()
            input_seqs, _ = parsers.parse_fasta(fasta_str)
            if len(input_seqs) != 1: 
                msg = f'More than one input_sequence found in {f}'
                if(args.raise_errors):
                    raise ValueError(msg)
                else:
                    logging.warning(msg)
            input_sequence = input_seqs[0]
            seqs[file_id] = input_sequence
        else:
            continue

        for name, seq in seqs.items():
            alignment_dir = os.path.join(args.output_dir, name)
            if(os.path.isdir(alignment_dir)):
                logging.info(f'{f} has already been processed. Skipping...')
                continue

            os.makedirs(alignment_dir)

            if(not is_fasta):
                fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
                with os.fdopen(fd, 'w') as fp:
                    fp.write(f'>query\n{seq}')

            alignment_runner.run(
                f if is_fasta else fasta_path, alignment_dir
            )

            if(not is_fasta):
                os.remove(fasta_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir", type=str,
        help="Path to directory containing mmCIF and/or FASTA files"
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
        "--cpus", type=int, default=4,
        help="Number of CPUs to use"
    )

    args = parser.parse_args()

    main(args)
