import argparse
import logging
import os
import tempfile

import openfold.data.mmcif_parsing as mmcif_parsing
from openfold.data.data_pipeline import AlignmentRunner
from openfold.np import protein, residue_constants

from utils import add_data_args

#python3 scripts/precompute_alignments.py mmcif_dir/ alignment_dir/     data/uniref90/uniref90.fasta     data/mgnify/mgy_clusters_2018_12.fa     data/pdb70/pdb70     data/pdb_mmcif/mmcif_files/     data/uniclust30/uniclust30_2018_08/uniclust30_2018_08     --bfd_database_path data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt     --cpus 16  --jackhmmer_binary_path /home/u00u98too4mkqFBu8M357/openfold/lib/conda/envs/openfold_venv/bin/jackhmmer  --hhblits_binary_path /home/u00u98too4mkqFBu8M357/openfold/lib/conda/envs/openfold_venv/bin/hhblits  --hhsearch_binary_path /home/u00u98too4mkqFBu8M357/openfold/lib/conda/envs/openfold_venv/bin/hhsearch  --kalign_binary_path /home/u00u98too4mkqFBu8M357/openfold/lib/conda/envs/openfold_venv/bin/kalign

logging.basicConfig(level=logging.DEBUG)


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
        file_id = os.path.splitext(f)[0]
        seqs = {}
        if(f.endswith('.cif')):
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
        elif(f.endswith('.fasta')):
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
        elif(f.endswith('.core')):
            with open(path, 'r') as fp:
                core_str = fp.read()
            core_prot = protein.from_proteinnet_string(core_str)
            aatype = core_prot.aatype
            seq = ''.join([
                residue_constants.restypes_with_x[aatype[i]] 
                for i in range(len(aatype))
            ])
            seqs[file_id] = seq
        else:
            continue

        for name, seq in seqs.items():
            alignment_dir = os.path.join(args.output_dir, name)
            if(os.path.isdir(alignment_dir)):
                logging.info(f'{f} has already been processed. Skipping...')
                continue

            os.makedirs(alignment_dir)

            fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
            with os.fdopen(fd, 'w') as fp:
                fp.write(f'>query\n{seq}')

            alignment_runner.run(
                fasta_path, alignment_dir
            )

            os.remove(fasta_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir", type=str,
        help="""Path to directory containing mmCIF, FASTA and/or ProteinNet
                .core files"""
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
