import argparse
import os
import pickle

from alphafold.data import pipeline, templates

from scripts.utils import add_data_args


def main(args):
    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=args.mmcif_dir,
        max_template_date=args.max_template_date,
        max_hits=20,
        kalign_binary_path=args.kalign_binary_path,
        release_dates_path=None,
        obsolete_pdbs_path=args.obsolete_pdbs_path,
    )

    data_pipeline = pipeline.DataPipeline(
        jackhmmer_binary_path=args.jackhmmer_binary_path,
        hhblits_binary_path=args.hhblits_binary_path,
        hhsearch_binary_path=args.hhsearch_binary_path,
        uniref90_database_path=args.uniref90_database_path,
        mgnify_database_path=args.mgnify_database_path,
        bfd_database_path=args.bfd_database_path,
        uniclust30_database_path=args.uniclust30_database_path,
        pdb70_database_path=args.pdb70_database_path,
        small_bfd_database_path=None,
        template_featurizer=template_featurizer,
        use_small_bfd=False,
    )

    feature_dict = data_pipeline.process(
        input_fasta_path=args.fasta_path,
        msa_output_dir=args.output_dir,
    )

    with open(os.path.join(args.output_dir, "feature_dict.pickle"), "wb") as fp:
        pickle.dump(feature_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_path", type=str)
    parser.add_argument("mmcif_dir", type=str)
    parser.add_argument("output_dir", type=str)
    add_data_args(parser)

    args = parser.parse_args()

    main(args)
