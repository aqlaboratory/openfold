import argparse
import copy
import errno
import logging
import multiprocessing
import os
import pathlib
import pickle
import shutil
from datetime import date
from multiprocessing import Pool
from typing import Optional

from openfold.data import templates
from openfold.data.data_pipeline import DataPipeline
from openfold.data.tools import hhsearch


class TemplatePipelineRunner:
    """Uses the MSA alignments to look up matching templates using
    hhsearch and creates template features using MMCIFs."""

    def __init__(
        self,
        hhsearch_binary_path: str,
        pdb70_database_path: str,
        template_mmcif_dir: str,
        max_template_date: str,
        max_templates: int,
        kalign_binary_path: str,
        release_dates_path: Optional[str],
        obsolete_pdbs_path: str,
        no_cpus: int,
    ):
        self.hhsearch_pdb70_runner = hhsearch.HHSearch(
            binary_path=hhsearch_binary_path,
            databases=[pdb70_database_path],
            n_cpu=no_cpus,
        )

        self.template_featurizer = templates.TemplateHitFeaturizer(
            mmcif_dir=template_mmcif_dir,
            max_template_date=max_template_date,
            max_hits=max_templates,
            kalign_binary_path=kalign_binary_path,
            release_dates_path=release_dates_path,
            obsolete_pdbs_path=obsolete_pdbs_path
        )

        self.data_pipeline = DataPipeline(
            template_featurizer=self.template_featurizer,
        )

    def run(
        self,
        a3m_dir: str,
        fasta_path: str,
    ):
        fasta_name = pathlib.Path(fasta_path).stem
        uniref90_msa_as_a3m = None
        for f in os.listdir(a3m_dir):
            # is_uniref90_a3m = f == "uniref90_hits.a3m"
            is_uniref90_a3m = f.endswith('.a3m')
            file_path = os.path.join(a3m_dir, f)
            if is_uniref90_a3m:
                with open(file_path, "r") as a3m:
                    uniref90_msa_as_a3m = a3m.read()
        if uniref90_msa_as_a3m is None:
            logging.warning(f'Failed to find uniref90_hits.a3m for fasta {fasta_name}')
            if args.raise_errors:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                        os.path.join(a3m_dir, "uniref90_hits.a3m"))
            else:
                return

        hhsearch_result = self.hhsearch_pdb70_runner.query(uniref90_msa_as_a3m)
        pdb70_out_path = os.path.join(a3m_dir, "pdb70_hits.hhr")
        with open(pdb70_out_path, "w") as f:
            f.write(hhsearch_result)

        feature_dict = self.data_pipeline.process_fasta(
            fasta_path=fasta_path,
            alignment_dir=a3m_dir,
        )

        return feature_dict





def main(args, template_pipeline_runner):
    output_dir_base = args.msa_dir
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)

    # Build the template pipeline
    # template_pipeline_runner = TemplatePipelineRunner(
    #     hhsearch_binary_path=args.hhsearch_binary_path,
    #     pdb70_database_path=args.pdb70_database_path,
    #     template_mmcif_dir=args.template_mmcif_dir,
    #     max_template_date=args.max_template_date,
    #     max_templates=args.max_templates,
    #     kalign_binary_path=args.kalign_binary_path,
    #     release_dates_path=args.release_dates_path,
    #     obsolete_pdbs_path=args.obsolete_pdbs_path,
    #     no_cpus=args.cpus,
    # )

    for f in os.listdir(args.fastas_dir):
        print(f'Processing {f}')
        is_fasta = f.endswith('.fa') or f.endswith('.fasta')
        if not is_fasta:
            continue
        fasta_file_path = os.path.join(args.fastas_dir, f)
        fasta_name = pathlib.Path(fasta_file_path).stem
        a3m_dir = os.path.join(args.msa_dir, fasta_name, "alignments")

        feature_dict = template_pipeline_runner.run(a3m_dir, fasta_file_path)
        if feature_dict is None:
            continue

        output_dir = os.path.join(output_dir_base, fasta_name)
        features_output_path = os.path.join(output_dir, 'features.pkl')
        with open(features_output_path, 'wb') as feat_path:
            pickle.dump(feature_dict, feat_path, protocol=4)

        print(f'Completed processing {f}')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "msa_dir", type=str,
        help="Path to the base directory containing the alignments in output-specific "
             "subdirectories, as outputted from the MSA alignment pipeline.",
    )
    parser.add_argument(
        "fastas_dir", type=str,
        help="Directory containing the FASTA files corresponding to the alignment files."
    )
    parser.add_argument(
        '--pdb70_database_path', type=str,
    )
    parser.add_argument(
        '--template_mmcif_dir', type=str,
    )
    parser.add_argument(
        '--hhsearch_binary_path', type=str, default='/usr/bin/hhsearch'
    )
    parser.add_argument(
        '--kalign_binary_path', type=str, default='/usr/bin/kalign'
    )
    parser.add_argument(
        '--max_template_date', type=str,
        default=date.today().strftime("%Y-%m-%d"),
    )
    parser.add_argument(
        '--max_templates', type=int, default=4
    )
    parser.add_argument(
        '--release_dates_path', type=str,
        default=None
    )
    parser.add_argument(
        '--obsolete_pdbs_path', type=str, default=None
    )
    parser.add_argument(
        '--raise_errors', type=bool, default=False,
        help="Whether to crash on errors"
    )
    parser.add_argument(
        '--cpus', type=int, default=2,
        help="Number of CPUs to use"
    )
    parser.add_argument(
        '--n_threads', type=int, default=1,
        help="Number of threads to use per node"
    )
    parser.add_argument(
        '--local_dir', type=str
    )

    args = parser.parse_args()

    template_pipeline_runner = TemplatePipelineRunner(
        hhsearch_binary_path=args.hhsearch_binary_path,
        pdb70_database_path=args.pdb70_database_path,
        template_mmcif_dir=args.template_mmcif_dir,
        max_template_date=args.max_template_date,
        max_templates=args.max_templates,
        kalign_binary_path=args.kalign_binary_path,
        release_dates_path=args.release_dates_path,
        obsolete_pdbs_path=args.obsolete_pdbs_path,
        no_cpus=args.cpus,
    )

    def walk_dirs(directory, batch_size):
        walk_dirs_generator = os.walk(directory)
        for dirname, subdirectories, filenames in walk_dirs_generator:
            for i in range(0, len(filenames), batch_size):
                # slice the subdirectories list
                yield [os.path.join(dirname, filename) for filename in filenames[i:i+batch_size] ]

    local_dir = args.local_dir
    # print('local_dir_temp', args.local_dir)
    fastas_dir = args.fastas_dir
    succesfully_processed = []
    already_processed = []
    if os.path.exists('succesfull.txt'):
        with open('succesfull.txt', 'r') as success_file:
            successes = success_file.readlines()
            successes = [entry.strip() for entry in successes]
            already_processed.extend(successes)
    print('Proteins already processed: ', len(already_processed))

    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    def process_fasta_batch(filenames_batch_for_thread):
        thread_id = multiprocessing.current_process().name
        args_batch = copy.deepcopy(args)
        local_dir_batch = os.path.join(local_dir, str(thread_id))

        template_pipeline_runner_for_thread = TemplatePipelineRunner(
            hhsearch_binary_path=args_batch.hhsearch_binary_path,
            pdb70_database_path=args_batch.pdb70_database_path,
            template_mmcif_dir=args_batch.template_mmcif_dir,
            max_template_date=args_batch.max_template_date,
            max_templates=args_batch.max_templates,
            kalign_binary_path=args_batch.kalign_binary_path,
            release_dates_path=args_batch.release_dates_path,
            obsolete_pdbs_path=args_batch.obsolete_pdbs_path,
            no_cpus=args_batch.cpus,
        )

        # Delete and recreate the local dir for this batch
        if os.path.exists(local_dir_batch):
            os.system("rm -rf " + local_dir_batch)
        os.makedirs(local_dir_batch)

        # Copy this subdirectory batch from fastas_dir to local_dir
        print('New batch')
        filenames_processed = []
        for filename in filenames_batch_for_thread:
            print(filename)
            fasta_filename = os.path.basename(filename)
            if fasta_filename in already_processed:
                print(f"Skipping as it's already processed: {fasta_filename}")
                logging.warning(f"Skipping as it's already processed {fasta_filename}")
                continue
            # ### TODO SK: To prevent dupes, change this
            # if os.path.exists(os.path.join(args_batch.msa_dir, pathlib.Path(filename).stem, 'alignments','pdb70_hits.hhr')):
            #     print(f"Skipping as it's already processed: {fasta_filename}")
            #     logging.warning(f"Skipping as it's already processed {fasta_filename}")
            #     filenames_processed.append(filename)
            #     continue

            src_dir_path = os.path.join(fastas_dir, fasta_filename)
            dest_dir_path = os.path.join(local_dir_batch, fasta_filename)
            shutil.copy(src_dir_path, dest_dir_path)
            filenames_processed.append(filename)

        args_batch.fastas_dir = local_dir_batch

        try:
            main(args_batch, template_pipeline_runner_for_thread)
        except Exception as e:
            # If any issues, we need to skip this batch entirely for now,
            logging.error('Error in template pipeline! Skipping this batch. Caught exception:')
            logging.exception(e)
            return []
        # rsync the outputs to exx6

        return filenames_processed


    for filenames_batch in walk_dirs(fastas_dir, args.n_threads*10):

        if os.path.exists(local_dir):
            os.system("rm -rf " + local_dir)
            os.makedirs(local_dir)

        filenames_batch_for_thread = [filenames_batch[i:i+10] for i in range(0, len(filenames_batch), 10)]
        pool = Pool(args.n_threads)
        succesfully_processed_in_batch = pool.map(process_fasta_batch, filenames_batch_for_thread)

        # for per_thread in succesfully_processed_in_batch:
        #     for filename in per_thread:
        #         fasta_filename = os.path.basename(filename)
        #         succesfully_processed.append(fasta_filename)

        with open('succesfull.txt', 'a') as f:
            for per_thread in succesfully_processed_in_batch:
                for filename in per_thread:
                    fasta_filename = os.path.basename(filename)
                    if fasta_filename not in already_processed:
                        f.write(fasta_filename+'\n')
            # for fasta_ in succesfully_processed:
            #     f.write(fasta_+'\n')

