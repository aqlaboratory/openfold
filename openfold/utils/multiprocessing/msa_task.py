# SPDX-FileCopyrightText: Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

import os
import logging
import shutil
from openfold.data.tools import hhsearch, hmmsearch
from openfold.data import data_pipeline
from scripts.precompute_embeddings import EmbeddingGenerator
from openfold.utils.multiprocessing.base_pipeline import BaseTask

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)


class MsaTask(BaseTask):
    def __init__(self, *args, **kwargs):
        super().__init__('msa', *args, **kwargs)

    def setup(self, params):
        self.params = params
        args = self.params['args']
        if args.use_precomputed_alignments is None:
            if "multimer" in args.config_preset:
                template_searcher = hmmsearch.Hmmsearch(
                    binary_path=args.hmmsearch_binary_path,
                    hmmbuild_binary_path=args.hmmbuild_binary_path,
                    database_path=args.pdb_seqres_database_path,
                )
            else:
                template_searcher = hhsearch.HHSearch(
                    binary_path=args.hhsearch_binary_path,
                    databases=[args.pdb70_database_path],
                )

            # In seqemb mode, use AlignmentRunner only to generate templates
            if args.use_single_seq_mode:
                self.alignment_runner = data_pipeline.AlignmentRunner(
                    jackhmmer_binary_path=args.jackhmmer_binary_path,
                    uniref90_database_path=args.uniref90_database_path,
                    template_searcher=template_searcher,
                    no_cpus=args.cpus,
                )
                self.embedding_generator = EmbeddingGenerator()
            else:
                self.alignment_runner = data_pipeline.AlignmentRunner(
                    jackhmmer_binary_path=args.jackhmmer_binary_path,
                    hhblits_binary_path=args.hhblits_binary_path,
                    uniref90_database_path=args.uniref90_database_path,
                    mgnify_database_path=args.mgnify_database_path,
                    bfd_database_path=args.bfd_database_path,
                    uniref30_database_path=args.uniref30_database_path,
                    uniclust30_database_path=args.uniclust30_database_path,
                    uniprot_database_path=args.uniprot_database_path,
                    template_searcher=template_searcher,
                    use_small_bfd=args.bfd_database_path is None,
                    no_cpus=args.cpus
                )

    def teardown(self):
        pass

    def process_one(self, item):
        tag = item['tag']
        tags = item['tags']
        seqs = item['seqs']
        args = self.params['args']
        alignment_dir = self.params['alignment_dir']
        for tag, seq in zip(tags, seqs):
            tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
            with open(tmp_fasta_path, "w") as fp:
                fp.write(f">{tag}\n{seq}")

            local_alignment_dir = os.path.join(alignment_dir, tag)

            if args.use_precomputed_alignments is None:
                logger.info(f"Generating alignments for {tag}...")

                os.makedirs(local_alignment_dir, exist_ok=True)
                
                # In seqemb mode, use AlignmentRunner only to generate templates
                if args.use_single_seq_mode:
                    self.embedding_generator.run(tmp_fasta_path, alignment_dir)

                self.alignment_runner.run(
                    tmp_fasta_path, local_alignment_dir
                )
            else:
                logger.info(
                    f"Using precomputed alignments for {tag} at {alignment_dir}..."
                )

            # Remove temporary FASTA file
            os.remove(tmp_fasta_path)

        return item
