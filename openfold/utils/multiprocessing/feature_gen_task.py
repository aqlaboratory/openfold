# SPDX-FileCopyrightText: Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

import os
import math
import torch
from openfold.data import feature_pipeline, templates, data_pipeline
from openfold.utils.trace_utils import pad_feature_dict_seq
from openfold.utils.multiprocessing.base_pipeline import BaseTask

TRACING_INTERVAL = 50

def round_up_seqlen(seqlen):
    return int(math.ceil(seqlen / TRACING_INTERVAL)) * TRACING_INTERVAL

class FeatureGenTask(BaseTask):
    def __init__(self):
        super().__init__('feature_gen')

    def setup(self, params):
        self.params = params
        args = self.params['args']
        config = self.params['config']
        is_multimer = "multimer" in args.config_preset
        is_custom_template = "use_custom_template" in args and args.use_custom_template
        if is_custom_template:
            template_featurizer = templates.CustomHitFeaturizer(
                mmcif_dir=args.template_mmcif_dir,
                max_template_date="9999-12-31", # just dummy, not used
                max_hits=-1, # just dummy, not used
                kalign_binary_path=args.kalign_binary_path
                )
        elif is_multimer:
            template_featurizer = templates.HmmsearchHitFeaturizer(
                mmcif_dir=args.template_mmcif_dir,
                max_template_date=args.max_template_date,
                max_hits=config.data.predict.max_templates,
                kalign_binary_path=args.kalign_binary_path,
                release_dates_path=args.release_dates_path,
                obsolete_pdbs_path=args.obsolete_pdbs_path
            )
        else:
            template_featurizer = templates.HhsearchHitFeaturizer(
                mmcif_dir=args.template_mmcif_dir,
                max_template_date=args.max_template_date,
                max_hits=config.data.predict.max_templates,
                kalign_binary_path=args.kalign_binary_path,
                release_dates_path=args.release_dates_path,
                obsolete_pdbs_path=args.obsolete_pdbs_path
            )
        self.data_processor = data_pipeline.DataPipeline(
            template_featurizer=template_featurizer,
        )
        if is_multimer:
            self.data_processor = data_pipeline.DataPipelineMultimer(
                monomer_data_pipeline=self.data_processor,
            )
        self.feature_processor = feature_pipeline.FeaturePipeline(config.data)
        self.feature_dicts = {}

    def teardown(self):
        pass

    def process_one(self, item): # item is a tuple of (tag, tags, seqs)
        tag = item['tag']
        tags = item['tags']
        seqs = item['seqs']
        alignment_dir = self.params['alignment_dir']
        args = self.params['args']
        is_multimer = "multimer" in args.config_preset

        feature_dict = self.feature_dicts.get(tag, None)
        if feature_dict is None:
            # Start of generate_feature_dict
            tmp_fasta_path = os.path.join(self.params['args'].output_dir, f"tmp_{os.getpid()}_{tag}.fasta")

            if "multimer" in args.config_preset:
                with open(tmp_fasta_path, "w") as fp:
                    fp.write(
                        '\n'.join([f">{tag}\n{seq}" for tag, seq in zip(tags, seqs)])
                    )
                feature_dict = self.data_processor.process_fasta(
                    fasta_path=tmp_fasta_path, alignment_dir=alignment_dir,
                )
            elif len(seqs) == 1:
                tag = tags[0]
                seq = seqs[0]
                with open(tmp_fasta_path, "w") as fp:
                    fp.write(f">{tag}\n{seq}")

                local_alignment_dir = os.path.join(alignment_dir, tag)
                feature_dict = self.data_processor.process_fasta(
                    fasta_path=tmp_fasta_path,
                    alignment_dir=local_alignment_dir,
                    seqemb_mode=args.use_single_seq_mode,
                )
            else:
                with open(tmp_fasta_path, "w") as fp:
                    fp.write(
                        '\n'.join([f">{tag}\n{seq}" for tag, seq in zip(tags, seqs)])
                    )
                feature_dict = self.data_processor.process_multiseq_fasta(
                    fasta_path=tmp_fasta_path, super_alignment_dir=alignment_dir,
                )
            os.remove(tmp_fasta_path)
            # End of generate_feature_dict

            if args.trace_model:
                n = feature_dict["aatype"].shape[-2]
                rounded_seqlen = round_up_seqlen(n)
                item['rounded_seqlen'] = rounded_seqlen
                feature_dict = pad_feature_dict_seq(
                    feature_dict, rounded_seqlen,
                )
            
            self.feature_dicts[tag] = feature_dict

        processed_feature_dict = self.feature_processor.process_features(
            feature_dict, mode='predict', is_multimer=is_multimer
        )

        processed_feature_dict = {
            k: torch.as_tensor(v, device=args.model_device)
            for k, v in processed_feature_dict.items()
        }
        
        # Move outputs to shared memory for zero-copy forwarding to inference process
        for k, v in processed_feature_dict.items():
            v.share_memory_()
        item['feature_dict'] = feature_dict
        item['processed_feature_dict'] = processed_feature_dict
        return item

