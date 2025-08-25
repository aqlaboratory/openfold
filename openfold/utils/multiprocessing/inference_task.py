# SPDX-FileCopyrightText: Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

import os
import time
import logging
import numpy as np
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.script_utils import (load_models_from_command_line, run_model, prep_output)
from openfold.utils.trace_utils import trace_model_
from openfold.np import protein
from openfold.utils.multiprocessing.base_pipeline import BaseTask

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

class InferenceTask(BaseTask):
    def __init__(self):
        super().__init__('infer')

    def setup(self, params):
        self.params = params
        args = params['args']
        self.model_generator = list(load_models_from_command_line(
            params['config'],
            args.model_device,
            args.openfold_checkpoint_path,
            args.jax_param_path,
            args.output_dir)
        )

    def process_one(self, item):
        tag = item['tag']
        args = self.params['args']
        config = self.params['config']
        feature_dict = item['feature_dict']
        processed_feature_dict = item['processed_feature_dict']
        item['unrelaxed_protein_list'] = []
        item['output_directory_list'] = []
        item['out'] = []
        for model, output_directory in self.model_generator:
            if args.trace_model:
                rounded_seqlen = item['rounded_seqlen']
                if rounded_seqlen > cur_tracing_interval:
                    logger.info(
                        f"Tracing model at {rounded_seqlen} residues..."
                    )
                    t = time.perf_counter()
                    trace_model_(model, processed_feature_dict)
                    tracing_time = time.perf_counter() - t
                    logger.info(
                        f"Tracing time: {tracing_time}"
                    )
                    cur_tracing_interval = rounded_seqlen

            out = run_model(model, processed_feature_dict, tag, args.output_dir)

            # Toss out the recycling dimensions --- we don't need them anymore
            processed_feature_dict = tensor_tree_map(
                lambda x: np.array(x[..., -1].cpu()),
                processed_feature_dict
            )
            out = tensor_tree_map(lambda x: np.array(x.cpu()), out)

            unrelaxed_protein = prep_output(
                out,
                processed_feature_dict,
                feature_dict,
                config,
                args.config_preset,
                args.multimer_ri_gap,
                args.subtract_plddt
            )
            unrelaxed_file_suffix = "_unrelaxed.pdb"
            if args.cif_output:
                unrelaxed_file_suffix = "_unrelaxed.cif"
            unrelaxed_output_path = os.path.join(
                output_directory, f'{item["output_name"]}{unrelaxed_file_suffix}'
            )

            with open(unrelaxed_output_path, 'w') as fp:
                if args.cif_output:
                    fp.write(protein.to_modelcif(unrelaxed_protein))
                else:
                    fp.write(protein.to_pdb(unrelaxed_protein))

            logger.info(f"Output written to {unrelaxed_output_path}...")

            item['unrelaxed_protein_list'].append(unrelaxed_protein)
            item['output_directory_list'].append(output_directory)
            if args.save_outputs:
                item['out'].append(out)
            else:
                item['out'].append(None) # Just to match counts with other lists in the item
        del item['feature_dict']
        del item['processed_feature_dict']
        return item