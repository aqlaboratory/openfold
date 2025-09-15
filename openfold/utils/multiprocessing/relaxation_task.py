# SPDX-FileCopyrightText: Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

import os
import pickle
import logging
from openfold.utils.script_utils import relax_protein
from openfold.utils.multiprocessing.base_pipeline import BaseTask

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

class RelaxationTask(BaseTask):
    def __init__(self):
        super().__init__('relax')

    def setup(self, params):
        self.params = params

    def teardown(self):
        pass

    def process_one(self, item):
        args = self.params['args']
        config = self.params['config']
        output_name = item['output_name']
        if (args.relaxation_device != 'skip') or (args.skip_relaxation):
            for unrelaxed_protein, output_directory, out in zip(item['unrelaxed_protein_list'], item['output_directory_list'], item['out']):
                relax_protein(config, 
                            args.relaxation_device, 
                            unrelaxed_protein, 
                            output_directory, 
                            output_name,
                            args.cif_output)

                if args.save_outputs:
                    output_dict_path = os.path.join(
                        output_directory, f'{output_name}_output_dict.pkl'
                    )
                    with open(output_dict_path, "wb") as fp:
                        pickle.dump(out, fp, protocol=pickle.HIGHEST_PROTOCOL)

                    logger.info(f"Model output written to {output_dict_path}...")

        return item
