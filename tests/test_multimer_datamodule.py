# Copyright 2021 AlQuraishi Laboratory
# Dingquan Yu @ EMBL-Hamburg Kosinski group
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from pathlib import Path
import shutil
import pickle
import torch
import torch.nn as nn
import numpy as np
from functools import partial
import unittest
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.config import model_config
from openfold.data.data_modules import OpenFoldMultimerDataModule,OpenFoldDataModule
from openfold.model.model import AlphaFold
from openfold.utils.loss import AlphaFoldMultimerLoss
from tests.config import consts
import logging
logger = logging.getLogger(__name__)
import os

class TestMultimerDataModule(unittest.TestCase):
    def setUp(self):
        """
        Set up model config

        use model_1_multimer_v3 for now
        """
        self.config = model_config(
        "model_1_multimer_v3", 
        train=True, 
        low_prec=True)
        self.data_module = OpenFoldMultimerDataModule(
        config=self.config.data, 
        batch_seed=42,
        train_epoch_len=100,
        template_mmcif_dir = "/g/alphafold/AlphaFold_DBs/2.3.0/pdb_mmcif/mmcif_files/",
        template_release_dates_cache_path=os.path.join(os.getcwd(),"tests/test_data/mmcif_cache.json"),
        max_template_date="2500-01-01",
        train_data_dir=os.path.join(os.getcwd(),"tests/test_data/mmcifs"),
        train_alignment_dir=os.path.join(os.getcwd(),"tests/test_data/alignments/"),
        kalign_binary_path=shutil.which('kalign'),
        train_mmcif_data_cache_path=os.path.join(os.getcwd(),
                                                 "tests/test_data/train_mmcifs_cache.json"),
        train_chain_data_cache_path=os.path.join(os.getcwd(),
                                                 "tests/test_data/train_chain_data_cache.json"),
    )
        # setup model
        self.c = model_config(consts.model, train=True)
        
        self.c.loss.masked_msa.num_classes = 22 # somehow need overwrite this part in multimer loss config
        self.c.model.evoformer_stack.no_blocks = 4  # no need to go overboard here
        self.c.model.evoformer_stack.blocks_per_ckpt = None  # don't want to set up
        # deepspeed for this test
        self.model = AlphaFold(self.c)
        self.multimer_loss = AlphaFoldMultimerLoss(self.c.loss)

    def testPrepareData(self):
        self.data_module.prepare_data()
        self.data_module.setup()
        train_dataset = self.data_module.train_dataset
        all_chain_features,ground_truth = train_dataset[1]
        add_batch_size_dimension = lambda t: (
            t.unsqueeze(0)
        )
        all_chain_features = tensor_tree_map(add_batch_size_dimension,all_chain_features)
        with torch.no_grad():
            out = self.model(all_chain_features)
            self.multimer_loss(out,(all_chain_features,ground_truth))