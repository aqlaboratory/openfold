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
import pickle
import torch
import torch.nn as nn
import numpy as np
import unittest
from openfold.config import model_config
from openfold.data.data_modules import OpenFoldDataModule
from openfold.utils.tensor_utils import tensor_tree_map
from tests.config import consts
import logging
logger = logging.getLogger(__name__)
import os
import io, contextlib
from tests.data_utils import (
    random_template_feats,
    random_extra_msa_feats,
    random_affines_vector, random_affines_4x4
)
from openfold.utils.rigid_utils import (
    Rotation,
    Rigid,
)

class TestPermutation(unittest.TestCase):
    def setUp(self):
        """
        Set up model config

        use model_1_multimer_v3 for now
        """
        self.config = model_config(
        "model_1_multimer_v3", 
        train=True, 
        low_prec=True)
        self.data_module = OpenFoldDataModule(
        config=self.config.data, 
        batch_seed=42,
        template_mmcif_dir: str,
        max_template_date: str,
        val_data_dir: Optional[str] = None,
        val_alignment_dir: Optional[str] = None,
        
    )