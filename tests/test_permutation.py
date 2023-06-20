# Copyright 2021 AlQuraishi Laboratory
#
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
from openfold.data import data_transforms
from openfold.model.model import AlphaFold
from openfold.utils.tensor_utils import tensor_tree_map
from tests.config import consts
from tests.data_utils import (
    random_template_feats,
    random_extra_msa_feats,
)
from tests.data_utils import load_labels
from openfold.data.data_transforms import make_msa_feat
import logging
logger = logging.getLogger(__name__)
import os

class TestPermutation(unittest.TestCase):
    def setUp(self):
        """
        Firstly setup model configs and model as in
        test_model.py

        In the test case, use PDB ID 1e4k as the label
        """
        self.multimer_feature_path=os.path.join(os.getcwd(),"tests/test_data/example_multimer_processed_feature.pkl")
        self.label_dir = os.path.join(os.getcwd(),"tests/test_data")

    def test_dry_run(self):
        c = model_config(consts.model, train=True)
        c.model.evoformer_stack.no_blocks = 4  # no need to go overboard here
        c.model.evoformer_stack.blocks_per_ckpt = None  # don't want to set up
        # deepspeed for this test

        model = AlphaFold(c)
        label_ids = ["1e4k_A","1e4k_B","1e4k_C"]
        sequence_ids = ["P01857","P01857","O75015"]
        features = pickle.load(open(self.multimer_feature_path,"rb"))

        #
        # I suppose between_segment_residues are always 0 ?
        # #
        num_res = features['aatype'].shape[0]
        protein = {'between_segment_residues': torch.tensor([0]*num_res,dtype=torch.int32),
                   'msa': torch.tensor(features['msa'], dtype=torch.int64),
                   'deletion_matrix': torch.tensor(features['deletion_matrix']),
                   'aatype': torch.tensor(features['aatype'],dtype=torch.int64)}
        protein = make_msa_feat.__wrapped__(protein)
        print(f"protein now is {type(protein)}")
        for k,v in protein.items():
            print(f"{k},{v.size()}")
        # if consts.is_multimer:
        #     #
        #     # Modify asym_id, entity_id and sym_id so that it encodes 
        #     # 2 chains
        #     # #
        #     asym_id = [1]*11 + [2]*11
        #     batch["asym_id"] = torch.tensor(asym_id,dtype=torch.float64)
        #     batch["entity_id"] = torch.randint(0, 1, size=(n_res,))
        #     batch["sym_id"] = torch.tensor(asym_id,dtype=torch.float64)
        #     batch["extra_deletion_matrix"] = torch.randint(0, 2, size=(n_extra_seq, n_res))
        # add_recycling_dims = lambda t: (
        #     t.unsqueeze(-1).expand(*t.shape, c.data.common.max_recycling_iters)
        # )
        # print(f"max_recycling_iters is {c.data.common.max_recycling_iters}")
        # batch = tensor_tree_map(add_recycling_dims, batch)

        # with torch.no_grad():
        #     out = model(batch)
        #     print("finished running multimer forward")
        #     print(f"out is {type(out)} and has keys {out.keys()}")
        #     print(f"final_atom_positions is {out['final_atom_positions'].shape}")