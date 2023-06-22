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
from .unifold_permutation import multi_chain_perm_align
import logging
logger = logging.getLogger(__name__)
import os
from tests.data_utils import (
    random_template_feats,
    random_extra_msa_feats,
)
class TestPermutation(unittest.TestCase):
    def setUp(self):
        """
        Firstly setup model configs and model as in
        test_model.py

        In the test case, use PDB ID 1e4k as the label
        """
        self.test_data_dir = os.path.join(os.getcwd(),"tests/test_data")
        self.label_ids = ['label_1','label_2']

    def test_dry_run(self):
        n_seq = consts.n_seq
        n_templ = consts.n_templ
        n_res = consts.n_res
        n_extra_seq = consts.n_extra

        c = model_config(consts.model, train=True)
        c.model.evoformer_stack.no_blocks = 4  # no need to go overboard here
        c.model.evoformer_stack.blocks_per_ckpt = None  # don't want to set up
        # deepspeed for this test

        model = AlphaFold(c)
        example_label = [pickle.load(open(os.path.join(self.test_data_dir,f"{i}.pkl"),'rb')) 
                         for i in self.label_ids]
        batch = {}
        tf = torch.randint(c.model.input_embedder.tf_dim - 1, size=(n_res,))
        batch["target_feat"] = nn.functional.one_hot(
            tf, c.model.input_embedder.tf_dim
        ).float()
        batch["aatype"] = torch.argmax(batch["target_feat"], dim=-1)
        print(f"target_feat shape is {batch['target_feat'].size()}")
        print(f"batch_dim is {batch['target_feat'].shape[:-2]}")
        batch["residue_index"] = torch.arange(n_res)

        batch["msa_feat"] = torch.rand((n_seq, n_res, c.model.input_embedder.msa_dim))
        t_feats = random_template_feats(n_templ, n_res)
        batch.update({k: torch.tensor(v) for k, v in t_feats.items()})
        extra_feats = random_extra_msa_feats(n_extra_seq, n_res)
        batch.update({k: torch.tensor(v) for k, v in extra_feats.items()})
        batch["msa_mask"] = torch.randint(
            low=0, high=2, size=(n_seq, n_res)
        ).float()
        batch["seq_mask"] = torch.randint(low=0, high=2, size=(n_res,)).float()
        batch.update(data_transforms.make_atom14_masks(batch))
        batch["no_recycling_iters"] = torch.tensor(2.)

        if consts.is_multimer:
            #
            # Modify asym_id, entity_id and sym_id so that it encodes 
            # 2 chains
            # #
            asym_id = [1]*9+[2]*13
            batch["asym_id"] = torch.tensor(asym_id,dtype=torch.float64)
            # batch["entity_id"] = torch.randint(0, 1, size=(n_res,))
            batch['entity_id'] = torch.tensor(asym_id,dtype=torch.float64)
            batch["sym_id"] = torch.tensor(asym_id,dtype=torch.float64)
            batch["num_sym"] = torch.tensor([2]*22,dtype=torch.int64) # currently there are just 2 chains
            batch["extra_deletion_matrix"] = torch.randint(0, 2, size=(n_extra_seq, n_res))
        add_recycling_dims = lambda t: (
            t.unsqueeze(-1).expand(*t.shape, c.data.common.max_recycling_iters)
        )
        print(f"max_recycling_iters is {c.data.common.max_recycling_iters}")
        input_batch = tensor_tree_map(add_recycling_dims, batch)

        with torch.no_grad():
            out = model(input_batch)
            print("finished running multimer forward")
            print(f"out is {type(out)} and has keys {out.keys()}")
            print(f"final_atom_positions is {out['final_atom_positions'].shape}")
            print(f"out itpm score is {out['iptm_score']}")
            multi_chain_perm_align(out,batch,example_label)