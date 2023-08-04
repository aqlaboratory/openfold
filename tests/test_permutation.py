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
from openfold.data import data_transforms
from openfold.model.model import AlphaFold
from openfold.utils.loss import AlphaFoldMultimerLoss
from openfold.utils.tensor_utils import tensor_tree_map
from tests.config import consts
import logging
logger = logging.getLogger(__name__)
import os
from tests.data_utils import (
    random_template_feats,
    random_extra_msa_feats,
    random_affines_vector
)
from openfold.utils.rigid_utils import (
    Rigid,
)

class TestPermutation(unittest.TestCase):
    def setUp(self):
        """
        Firstly setup model configs and model as in
        test_model.py

        In the test case, use PDB ID 1e4k as the label
        """
        self.test_data_dir = os.path.join(os.getcwd(),"tests/test_data")
        self.label_ids = ['label_1','label_1','label_2','label_2','label_2']
        self.asym_id = [0]*9+[1]*9+[2]*13+[3]*13 + [4]*13

    def affine_vector_to_4x4(self,affine):
        r = Rigid.from_tensor_7(affine)
        return r.to_tensor_4x4()
    
    def test_dry_run(self):
        os.environ["CUDA_VISIBLE_DEVICES"] = "0"
        os.environ["CUDA_LAUNCH_BLOCKING"] = "1"
        n_seq = consts.n_seq
        n_templ = consts.n_templ
        n_res = len(self.asym_id)
        n_extra_seq = consts.n_extra

        c = model_config(consts.model, train=True)
        
        c.loss.masked_msa.num_classes = 22 # somehow need overwrite this part in multimer loss config
        c.model.evoformer_stack.no_blocks = 4  # no need to go overboard here
        c.model.evoformer_stack.blocks_per_ckpt = None  # don't want to set up
        # deepspeed for this test

        model = AlphaFold(c)
        multimer_loss = AlphaFoldMultimerLoss(c.loss)
        example_label = [pickle.load(open(os.path.join(self.test_data_dir,f"{i}.pkl"),'rb')) 
                         for i in self.label_ids]
        batch = {}
        tf = torch.randint(c.model.input_embedder.tf_dim - 1, size=(n_res,))
        batch["target_feat"] = nn.functional.one_hot(
            tf, c.model.input_embedder.tf_dim
        ).float()
        batch["aatype"] = torch.argmax(batch["target_feat"], dim=-1)
        batch["residue_index"] = torch.arange(n_res)

        backbone_dict ={
            "backbone_affine_tensor": torch.tensor(random_affines_vector((n_res,))),
            "backbone_affine_mask": torch.from_numpy(np.random.randint(0, 2, (n_res,)).astype(
                np.float32
            )),
            "use_clamped_fape": torch.from_numpy(np.array(0.0)),
        }
        batch['backbone_rigid_tensor'] = self.affine_vector_to_4x4(backbone_dict['backbone_affine_tensor'])
        batch['backbone_rigid_mask'] = backbone_dict['backbone_affine_mask']
        
        true_msa_dict ={
            "true_msa": torch.tensor(np.random.randint(0, 21, (n_seq,n_res))),
            "bert_mask": torch.tensor(np.random.randint(0, 2, (n_seq,n_res)).astype(
                np.float32)
            )
        }
        
        batch.update(true_msa_dict)

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
        batch["seq_length"] = torch.from_numpy(np.array([n_res] * n_res))

        if consts.is_multimer:
            #
            # Modify asym_id, entity_id and sym_id so that it encodes 
            # 2 chains
            # #
            asym_id = self.asym_id
            batch["asym_id"] = torch.tensor(asym_id,dtype=torch.float64)
            # batch["entity_id"] = torch.randint(0, 1, size=(n_res,))
            batch['entity_id'] = torch.tensor([0]*18+[1]*39,dtype=torch.float64)
            batch["sym_id"] = torch.tensor(asym_id,dtype=torch.float64)
            # batch["num_sym"] = torch.tensor([1]*18+[2]*13,dtype=torch.int64) # currently there are just 2 chains
            batch["extra_deletion_matrix"] = torch.randint(0, 2, size=(n_extra_seq, n_res))
        add_recycling_dims = lambda t: (
            t.unsqueeze(-1).expand(*t.shape, c.data.common.max_recycling_iters)
        )
        add_batch_size_dimension = lambda t: (
            t.unsqueeze(0)
        )
        batch = tensor_tree_map(add_recycling_dims, batch)
        batch = tensor_tree_map(add_batch_size_dimension, batch)
     
        with torch.no_grad():
            out = model(batch)
            print(f"finished foward on batch with batch_size dim")
            multimer_loss(out,(batch,example_label))
            # print(f"permuated_labels is {type(permutated_labels)} and keys are:\n {permutated_labels.keys()}")
