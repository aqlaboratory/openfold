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

import torch

import unittest
from openfold.utils.loss import AlphaFoldMultimerLoss
from openfold.utils.loss import get_least_asym_entity_or_longest_length,merge_labels,pad_features
from openfold.utils.tensor_utils import tensor_tree_map
import math

class TestPermutation(unittest.TestCase):
    def setUp(self):
        """
        create fake input structure features 
        and rotation matrices 
        """

        theta = math.pi/4
        self.rotation_matrix_z = torch.tensor([
            [math.cos(theta),-math.sin(theta),0],
            [math.sin(theta),math.cos(theta),0],
            [0,0,1]
        ],device='cuda')
        self.rotation_matrix_x = torch.tensor([
            [1,0,0],
            [0,math.cos(theta),-math.sin(theta)],
            [0,math.sin(theta),math.cos(theta)],
        ],device='cuda')
        self.rotation_matrix_y = torch.tensor([
            [math.cos(theta),0,math.sin(theta)],
            [0,1,0],
            [-math.sin(theta),1,math.cos(theta)],
        ],device='cuda')
        self.chain_a_num_res=9
        self.chain_b_num_res=13
        # below create default fake ground truth structures for a hetero-pentamer A2B3
        self.residue_index=list(range(self.chain_a_num_res))*2 + list(range(self.chain_b_num_res))*3
        self.num_res = self.chain_a_num_res*2 + self.chain_b_num_res*3
        self.asym_id = torch.tensor([[1]*self.chain_a_num_res+[2]*self.chain_a_num_res+[3]*self.chain_b_num_res+[4]*self.chain_b_num_res+[5]*self.chain_b_num_res],device='cuda')
        self.sym_id = self.asym_id
        self.entity_id = torch.tensor([[1]*(self.chain_a_num_res*2)+[2]*(self.chain_b_num_res*3)],device='cuda')

    def test_1_selecting_anchors(self):
        self.batch = {
            'asym_id':self.asym_id,
            'sym_id':self.sym_id,
            'entity_id':self.entity_id,
            'seq_length':torch.tensor([57])
        }
        anchor_gt_asym, anchor_pred_asym=get_least_asym_entity_or_longest_length(self.batch)
        self.assertIn(int(anchor_gt_asym),[1,2])
        self.assertNotIn(int(anchor_gt_asym),[3,4,5])
        self.assertIn(int(anchor_pred_asym),[1,2])
        self.assertNotIn(int(anchor_pred_asym),[3,4,5])

    def test_2_permutation_pentamer(self):
        batch = {
            'asym_id':self.asym_id,
            'sym_id':self.sym_id,
            'entity_id':self.entity_id,
            'seq_length':torch.tensor([57]),
            'aatype':torch.randint(21,size=(1,57))
        }
        batch['asym_id'] = batch['asym_id'].reshape(1,self.num_res)
        batch["residue_index"] = torch.tensor([self.residue_index],device='cuda')
        # create fake ground truth atom positions        
        chain_a1_pos = torch.randint(15,(self.chain_a_num_res,3*37),
                                     device='cuda',dtype=torch.float).reshape(1,self.chain_a_num_res,37,3)
        chain_a2_pos = torch.matmul(chain_a1_pos,self.rotation_matrix_x)+10

        chain_b1_pos = torch.randint(low=15,high=30,size=(self.chain_b_num_res,3*37),
                                     device='cuda',dtype=torch.float).reshape(1,self.chain_b_num_res,37,3)
        chain_b2_pos = torch.matmul(chain_b1_pos,self.rotation_matrix_y)+10
        chain_b3_pos =  torch.matmul(torch.matmul(chain_b1_pos,self.rotation_matrix_z),self.rotation_matrix_x)+30
        # Below permutate predicted chain positions 
        pred_atom_position = torch.cat((chain_a2_pos,chain_a1_pos,chain_b2_pos,chain_b3_pos,chain_b1_pos),dim=1)
        pred_atom_mask = torch.ones((1,self.num_res,37),device='cuda')
        out = {
            'final_atom_positions':pred_atom_position,
            'final_atom_mask':pred_atom_mask
        }

        true_atom_position = torch.cat((chain_a1_pos,chain_a2_pos,chain_b1_pos,chain_b2_pos,chain_b3_pos),dim=1)
        true_atom_mask = torch.cat((torch.ones((1,self.chain_a_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_a_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_b_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_b_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_b_num_res,37),device='cuda')),dim=1)
        batch['all_atom_positions'] = true_atom_position
        batch['all_atom_mask'] = true_atom_mask
        
        dim_dict = AlphaFoldMultimerLoss.determine_split_dim(batch) 
        aligns,_ = AlphaFoldMultimerLoss.multi_chain_perm_align(out,batch,
                                                                dim_dict,
                                                                permutate_chains=True)
        possible_outcome = [[(0,1),(1,0),(2,3),(3,4),(4,2)],[(0,0),(1,1),(2,3),(3,4),(4,2)]]
        wrong_outcome = [[(0,1),(1,0),(2,4),(3,2),(4,3)],[(0,0),(1,1),(2,2),(3,3),(4,4)]]
        self.assertIn(aligns,possible_outcome)
        self.assertNotIn(aligns,wrong_outcome)

    def test_3_merge_labels(self):
        nres_pad = 325 - 57 # suppose the cropping size is 325
        batch = {
            'asym_id':pad_features(self.asym_id,nres_pad,pad_dim=1),
            'sym_id':pad_features(self.sym_id,nres_pad,pad_dim=1),
            'entity_id':pad_features(self.entity_id,nres_pad,pad_dim=1),
            'aatype':torch.randint(21,size=(1,325)),
            'seq_length':torch.tensor([57])
        }
        batch['asym_id'] = batch['asym_id'].reshape(1,325)
        batch["residue_index"] = pad_features(torch.tensor(self.residue_index).reshape(1,57),nres_pad,pad_dim=1)
        # create fake ground truth atom positions        
        chain_a1_pos = torch.randint(15,(self.chain_a_num_res,3*37),
                                     device='cuda',dtype=torch.float).reshape(1,self.chain_a_num_res,37,3)
        chain_a2_pos = torch.matmul(chain_a1_pos,self.rotation_matrix_x)+10

        chain_b1_pos = torch.randint(low=15,high=30,size=(self.chain_b_num_res,3*37),
                                     device='cuda',dtype=torch.float).reshape(1,self.chain_b_num_res,37,3)
        chain_b2_pos = torch.matmul(chain_b1_pos,self.rotation_matrix_y)+10
        chain_b3_pos =  torch.matmul(torch.matmul(chain_b1_pos,self.rotation_matrix_z),self.rotation_matrix_x)+30
        # Below permutate predicted chain positions 
        pred_atom_position = torch.cat((chain_a2_pos,chain_a1_pos,chain_b2_pos,chain_b3_pos,chain_b1_pos),dim=1)
        pred_atom_mask = torch.ones((1,self.num_res,37),device='cuda')
        pred_atom_position = pad_features(pred_atom_position,nres_pad,pad_dim=1)
        pred_atom_mask = pad_features(pred_atom_mask,nres_pad,pad_dim=1)
        out = {
            'final_atom_positions':pred_atom_position,
            'final_atom_mask':pred_atom_mask
        }
        true_atom_position = torch.cat((chain_a1_pos,chain_a2_pos,chain_b1_pos,chain_b2_pos,chain_b3_pos),dim=1)
        true_atom_mask = torch.cat((torch.ones((1,self.chain_a_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_a_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_b_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_b_num_res,37),device='cuda'),
                                    torch.ones((1,self.chain_b_num_res,37),device='cuda')),dim=1)
        batch['all_atom_positions'] = pad_features(true_atom_position,nres_pad,pad_dim=1)
        batch['all_atom_mask'] = pad_features(true_atom_mask,nres_pad=nres_pad,pad_dim=1)
        
        tensor_to_cuda = lambda t: t.to('cuda')
        batch = tensor_tree_map(tensor_to_cuda,batch)
        dim_dict = AlphaFoldMultimerLoss.determine_split_dim(batch) 
        aligns,per_asym_residue_index = AlphaFoldMultimerLoss.multi_chain_perm_align(out,
                                                                batch,
                                                                dim_dict,
                                                                permutate_chains=True)
        
        labels = AlphaFoldMultimerLoss.split_ground_truth_labels(batch,dim_dict=dim_dict,
                                                                 REQUIRED_FEATURES=[i for i in batch.keys() if i in dim_dict])
        
        labels = merge_labels(per_asym_residue_index,labels,aligns,
                              original_nres=batch['aatype'].shape[-1])

        self.assertTrue(torch.equal(labels['residue_index'],batch['residue_index']))
        
        expected_permutated_gt_pos = torch.cat((chain_a2_pos,chain_a1_pos,chain_b2_pos,chain_b3_pos,chain_b1_pos),dim=1)
        expected_permutated_gt_pos = pad_features(expected_permutated_gt_pos,nres_pad,pad_dim=1)
        self.assertTrue(torch.equal(labels['all_atom_positions'],expected_permutated_gt_pos))