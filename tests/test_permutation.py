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

import math
import torch

import unittest
from openfold.utils.multi_chain_permutation import (pad_features, get_least_asym_entity_or_longest_length,
                                                    compute_permutation_alignment, split_ground_truth_labels,
                                                    merge_labels)


class TestPermutation(unittest.TestCase):
    def setUp(self):
        """
        create fake input structure features 
        and rotation matrices 
        """

        theta = math.pi / 4
        device = 'cpu'
        self.rotation_matrix_z = torch.tensor([
            [math.cos(theta), -math.sin(theta), 0],
            [math.sin(theta), math.cos(theta), 0],
            [0, 0, 1]
        ], device=device)
        self.rotation_matrix_x = torch.tensor([
            [1, 0, 0],
            [0, math.cos(theta), -math.sin(theta)],
            [0, math.sin(theta), math.cos(theta)],
        ], device=device)
        self.rotation_matrix_y = torch.tensor([
            [math.cos(theta), 0, math.sin(theta)],
            [0, 1, 0],
            [-math.sin(theta), 1, math.cos(theta)],
        ], device=device)
        self.chain_a_num_res = 9
        self.chain_b_num_res = 13
        # below create default fake ground truth structures for a hetero-pentamer A2B3
        self.residue_index = list(
            range(self.chain_a_num_res)) * 2 + list(range(self.chain_b_num_res)) * 3
        self.num_res = self.chain_a_num_res * 2 + self.chain_b_num_res * 3
        self.asym_id = torch.tensor([[1] * self.chain_a_num_res + [2] * self.chain_a_num_res + [
            3] * self.chain_b_num_res + [4] * self.chain_b_num_res + [5] * self.chain_b_num_res], device=device)
        self.sym_id = self.asym_id
        self.entity_id = torch.tensor([[1] * (self.chain_a_num_res * 2) + [2] * (self.chain_b_num_res * 3)],
                                      device=device)

    def test_1_selecting_anchors(self):
        batch = {
            'asym_id': self.asym_id,
            'sym_id': self.sym_id,
            'entity_id': self.entity_id,
            'seq_length': torch.tensor([57])
        }
        anchor_gt_asym, anchor_pred_asym = get_least_asym_entity_or_longest_length(
            batch, batch['asym_id'])
        anchor_gt_asym = int(anchor_gt_asym)
        anchor_pred_asym = {int(i) for i in anchor_pred_asym}
        expected_anchors = {1, 2}
        expected_non_anchors = {3, 4, 5}

        self.assertIn(anchor_gt_asym, expected_anchors)
        self.assertNotIn(anchor_gt_asym, expected_non_anchors)
        # Check that predicted anchors are within expected anchor set
        self.assertEqual(anchor_pred_asym, expected_anchors & anchor_pred_asym)
        self.assertEqual(set(), anchor_pred_asym & expected_non_anchors)

    def test_2_permutation_pentamer(self):
        """
        Test the permutation results on a pentamer A2B3, in which protein A has 9 residues
        and protein B has 13 residues. 

        Expected outputs:
            Only protein A should be selected as an anchor thus, in the output list, either [(0,1), (1,0)] or [(0,0), (1,1)] are allowed 
            The 3 chains from protein B should ALWAYS be aligned in a way that predicted b1 to be aligned with ground truth b1, pred b2 to ground truth b2 
            as shown below:

            predicted structure:              a2 - a1 - b2 - b3 - b1
            indexes in the predicted list:    0    1    2    3    4

            ground truth structure:           a1 - a2 - b1 - b2 - b3
            indexes in the ground truth list: 0    1    2    3    4

            then the 2 protein A chains are free to be aligned by either order, thus either [(0,1),(1,0)] or [(0,0),(1,1)] is valid.

            However, the 3 protein B chains should be strictly aligned in the following order: 
            [(2,3), (3,4), (4,2)], regardless of how protein A chains are aligned. 

            Therefore, the only 2 correct permutations are :
            [(0, 1), (1, 0), (2, 3), (3, 4), (4, 2)] and 
            [(0, 0), (1, 1), (2, 3), (3, 4), (4, 2)]
        """
        batch = {
            'asym_id': self.asym_id,
            'sym_id': self.sym_id,
            'entity_id': self.entity_id,
            'seq_length': torch.tensor([57]),
            'aatype': torch.randint(21, size=(1, 57))
        }
        batch['asym_id'] = batch['asym_id'].reshape(1, self.num_res)
        batch["residue_index"] = torch.tensor([self.residue_index])
        # create fake ground truth atom positions
        chain_a1_pos = torch.randint(15, (self.chain_a_num_res, 3 * 37),
                                     dtype=torch.float).reshape(1, self.chain_a_num_res, 37, 3)
        chain_a2_pos = torch.matmul(chain_a1_pos, self.rotation_matrix_x) + 10

        chain_b1_pos = torch.randint(low=15, high=30, size=(self.chain_b_num_res, 3 * 37),
                                     dtype=torch.float).reshape(1, self.chain_b_num_res, 37, 3)
        chain_b2_pos = torch.matmul(chain_b1_pos, self.rotation_matrix_y) + 10
        chain_b3_pos = torch.matmul(torch.matmul(
            chain_b1_pos, self.rotation_matrix_z), self.rotation_matrix_x) + 30
        # Below permutate predicted chain positions
        # here the b2 chain from the ground truth is deliberately put in b1 chain's position, and predicted b3 chain to b2's position
        # and predicted b1 chain to b3's position
        pred_atom_position = torch.cat(
            (chain_a2_pos, chain_a1_pos, chain_b2_pos, chain_b3_pos, chain_b1_pos), dim=1)

        pred_atom_mask = torch.ones((1, self.num_res, 37))
        out = {
            'final_atom_positions': pred_atom_position,
            'final_atom_mask': pred_atom_mask
        }

        true_atom_position = torch.cat(
            (chain_a1_pos, chain_a2_pos, chain_b1_pos, chain_b2_pos, chain_b3_pos), dim=1)
        true_atom_mask = torch.cat((torch.ones((1, self.chain_a_num_res, 37)),
                                    torch.ones((1, self.chain_a_num_res, 37)),
                                    torch.ones((1, self.chain_b_num_res, 37)),
                                    torch.ones((1, self.chain_b_num_res, 37)),
                                    torch.ones((1, self.chain_b_num_res, 37))), dim=1)
        batch['all_atom_positions'] = true_atom_position
        batch['all_atom_mask'] = true_atom_mask

        aligns, per_asym_residue_index = compute_permutation_alignment(out, batch,
                                                                       batch)

        expected_asym_residue_index = {
            1: torch.tensor(list(range(self.chain_a_num_res))),
            2: torch.tensor(list(range(self.chain_a_num_res))),
            3: torch.tensor(list(range(self.chain_b_num_res))),
            4: torch.tensor(list(range(self.chain_b_num_res))),
            5: torch.tensor(list(range(self.chain_b_num_res)))
        }
        chain_a_permutated_chain_b_permutated = [
            (0, 1), (1, 0), (2, 3), (3, 4), (4, 2)]
        chain_a_not_permutated_chain_b_permutated = [
            (0, 0), (1, 1), (2, 3), (3, 4), (4, 2)]
        chain_a_permutated_chain_b_not_permuated = [
            (0, 1), (1, 0), (2, 2), (3, 3), (4, 4)]
        chain_a_not_permutated_chain_b_not_permuated = [
            (0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
        
        # test on the permutation alignments 
        self.assertIn(aligns, [chain_a_permutated_chain_b_permutated,
                      chain_a_not_permutated_chain_b_permutated])
        self.assertNotIn(aligns, [chain_a_permutated_chain_b_not_permuated,
                         chain_a_not_permutated_chain_b_not_permuated])

        # test on the per_aysm_residue_index 
        for k, v in expected_asym_residue_index.items():
            self.assertTrue(torch.equal(v, per_asym_residue_index[k]))

    def test_3_merge_labels(self):
        nres_pad = 325 - 57  # suppose the cropping size is 325
        batch = {
            'asym_id': self.asym_id,
            'sym_id': self.sym_id,
            'entity_id': self.entity_id,
            'aatype': torch.randint(21, size=(1, 57)),
            'seq_length': torch.tensor([57])
        }
        batch['asym_id'] = batch['asym_id'].reshape(1, 57)
        batch["residue_index"] = torch.tensor([self.residue_index])
        # create fake ground truth atom positions
        chain_a1_pos = torch.randint(15, (self.chain_a_num_res, 3 * 37),
                                     dtype=torch.float).reshape(1, self.chain_a_num_res, 37, 3)
        chain_a2_pos = torch.matmul(chain_a1_pos, self.rotation_matrix_x) + 10

        chain_b1_pos = torch.randint(low=15, high=30, size=(self.chain_b_num_res, 3 * 37),
                                     dtype=torch.float).reshape(1, self.chain_b_num_res, 37, 3)
        chain_b2_pos = torch.matmul(chain_b1_pos, self.rotation_matrix_y) + 10
        chain_b3_pos = torch.matmul(torch.matmul(
            chain_b1_pos, self.rotation_matrix_z), self.rotation_matrix_x) + 30
        # Below permutate predicted chain positions
        pred_atom_position = torch.cat(
            (chain_a2_pos, chain_a1_pos, chain_b2_pos, chain_b3_pos, chain_b1_pos), dim=1)
        pred_atom_mask = torch.ones((1, self.num_res, 37))
        pred_atom_position = pad_features(
            pred_atom_position, nres_pad, pad_dim=1)
        pred_atom_mask = pad_features(pred_atom_mask, nres_pad, pad_dim=1)
        out = {
            'final_atom_positions': pred_atom_position,
            'final_atom_mask': pred_atom_mask
        }
        true_atom_position = torch.cat(
            (chain_a1_pos, chain_a2_pos, chain_b1_pos, chain_b2_pos, chain_b3_pos), dim=1)
        true_atom_mask = torch.cat((torch.ones((1, self.chain_a_num_res, 37)),
                                    torch.ones((1, self.chain_a_num_res, 37)),
                                    torch.ones((1, self.chain_b_num_res, 37)),
                                    torch.ones((1, self.chain_b_num_res, 37)),
                                    torch.ones((1, self.chain_b_num_res, 37))), dim=1)

        batch['all_atom_positions'] = true_atom_position
        batch['all_atom_mask'] = true_atom_mask

        # Below create a fake_input_features
        fake_input_features = {
            'asym_id': pad_features(self.asym_id, nres_pad, pad_dim=1),
            'sym_id': pad_features(self.sym_id, nres_pad, pad_dim=1),
            'entity_id': pad_features(self.entity_id, nres_pad, pad_dim=1),
            'aatype': torch.randint(21, size=(1, 325)),
            'seq_length': torch.tensor([57])
        }
        fake_input_features['asym_id'] = fake_input_features['asym_id'].reshape(
            1, 325)
        fake_input_features["residue_index"] = pad_features(
            torch.tensor(self.residue_index).reshape(1, 57), nres_pad, pad_dim=1)
        fake_input_features['all_atom_positions'] = pad_features(
            true_atom_position, nres_pad, pad_dim=1)
        fake_input_features['all_atom_mask'] = pad_features(
            true_atom_mask, nres_pad=nres_pad, pad_dim=1)

        # NOTE
        # batch: simulates ground_truth features
        # fake_input_features: simulates the data that are going be used as input for model.forward(fake_input_features)
        # out: simulates the output of model.forward(fake_input_features)
        aligns, per_asym_residue_index = compute_permutation_alignment(out,
                                                                       fake_input_features,
                                                                       batch)
        labels = split_ground_truth_labels(batch)

        labels = merge_labels(per_asym_residue_index, labels, aligns,
                              original_nres=batch['aatype'].shape[-1])

        self.assertTrue(torch.equal(
            labels['residue_index'], batch['residue_index']))

        expected_permutated_gt_pos = torch.cat((chain_a2_pos, chain_a1_pos, chain_b2_pos, chain_b3_pos, chain_b1_pos),
                                               dim=1)

        self.assertTrue(torch.equal(
            labels['all_atom_positions'], expected_permutated_gt_pos))
