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

import math
import torch
import numpy as np
import unittest

from alphafold.utils.loss import *
from alphafold.utils.utils import T


class TestLoss(unittest.TestCase):
    def test_run_torsion_angle_loss(self):
        batch_size = 2
        n = 5

        a = torch.rand((batch_size, n, 7, 2))
        a_gt = torch.rand((batch_size, n, 7, 2))
        a_alt_gt = torch.rand((batch_size, n, 7, 2))

        loss = torsion_angle_loss(a, a_gt, a_alt_gt)

    def test_run_fape(self):
        batch_size = 2
        n_frames = 7
        n_atoms = 5

        x = torch.rand((batch_size, n_atoms, 3))
        x_gt = torch.rand((batch_size, n_atoms, 3))
        rots = torch.rand((batch_size, n_frames, 3, 3))
        rots_gt = torch.rand((batch_size, n_frames, 3, 3))
        trans = torch.rand((batch_size, n_frames, 3))
        trans_gt = torch.rand((batch_size, n_frames, 3))
        t = T(rots, trans)
        t_gt = T(rots_gt, trans_gt)

        loss = compute_fape(t, x, t_gt, x_gt)

    def test_between_residue_bond_loss(self):
        bs = 2
        n = 10
        pred_pos = torch.rand(bs, n, 14, 3)
        pred_atom_mask = torch.randint(0, 2, (bs, n, 14))
        residue_index = torch.arange(n).unsqueeze(0)
        aatype = torch.randint(0, 22, (bs, n,))
        
        between_residue_bond_loss(
            pred_pos,
            pred_atom_mask,
            residue_index, 
            aatype,
        )

    def test_between_residue_clash_loss(self):
        bs = 2
        n = 10
        pred_pos = torch.rand(bs, n, 14, 3)
        pred_atom_mask = torch.randint(0, 2, (bs, n, 14))
        atom14_atom_radius = torch.rand(bs, n, 14)
        residue_index = torch.arange(n).unsqueeze(0)

        loss = between_residue_clash_loss(
            pred_pos,
            pred_atom_mask,
            atom14_atom_radius, 
            residue_index,
        )

    def test_find_structural_violations(self):
        n = 10

        batch = {
            "atom14_atom_exists": torch.randint(0, 2, (n, 14)),
            "residue_index": torch.arange(n),
            "aatype": torch.randint(0, 21, (n,)),
            "residx_atom14_to_atom37": torch.randint(0, 37, (n, 14)).long(),
        }

        pred_pos = torch.rand(n, 14, 3)
        
        config = ml_collections.ConfigDict({
            "clash_overlap_tolerance": 1.5,
            "violation_tolerance_factor": 12.0,
        })

        find_structural_violations(batch, pred_pos, config)


if __name__ == "__main__":
    unittest.main()

