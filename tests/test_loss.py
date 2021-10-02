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

from openfold.utils.loss import (
    torsion_angle_loss,
    compute_fape,
    between_residue_bond_loss,
    between_residue_clash_loss,
    find_structural_violations,
)
from openfold.utils.affine_utils import T
from openfold.utils.tensor_utils import tensor_tree_map
import tests.compare_utils as compare_utils
from tests.config import consts

if(compare_utils.alphafold_is_installed()):
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestLoss(unittest.TestCase):
    def test_run_torsion_angle_loss(self):
        batch_size = consts.batch_size
        n_res = consts.n_res

        a = torch.rand((batch_size, n_res, 7, 2))
        a_gt = torch.rand((batch_size, n_res, 7, 2))
        a_alt_gt = torch.rand((batch_size, n_res, 7, 2))

        loss = torsion_angle_loss(a, a_gt, a_alt_gt)

    def test_run_fape(self):
        batch_size = consts.batch_size
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
        frames_mask = torch.randint(0, 2, (batch_size, n_frames)).float()
        positions_mask = torch.randint(0, 2, (batch_size, n_atoms)).float()
        length_scale = 10

        loss = compute_fape(
            pred_frames=t,
            target_frames=t_gt,
            frames_mask=frames_mask,
            pred_positions=x,
            target_positions=x_gt,
            positions_mask=positions_mask,
            length_scale=length_scale,
        )

    def test_run_between_residue_bond_loss(self):
        bs = consts.batch_size
        n = consts.n_res
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

    @compare_utils.skip_unless_alphafold_installed()
    def test_between_residue_bond_loss_compare(self):
        def run_brbl(pred_pos, pred_atom_mask, residue_index, aatype):
            return alphafold.model.all_atom.between_residue_bond_loss(
                pred_pos,
                pred_atom_mask,
                residue_index,
                aatype,
            )
    
        f = hk.transform(run_brbl)
    
        n_res = consts.n_res 
        pred_pos = np.random.rand(n_res, 14, 3).astype(np.float32)
        pred_atom_mask = np.random.randint(0, 2, (n_res, 14)).astype(np.float32)
        residue_index = np.arange(n_res)
        aatype = np.random.randint(0, 22, (n_res,))
        
        out_gt = f.apply(
            {}, None, 
            pred_pos, 
            pred_atom_mask, 
            residue_index,
            aatype,
        )
        out_gt = jax.tree_map(lambda x: x.block_until_ready(), out_gt)
        out_gt = jax.tree_map(lambda x: torch.tensor(np.copy(x)), out_gt)
    
        out_repro = between_residue_bond_loss(
            torch.tensor(pred_pos).cuda(),
            torch.tensor(pred_atom_mask).cuda(),
            torch.tensor(residue_index).cuda(),
            torch.tensor(aatype).cuda(),
        )
        out_repro = tensor_tree_map(lambda x: x.cpu(), out_repro)
    
        for k in out_gt.keys():
            self.assertTrue(
                torch.max(torch.abs(out_gt[k] - out_repro[k])) < consts.eps
            )


    def test_between_residue_clash_loss(self):
        bs = consts.batch_size
        n = consts.n_res

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
        n = consts.n_res

        batch = {
            "atom14_atom_exists": torch.randint(0, 2, (n, 14)),
            "residue_index": torch.arange(n),
            "aatype": torch.randint(0, 21, (n,)),
            "residx_atom14_to_atom37": torch.randint(0, 37, (n, 14)).long(),
        }

        pred_pos = torch.rand(n, 14, 3)
        
        config = {
            "clash_overlap_tolerance": 1.5,
            "violation_tolerance_factor": 12.0,
        }

        find_structural_violations(batch, pred_pos, **config)


if __name__ == "__main__":
    unittest.main()

