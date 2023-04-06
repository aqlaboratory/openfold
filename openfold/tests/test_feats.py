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

import torch
import numpy as np
import unittest

import openfold.data.data_transforms as data_transforms
from openfold.np.residue_constants import (
    restype_rigid_group_default_frame,
    restype_atom14_to_rigid_group,
    restype_atom14_mask,
    restype_atom14_rigid_group_positions,
)
import openfold.utils.feats as feats
from openfold.utils.rigid_utils import Rotation, Rigid
from openfold.utils.tensor_utils import (
    tree_map,
    tensor_tree_map,
)
import tests.compare_utils as compare_utils
from tests.config import consts
from tests.data_utils import random_affines_4x4

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestFeats(unittest.TestCase):
    @compare_utils.skip_unless_alphafold_installed()
    def test_pseudo_beta_fn_compare(self):
        def test_pbf(aatype, all_atom_pos, all_atom_mask):
            return alphafold.model.modules.pseudo_beta_fn(
                aatype,
                all_atom_pos,
                all_atom_mask,
            )

        f = hk.transform(test_pbf)

        n_res = consts.n_res

        aatype = np.random.randint(0, 22, (n_res,))
        all_atom_pos = np.random.rand(n_res, 37, 3).astype(np.float32)
        all_atom_mask = np.random.randint(0, 2, (n_res, 37))

        out_gt_pos, out_gt_mask = f.apply(
            {}, None, aatype, all_atom_pos, all_atom_mask
        )
        out_gt_pos = torch.tensor(np.array(out_gt_pos.block_until_ready()))
        out_gt_mask = torch.tensor(np.array(out_gt_mask.block_until_ready()))

        out_repro_pos, out_repro_mask = feats.pseudo_beta_fn(
            torch.tensor(aatype).cuda(),
            torch.tensor(all_atom_pos).cuda(),
            torch.tensor(all_atom_mask).cuda(),
        )
        out_repro_pos = out_repro_pos.cpu()
        out_repro_mask = out_repro_mask.cpu()

        self.assertTrue(
            torch.max(torch.abs(out_gt_pos - out_repro_pos)) < consts.eps
        )
        self.assertTrue(
            torch.max(torch.abs(out_gt_mask - out_repro_mask)) < consts.eps
        )

    @compare_utils.skip_unless_alphafold_installed()
    def test_atom37_to_torsion_angles_compare(self):
        def run_test(aatype, all_atom_pos, all_atom_mask):
            return alphafold.model.all_atom.atom37_to_torsion_angles(
                aatype,
                all_atom_pos,
                all_atom_mask,
                placeholder_for_undefined=False,
            )

        f = hk.transform(run_test)

        n_templ = 7
        n_res = 13

        aatype = np.random.randint(0, 22, (n_templ, n_res)).astype(np.int64)
        all_atom_pos = np.random.rand(n_templ, n_res, 37, 3).astype(np.float32)
        all_atom_mask = np.random.randint(0, 2, (n_templ, n_res, 37)).astype(
            np.float32
        )

        out_gt = f.apply({}, None, aatype, all_atom_pos, all_atom_mask)
        out_gt = jax.tree_map(lambda x: torch.as_tensor(np.array(x)), out_gt)

        out_repro = data_transforms.atom37_to_torsion_angles()(
            {
                "aatype": torch.as_tensor(aatype).cuda(),
                "all_atom_positions": torch.as_tensor(all_atom_pos).cuda(),
                "all_atom_mask": torch.as_tensor(all_atom_mask).cuda(),
            },
        )
        tasc = out_repro["torsion_angles_sin_cos"].cpu()
        atasc = out_repro["alt_torsion_angles_sin_cos"].cpu()
        tam = out_repro["torsion_angles_mask"].cpu()

        # This function is extremely sensitive to floating point imprecisions,
        # so it is given much greater latitude in comparison tests.
        self.assertTrue(
            torch.mean(torch.abs(out_gt["torsion_angles_sin_cos"] - tasc))
            < 0.01
        )
        self.assertTrue(
            torch.mean(torch.abs(out_gt["alt_torsion_angles_sin_cos"] - atasc))
            < 0.01
        )
        self.assertTrue(
            torch.max(torch.abs(out_gt["torsion_angles_mask"] - tam))
            < consts.eps
        )

    @compare_utils.skip_unless_alphafold_installed()
    def test_atom37_to_frames_compare(self):
        def run_atom37_to_frames(aatype, all_atom_positions, all_atom_mask):
            return alphafold.model.all_atom.atom37_to_frames(
                aatype, all_atom_positions, all_atom_mask
            )

        f = hk.transform(run_atom37_to_frames)

        n_res = consts.n_res

        batch = {
            "aatype": np.random.randint(0, 21, (n_res,)),
            "all_atom_positions": np.random.rand(n_res, 37, 3).astype(
                np.float32
            ),
            "all_atom_mask": np.random.randint(0, 2, (n_res, 37)).astype(
                np.float32
            ),
        }

        out_gt = f.apply({}, None, **batch)
        to_tensor = lambda t: torch.tensor(np.array(t))
        out_gt = {k: to_tensor(v) for k, v in out_gt.items()}

        def flat12_to_4x4(flat12):
            rot = flat12[..., :9].view(*flat12.shape[:-1], 3, 3)
            trans = flat12[..., 9:]

            four_by_four = torch.zeros(*flat12.shape[:-1], 4, 4)
            four_by_four[..., :3, :3] = rot
            four_by_four[..., :3, 3] = trans
            four_by_four[..., 3, 3] = 1

            return four_by_four

        out_gt["rigidgroups_gt_frames"] = flat12_to_4x4(
            out_gt["rigidgroups_gt_frames"]
        )
        out_gt["rigidgroups_alt_gt_frames"] = flat12_to_4x4(
            out_gt["rigidgroups_alt_gt_frames"]
        )

        to_tensor = lambda t: torch.tensor(np.array(t)).cuda()
        batch = tree_map(to_tensor, batch, np.ndarray)

        out_repro = data_transforms.atom37_to_frames(batch)
        out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)

        for k, v in out_gt.items():
            self.assertTrue(
                torch.max(torch.abs(out_gt[k] - out_repro[k])) < consts.eps
            )

    def test_torsion_angles_to_frames_shape(self):
        batch_size = 2
        n = 5
        rots = torch.rand((batch_size, n, 3, 3))
        trans = torch.rand((batch_size, n, 3))
        ts = Rigid(Rotation(rot_mats=rots), trans)

        angles = torch.rand((batch_size, n, 7, 2))

        aas = torch.tensor([i % 2 for i in range(n)])
        aas = torch.stack([aas for _ in range(batch_size)])

        frames = feats.torsion_angles_to_frames(
            ts,
            angles,
            aas,
            torch.tensor(restype_rigid_group_default_frame),
        )

        self.assertTrue(frames.shape == (batch_size, n, 8))

    @compare_utils.skip_unless_alphafold_installed()
    def test_torsion_angles_to_frames_compare(self):
        def run_torsion_angles_to_frames(
            aatype, backb_to_global, torsion_angles_sin_cos
        ):
            return alphafold.model.all_atom.torsion_angles_to_frames(
                aatype,
                backb_to_global,
                torsion_angles_sin_cos,
            )

        f = hk.transform(run_torsion_angles_to_frames)

        n_res = consts.n_res

        aatype = np.random.randint(0, 21, size=(n_res,))

        affines = random_affines_4x4((n_res,))
        rigids = alphafold.model.r3.rigids_from_tensor4x4(affines)
        transformations = Rigid.from_tensor_4x4(
            torch.as_tensor(affines).float()
        )

        torsion_angles_sin_cos = np.random.rand(n_res, 7, 2)

        out_gt = f.apply({}, None, aatype, rigids, torsion_angles_sin_cos)

        jax.tree_map(lambda x: x.block_until_ready(), out_gt)

        out = feats.torsion_angles_to_frames(
            transformations.cuda(),
            torch.as_tensor(torsion_angles_sin_cos).cuda(),
            torch.as_tensor(aatype).cuda(),
            torch.tensor(restype_rigid_group_default_frame).cuda(),
        )

        # Convert the Rigids to 4x4 transformation tensors
        rots_gt = list(map(lambda x: torch.as_tensor(np.array(x)), out_gt.rot))
        trans_gt = list(
            map(lambda x: torch.as_tensor(np.array(x)), out_gt.trans)
        )
        rots_gt = torch.cat([x.unsqueeze(-1) for x in rots_gt], dim=-1)
        rots_gt = rots_gt.view(*rots_gt.shape[:-1], 3, 3)
        trans_gt = torch.cat([x.unsqueeze(-1) for x in trans_gt], dim=-1)
        transforms_gt = torch.cat([rots_gt, trans_gt.unsqueeze(-1)], dim=-1)
        bottom_row = torch.zeros((*rots_gt.shape[:-2], 1, 4))
        bottom_row[..., 3] = 1
        transforms_gt = torch.cat([transforms_gt, bottom_row], dim=-2)

        transforms_repro = out.to_tensor_4x4().cpu()

        self.assertTrue(
            torch.max(torch.abs(transforms_gt - transforms_repro) < consts.eps)
        )

    def test_frames_and_literature_positions_to_atom14_pos_shape(self):
        batch_size = consts.batch_size
        n_res = consts.n_res

        rots = torch.rand((batch_size, n_res, 8, 3, 3))
        trans = torch.rand((batch_size, n_res, 8, 3))
        ts = Rigid(Rotation(rot_mats=rots), trans)

        f = torch.randint(low=0, high=21, size=(batch_size, n_res)).long()

        xyz = feats.frames_and_literature_positions_to_atom14_pos(
            ts,
            f,
            torch.tensor(restype_rigid_group_default_frame),
            torch.tensor(restype_atom14_to_rigid_group),
            torch.tensor(restype_atom14_mask),
            torch.tensor(restype_atom14_rigid_group_positions),
        )

        self.assertTrue(xyz.shape == (batch_size, n_res, 14, 3))

    @compare_utils.skip_unless_alphafold_installed()
    def test_frames_and_literature_positions_to_atom14_pos_compare(self):
        def run_f(aatype, affines):
            am = alphafold.model
            return am.all_atom.frames_and_literature_positions_to_atom14_pos(
                aatype, affines
            )

        f = hk.transform(run_f)

        n_res = consts.n_res

        aatype = np.random.randint(0, 21, size=(n_res,))

        affines = random_affines_4x4((n_res, 8))
        rigids = alphafold.model.r3.rigids_from_tensor4x4(affines)
        transformations = Rigid.from_tensor_4x4(
            torch.as_tensor(affines).float()
        )

        out_gt = f.apply({}, None, aatype, rigids)
        jax.tree_map(lambda x: x.block_until_ready(), out_gt)
        out_gt = torch.stack(
            [torch.as_tensor(np.array(x)) for x in out_gt], dim=-1
        )

        out_repro = feats.frames_and_literature_positions_to_atom14_pos(
            transformations.cuda(),
            torch.as_tensor(aatype).cuda(),
            torch.tensor(restype_rigid_group_default_frame).cuda(),
            torch.tensor(restype_atom14_to_rigid_group).cuda(),
            torch.tensor(restype_atom14_mask).cuda(),
            torch.tensor(restype_atom14_rigid_group_positions).cuda(),
        ).cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro) < consts.eps))


if __name__ == "__main__":
    unittest.main()
