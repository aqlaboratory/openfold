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

import os
import torch
import numpy as np
from pathlib import Path
import unittest
import ml_collections as mlc

from openfold.data import data_transforms
from openfold.np import residue_constants
from openfold.utils.rigid_utils import (
    Rotation,
    Rigid,
)
from openfold.utils.loss import (
    torsion_angle_loss,
    compute_fape,
    between_residue_bond_loss,
    between_residue_clash_loss,
    find_structural_violations,
    compute_renamed_ground_truth,
    masked_msa_loss,
    distogram_loss,
    experimentally_resolved_loss,
    violation_loss,
    fape_loss,
    lddt_loss,
    supervised_chi_loss,
    backbone_loss,
    sidechain_loss,
    tm_loss,
    compute_plddt,
    compute_tm,
    chain_center_of_mass_loss
)
from openfold.utils.tensor_utils import (
    tree_map,
    tensor_tree_map,
    dict_multimap,
)
import tests.compare_utils as compare_utils
from tests.config import consts
from tests.data_utils import random_affines_vector, random_affines_4x4, random_asym_ids

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


def affine_vector_to_4x4(affine):
    r = Rigid.from_tensor_7(affine)
    return r.to_tensor_4x4()


def affine_vector_to_rigid(am_rigid, affine):
    rigid_flat = np.split(affine, 7, axis=-1)
    rigid_flat = [r.squeeze(-1) for r in rigid_flat]
    qw, qx, qy, qz = rigid_flat[:4]
    trans = rigid_flat[4:]
    rotation = am_rigid.Rot3Array.from_quaternion(qw, qx, qy, qz, normalize=True)
    translation = am_rigid.Vec3Array(*trans)
    return am_rigid.Rigid3Array(rotation, translation)


class TestLoss(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if compare_utils.alphafold_is_installed():
            if consts.is_multimer:
                cls.am_atom = alphafold.model.all_atom_multimer
                cls.am_fold = alphafold.model.folding_multimer
                cls.am_modules = alphafold.model.modules_multimer
                cls.am_rigid = alphafold.model.geometry
            else:
                cls.am_atom = alphafold.model.all_atom
                cls.am_fold = alphafold.model.folding
                cls.am_modules = alphafold.model.modules
                cls.am_rigid = alphafold.model.r3

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
        t = Rigid(Rotation(rot_mats=rots), trans)
        t_gt = Rigid(Rotation(rot_mats=rots_gt), trans_gt)
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
        aatype = torch.randint(
            0,
            22,
            (
                bs,
                n,
            ),
        )

        between_residue_bond_loss(
            pred_pos,
            pred_atom_mask,
            residue_index,
            aatype,
        )

    @compare_utils.skip_unless_alphafold_installed()
    def test_between_residue_bond_loss_compare(self):
        def run_brbl(pred_pos, pred_atom_mask, residue_index, aatype):
            if consts.is_multimer:
                pred_pos = self.am_rigid.Vec3Array.from_array(pred_pos)

            return self.am_atom.between_residue_bond_loss(
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
            {},
            None,
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

    def test_run_between_residue_clash_loss(self):
        bs = consts.batch_size
        n = consts.n_res

        pred_pos = torch.rand(bs, n, 14, 3)
        pred_atom_mask = torch.randint(0, 2, (bs, n, 14)).float()
        atom14_atom_radius = torch.rand(bs, n, 14)
        residue_index = torch.arange(n).unsqueeze(0)

        loss = between_residue_clash_loss(
            pred_pos,
            pred_atom_mask,
            atom14_atom_radius,
            residue_index,
        )

    @compare_utils.skip_unless_alphafold_installed()
    def test_between_residue_clash_loss_compare(self):
        def run_brcl(pred_pos, atom_exists, atom_radius, res_ind, asym_id):
            if consts.is_multimer:
                pred_pos = self.am_rigid.Vec3Array.from_array(pred_pos)
                return self.am_atom.between_residue_clash_loss(
                    pred_pos,
                    atom_exists,
                    atom_radius,
                    res_ind,
                    asym_id
                )

            return self.am_atom.between_residue_clash_loss(
                pred_pos,
                atom_exists,
                atom_radius,
                res_ind
            )

        f = hk.transform(run_brcl)

        n_res = consts.n_res

        pred_pos = np.random.rand(n_res, 14, 3).astype(np.float32)
        atom_exists = np.random.randint(0, 2, (n_res, 14)).astype(np.float32)
        res_ind = np.arange(
            n_res,
        )
        residx_atom14_to_atom37 = np.random.randint(0, 37, (n_res, 14)).astype(np.int64)

        atomtype_radius = [
            residue_constants.van_der_waals_radius[name[0]]
            for name in residue_constants.atom_types
        ]
        atomtype_radius = np.array(atomtype_radius).astype(np.float32)
        atom_radius = (
                atom_exists
                * atomtype_radius[residx_atom14_to_atom37]
        )

        asym_id = None
        if consts.is_multimer:
            asym_id = random_asym_ids(n_res)

        out_gt = f.apply(
            {},
            None,
            pred_pos,
            atom_exists,
            atom_radius,
            res_ind,
            asym_id
        )
        out_gt = jax.tree_map(lambda x: x.block_until_ready(), out_gt)
        out_gt = jax.tree_map(lambda x: torch.tensor(np.copy(x)), out_gt)

        out_repro = between_residue_clash_loss(
            torch.tensor(pred_pos).cuda(),
            torch.tensor(atom_exists).cuda(),
            torch.tensor(atom_radius).cuda(),
            torch.tensor(res_ind).cuda(),
            torch.tensor(asym_id).cuda() if asym_id is not None else None,
        )
        out_repro = tensor_tree_map(lambda x: x.cpu(), out_repro)

        for k in out_gt.keys():
            self.assertTrue(
                torch.max(torch.abs(out_gt[k] - out_repro[k])) < consts.eps
            )

    @compare_utils.skip_unless_alphafold_installed()
    def test_compute_plddt_compare(self):
        n_res = consts.n_res

        logits = np.random.rand(n_res, 50)

        out_gt = alphafold.common.confidence.compute_plddt(logits)
        out_gt = torch.tensor(out_gt)
        logits_t = torch.tensor(logits)
        out_repro = compute_plddt(logits_t)

        self.assertTrue(
            torch.max(torch.abs(out_gt - out_repro)) < consts.eps
        )

    @compare_utils.skip_unless_alphafold_installed()
    def test_compute_ptm_compare(self):
        n_res = consts.n_res
        max_bin = 31
        no_bins = 64

        logits = np.random.rand(n_res, n_res, no_bins)
        boundaries = np.linspace(0, max_bin, num=(no_bins - 1))

        ptm_gt = alphafold.common.confidence.predicted_tm_score(logits, boundaries)
        ptm_gt = torch.tensor(ptm_gt)
        logits_t = torch.tensor(logits)
        ptm_repro = compute_tm(logits_t, no_bins=no_bins, max_bin=max_bin)

        self.assertTrue(
            torch.max(torch.abs(ptm_gt - ptm_repro)) < consts.eps
        )

        if consts.is_multimer:
            asym_id = random_asym_ids(n_res)
            iptm_gt = alphafold.common.confidence.predicted_tm_score(logits, boundaries,
                                                                     asym_id=asym_id, interface=True)
            iptm_gt = torch.tensor(iptm_gt)
            iptm_repro = compute_tm(logits_t, no_bins=no_bins, max_bin=max_bin,
                                    asym_id=torch.tensor(asym_id), interface=True)

            self.assertTrue(
                torch.max(torch.abs(iptm_gt - iptm_repro)) < consts.eps
            )

    def test_find_structural_violations(self):
        n = consts.n_res

        batch = {
            "atom14_atom_exists": torch.randint(0, 2, (n, 14)),
            "residue_index": torch.arange(n),
            "aatype": torch.randint(0, 20, (n,)),
            "residx_atom14_to_atom37": torch.randint(0, 37, (n, 14)).long(),
        }

        pred_pos = torch.rand(n, 14, 3)

        config = {
            "clash_overlap_tolerance": 1.5,
            "violation_tolerance_factor": 12.0,
        }

        find_structural_violations(batch, pred_pos, **config)

    @compare_utils.skip_unless_alphafold_installed()
    def test_find_structural_violations_compare(self):
        def run_fsv(batch, pos, config):
            cwd = os.getcwd()
            fpath = Path(__file__).parent.resolve() / "test_data"
            os.chdir(str(fpath))

            if consts.is_multimer:
                atom14_pred_pos = self.am_rigid.Vec3Array.from_array(pos)
                return self.am_fold.find_structural_violations(
                    batch['aatype'],
                    batch['residue_index'],
                    batch['atom14_atom_exists'],
                    atom14_pred_pos,
                    config,
                    batch['asym_id']
                )

            loss = self.am_fold.find_structural_violations(
                batch,
                pos,
                config,
            )
            os.chdir(cwd)
            return loss

        f = hk.transform(run_fsv)

        n_res = consts.n_res

        batch = {
            "atom14_atom_exists": np.random.randint(0, 2, (n_res, 14)),
            "residue_index": np.arange(n_res),
            "aatype": np.random.randint(0, 20, (n_res,)),
            "residx_atom14_to_atom37": np.random.randint(
                0, 37, (n_res, 14)
            ).astype(np.int64),
        }

        if consts.is_multimer:
            batch["asym_id"] = random_asym_ids(n_res)

        pred_pos = np.random.rand(n_res, 14, 3)

        config = mlc.ConfigDict(
            {
                "clash_overlap_tolerance": 1.5,
                "violation_tolerance_factor": 12.0,
            }
        )

        out_gt = f.apply({}, None, batch, pred_pos, config)
        out_gt = jax.tree_map(lambda x: x.block_until_ready(), out_gt)
        out_gt = jax.tree_map(lambda x: torch.tensor(np.copy(x)), out_gt)

        batch = tree_map(lambda x: torch.tensor(x).cuda(), batch, np.ndarray)
        out_repro = find_structural_violations(
            batch,
            torch.tensor(pred_pos).cuda(),
            **config,
        )
        out_repro = tensor_tree_map(lambda x: x.cpu(), out_repro)

        def compare(out):
            gt, repro = out
            assert torch.max(torch.abs(gt - repro)) < consts.eps

        dict_multimap(compare, [out_gt, out_repro])

    @compare_utils.skip_unless_alphafold_installed()
    def test_compute_renamed_ground_truth_compare(self):
        def run_crgt(batch, atom14_pred_pos):
            return alphafold.model.folding.compute_renamed_ground_truth(
                batch,
                atom14_pred_pos,
            )

        f = hk.transform(run_crgt)

        n_res = consts.n_res

        batch = {
            "seq_mask": np.random.randint(0, 2, (n_res,)).astype(np.float32),
            "aatype": np.random.randint(0, 20, (n_res,)),
            "atom14_gt_positions": np.random.rand(n_res, 14, 3),
            "atom14_gt_exists": np.random.randint(0, 2, (n_res, 14)).astype(
                np.float32
            ),
            "all_atom_mask": np.random.randint(0, 2, (n_res, 37)).astype(
                np.float32
            ),
            "all_atom_positions": np.random.rand(n_res, 37, 3).astype(
                np.float32
            ),
        }

        def _build_extra_feats_np():
            b = tree_map(lambda n: torch.tensor(n), batch, np.ndarray)
            b = data_transforms.make_atom14_masks(b)
            b = data_transforms.make_atom14_positions(b)
            return tensor_tree_map(lambda t: np.array(t), b)

        batch = _build_extra_feats_np()

        atom14_pred_pos = np.random.rand(n_res, 14, 3)

        out_gt = f.apply({}, None, batch, atom14_pred_pos)
        out_gt = jax.tree_map(lambda x: torch.tensor(np.array(x)), out_gt)

        batch = tree_map(lambda x: torch.tensor(x).cuda(), batch, np.ndarray)
        atom14_pred_pos = torch.tensor(atom14_pred_pos).cuda()

        out_repro = compute_renamed_ground_truth(batch, atom14_pred_pos)
        out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)

        for k in out_repro:
            self.assertTrue(
                torch.max(torch.abs(out_gt[k] - out_repro[k])) < consts.eps
            )

    @compare_utils.skip_unless_alphafold_installed()
    def test_msa_loss_compare(self):
        def run_msa_loss(value, batch):
            config = compare_utils.get_alphafold_config()
            msa_head = alphafold.model.modules.MaskedMsaHead(
                config.model.heads.masked_msa, config.model.global_config
            )
            return msa_head.loss(value, batch)

        f = hk.transform(run_msa_loss)

        n_res = consts.n_res
        n_seq = consts.n_seq

        value = {
            "logits": np.random.rand(n_res, n_seq, consts.msa_logits).astype(np.float32),
        }

        batch = {
            "true_msa": np.random.randint(0, 21, (n_res, n_seq)),
            "bert_mask": np.random.randint(0, 2, (n_res, n_seq)).astype(
                np.float32
            )
        }

        out_gt = f.apply({}, None, value, batch)["loss"]
        out_gt = torch.tensor(np.array(out_gt))

        value = tree_map(lambda x: torch.tensor(x).cuda(), value, np.ndarray)
        batch = tree_map(lambda x: torch.tensor(x).cuda(), batch, np.ndarray)

        with torch.no_grad():
            out_repro = masked_msa_loss(
                value["logits"],
                batch["true_msa"],
                batch["bert_mask"],
                consts.msa_logits
            )
        out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_distogram_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_distogram = config.model.heads.distogram

        def run_distogram_loss(value, batch):
            dist_head = alphafold.model.modules.DistogramHead(
                c_distogram, config.model.global_config
            )
            return dist_head.loss(value, batch)

        f = hk.transform(run_distogram_loss)

        n_res = consts.n_res

        value = {
            "logits": np.random.rand(n_res, n_res, c_distogram.num_bins).astype(
                np.float32
            ),
            "bin_edges": np.linspace(
                c_distogram.first_break,
                c_distogram.last_break,
                c_distogram.num_bins,
            ),
        }

        batch = {
            "pseudo_beta": np.random.rand(n_res, 3).astype(np.float32),
            "pseudo_beta_mask": np.random.randint(0, 2, (n_res,)),
        }

        out_gt = f.apply({}, None, value, batch)["loss"]
        out_gt = torch.tensor(np.array(out_gt))

        value = tree_map(lambda x: torch.tensor(x).cuda(), value, np.ndarray)

        batch = tree_map(lambda x: torch.tensor(x).cuda(), batch, np.ndarray)

        with torch.no_grad():
            out_repro = distogram_loss(
                logits=value["logits"],
                min_bin=c_distogram.first_break,
                max_bin=c_distogram.last_break,
                no_bins=c_distogram.num_bins,
                **batch,
            )

        out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_experimentally_resolved_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_experimentally_resolved = config.model.heads.experimentally_resolved

        def run_experimentally_resolved_loss(value, batch):
            er_head = alphafold.model.modules.ExperimentallyResolvedHead(
                c_experimentally_resolved, config.model.global_config
            )
            return er_head.loss(value, batch)

        f = hk.transform(run_experimentally_resolved_loss)

        n_res = consts.n_res

        value = {
            "logits": np.random.rand(n_res, 37).astype(np.float32),
        }

        batch = {
            "all_atom_mask": np.random.randint(0, 2, (n_res, 37)),
            "atom37_atom_exists": np.random.randint(0, 2, (n_res, 37)),
            "resolution": np.array(1.0),
        }

        out_gt = f.apply({}, None, value, batch)["loss"]
        out_gt = torch.tensor(np.array(out_gt))

        value = tree_map(lambda x: torch.tensor(x).cuda(), value, np.ndarray)

        batch = tree_map(lambda x: torch.tensor(x).cuda(), batch, np.ndarray)

        with torch.no_grad():
            out_repro = experimentally_resolved_loss(
                logits=value["logits"],
                min_resolution=c_experimentally_resolved.min_resolution,
                max_resolution=c_experimentally_resolved.max_resolution,
                **batch,
            )

        out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_supervised_chi_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_chi_loss = config.model.heads.structure_module

        def run_supervised_chi_loss(value, batch):
            if consts.is_multimer:
                pred_angles = np.reshape(
                    value['sidechains']['angles_sin_cos'], [-1, consts.n_res, 7, 2])

                unnormed_angles = np.reshape(
                    value['sidechains']['unnormalized_angles_sin_cos'], [-1, consts.n_res, 7, 2])

                chi_loss, _, _ = self.am_fold.supervised_chi_loss(
                    batch['seq_mask'],
                    batch['chi_mask'],
                    batch['aatype'],
                    batch['chi_angles'],
                    pred_angles,
                    unnormed_angles,
                    c_chi_loss
                )
                return chi_loss

            ret = {
                "loss": jax.numpy.array(0.0),
            }
            self.am_fold.supervised_chi_loss(
                ret, batch, value, c_chi_loss
            )
            return ret["loss"]

        f = hk.transform(run_supervised_chi_loss)

        n_res = consts.n_res

        value = {
            "sidechains": {
                "angles_sin_cos": np.random.rand(8, n_res, 7, 2).astype(
                    np.float32
                ),
                "unnormalized_angles_sin_cos": np.random.rand(
                    8, n_res, 7, 2
                ).astype(np.float32),
            }
        }

        batch = {
            "aatype": np.random.randint(0, 21, (n_res,)),
            "seq_mask": np.random.randint(0, 2, (n_res,)),
            "chi_mask": np.random.randint(0, 2, (n_res, 4)),
            "chi_angles": np.random.rand(n_res, 4).astype(np.float32),
        }

        out_gt = f.apply({}, None, value, batch)
        out_gt = torch.tensor(np.array(out_gt.block_until_ready()))
        value = tree_map(lambda x: torch.tensor(x).cuda(), value, np.ndarray)

        batch = tree_map(lambda x: torch.tensor(x).cuda(), batch, np.ndarray)

        batch["chi_angles_sin_cos"] = torch.stack(
            [
                torch.sin(batch["chi_angles"]),
                torch.cos(batch["chi_angles"]),
            ],
            dim=-1,
        )

        with torch.no_grad():
            out_repro = supervised_chi_loss(
                chi_weight=c_chi_loss.chi_weight,
                angle_norm_weight=c_chi_loss.angle_norm_weight,
                **{**batch, **value["sidechains"]},
            )

        out_repro = tensor_tree_map(lambda t: t.cpu(), out_repro)

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_violation_loss(self):
        config = compare_utils.get_alphafold_config()
        c_viol = config.model.heads.structure_module
        n_res = consts.n_res

        batch = {
            "seq_mask": np.random.randint(0, 2, (n_res,)).astype(np.float32),
            "residue_index": np.arange(n_res),
            "aatype": np.random.randint(0, 21, (n_res,)),
        }

        if consts.is_multimer:
            batch["asym_id"] = random_asym_ids(n_res)

        batch = tree_map(lambda n: torch.tensor(n).cuda(), batch, np.ndarray)

        atom14_pred_pos = np.random.rand(n_res, 14, 3).astype(np.float32)
        atom14_pred_pos = torch.tensor(atom14_pred_pos).cuda()

        batch = data_transforms.make_atom14_masks(batch)

        loss_sum_clash = violation_loss(
            find_structural_violations(batch, atom14_pred_pos, **c_viol),
            average_clashes=False, **batch
        )
        loss_sum_clash = loss_sum_clash.cpu()

        loss_avg_clash = violation_loss(
            find_structural_violations(batch, atom14_pred_pos, **c_viol),
            average_clashes=True, **batch
        )
        loss_avg_clash = loss_avg_clash.cpu()

    @compare_utils.skip_unless_alphafold_installed()
    def test_violation_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_viol = config.model.heads.structure_module

        def run_viol_loss(batch, atom14_pred_pos):
            ret = {
                "loss": np.array(0.0).astype(np.float32),
            }

            if consts.is_multimer:
                atom14_pred_pos = self.am_rigid.Vec3Array.from_array(atom14_pred_pos)
                viol = self.am_fold.find_structural_violations(
                    batch['aatype'],
                    batch['residue_index'],
                    batch['atom14_atom_exists'],
                    atom14_pred_pos,
                    c_viol,
                    batch['asym_id']
                )
                return self.am_fold.structural_violation_loss(mask=batch['atom14_atom_exists'],
                                                              violations=viol,
                                                              config=c_viol)

            value = {}
            value[
                "violations"
            ] = self.am_fold.find_structural_violations(
                batch,
                atom14_pred_pos,
                c_viol,
            )

            self.am_fold.structural_violation_loss(
                ret,
                batch,
                value,
                c_viol,
            )
            return ret["loss"]

        f = hk.transform(run_viol_loss)

        n_res = consts.n_res

        batch = {
            "seq_mask": np.random.randint(0, 2, (n_res,)).astype(np.float32),
            "residue_index": np.arange(n_res),
            "aatype": np.random.randint(0, 21, (n_res,))
        }

        if consts.is_multimer:
            batch["asym_id"] = random_asym_ids(n_res)

        atom14_pred_pos = np.random.rand(n_res, 14, 3).astype(np.float32)

        alphafold.model.tf.data_transforms.make_atom14_masks(batch)
        batch = {k: np.array(v) for k, v in batch.items()}

        out_gt = f.apply({}, None, batch, atom14_pred_pos)
        out_gt = torch.tensor(np.array(out_gt.block_until_ready()))

        batch = tree_map(lambda n: torch.tensor(n).cuda(), batch, np.ndarray)
        atom14_pred_pos = torch.tensor(atom14_pred_pos).cuda()

        batch = data_transforms.make_atom14_masks(batch)

        out_repro = violation_loss(
            find_structural_violations(batch, atom14_pred_pos, **c_viol),
            **batch,
        )
        out_repro = out_repro.cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_lddt_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_plddt = config.model.heads.predicted_lddt

        def run_plddt_loss(value, batch):
            head = alphafold.model.modules.PredictedLDDTHead(
                c_plddt, config.model.global_config
            )
            return head.loss(value, batch)

        f = hk.transform(run_plddt_loss)

        n_res = consts.n_res

        value = {
            "predicted_lddt": {
                "logits": np.random.rand(n_res, c_plddt.num_bins).astype(
                    np.float32
                ),
            },
            "structure_module": {
                "final_atom_positions": np.random.rand(n_res, 37, 3).astype(
                    np.float32
                ),
            },
        }

        batch = {
            "all_atom_positions": np.random.rand(n_res, 37, 3).astype(
                np.float32
            ),
            "all_atom_mask": np.random.randint(0, 2, (n_res, 37)).astype(
                np.float32
            ),
            "resolution": np.array(1.0).astype(np.float32),
        }

        out_gt = f.apply({}, None, value, batch)
        out_gt = torch.tensor(np.array(out_gt["loss"]))

        to_tensor = lambda t: torch.tensor(t).cuda()
        value = tree_map(to_tensor, value, np.ndarray)
        batch = tree_map(to_tensor, batch, np.ndarray)

        out_repro = lddt_loss(
            logits=value["predicted_lddt"]["logits"],
            all_atom_pred_pos=value["structure_module"]["final_atom_positions"],
            **{**batch, **c_plddt},
        )
        out_repro = out_repro.cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_backbone_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_sm = config.model.heads.structure_module

        def run_bb_loss(batch, value):
            if consts.is_multimer:
                intra_chain_mask = (batch["asym_id"][..., None] == batch["asym_id"][..., None, :]).astype(np.float32)
                gt_rigid = affine_vector_to_rigid(self.am_rigid, batch["backbone_affine_tensor"])
                target_rigid = affine_vector_to_rigid(self.am_rigid, value['traj'])
                intra_chain_bb_loss, intra_chain_fape = self.am_fold.backbone_loss(
                    gt_rigid=gt_rigid,
                    gt_frames_mask=batch["backbone_affine_mask"],
                    gt_positions_mask=batch["backbone_affine_mask"],
                    target_rigid=target_rigid,
                    config=c_sm.intra_chain_fape,
                    pair_mask=intra_chain_mask)
                interface_bb_loss, interface_fape = self.am_fold.backbone_loss(
                    gt_rigid=gt_rigid,
                    gt_frames_mask=batch["backbone_affine_mask"],
                    gt_positions_mask=batch["backbone_affine_mask"],
                    target_rigid=target_rigid,
                    config=c_sm.interface_fape,
                    pair_mask=1. - intra_chain_mask)

                return intra_chain_bb_loss + interface_bb_loss

            ret = {
                "loss": np.array(0.0),
            }
            self.am_fold.backbone_loss(ret, batch, value, c_sm)
            return ret["loss"]

        f = hk.transform(run_bb_loss)

        n_res = consts.n_res

        batch = {
            "backbone_affine_tensor": random_affines_vector((n_res,)),
            "backbone_affine_mask": np.random.randint(0, 2, (n_res,)).astype(
                np.float32
            ),
            "use_clamped_fape": np.array(0.0)
        }

        value = {
            "traj": random_affines_vector(
                (
                    c_sm.num_layer,
                    n_res,
                )
            ),
        }

        if consts.is_multimer:
            batch["asym_id"] = random_asym_ids(n_res)

        out_gt = f.apply({}, None, batch, value)
        out_gt = torch.tensor(np.array(out_gt.block_until_ready()))

        to_tensor = lambda t: torch.tensor(t).cuda()
        batch = tree_map(to_tensor, batch, np.ndarray)
        value = tree_map(to_tensor, value, np.ndarray)

        batch["backbone_rigid_tensor"] = affine_vector_to_4x4(
            batch["backbone_affine_tensor"]
        )
        batch["backbone_rigid_mask"] = batch["backbone_affine_mask"]
        
        if consts.is_multimer:
            intra_chain_mask = (batch["asym_id"][..., None]
                                == batch["asym_id"][..., None, :]).to(dtype=value["traj"].dtype)
            intra_chain_out = backbone_loss(traj=value["traj"], pair_mask=intra_chain_mask,
                                            **{**batch, **c_sm.intra_chain_fape})
            interface_out = backbone_loss(traj=value["traj"], pair_mask=1. - intra_chain_mask,
                                          **{**batch, **c_sm.interface_fape})
            out_repro = intra_chain_out + interface_out
            out_repro = out_repro.cpu()
        else:
            out_repro = backbone_loss(traj=value["traj"], **{**batch, **c_sm})
            out_repro = out_repro.cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_sidechain_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_sm = config.model.heads.structure_module

        def run_sidechain_loss(batch, value, atom14_pred_positions):
            if consts.is_multimer:
                atom14_pred_positions = self.am_rigid.Vec3Array.from_array(atom14_pred_positions)
                all_atom_positions = self.am_rigid.Vec3Array.from_array(batch["all_atom_positions"])
                gt_positions, gt_mask, alt_naming_is_better = self.am_fold.compute_atom14_gt(
                    aatype=batch["aatype"], all_atom_positions=all_atom_positions,
                    all_atom_mask=batch["all_atom_mask"], pred_pos=atom14_pred_positions)
                pred_frames = self.am_rigid.Rigid3Array.from_array4x4(value["sidechains"]["frames"])
                pred_positions = self.am_rigid.Vec3Array.from_array(value["sidechains"]["atom_pos"])
                gt_sc_frames, gt_sc_frames_mask = self.am_fold.compute_frames(
                    aatype=batch["aatype"],
                    all_atom_positions=all_atom_positions,
                    all_atom_mask=batch["all_atom_mask"],
                    use_alt=alt_naming_is_better)
                return self.am_fold.sidechain_loss(gt_sc_frames,
                                                   gt_sc_frames_mask,
                                                   gt_positions,
                                                   gt_mask,
                                                   pred_frames,
                                                   pred_positions,
                                                   c_sm)['loss']
            batch = {
                **batch,
                **self.am_atom.atom37_to_frames(
                    batch["aatype"],
                    batch["all_atom_positions"],
                    batch["all_atom_mask"],
                ),
            }
            v = {}
            v["sidechains"] = {}
            v["sidechains"][
                "frames"
            ] = self.am_rigid.rigids_from_tensor4x4(
                value["sidechains"]["frames"]
            )
            v["sidechains"]["atom_pos"] = self.am_rigid.vecs_from_tensor(
                value["sidechains"]["atom_pos"]
            )
            v.update(
                self.am_fold.compute_renamed_ground_truth(
                    batch,
                    atom14_pred_positions,
                )
            )
            value = v

            ret = self.am_fold.sidechain_loss(batch, value, c_sm)
            return ret["loss"]

        f = hk.transform(run_sidechain_loss)

        n_res = consts.n_res

        batch = {
            "seq_mask": np.random.randint(0, 2, (n_res,)).astype(np.float32),
            "aatype": np.random.randint(0, 20, (n_res,)),
            "atom14_gt_positions": np.random.rand(n_res, 14, 3).astype(
                np.float32
            ),
            "atom14_gt_exists": np.random.randint(0, 2, (n_res, 14)).astype(
                np.float32
            ),
            "all_atom_positions": np.random.rand(n_res, 37, 3).astype(
                np.float32
            ),
            "all_atom_mask": np.random.randint(0, 2, (n_res, 37)).astype(
                np.float32
            ),
        }

        def _build_extra_feats_np():
            b = tree_map(lambda n: torch.tensor(n), batch, np.ndarray)
            b = data_transforms.make_atom14_masks(b)
            b = data_transforms.make_atom14_positions(b)
            return tensor_tree_map(lambda t: np.array(t), b)

        batch = _build_extra_feats_np()

        value = {
            "sidechains": {
                "frames": random_affines_4x4((c_sm.num_layer, n_res, 8)),
                "atom_pos": np.random.rand(c_sm.num_layer, n_res, 14, 3).astype(
                    np.float32
                ),
            }
        }

        atom14_pred_pos = np.random.rand(n_res, 14, 3).astype(np.float32)

        out_gt = f.apply({}, None, batch, value, atom14_pred_pos)
        out_gt = torch.tensor(np.array(out_gt.block_until_ready()))

        to_tensor = lambda t: torch.tensor(t).cuda()
        batch = tree_map(to_tensor, batch, np.ndarray)
        value = tree_map(to_tensor, value, np.ndarray)
        atom14_pred_pos = to_tensor(atom14_pred_pos)

        batch = data_transforms.atom37_to_frames(batch)
        batch.update(compute_renamed_ground_truth(batch, atom14_pred_pos))

        out_repro = sidechain_loss(
            sidechain_frames=value["sidechains"]["frames"],
            sidechain_atom_pos=value["sidechains"]["atom_pos"],
            **{**batch, **c_sm},
        )
        out_repro = out_repro.cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    @unittest.skipIf(not consts.is_multimer and "ptm" not in consts.model, "Not enabled for non-ptm models.")
    def test_tm_loss_compare(self):
        config = compare_utils.get_alphafold_config()
        c_tm = config.model.heads.predicted_aligned_error

        def run_tm_loss(representations, batch, value):
            head = alphafold.model.modules.PredictedAlignedErrorHead(
                c_tm, config.model.global_config
            )
            v = {}
            v.update(value)
            v["predicted_aligned_error"] = head(representations, batch, False)
            return head.loss(v, batch)["loss"]

        f = hk.transform(run_tm_loss)

        np.random.seed(42)

        n_res = consts.n_res

        representations = {
            "pair": np.random.rand(n_res, n_res, consts.c_z).astype(np.float32),
        }

        batch = {
            "backbone_affine_tensor": random_affines_vector((n_res,)),
            "backbone_affine_mask": np.random.randint(0, 2, (n_res,)).astype(
                np.float32
            ),
            "resolution": np.array(1.0).astype(np.float32),
        }

        value = {
            "structure_module": {
                "final_affines": random_affines_vector((n_res,)),
            }
        }

        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/predicted_aligned_error_head"
        )

        out_gt = f.apply(params, None, representations, batch, value)
        out_gt = torch.tensor(np.array(out_gt.block_until_ready()))

        to_tensor = lambda n: torch.tensor(n).cuda()
        representations = tree_map(to_tensor, representations, np.ndarray)
        batch = tree_map(to_tensor, batch, np.ndarray)
        value = tree_map(to_tensor, value, np.ndarray)

        batch["backbone_rigid_tensor"] = affine_vector_to_4x4(
            batch["backbone_affine_tensor"]
        )
        batch["backbone_rigid_mask"] = batch["backbone_affine_mask"]

        model = compare_utils.get_global_pretrained_openfold()
        logits = model.aux_heads.tm(representations["pair"])

        out_repro = tm_loss(
            logits=logits,
            final_affine_tensor=value["structure_module"]["final_affines"],
            **{**batch, **c_tm},
        )
        out_repro = out_repro.cpu()

        self.assertTrue(torch.max(torch.abs(out_gt - out_repro)) < consts.eps)

    @compare_utils.skip_unless_alphafold_installed()
    def test_chain_center_of_mass_loss(self):
        batch_size = consts.batch_size
        n_res = consts.n_res

        batch = {
            "all_atom_positions": np.random.rand(batch_size, n_res, 37, 3).astype(np.float32) * 10.0,
            "all_atom_mask": np.random.randint(0, 2, (batch_size, n_res, 37)).astype(np.float32),
            "asym_id": np.stack([random_asym_ids(n_res) for _ in range(batch_size)])
        }

        config = {
            "weight": 0.05,
            "clamp_distance": -4.0,
        }

        final_atom_positions = torch.rand(batch_size, n_res, 37, 3).cuda()

        to_tensor = lambda t: torch.tensor(t).cuda()
        batch = tree_map(to_tensor, batch, np.ndarray)

        out_repro = chain_center_of_mass_loss(
            all_atom_pred_pos=final_atom_positions,
            **{**batch, **config},
        )
        out_repro = out_repro.cpu()


if __name__ == "__main__":
    unittest.main()
