# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
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

import ml_collections
import numpy as np
import torch
import torch.nn as nn
from typing import Dict, Optional

from alphafold.np import residue_constants
from alphafold.model.primitives import Linear
from alphafold.utils.affine_utils import T
from alphafold.utils.tensor_utils import (
    tree_map, 
    tensor_tree_map, 
    masked_mean,
)


def softmax_cross_entropy(logits, labels):
    loss = -1 * torch.sum(
        labels * torch.nn.functional.log_softmax(logits),
        dim=-1,
    )
    return loss


def torsion_angle_loss(
    a,           # [*, N, 7, 2]
    a_gt,        # [*, N, 7, 2]
    a_alt_gt,    # [*, N, 7, 2]
):
    # [*, N, 7]
    norm = torch.norm(a, dim=-1)

    # [*, N, 7, 2]
    a = a / norm.unsqueeze(-1)

    # [*, N, 7]
    diff_norm_gt = torch.norm(a - a_gt, dim=-1)
    diff_norm_alt_gt = torch.norm(a - a_alt_gt, dim=-1)
    min_diff = torch.minimum(diff_norm_gt ** 2, diff_norm_alt_gt ** 2)

    # [*]
    l_torsion = torch.mean(min_diff, dim=(-1, -2))
    l_angle_norm = torch.mean(torch.abs(norm - 1), dim=(-1, -2))

    an_weight = 0.02
    return l_torsion + an_weight * l_angle_norm


def compute_fape(
    pred_frames: T,
    target_frames: T,
    frames_mask: torch.Tensor,
    pred_positions: torch.Tensor,
    target_positions: torch.Tensor,
    positions_mask: torch.Tensor,
    length_scale: float,
    l1_clamp_distance: Optional[float] = None,
    eps=1e-4
) -> torch.Tensor:
    # [*, N_frames, N_pts, 3]
    local_pred_pos = pred_frames.invert()[..., None].apply(
        pred_positions[..., None, :, :],
    )
    local_target_pos = target_frames.invert()[..., None].apply(
        target_positions[..., None, :, :],
    )
    error_dist = torch.sqrt(
        (pred_positions - target_positions)**2 + eps
    )

    if(l1_clamp_distance is not None):
        error_dist = torch.clamp(error_dist, min=0, max=l1_clamp_distance)

    normed_error = error_dist / length_scale
    normed_error *= frames_mask.unsqueeze(-1)
    normed_error *= positions_mask.unsqueeze(-2)

    norm_factor = (
        torch.sum(frames_mask, dim=-1) *
        torch.sum(positions_mask, dim=-1)
    )

    normed_error = torch.sum(normed_error, dim=(-1, -2)) / (eps + norm_factor)

    return normed_error


def backbone_loss(
    batch: Dict[str, torch.Tensor],
    pred_aff: T,
    clamp_distance: float = 10.,
    loss_unit_distance: float = 10.,
) -> torch.Tensor:
    gt_aff = T.from_tensor(batch['backbone_affine_tensor'])
    backbone_mask = batch['backbone_affine_mask']
    
    fape_loss = compute_fape(
        pred_aff,
        gt_aff,
        backbone_mask,
        pred_aff.get_trans(),
        gt_aff.get_trans(),
        backbone_mask,
        l1_clamp_distance=clamp_distance,
        length_scale=loss_unit_distance,
    )

    if('use_clamped_fape' in batch):
        use_clamped_fape = batch["use_clamped_fape"]
        unclamped_fape_loss = compute_fape(
            pred_aff,
            gt_aff,
            backbone_mask,
            pred_aff.get_trans(),
            gt_aff.get_trans(),
            backbone_mask,
            l1_clamp_distance=None,
            length_scale=loss_unit_distance,
        )

        fape_loss = (
            fape_loss * use_clamped_fape +
            fape_loss_unclamped * (1 - use_clamped_fape)
        )

    return torch.mean(fape_loss, dim=backbone_mask.shape[:-1])


def sidechain_loss(
    sidechain_frames,
    sidechain_atom_pos,
    gt_frames,
    alt_gt_frames,
    gt_exists,
    renamed_atom14_gt_positions,
    renamed_atom14_gt_exists,
    alt_naming_is_better,
    clamp_distance=10.,
    length_scale=10.,
):
    renamed_gt_frames = (
        (1. - alt_naming_is_better[..., None, None, None, None]) *
        gt_frames +
        alt_naming_is_better[..., None, None, None, None] *
        alt_gt_frames
    )

    renamed_gt_frames = T.from_4x4(renamed_gt_frames) 

    fape = compute_fape(
        sidechain_frames,
        renamed_gt_frames,
        gt_exists,
        sidechain_atom_pos,
        renamed_atom14_gt_positions,
        renamed_atom14_gt_exists,
        l1_clamp_distance=clamp_distance,
        length_scale=length_scale,
    )

    return fape
    

def compute_plddt(logits: torch.Tensor) -> torch.Tensor:
    num_bins = logits.shape[-1]
    bin_width = 1. / num_bins
    bounds = torch.arange(
        start=0.5 * bin_width, end=1.0, step=bin_width, device=logits.device
    )
    probs = torch.nn.functional.softmax(logits, dim=-1)
    pred_lddt_ca = torch.sum(
        probs * 
        bounds.view(*((1,) * len(probs.shape[:-1])), *bounds.shape),
        dim=-1,
    )
    return pred_lddt_ca * 100


def lddt_loss(
    batch: Dict[str, torch.Tensor],
    cutoff: float = 15.,
    num_bins: int = 50,
    min_resolution: float = 0.1,
    max_resolution: float = 3.0,
    eps: float = 1e-10,
) -> torch.Tensor:
    all_atom_pred_pos = batch["sm"]["pred_pos"][-1]
    all_atom_true_pos = batch["all_atom_positions"]
    all_atom_mask = batch["all_atom_mask"]
    logits = batch["predicted_lddt_logits"]

    n = all_atom_mask.shape[-1]
   
    ca_pos = residue_constants.atom_order['CA']
    all_atom_pred_pos = all_atom_pred_pos[..., :, ca_pos, :]
    all_atom_true_pos = all_atom_true_pos[..., :, ca_pos, :]
    all_atom_mask = all_atom_mask[..., :, ca_pos:(ca_pos + 1)] # keep dim

    dmat_true = torch.sqrt(
        eps +
        torch.sum(
            (
                all_atom_true_pos[..., None] - 
                all_atom_true_pos[..., None, :]
            )**2,
            dim=-1,
        )
    )

    dmat_pred = torch.sqrt(
        eps +
        torch.sum(
            (
                all_atom_pred_pos[..., None] - 
                all_atom_pred_pos[..., None, :]

            )**2,
            dim=-1,
        )
    )

    dists_to_score = (
        (dmat_true < cutoff) * all_atom_mask *
        permute_final_dims(all_atom_mask, 1, 0) *
        (1. - torch.eye(n, device=all_atom_mask.device))
    )

    dist_l1 = torch.abs(dmat_true - dmat_pred)

    score = (
        (dist_l1 < 0.5) + 
        (dist_l1 < 1.0) +
        (dist_l1 < 2.0) +
        (dist_l1 < 4.0)
    )
    score *= 0.25

    norm = 1. / (eps + torch.sum(dists_to_score, dim=-1))
    score = norm * (eps + torch.sum(dists_to_score * score, dim=-1))

    # TODO: this feels a bit weird, but it's in the source
    score = score.detach() 

    bin_index = torch.floor(lddt_ca * num_bins)

    lddt_ca_one_hot = torch.nn.functional.one_hot(
        bin_index, num_classes=num_bins
    )

    errors = softmax_cross_entropy(logits, lddt_ca_one_hot)

    loss = torch.sum(errors * all_atom_mask) / (torch.sum(mask_ca) + eps)

    loss *= (
        (batch["resolution"] >= min_resolution) &
        (batch["resolution"] <= max_resolution)
    )

    return loss


def distogram_loss(
    pred_distr, 
    gt, 
    mask, 
    min_bin=2.3125, max_bin=21.6875, no_bins=64, eps=1e-6
):
    boundaries = torch.linspace(
        min_bin, max_bin, no_bins - 1, device=pred_distr.device,
    )
    boundaries = boundaries ** 2

    dists = torch.sum(
        (gt[..., None, :] - gt[..., None, :, :]) ** 2, dim=-1, keepdims=True
    )

    true_bins = torch.sum(dists > sq_breaks, dim=-1)

    errors = softmax_cross_entropy(
        pred_distr,
        torch.nn.functional.one_hot(true_bins, num_bins),
    )

    square_mask = mask[..., None] * mask[..., None, :]

    mean = (
        torch.sum(errors * square_mask, dim=(-1, -2)) /
        (eps + torch.sum(square_mask, dim=(-1, -2)))
    )

    return mean


def tm_score(
    logits,
    t_pred, 
    t_gt, 
    mask, 
    resolution,
    max_bin=31, 
    no_bins=64, 
    min_resolution: float = 0.1,
    max_resolution: float = 3.0,
    eps=1e-8
):
    boundaries = torch.linspace(
        min=0, 
        max=max_bin, 
        steps=(no_bins - 1), 
        device=logits.device
    )
    boundaries = boundaries ** 2

    def _points(affine):
        pts = affine.trans.unsqueeze(-3)
        return affine.invert().apply(pts, addl_dims=1)

    sq_diff = torch.sum((_points(t_pred) - _points(t_gt)) ** 2, dim=-1)
    sq_diff = sq_diff.detach()

    true_bins = torch.sum(
        sq_diff[..., None] > boundaries
    ).float()

    errors = softmax_cross_entropy(
        logits,
        torch.nn.functional.one_hot(true_bins, no_bins)
    )

    square_mask = mask[..., None] * mask[..., None, :]

    loss = (
        torch.sum(loss, dim=(-1, -2)) /
        (eps + torch.sum(square_mask, dim=(-1, -2)))
    )

    loss *= (
        (resolution >= min_resolution) &
        (resolution <= max_resolution)
    )

    return loss


def between_residue_bond_loss(
    pred_atom_positions: torch.Tensor,  # (N, 37(14), 3)
    pred_atom_mask: torch.Tensor,  # (N, 37(14))
    residue_index: torch.Tensor,  # (N)
    aatype: torch.Tensor,  # (N)
    tolerance_factor_soft=12.0,
    tolerance_factor_hard=12.0,
    eps=1e-6,
) -> Dict[str, torch.Tensor]:
    """Flat-bottom loss to penalize structural violations between residues.
  
    This is a loss penalizing any violation of the geometry around the peptide
    bond between consecutive amino acids. This loss corresponds to
    Jumper et al. (2021) Suppl. Sec. 1.9.11, eq 44, 45.
  
    Args:
      pred_atom_positions: Atom positions in atom37/14 representation
      pred_atom_mask: Atom mask in atom37/14 representation
      residue_index: Residue index for given amino acid, this is assumed to be
        monotonically increasing.
      aatype: Amino acid type of given residue
      tolerance_factor_soft: soft tolerance factor measured in standard deviations
        of pdb distributions
      tolerance_factor_hard: hard tolerance factor measured in standard deviations
        of pdb distributions
  
    Returns:
      Dict containing:
        * 'c_n_loss_mean': Loss for peptide bond length violations
        * 'ca_c_n_loss_mean': Loss for violations of bond angle around C spanned
            by CA, C, N
        * 'c_n_ca_loss_mean': Loss for violations of bond angle around N spanned
            by C, N, CA
        * 'per_residue_loss_sum': sum of all losses for each residue
        * 'per_residue_violation_mask': mask denoting all residues with violation
            present.
    """
    # Get the positions of the relevant backbone atoms.
    this_ca_pos = pred_atom_positions[..., :-1, 1, :]
    this_ca_mask = pred_atom_mask[..., :-1, 1]
    this_c_pos = pred_atom_positions[..., :-1, 2, :]
    this_c_mask = pred_atom_mask[..., :-1, 2]
    next_n_pos = pred_atom_positions[..., 1:, 0, :]
    next_n_mask = pred_atom_mask[..., 1:, 0]
    next_ca_pos = pred_atom_positions[..., 1:, 1, :]
    next_ca_mask = pred_atom_mask[..., 1:, 1]
    has_no_gap_mask = (
        (residue_index[..., 1:] - residue_index[..., :-1]) == 1.0
    )
  
    # Compute loss for the C--N bond.
    c_n_bond_length = torch.sqrt(
        eps + 
        torch.sum(
            (this_c_pos - next_n_pos)**2, dim=-1
        )
    )
  
    # The C-N bond to proline has slightly different length because of the ring.
    next_is_proline = (
        aatype[..., 1:] == residue_constants.resname_to_idx['PRO']
    )
    gt_length = (
        (~next_is_proline) * residue_constants.between_res_bond_length_c_n[0]
        + next_is_proline * residue_constants.between_res_bond_length_c_n[1]
    )
    gt_stddev = (
        (~next_is_proline) *
        residue_constants.between_res_bond_length_stddev_c_n[0] +
        next_is_proline * 
        residue_constants.between_res_bond_length_stddev_c_n[1]
    )
    c_n_bond_length_error = torch.sqrt(
        eps + (c_n_bond_length - gt_length)**2
    )
    c_n_loss_per_residue = torch.nn.functional.relu(
        c_n_bond_length_error - tolerance_factor_soft * gt_stddev
    )
    mask = this_c_mask * next_n_mask * has_no_gap_mask
    c_n_loss = torch.sum(mask * c_n_loss_per_residue) / (torch.sum(mask) + eps)
    c_n_violation_mask = mask * (
        c_n_bond_length_error > (tolerance_factor_hard * gt_stddev)
    )
  
    # Compute loss for the angles.
    ca_c_bond_length = torch.sqrt(
        eps + torch.sum((this_ca_pos - this_c_pos)**2, dim=-1)
    )
    n_ca_bond_length = torch.sqrt(
        eps + torch.sum((next_n_pos - next_ca_pos)**2, dim=-1)
    )
  
    c_ca_unit_vec = (this_ca_pos - this_c_pos) / ca_c_bond_length[..., None]
    c_n_unit_vec = (next_n_pos - this_c_pos) / c_n_bond_length[..., None]
    n_ca_unit_vec = (next_ca_pos - next_n_pos) / n_ca_bond_length[..., None]
  
    ca_c_n_cos_angle = torch.sum(c_ca_unit_vec * c_n_unit_vec, dim=-1)
    gt_angle = residue_constants.between_res_cos_angles_ca_c_n[0]
    gt_stddev = residue_constants.between_res_bond_length_stddev_c_n[0]
    ca_c_n_cos_angle_error = torch.sqrt(
        eps + (ca_c_n_cos_angle - gt_angle)**2
    )
    ca_c_n_loss_per_residue = torch.nn.functional.relu(
        ca_c_n_cos_angle_error - tolerance_factor_soft * gt_stddev
    )
    mask = this_ca_mask * this_c_mask * next_n_mask * has_no_gap_mask
    ca_c_n_loss = (
        torch.sum(mask * ca_c_n_loss_per_residue) / (torch.sum(mask) + eps)
    )
    ca_c_n_violation_mask = mask * (ca_c_n_cos_angle_error >
                                    (tolerance_factor_hard * gt_stddev))
  
    c_n_ca_cos_angle = torch.sum((-c_n_unit_vec) * n_ca_unit_vec, dim=-1)
    gt_angle = residue_constants.between_res_cos_angles_c_n_ca[0]
    gt_stddev = residue_constants.between_res_cos_angles_c_n_ca[1]
    c_n_ca_cos_angle_error = torch.sqrt(
        eps + torch.square(c_n_ca_cos_angle - gt_angle))
    c_n_ca_loss_per_residue = torch.nn.functional.relu(
        c_n_ca_cos_angle_error - tolerance_factor_soft * gt_stddev
    )
    mask = this_c_mask * next_n_mask * next_ca_mask * has_no_gap_mask
    c_n_ca_loss = (
        torch.sum(mask * c_n_ca_loss_per_residue) / (torch.sum(mask) + eps)
    )
    c_n_ca_violation_mask = mask * (
        c_n_ca_cos_angle_error > (tolerance_factor_hard * gt_stddev)
    )
  
    # Compute a per residue loss (equally distribute the loss to both
    # neighbouring residues).
    per_residue_loss_sum = (c_n_loss_per_residue +
                            ca_c_n_loss_per_residue +
                            c_n_ca_loss_per_residue)
    per_residue_loss_sum = 0.5 * (
        torch.nn.functional.pad(per_residue_loss_sum, (0, 1)) +
        torch.nn.functional.pad(per_residue_loss_sum, (1, 0))
    )
 
    # Compute hard violations.
    violation_mask = torch.max(
        torch.stack(
            [
                c_n_violation_mask,
                ca_c_n_violation_mask,
                c_n_ca_violation_mask
            ]
        ), 
        dim=-2
    )[0]
    violation_mask = torch.maximum(
        torch.nn.functional.pad(violation_mask, (0, 1)),
        torch.nn.functional.pad(violation_mask, (1, 0))
    )
  
    return {
        'c_n_loss_mean': c_n_loss,
        'ca_c_n_loss_mean': ca_c_n_loss,
        'c_n_ca_loss_mean': c_n_ca_loss,
        'per_residue_loss_sum': per_residue_loss_sum,
        'per_residue_violation_mask': violation_mask
    }


def between_residue_clash_loss(
    atom14_pred_positions: torch.Tensor,
    atom14_atom_exists: torch.Tensor,
    atom14_atom_radius: torch.Tensor,
    residue_index: torch.Tensor,
    overlap_tolerance_soft=1.5,
    overlap_tolerance_hard=1.5,
    eps=1e-10,
) -> Dict[str, torch.Tensor]:
    """Loss to penalize steric clashes between residues.
  
    This is a loss penalizing any steric clashes due to non bonded atoms in
    different peptides coming too close. This loss corresponds to the part with
    different residues of
    Jumper et al. (2021) Suppl. Sec. 1.9.11, eq 46.
  
    Args:
      atom14_pred_positions: Predicted positions of atoms in
        global prediction frame
      atom14_atom_exists: Mask denoting whether atom at positions exists for given
        amino acid type
      atom14_atom_radius: Van der Waals radius for each atom.
      residue_index: Residue index for given amino acid.
      overlap_tolerance_soft: Soft tolerance factor.
      overlap_tolerance_hard: Hard tolerance factor.
  
    Returns:
      Dict containing:
        * 'mean_loss': average clash loss
        * 'per_atom_loss_sum': sum of all clash losses per atom, shape (N, 14)
        * 'per_atom_clash_mask': mask whether atom clashes with any other atom
            shape (N, 14)
    """
    fp_type = atom14_pred_positions.dtype
    
    # Create the distance matrix.
    # (N, N, 14, 14)
    dists = torch.sqrt(
        eps + 
        torch.sum(
            (
                atom14_pred_positions[..., :, None, :, None, :] -
                atom14_pred_positions[..., None, :, None, :, :]
            )**2,
            dim=-1)
    )
  
    # Create the mask for valid distances.
    # shape (N, N, 14, 14)
    dists_mask = (
        atom14_atom_exists[..., :, None, :, None] *
        atom14_atom_exists[..., None, :, None, :]
    ).type(fp_type)
  
    # Mask out all the duplicate entries in the lower triangular matrix.
    # Also mask out the diagonal (atom-pairs from the same residue) -- these atoms
    # are handled separately.
    dists_mask *= (
        residue_index[..., :, None, None, None] < residue_index[..., None, :, None, None]
    )
  
    # Backbone C--N bond between subsequent residues is no clash.
    c_one_hot = torch.nn.functional.one_hot(
        residue_index.new_tensor(2), num_classes=14
    )
    c_one_hot = c_one_hot.reshape(
        *((1,) * len(residue_index.shape[:-1])), *c_one_hot.shape
    )
    c_one_hot = c_one_hot.type(fp_type)
    n_one_hot = torch.nn.functional.one_hot(
        residue_index.new_tensor(0), num_classes=14
    )
    n_one_hot = n_one_hot.reshape(
        *((1,) * len(residue_index.shape[:-1])), *n_one_hot.shape
    )
    n_one_hot = n_one_hot.type(fp_type)

    neighbour_mask = (
        (residue_index[..., :, None, None, None] + 1) == 
        residue_index[..., None, :, None, None]
    )
    c_n_bonds = (
        neighbour_mask * 
        c_one_hot[..., None, None, :, None] * 
        n_one_hot[..., None, None, None, :]
    )
    dists_mask *= (1. - c_n_bonds)
  
    # Disulfide bridge between two cysteines is no clash.
    cys = residue_constants.restype_name_to_atom14_names['CYS']
    cys_sg_idx = cys.index('SG')
    cys_sg_idx = residue_index.new_tensor(cys_sg_idx)
    cys_sg_idx = cys_sg_idx.reshape(
        *((1,) * len(residue_index.shape[:-1])), 1 
    ).squeeze(-1)
    cys_sg_one_hot = torch.nn.functional.one_hot(
        cys_sg_idx, num_classes=14
    )
    disulfide_bonds = (
        cys_sg_one_hot[..., None, None, :, None] *
        cys_sg_one_hot[..., None, None, None, :])
    dists_mask *= (1. - disulfide_bonds)
  
    # Compute the lower bound for the allowed distances.
    # shape (N, N, 14, 14)
    dists_lower_bound = dists_mask * (
        atom14_atom_radius[..., :, None, :, None] +
        atom14_atom_radius[..., None, :, None, :]
    )
  
    # Compute the error.
    # shape (N, N, 14, 14)
    dists_to_low_error = dists_mask * torch.nn.functional.relu(
        dists_lower_bound - overlap_tolerance_soft - dists
    )
  
    # Compute the mean loss.
    # shape ()
    mean_loss = (
        torch.sum(dists_to_low_error) / (1e-6 + torch.sum(dists_mask))
    )
  
    # Compute the per atom loss sum.
    # shape (N, 14)
    per_atom_loss_sum = (
        torch.sum(dists_to_low_error, dim=(-4, -2)) +
        torch.sum(dists_to_low_error, axis=(-3, -1))
    )
  
    # Compute the hard clash mask.
    # shape (N, N, 14, 14)
    clash_mask = dists_mask * (
        dists < (dists_lower_bound - overlap_tolerance_hard)
    )
  
    # Compute the per atom clash.
    # shape (N, 14)
    per_atom_clash_mask = torch.maximum(
        torch.amax(clash_mask, axis=(-4, -2)),
        torch.amax(clash_mask, axis=(-3, -1)),
    )
  
    return {
        'mean_loss': mean_loss,  # shape ()
        'per_atom_loss_sum': per_atom_loss_sum,  # shape (N, 14)
        'per_atom_clash_mask': per_atom_clash_mask  # shape (N, 14)
    }


def within_residue_violations(
    atom14_pred_positions: torch.Tensor,
    atom14_atom_exists: torch.Tensor,
    atom14_dists_lower_bound: torch.Tensor,
    atom14_dists_upper_bound: torch.Tensor,
    tighten_bounds_for_loss=0.0,
    eps=1e-10,
) -> Dict[str, torch.Tensor]:
    """Loss to penalize steric clashes within residues.
  
    This is a loss penalizing any steric violations or clashes of non-bonded atoms
    in a given peptide. This loss corresponds to the part with
    the same residues of
    Jumper et al. (2021) Suppl. Sec. 1.9.11, eq 46.
  
    Args:
        atom14_pred_positions ([*, N, 14, 3]): 
            Predicted positions of atoms in global prediction frame.
        atom14_atom_exists ([*, N, 14]): 
            Mask denoting whether atom at positions exists for given
            amino acid type
        atom14_dists_lower_bound ([*, N, 14]): 
            Lower bound on allowed distances.
        atom14_dists_upper_bound ([*, N, 14]): 
            Upper bound on allowed distances
        tighten_bounds_for_loss ([*, N]): 
            Extra factor to tighten loss
  
    Returns:
      Dict containing:
        * 'per_atom_loss_sum' ([*, N, 14]): 
              sum of all clash losses per atom, shape
        * 'per_atom_clash_mask' ([*, N, 14]): 
              mask whether atom clashes with any other atom shape 
    """  
    # Compute the mask for each residue.
    dists_masks = (
        1. - torch.eye(14, device=atom14_atom_exists.device)[None]
    )
    dists_masks = dists_masks.reshape(
        *((1,) * len(atom14_atom_exists.shape[:-2])), *dists_masks.shape
    )
    dists_masks = (
        atom14_atom_exists[..., :, :, None] *
        atom14_atom_exists[..., :, None, :] *
        dists_masks
    )
  
    # Distance matrix
    dists = torch.sqrt(
        eps + 
        torch.sum(
            (
                atom14_pred_positions[..., :, :, None, :] -
                atom14_pred_positions[..., :, None, :, :]
            )**2,
            dim=-1
        )
    )

    # Compute the loss.
    dists_to_low_error = torch.nn.functional.relu(
        atom14_dists_lower_bound + tighten_bounds_for_loss - dists
    )
    dists_to_high_error = torch.nn.functional.relu(
        dists - (atom14_dists_upper_bound - tighten_bounds_for_loss)
    )
    loss = dists_masks * (dists_to_low_error + dists_to_high_error)
  
    # Compute the per atom loss sum.
    per_atom_loss_sum = (
        torch.sum(loss, dim=-2) +
        torch.sum(loss, dim=-1)
    )
  
    # Compute the violations mask.
    violations = (
        dists_masks * 
        (
            (dists < atom14_dists_lower_bound) |
            (dists > atom14_dists_upper_bound)
        )
    )
  
    # Compute the per atom violations.
    per_atom_violations = torch.maximum(
        torch.max(violations, dim=-2)[0], torch.max(violations, axis=-1)[0]
    )
  
    return {
        'per_atom_loss_sum': per_atom_loss_sum,
        'per_atom_violations': per_atom_violations
    }



def find_structural_violations(
    batch: Dict[str, torch.Tensor],
    atom14_pred_positions: torch.Tensor,
    config: ml_collections.ConfigDict
) -> Dict[str, torch.Tensor]:
    """Computes several checks for structural violations."""
  
    # Compute between residue backbone violations of bonds and angles.
    connection_violations = between_residue_bond_loss(
        pred_atom_positions=atom14_pred_positions,
        pred_atom_mask=batch['atom14_atom_exists'],
        residue_index=batch['residue_index'],
        aatype=batch['aatype'],
        tolerance_factor_soft=config.violation_tolerance_factor,
        tolerance_factor_hard=config.violation_tolerance_factor
    )
  
    # Compute the Van der Waals radius for every atom
    # (the first letter of the atom name is the element type).
    # Shape: (N, 14).
    atomtype_radius = [
        residue_constants.van_der_waals_radius[name[0]]
        for name in residue_constants.atom_types
    ]
    atomtype_radius = atom14_pred_positions.new_tensor(
        atomtype_radius 
    )
    atom14_atom_radius = (
        batch['atom14_atom_exists'] *
        atomtype_radius[batch['residx_atom14_to_atom37']]
    )
  
    # Compute the between residue clash loss.
    between_residue_clashes = between_residue_clash_loss(
        atom14_pred_positions=atom14_pred_positions,
        atom14_atom_exists=batch['atom14_atom_exists'],
        atom14_atom_radius=atom14_atom_radius,
        residue_index=batch['residue_index'],
        overlap_tolerance_soft=config.clash_overlap_tolerance,
        overlap_tolerance_hard=config.clash_overlap_tolerance
    )
  
    # Compute all within-residue violations (clashes,
    # bond length and angle violations).
    restype_atom14_bounds = residue_constants.make_atom14_dists_bounds(
        overlap_tolerance=config.clash_overlap_tolerance,
        bond_length_tolerance_factor=config.violation_tolerance_factor
    )
    atom14_dists_lower_bound = restype_atom14_bounds['lower_bound'][
        batch['aatype']
    ]
    atom14_dists_upper_bound = restype_atom14_bounds['upper_bound'][
        batch['aatype']
    ]
    atom14_dists_lower_bound = atom14_pred_positions.new_tensor(
        atom14_dists_lower_bound
    )
    atom14_dists_upper_bound = atom14_pred_positions.new_tensor(
        atom14_dists_upper_bound
    )
    residue_violations = within_residue_violations(
        atom14_pred_positions=atom14_pred_positions,
        atom14_atom_exists=batch['atom14_atom_exists'],
        atom14_dists_lower_bound=atom14_dists_lower_bound,
        atom14_dists_upper_bound=atom14_dists_upper_bound,
        tighten_bounds_for_loss=0.0
    )
  
    # Combine them to a single per-residue violation mask (used later for LDDT).
    per_residue_violations_mask = torch.max(
        torch.stack(
            [
                connection_violations['per_residue_violation_mask'],
                torch.max(
                    between_residue_clashes['per_atom_clash_mask'], dim=-1
                )[0],
                torch.max(
                    residue_violations['per_atom_violations'], dim=-1
                )[0],
            ], 
            dim=-1,
        ), 
        dim=-1,
    )[0]

    return {
        'between_residues': {
            'bonds_c_n_loss_mean':
                connection_violations['c_n_loss_mean'],  # ()
            'angles_ca_c_n_loss_mean':
                connection_violations['ca_c_n_loss_mean'],  # ()
            'angles_c_n_ca_loss_mean':
                connection_violations['c_n_ca_loss_mean'],  # ()
            'connections_per_residue_loss_sum':
                connection_violations['per_residue_loss_sum'],  # (N)
            'connections_per_residue_violation_mask':
                connection_violations['per_residue_violation_mask'],  # (N)
            'clashes_mean_loss':
                between_residue_clashes['mean_loss'],  # ()
            'clashes_per_atom_loss_sum':
                between_residue_clashes['per_atom_loss_sum'],  # (N, 14)
            'clashes_per_atom_clash_mask':
                between_residue_clashes['per_atom_clash_mask'],  # (N, 14)
        },
        'within_residues': {
            'per_atom_loss_sum':
                residue_violations['per_atom_loss_sum'],  # (N, 14)
            'per_atom_violations':
                residue_violations['per_atom_violations'],  # (N, 14),
        },
        'total_per_residue_violations_mask':
            per_residue_violations_mask,  # (N)
    }


def find_structural_violations_np(
    batch: Dict[str, np.ndarray],
    atom14_pred_positions: np.ndarray,
    config: ml_collections.ConfigDict
) -> Dict[str, np.ndarray]:
    to_tensor = lambda x: torch.tensor(x, requires_grad=False)
    batch = tree_map(to_tensor, batch, np.ndarray)
    atom14_pred_positions = to_tensor(atom14_pred_positions)

    out = find_structural_violations(batch, atom14_pred_positions, config)

    to_np = lambda x: np.array(x)
    np_out = tensor_tree_map(to_np, out)

    return np_out


def extreme_ca_ca_distance_violations(
      pred_atom_positions: torch.Tensor,  # (N, 37(14), 3)
      pred_atom_mask: torch.Tensor,  # (N, 37(14))
      residue_index: torch.Tensor,  # (N)
      max_angstrom_tolerance=1.5,
      eps=1e-6,
) -> torch.Tensor:
    """Counts residues whose Ca is a large distance from its neighbour.
  
    Measures the fraction of CA-CA pairs between consecutive amino acids that are
    more than 'max_angstrom_tolerance' apart.
  
    Args:
      pred_atom_positions: Atom positions in atom37/14 representation
      pred_atom_mask: Atom mask in atom37/14 representation
      residue_index: Residue index for given amino acid, this is assumed to be
        monotonically increasing.
      max_angstrom_tolerance: Maximum distance allowed to not count as violation.
    Returns:
      Fraction of consecutive CA-CA pairs with violation.
    """
    this_ca_pos = pred_atom_positions[..., :-1, 1, :]
    this_ca_mask = pred_atom_mask[..., :-1, 1]
    next_ca_pos = pred_atom_positions[..., 1:, 1, :]
    next_ca_mask = pred_atom_mask[..., 1:, 1]
    has_no_gap_mask = ((residue_index[..., 1:] - residue_index[..., :-1]) == 1.0)
    ca_ca_distance = torch.sqrt(
        eps + torch.sum((this_ca_pos - next_ca_pos)**2, dim=-1)
    )
    violations = (
        (ca_ca_distance - residue_constants.ca_ca) > max_angstrom_tolerance
    )
    mask = this_ca_mask * next_ca_mask * has_no_gap_mask
    mean = masked_mean(mask, violations, -1)
    return mean


def compute_violation_metrics(
    batch: Dict[str, torch.Tensor],
    atom14_pred_positions: torch.Tensor,  # (N, 14, 3)
    violations: Dict[str, torch.Tensor],
) -> Dict[str, torch.Tensor]:
    """Compute several metrics to assess the structural violations.""" 
    ret = {}
    extreme_ca_ca_violations = extreme_ca_ca_distance_violations(
        pred_atom_positions=atom14_pred_positions,
        pred_atom_mask=batch['atom14_atom_exists'],
        residue_index=batch['residue_index']
    )
    ret['violations_extreme_ca_ca_distance'] = extreme_ca_ca_violations
    ret['violations_between_residue_bond'] = masked_mean(
        batch['seq_mask'],
        violations['between_residues'][
            'connections_per_residue_violation_mask'
        ],
        dim=-1,
    )
    ret['violations_between_residue_clash'] = masked_mean(
        mask=batch['seq_mask'],
        value=torch.max(
            violations['between_residues']['clashes_per_atom_clash_mask'],
            dim=-1
        )[0],
        dim=-1,
    )
    ret['violations_within_residue'] = masked_mean(
        mask=batch['seq_mask'],
        value=torch.max(
            violations['within_residues']['per_atom_violations'], dim=-1
        )[0],
        dim=-1,
    )
    ret['violations_per_residue'] = masked_mean(
        mask=batch['seq_mask'],
        value=violations['total_per_residue_violations_mask'],
        dim=-1,
    )
    return ret


def compute_violation_metrics_np(
    batch: Dict[str, np.ndarray],
    atom14_pred_positions: np.ndarray,
    violations: Dict[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    to_tensor = lambda x: torch.tensor(x, requires_grad=False)
    batch = tree_map(to_tensor, batch, np.ndarray)
    atom14_pred_positions = to_tensor(atom14_pred_positions)
    violations = tree_map(to_tensor, violations, np.ndarray)


    out = compute_violation_metrics(batch, atom14_pred_positions, violations)

    to_np = lambda x: np.array(x)
    return tree_map(to_np, out, torch.Tensor)


def compute_renamed_ground_truth(
    batch: Dict[str, torch.Tensor],
    atom14_pred_positions: torch.Tensor,
    eps=1e-10,
) -> Dict[str, torch.Tensor]:
    """
    Find optimal renaming of ground truth based on the predicted positions.
  
    Alg. 26 "renameSymmetricGroundTruthAtoms"
  
    This renamed ground truth is then used for all losses,
    such that each loss moves the atoms in the same direction.
  
    Args:
      batch: Dictionary containing:
        * atom14_gt_positions: Ground truth positions.
        * atom14_alt_gt_positions: Ground truth positions with renaming swaps.
        * atom14_atom_is_ambiguous: 1.0 for atoms that are affected by
            renaming swaps.
        * atom14_gt_exists: Mask for which atoms exist in ground truth.
        * atom14_alt_gt_exists: Mask for which atoms exist in ground truth
            after renaming.
        * atom14_atom_exists: Mask for whether each atom is part of the given
            amino acid type.
      atom14_pred_positions: Array of atom positions in global frame with shape
    Returns:
      Dictionary containing:
        alt_naming_is_better: Array with 1.0 where alternative swap is better.
        renamed_atom14_gt_positions: Array of optimal ground truth positions
          after renaming swaps are performed.
        renamed_atom14_gt_exists: Mask after renaming swap is performed.
    """

    pred_dists = torch.sqrt(
        eps + 
        torch.sum(
            (
                atom14_pred_positions[...,    None, :, None, :] -
                atom14_pred_positions[..., None, :, None, :, :]
            )**2,
            dim=-1,
        )
    )
    
    atom14_gt_positions = batch['atom14_gt_positions']
    gt_dists = torch.sqrt(
        eps + 
        torch.sum(
            (
                atom14_gt_positions[...,    None, :, None, :] -
                atom14_gt_positions[..., None, :, None, :, :]
            )**2,
            dim=-1,
        )
    )

    atom14_alt_gt_positions = batch['atom14_alt_gt_positions']
    alt_gt_dists = torch.sqrt(
        eps + 
        torch.sum(
            (
                atom14_alt_gt_positions[...,    None, :, None, :] -
                atom14_alt_gt_positions[..., None, :, None, :, :]
            )**2,
            dim=-1,
        )
    )

    lddt = torch.sqrt(eps + (pred_dists - gt_dists)**2)
    alt_lddt = torch.sqrt(eps + (pred_dists - alt_gt_dists)**2)

    atom14_gt_exists = batch['atom14_gt_exists']
    atom14_atom_is_ambiguous = batch['atom14_atom_is_ambiguous']
    mask = (
        atom14_gt_exists[..., None, :, None] *
        atom14_atom_is_ambiguous[..., None, :, None] *
        atom14_gt_exists[..., None, :, None, :] *
        (1. - atom14_atom_is_ambiguous[..., None, :, None, :])
    )

    per_res_lddt = torch.sum(mask * lddt, dim=(-1, -2, -3))
    alt_per_res_lddt = torch.sum(mask * alt_lddt, dim=(-1, -2, -3))

    fp_type = atom14_pred_positions.dtype
    alt_naming_is_better = (alt_per_res_lddt < per_res_lddt).type(fp_type)

    renamed_atom14_gt_positions = (
        (1. - alt_naming_is_better[..., None, None]) *
        atom14_gt_positions +
        alt_naming_is_better[..., None, None] *
        atom14_alt_gt_positions
    )

    renamed_atom14_gt_mask = (
        (1. - alt_naming_is_better[..., None]) * atom14_gt_exists +
        alt_naming_is_better[..., None] * batch['atom14_alt_gt_exists']
    )

    return {
        'alt_naming_is_better': alt_naming_is_better,
        'renamed_atom14_gt_positions': renamed_atom14_gt_positions,
        'renamed_atom14_gt_exists': renamed_atom14_gt_mask,
    }


def experimentally_resolved_loss(
    logits: torch.Tensor,
    atom37_atom_exists: torch.Tensor,
    all_atom_mask: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    errors = sigmoid_cross_entropy(logits, all_atom_mask)
    loss_num = torch.sum(errors * atom37_atom_exists, dim=(-1, -2))
    loss = loss_num / (eps + torch.sum(atom37_atom_exists, dim=(-1, -2)))
    return loss
