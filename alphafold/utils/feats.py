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

import numpy as np
import torch
import torch.nn as nn
from typing import Dict

import alphafold.np.residue_constants as residue_constants
from alphafold.utils.affine_utils import T
from alphafold.utils.tensor_utils import (
    batched_gather, 
    one_hot,
)


def pseudo_beta_fn(aatype, all_atom_positions, all_atom_masks):
    is_gly = (aatype == residue_constants.restype_order['G'])
    ca_idx = residue_constants.atom_order['CA']
    cb_idx = residue_constants.atom_order['CB']
    pseudo_beta = torch.where(
        is_gly[..., None].expand(*((-1,) * len(is_gly.shape)), 3),
        all_atom_positions[..., ca_idx, :],
        all_atom_positions[..., cb_idx, :]
    )

    if(all_atom_masks is not None):
        pseudo_beta_mask = torch.where(
            is_gly,
            all_atom_masks[..., ca_idx],
            all_atom_masks[..., cb_idx],
        )
        return pseudo_beta, pseudo_beta_mask
    else:
        return pseudo_beta


def get_chi_atom_indices():
  """Returns atom indices needed to compute chi angles for all residue types.

  Returns:
    A tensor of shape [residue_types=21, chis=4, atoms=4]. The residue types are
    in the order specified in residue_constants.restypes + unknown residue type
    at the end. For chi angles which are not defined on the residue, the
    positions indices are by default set to 0.
  """
  chi_atom_indices = []
  for residue_name in residue_constants.restypes:
    residue_name = residue_constants.restype_1to3[residue_name]
    residue_chi_angles = residue_constants.chi_angles_atoms[residue_name]
    atom_indices = []
    for chi_angle in residue_chi_angles:
      atom_indices.append(
          [residue_constants.atom_order[atom] for atom in chi_angle])
    for _ in range(4 - len(atom_indices)):
      atom_indices.append([0, 0, 0, 0])  # For chi angles not defined on the AA.
    chi_atom_indices.append(atom_indices)

  chi_atom_indices.append([[0, 0, 0, 0]] * 4)  # For UNKNOWN residue.

  return chi_atom_indices


def compute_residx(batch):
    aatype = batch["aatype"]

    restype_atom14_to_atom37 = []  # mapping (restype, atom37) --> atom14
    restype_atom37_to_atom14 = []  # mapping (restype, atom37) --> atom14
    restype_atom14_mask = []
    
    for rt in residue_constants.restypes:
        atom_names = residue_constants.restype_name_to_atom14_names[
          residue_constants.restype_1to3[rt]]
 
        restype_atom14_to_atom37.append([
            (residue_constants.atom_order[name] if name else 0)
            for name in atom_names
        ])

        atom_name_to_idx14 = {name: i for i, name in enumerate(atom_names)}
        restype_atom37_to_atom14.append([
            (atom_name_to_idx14[name] if name in atom_name_to_idx14 else 0)
            for name in residue_constants.atom_types
        ])

        restype_atom14_mask.append(
            [(1. if name else 0.) for name in atom_names]
        )
 
    # Add dummy mapping for restype 'UNK'
    restype_atom14_to_atom37.append([0] * 14)
    restype_atom37_to_atom14.append([0] * 37)
    restype_atom14_mask.append([0.] * 14)

    restype_atom14_to_atom37 = np.array(restype_atom14_to_atom37, dtype=np.int32)
    restype_atom37_to_atom14 = np.array(restype_atom37_to_atom14, dtype=np.int32)
    restype_atom14_mask = np.array(restype_atom14_mask, dtype=np.float32)
 
    residx_atom14_to_atom37 = np.take_along_axis(
        restype_atom14_to_atom37,
        aatype[..., None],
        axis=0
    )
    residx_atom14_mask = np.take_along_axis(
        restype_atom14_mask,
        aatype[..., None],
        axis=0,
    )

    batch['atom14_atom_exists'] = residx_atom14_mask
    batch['residx_atom14_to_atom37'] = residx_atom14_to_atom37.long()

    # create the gather indices for mapping back
    residx_atom37_to_atom14 = np.take_along_axis(
        restype_atom37_to_atom14,
        aatype[..., None],
        axis=0,
    )
    batch['residx_atom37_to_atom14'] = residx_atom37_to_atom14.long()

    # create the corresponding mask
    restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
    for restype, restype_letter in enumerate(residue_constants.restypes):
      restype_name = residue_constants.restype_1to3[restype_letter]
      atom_names = residue_constants.residue_atoms[restype_name]
      for atom_name in atom_names:
        atom_type = residue_constants.atom_order[atom_name]
        restype_atom37_mask[restype, atom_type] = 1

    residx_atom37_mask = np.take_along_axis(
        restype_atom37_mask,
        aatype[..., None],
        axis=0,
    )
    batch['atom37_atom_exists'] = residx_atom37_mask


def atom14_to_atom37(atom14, batch):
    atom37_data = batched_gather(
        atom14,
        batch["residx_atom37_to_atom14"],
        dim=-2,
        no_batch_dims=len(atom14.shape[:-2]),
    )

    atom37_data *= batch["atom37_atom_exists"][..., None]

    return atom37_data


def atom37_to_torsion_angles(
    aatype: torch.Tensor, 
    all_atom_pos: torch.Tensor, 
    all_atom_mask: torch.Tensor, 
    eps: float = 1e-8,
) -> Dict[str, torch.Tensor]:
    """
        Args:
            aatype:
                [*, N_res] residue indices
            all_atom_pos:
                [*, N_res, 37, 3] atom positions (in atom37 
                format)
            all_atom_mask:
                [*, N_res, 37] atom position mask
        Returns:
            Dictionary of the following features:
            
            "torsion_angles_sin_cos" ([*, N_res, 7, 2])
                Torsion angles
            "alt_torsion_angles_sin_cos" ([*, N_res, 7, 2])
                Alternate torsion angles (accounting for 180-degree symmetry)
            "torsion_angles_mask" ([*, N_res, 7])
                Torsion angles mask
    """
    aatype = torch.clamp(aatype, max=20)
    
    pad = all_atom_pos.new_zeros([*all_atom_pos.shape[:-3], 1, 37, 3])
    prev_all_atom_pos = torch.cat([pad, all_atom_pos[..., :-1, :, :]], dim=-3)

    pad = all_atom_mask.new_zeros([*all_atom_mask.shape[:-2], 1, 37])
    prev_all_atom_mask = torch.cat([pad, all_atom_mask[..., :-1, :]], dim=-2)

    pre_omega_atom_pos = torch.cat(
        [
            prev_all_atom_pos[..., 1:3, :],
            all_atom_pos[..., :2, :]
        ], dim=-2
    )
    phi_atom_pos = torch.cat(
        [
            prev_all_atom_pos[..., 2:3, :],
            all_atom_pos[..., :3, :]
        ], dim=-2
    )
    psi_atom_pos = torch.cat(
        [
            all_atom_pos[..., :3, :],
            all_atom_pos[..., 4:5, :]
        ], dim=-2
    )

    pre_omega_mask = (
        torch.prod(prev_all_atom_mask[..., 1:3], dim=-1) *
        torch.prod(all_atom_mask[..., :2], dim=-1)
    )
    phi_mask = (
        prev_all_atom_mask[..., 2] *
        torch.prod(all_atom_mask[..., :3], dim=-1)
    )
    psi_mask = (
        torch.prod(all_atom_mask[..., :3], dim=-1) *
        all_atom_mask[..., 4]
    )

    chi_atom_indices = torch.as_tensor(
        get_chi_atom_indices(), device=aatype.device
    )

    atom_indices = chi_atom_indices[..., aatype, :, :]
    chis_atom_pos = batched_gather(
        all_atom_pos, atom_indices, -2, len(atom_indices.shape[:-2])
    )

    chi_angles_mask = list(residue_constants.chi_angles_mask)
    chi_angles_mask.append([0., 0., 0., 0.])
    chi_angles_mask = all_atom_pos.new_tensor(chi_angles_mask)
    
    chis_mask = chi_angles_mask[aatype, :]

    chi_angle_atoms_mask = batched_gather(
        all_atom_mask, 
        atom_indices, 
        dim=-1, 
        no_batch_dims=len(atom_indices.shape[:-2])
    )
    chi_angle_atoms_mask = torch.prod(chi_angle_atoms_mask, dim=-1)
    chis_mask = chis_mask * chi_angle_atoms_mask

    torsions_atom_pos = torch.cat(
        [
            pre_omega_atom_pos[..., None, :, :],
            phi_atom_pos[..., None, :, :],
            psi_atom_pos[..., None, :, :],
            chis_atom_pos,
        ], dim=-3
    )

    torsion_angles_mask = torch.cat(
        [
            pre_omega_mask[..., None],
            phi_mask[..., None],
            psi_mask[..., None],
            chis_mask,
        ], dim=-1
    )

    torsion_frames = T.from_3_points(
        torsions_atom_pos[..., 1, :],
        torsions_atom_pos[..., 2, :],
        torsions_atom_pos[..., 0, :],
    )

    fourth_atom_rel_pos = torsion_frames.invert().apply(
        torsions_atom_pos[..., 3, :]
    )

    torsion_angles_sin_cos = torch.stack(
        [fourth_atom_rel_pos[..., 2], fourth_atom_rel_pos[..., 1]], dim=-1)
    denom = torch.sqrt(
        torch.sum(
            torch.square(torsion_angles_sin_cos), dim=-1, keepdims=True
        ) + eps
    )
    torsion_angles_sin_cos /= denom

    torsion_angles_sin_cos *= torch.tensor(
        [1., 1., -1., 1., 1., 1., 1.], device=aatype.device,
    )[((None,) * len(torsion_angles_sin_cos.shape[:-2])) + (slice(None), None)]

    chi_is_ambiguous = torsion_angles_sin_cos.new_tensor(
        residue_constants.chi_pi_periodic,
    )[aatype, ...]

    mirror_torsion_angles = torch.cat(
        [
            aatype.new_ones(*aatype.shape, 3),
            1. - 2. * chi_is_ambiguous
        ], dim=-1
    )

    alt_torsion_angles_sin_cos = (
        torsion_angles_sin_cos * mirror_torsion_angles[..., None]
    )

    return {
        "torsion_angles_sin_cos": torsion_angles_sin_cos,
        "alt_torsion_angles_sin_cos": alt_torsion_angles_sin_cos,
        "torsion_angles_mask": torsion_angles_mask,
    }


def atom37_to_frames(
    aatype: torch.Tensor,
    all_atom_positions: torch.Tensor,
    all_atom_mask: torch.Tensor,
) -> Dict[str, torch.Tensor]:
    batch_dims = len(aatype.shape[:-1])

    restype_rigidgroup_base_atom_names = np.full([21, 8, 3], '', dtype=object)
    restype_rigidgroup_base_atom_names[:, 0, :] = ['C', 'CA', 'N']
    restype_rigidgroup_base_atom_names[:, 3, :] = ['CA', 'C', 'O']
    
    for restype, restype_letter in enumerate(residue_constants.restypes):
        resname = residue_constants.restype_1to3[restype_letter]
        for chi_idx in range(4):
            if(residue_constants.chi_angles_mask[restype][chi_idx]):
                names = residue_constants.chi_angles_atoms[resname][chi_idx]
                restype_rigidgroup_base_atom_names[
                    restype, chi_idx + 4, :] = atom_names[1:]

    restype_rigidgroup_mask = torch.zeros(
        (*aatype.shape[:-1], 21, 8), 
        dtype=torch.float, 
        device=aatype.device, 
        requires_grad=False
    )
    restype_rigidgroup_mask[:, 0] = 1
    restype_rigidgroup_mask[:, 3] = 1
    restype_rigidgroup_mask[:20, 4:] = residue_constants.chi_angles_mask

    lookuptable = residue_constants.atom_order.copy()
    lookuptable[''] = 0
    lookup = np.vectorize(lambda x: lookuptable[x])
    restype_rigidgroup_base_atom37_idx = lookup(
        restype_rigidgroup_base_atom_names,
    )
    restype_rigidgroup_base_atom37_idx = aatype.new_tensor(
        restype_rigidgroup_base_atom37_idx,
    )
    restype_rigidgroup_base_atom37_idx = (
        restype_rigidgroup_base_atom37_idx.view(
            *((1,) * batch_dims), 
            *restype_rigidgroup_base_atom37_idx.shape
        )
    )

    residx_rigidgroup_base_atom37_idx = batched_gather(
        residx_rigidgroup_base_atom37_idx,
        aatype,
        dim=-3,
        no_batch_dims=batch_dims,
    )
    
    base_atom_pos = batched_gather(
        all_atom_positions,
        residx_rigidgroup_base_atom37_idx,
        dim=-2,
        no_batch_dims=len(all_atom_positions.shape[:-2]),
    )

    gt_frames = T.from_3_points(
        point_on_neg_x_axis=base_atom_pos[..., 0, :],
        origin=base_atom_pos[..., 1, :],
        point_on_xy_plane=base_atom_pos[..., 2, :],
    )

    group_exists = batched_gather(
        restype_rigidgroup_mask, 
        aatype, 
        dim=-2, 
        no_batch_dims=batch_dims,
    )

    gt_atoms_exist = batched_gather(
        all_atom_mask.float(),
        residx_rigidgroup_base_atom37_idx,
        dim=-1,
        no_batch_dims=len(all_atom_mask.shape[:-1])
    )
    gt_exists = torch.min(gt_atoms_exist, dim=-1) * group_exists

    rots = torch.eye(3, device=aatype.device, requires_grad=False)
    rots = rots.view(*((1,) * batch_dims), 1, 3, 3)
    rots = rots.expand(*((-1,) * batch_dims), 8, -1, -1)
    rots[..., 0, 0, 0] = -1
    rots[..., 0, 2, 2] = -1
    gt_frames = gt_frames.compose(T(rots, None)) 

    restype_rigidgroup_is_ambiguous = all_atom_mask.new_zeros(
        *((1,) * batch_dims), 21, 8
    )
    restype_rigidgroup_rots = torch.eye(
        3, device=aatype.device, requires_grad=False
    )
    restype_rigidgroup_rots = restype_rigidgroup_rots.view(
        *((1,) * batch_dims), 1, 1, 3, 3
    )
    restype_rigidgroup_rots = restype_rigidgroup_rots.expand(
        *((-1,) * batch_dims), 21, 8, 3, 3
    )

    for resname, _ in residue_constants.residue_atom_renaming_swaps.items():
        restype = residue_constants.restype_order[
            residue_constants.restype3to1[resname]
        ]
        chi_idx = int(sum(residue_constants.chi_angles_mask[restype]) - 1)
        restype_rigidgroup_is_ambiguous[..., restype, chi_idx + 4] = 1
        restype_rigidgroup_rots[..., restype, chi_idx + 4,  1, 1] = -1
        restype_rigidgroup_rots[..., restype, chi_idx + 4, 2, 2] = -1

    residx_rigidgroup_is_ambiguous = batched_gather(
        restype_rigidgroup_is_ambiguous,
        aatype,
        dim=-2,
        no_batch_dims=batch_dims,
    )

    residx_rigidgroup_ambiguity_rot = utils.batched_gather(
        restype_rigidgroup_rots,
        aatype,
        dim=-4,
        no_batch_dims=batch_dims,
    )

    alt_gt_frames = gt_frames.apply(T(residx_rigidgroup_ambiguity_rot, None))

    # TODO: Verify that I can get away with skipping the flat12 format
    gt_frames_tensor = gt_frames.to_tensor()
    alt_gt_frames_tensor = alt_gt_frames.to_tensor()

    return {
        'rigidgroups_gt_frames': gt_frames_tensor,
        'rigidgroups_gt_exists': gt_exists,
        'rigidgroups_group_exists': group_exists,
        'rigidgroups_group_is_ambiguous': residx_rigidgroup_is_ambiguous,
        'rigidgroups_alt_gt_frames': alt_gt_frames_tensor,
    }


def build_template_angle_feat(angle_feats, template_aatype):
    torsion_angles_sin_cos = angle_feats["torsion_angles_sin_cos"]
    alt_torsion_angles_sin_cos = angle_feats["alt_torsion_angles_sin_cos"]
    torsion_angles_mask = angle_feats["torsion_angles_mask"]
    template_angle_feat = torch.cat(
        [
            nn.functional.one_hot(template_aatype, 22),
            torsion_angles_sin_cos.reshape(
                *torsion_angles_sin_cos.shape[:-2], 14
            ),
            alt_torsion_angles_sin_cos.reshape(
                *alt_torsion_angles_sin_cos.shape[:-2], 14
            ),
            torsion_angles_mask,
        ], 
        dim=-1,
    )
    
    return template_angle_feat


def build_template_pair_feat(batch, min_bin, max_bin, no_bins, eps=1e-6, inf=1e8):
    template_mask = batch["template_pseudo_beta_mask"]
    template_mask_2d = template_mask[..., None] * template_mask[..., None, :]

    # Compute distogram (this seems to differ slightly from Alg. 5)
    tpb = batch["template_pseudo_beta"]
    dgram = torch.sum(
        (tpb[..., None, :] - tpb[..., None, :, :]) ** 2, dim=-1, keepdim=True)
    lower = torch.linspace(min_bin, max_bin, no_bins, device=tpb.device) ** 2
    upper = torch.cat([lower[:-1], lower.new_tensor([inf])], dim=-1)
    dgram = ((dgram > lower) * (dgram < upper)).type(dgram.dtype)

    to_concat = [dgram, template_mask_2d[..., None]]

    aatype_one_hot = nn.functional.one_hot(
        batch["template_aatype"], batch["target_feat"].shape[-1]
    )

    n_res = batch["template_aatype"].shape[-1]
    to_concat.append(
        aatype_one_hot[..., None, :, :].expand(
            *aatype_one_hot.shape[:-2], n_res, -1, -1
        )
    )
    to_concat.append(
        aatype_one_hot[..., None, :].expand(
            *aatype_one_hot.shape[:-2], -1, n_res, -1
        )
    )

    n, ca, c = [residue_constants.atom_order[a] for a in ['N', 'CA', 'C']]
    #t_aa_pos = batch["template_all_atom_positions"]
    #affines = T.make_transform_from_reference(
    #    n_xyz=t_aa_pos[..., n],
    #    ca_xyz=t_aa_pos[..., ca],
    #    c_xyz=t_aa_pos[..., c],
    #)
    #rots = affines.rots
    #trans = affines.trans
    #affine_vec = rot_mul_vec(
    #    rots.transpose(-1, -2), 
    #    trans[..., None, :, :] - trans[..., None, :],
    #)
    #inverted_dists = torch.rsqrt(eps + torch.sum(inverted_dists**2, dim=-1))

    t_aa_masks = batch["template_all_atom_masks"]
    template_mask = (
        t_aa_masks[..., n] * t_aa_masks[..., ca] * t_aa_masks[..., c]
    )
    template_mask_2d = template_mask[..., None] * template_mask[..., None, :]

    #inverted_dists *= template_mask_2d

    #unit_vector = affine_vec * inverted_dists.unsqueeze(-1)
    #unit_vector = unit_vector.unsqueeze(-2)
    unit_vector = template_mask_2d.new_zeros(*template_mask_2d.shape, 3)
    to_concat.append(unit_vector)
    to_concat.append(template_mask_2d[..., None])
    
    act = torch.cat(to_concat, dim=-1)

    act *= template_mask_2d[..., None]

    return act


def build_extra_msa_feat(batch):
    msa_1hot = nn.functional.one_hot(batch["extra_msa"], 23)
    msa_feat = [
        msa_1hot,
        batch["extra_has_deletion"].unsqueeze(-1),
        batch["extra_deletion_value"].unsqueeze(-1),
    ]
    return torch.cat(msa_feat, dim=-1)


# adapted from model/tf/data_transforms.py
def build_msa_feat(protein):
  """Create and concatenate MSA features."""
  # Whether there is a domain break. Always zero for chains, but keeping
  # for compatibility with domain datasets.
  has_break = batch["between_segment_residues"] 
  aatype_1hot = nn.functional.one_hot(batch['aatype'], num_classes=21)

  target_feat = [
      has_break.unsqueeze(-1),
      aatype_1hot,  # Everyone gets the original sequence.
  ]

  msa_1hot = nn.functional.one_hot(batch['msa'], num_classes=23)
  has_deletion = batch["deletion_matrix"]
  deletion_value = torch.atan(batch['deletion_matrix'] / 3.) * (2. / math.pi)

  msa_feat = [
      msa_1hot,
      has_deletion.unsqueeze(-1),
      deletion_value.unsqueeze(-1),
  ]

  if 'cluster_profile' in protein:
    deletion_mean_value = (
        tf.atan(batch['cluster_deletion_mean'] / 3.) * (2. / np.pi))
    msa_feat.extend([
        batch['cluster_profile'],
        tf.expand_dims(deletion_mean_value, axis=-1),
    ])

  if 'extra_deletion_matrix' in protein:
    batch['extra_has_deletion'] = tf.clip_by_value(
        batch['extra_deletion_matrix'], 0., 1.)
    batch['extra_deletion_value'] = tf.atan(
        batch['extra_deletion_matrix'] / 3.) * (2. / np.pi)

  batch['msa_feat'] = torch.cat(msa_feat, dim=-1)
  batch['target_feat'] = torch.cat(target_feat, dim=-1)
  return protein
