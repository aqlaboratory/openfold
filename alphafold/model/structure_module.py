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

import math
import torch
import torch.nn as nn
from typing import Optional

from alphafold.model.primitives import Linear, ipa_point_weights_init_
from alphafold.np.residue_constants import (
    restype_rigid_group_default_frame,
    restype_atom14_to_rigid_group,
    restype_atom14_mask,
    restype_atom14_rigid_group_positions,
)
from alphafold.utils.affine_utils import T, quat_to_rot 
from alphafold.utils.tensor_utils import (
    stack_tensor_dicts, 
    permute_final_dims, 
    flatten_final_dims,
)


class AngleResnetBlock(nn.Module):
    def __init__(self, c_hidden):
        """
            Args:
                c_hidden:
                    Hidden channel dimension
        """
        super(AngleResnetBlock, self).__init__()

        self.c_hidden = c_hidden

        self.linear_1 = Linear(self.c_hidden, self.c_hidden, init="relu")
        self.linear_2 = Linear(self.c_hidden, self.c_hidden, init="final")

        self.relu = nn.ReLU()

    def forward(self, a):

        s_initial = a

        a = self.relu(a)
        a = self.linear_1(a)
        a = self.relu(a)
        a = self.linear_2(a)

        return a + s_initial


class AngleResnet(nn.Module):
    """
        Implements Algorithm 20, lines 11-14
    """
    def __init__(self, c_in, c_hidden, no_blocks, no_angles, epsilon):
        """
            Args:
                c_in:
                    Input channel dimension
                c_hidden:
                    Hidden channel dimension
                no_blocks:
                    Number of resnet blocks
                no_angles:
                    Number of torsion angles to generate
                epsilon:
                    Small constant for normalization
        """
        super(AngleResnet, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_blocks = no_blocks
        self.no_angles = no_angles
        self.epsilon = epsilon

        self.linear_in = Linear(self.c_in, self.c_hidden)
        self.linear_initial = Linear(self.c_in, self.c_hidden)

        self.layers = nn.ModuleList()
        for _ in range(self.no_blocks):
            layer = AngleResnetBlock(c_hidden=self.c_hidden)
            self.layers.append(layer)

        self.linear_out = Linear(self.c_hidden, self.no_angles * 2)

        self.relu = nn.ReLU()

    def forward(self, s, s_initial):
        """
            Args:
                s:
                    [*, C_hidden] single embedding
                s_initial:
                    [*, C_hidden] single embedding as of the start of the 
                    StructureModule
            Returns:
                [*, no_angles, 2] predicted angles
        """       
        # NOTE: The ReLU's applied to the inputs are absent from the supplement 
        # pseudocode but present in the source. For maximal compatibility with 
        # the pretrained weights, I'm going with the source.

        # [*, C_hidden]
        s_initial = self.relu(s_initial)
        s_initial = self.linear_initial(s_initial)
        s = self.relu(s)
        s = self.linear_in(s)
        s = s + s_initial

        for l in self.layers:
            s = l(s)

        s = self.relu(s)

        # [*, no_angles * 2]
        s = self.linear_out(s)

        # [*, no_angles, 2]
        s = s.view(*s.shape[:-1], -1, 2)
        norm_denom = torch.sqrt(
            torch.clamp(
                torch.sum(s ** 2, dim=-1, keepdims=True),
                min=self.epsilon,
            )
        )
        s = s / norm_denom

        return s


class InvariantPointAttention(nn.Module):
    """
        Implements Algorithm 22.
    """
    def __init__(self,
        c_s,
        c_z,
        c_hidden,
        no_heads,
        no_qk_points,
        no_v_points,
        inf=1e5,
        eps=1e-8,
    ):
        """
            Args:
                c_s:
                    Single representation channel dimension
                c_z:
                    Pair representation channel dimension
                c_hidden:
                    Hidden channel dimension
                no_heads:
                    Number of attention heads
                no_qk_points:
                    Number of query/key points to generate
                no_v_points:
                    Number of value points to generate
        """
        super(InvariantPointAttention, self).__init__()

        self.c_s = c_s
        self.c_z = c_z
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.no_qk_points = no_qk_points
        self.no_v_points = no_v_points
        self.inf = inf
        self.eps = eps

        # These linear layers differ from their specifications in the
        # supplement. There, they lack bias and use Glorot initialization.
        # Here as in the official source, they have bias and use the default
        # Lecun initialization.
        hc = self.c_hidden * self.no_heads
        self.linear_q = Linear(self.c_s, hc)
        self.linear_kv = Linear(self.c_s, 2 * hc)

        hpq = self.no_heads * self.no_qk_points * 3
        self.linear_q_points = Linear(self.c_s, hpq)

        hpkv = self.no_heads * (self.no_qk_points + self.no_v_points) * 3
        self.linear_kv_points = Linear(self.c_s, hpkv)

        hpv = self.no_heads * self.no_v_points * 3

        self.linear_b = Linear(self.c_z, self.no_heads)

        self.head_weights = nn.Parameter(torch.zeros((no_heads)))
        ipa_point_weights_init_(self.head_weights)

        concat_out_dim = self.no_heads * (self.c_z
                                          + self.c_hidden
                                          + self.no_v_points * 4)
        self.linear_out = Linear(concat_out_dim, self.c_s, init="final")

        self.softmax = nn.Softmax(dim=-1)
        self.softplus = nn.Softplus()

    def forward(self, 
        s: torch.Tensor, 
        z: torch.Tensor, 
        t: T,
        mask: torch.Tensor,
    ):
        """
            Args:
                s:
                    [*, N_res, C_s] single representation
                z:
                    [*, N_res, N_res, C_z] pair representation
                t:
                    [*, N_res] affine transformation object
                mask:
                    [*, N_res] mask
            Returns:
                [*, N_res, C_s] single representation update
        """
        #######################################
        # Generate scalar and point activations
        #######################################

        # [*, N_res, H * C_hidden]
        q = self.linear_q(s)
        kv = self.linear_kv(s)

        # [*, N_res, H, C_hidden]
        q = q.view(*q.shape[:-1], self.no_heads, -1)

        # [*, N_res, H, 2 * C_hidden]
        kv = kv.view(*kv.shape[:-1], self.no_heads, -1)

        # [*, N_res, H, C_hidden]
        k, v = torch.split(kv, self.c_hidden, dim=-1)

        # [*, N_res, H * P_q * 3]
        q_pts = self.linear_q_points(s)

        # This is kind of clunky, but it's how the original does it
        # [*, N_res, H * P_q, 3]
        q_pts = torch.split(q_pts, q_pts.shape[-1] // 3, dim=-1)
        q_pts = torch.stack(q_pts, dim=-1) 
        q_pts = t[..., None].apply(q_pts)

        # [*, N_res, H, P_q, 3]
        q_pts = q_pts.view(
             *q_pts.shape[:-2], self.no_heads, self.no_qk_points, 3
        )

        # [*, N_res, H * (P_q + P_v) * 3]
        kv_pts = self.linear_kv_points(s)

        # [*, N_res, H * (P_q + P_v), 3]
        kv_pts = torch.split(kv_pts, kv_pts.shape[-1] // 3, dim=-1)
        kv_pts = torch.stack(kv_pts, dim=-1)
        kv_pts = t[..., None].apply(kv_pts)

        # [*, N_res, H, (P_q + P_v), 3]
        kv_pts = kv_pts.view(
            *kv_pts.shape[:-2], self.no_heads, -1, 3
        )

        # [*, N_res, H, P_q/P_v, 3]
        k_pts, v_pts = torch.split(
            kv_pts, 
            [self.no_qk_points, self.no_v_points], 
            dim=-2
        )

        ##########################
        # Compute attention scores
        ##########################

        # [*, N_res, N_res, H]
        b = self.linear_b(z)

        # [*, H, N_res, N_res]
        a = torch.matmul(
            permute_final_dims(q, 1, 0, 2), # [*, H, N_res, C_hidden]
            permute_final_dims(k, 1, 2, 0), # [*, H, C_hidden, N_res]
        )
        a *= math.sqrt(1. / (3 * self.c_hidden))
        a += math.sqrt(1. / 3) * permute_final_dims(b, 2, 0, 1)
        
        # [*, N_res, N_res, H, P_q, 3]
        pt_att = q_pts.unsqueeze(-4) - k_pts.unsqueeze(-5)
        pt_att = pt_att ** 2

        # [*, N_res, N_res, H, P_q]
        pt_att = torch.sum(pt_att, dim=-1)
        head_weights = self.softplus(self.head_weights).view(
            *((1,) * len(pt_att.shape[:-2]) + (-1, 1))
        ) 
        head_weights *= math.sqrt(1. / (3 * (self.no_qk_points * 9. / 2)))
        pt_att = pt_att * head_weights 
        
        # [*, N_res, N_res, H]
        pt_att = torch.sum(pt_att, dim=-1) * (-0.5)

        # [*, N_res, N_res]
        square_mask = mask.unsqueeze(-1) * mask.unsqueeze(-2)
        square_mask = self.inf * (square_mask - 1)

        # [*, H, N_res, N_res]
        pt_att = permute_final_dims(pt_att, 2, 0, 1)
        a += pt_att
        a += square_mask.unsqueeze(-3)
        a = self.softmax(a)

        ################
        # Compute output
        ################

        # [*, N_res, H, C_hidden]
        o = torch.matmul(a, v.transpose(-2, -3)).transpose(-2, -3)

        # [*, N_res, H * C_hidden]
        o = flatten_final_dims(o, 2)

        # [*, H, 3, N_res, P_v]
        o_pt = torch.matmul(
            a.unsqueeze(-3),                       # [*, H, 1, N_res, N_res]
            permute_final_dims(v_pts, 1, 3, 0, 2), # [*, H, 3, N_res, P_v]
        )

        # [*, N_res, H, P_v, 3]
        o_pt = permute_final_dims(o_pt, 2, 0, 3, 1)
        o_pt = t[..., None, None].invert_apply(o_pt)

        # [*, N_res, H * P_v]
        o_pt_norm = flatten_final_dims(
            torch.sqrt(torch.sum(o_pt ** 2, dim=-1) + self.eps), 
            2
        )

        # [*, N_res, H * P_v, 3]
        o_pt = o_pt.view(*o_pt.shape[:-3], -1, 3)

        # [*, N_res, H, C_z]
        o_pair = torch.matmul(a.transpose(-2, -3), z)

        # [*, N_res, H * C_z]
        o_pair = flatten_final_dims(o_pair, 2)

        # [*, N_res, C_s]
        s = self.linear_out(
                torch.cat((
                    o, 
                    *torch.unbind(o_pt, dim=-1), 
                    o_pt_norm, 
                    o_pair
                ), dim=-1)
            )

        return s


class BackboneUpdate(nn.Module):
    """
        Implements Algorithm 23.
    """
    def __init__(self, c_s):
        """
            Args:
                c_s:
                    Single representation channel dimension
        """
        super(BackboneUpdate, self).__init__()

        self.c_s = c_s

        self.linear = Linear(self.c_s, 6, init="final")

    def forward(self, s):
        """
            Args:
                [*, N_res, C_s] single representation
            Returns:
                [*, N_res] affine transformation object
        """
        # [*, 6]
        params = self.linear(s)

        # [*, 3]
        quats, trans = params[...,:3], params[...,3:]

        # [*]
        norm_denom = torch.sqrt(torch.sum(quats ** 2, dim=-1) + 1)

        # As many ones as there are dimensions in quats
        ones = s.new_ones((1,) * len(quats.shape))

        # [*, 4]
        quats = torch.cat((ones.expand(*quats.shape[:-1], 1), quats), dim=-1)
        quats = quats / norm_denom.unsqueeze(-1)

        # [*, 3, 3]
        rots = quat_to_rot(quats)

        return T(rots, trans)


def _torsion_angles_to_frames(t, alpha, f, rrgdf):
    # [*, N, 8, 4, 4]
    default_4x4 = rrgdf[f,...]
    
    # [*, N, 8] transformations, i.e.
    #   One [*, N, 8, 3, 3] rotation matrix and
    #   One [*, N, 8, 3]    translation matrix
    default_t = T.from_4x4(default_4x4)

    bb_rot = alpha.new_zeros((*((1,) * len(alpha.shape[:-1])), 2))
    bb_rot[..., 1] = 1
    
    # [*, N, 8, 2]
    alpha = torch.cat(
        [bb_rot.expand(*alpha.shape[:-2], -1, -1), alpha], 
        dim=-2
    )

    # [*, N, 8, 3, 3]
    # Produces rotation matrices of the form:
    # [
    #   [1, 0  , 0  ],
    #   [0, a_2,-a_1],
    #   [0, a_1, a_2]
    # ]
    # This follows the original code rather than the supplement, which uses
    # different indices.
        
    all_rots = alpha.new_zeros(default_t.rots.shape)
    all_rots[..., 0, 0] = 1
    all_rots[..., 1, 1] = alpha[..., 1]
    all_rots[..., 1, 2] = -alpha[..., 0]
    all_rots[..., 2, 1:] = alpha

    all_rots = T(all_rots, None)

    all_frames = default_t.compose(all_rots)

    chi2_frame_to_frame = all_frames[..., 5]
    chi3_frame_to_frame = all_frames[..., 6]
    chi4_frame_to_frame = all_frames[..., 7]

    chi1_frame_to_bb = all_frames[..., 4]
    chi2_frame_to_bb = chi1_frame_to_bb.compose(chi2_frame_to_frame)
    chi3_frame_to_bb = chi2_frame_to_bb.compose(chi3_frame_to_frame)
    chi4_frame_to_bb = chi3_frame_to_bb.compose(chi4_frame_to_frame)

    all_frames_to_bb = T.concat([
            all_frames[..., :5],
            chi2_frame_to_bb.unsqueeze(-1),
            chi3_frame_to_bb.unsqueeze(-1),
            chi4_frame_to_bb.unsqueeze(-1),
        ], dim=-1,
    )

    all_frames_to_global = t[..., None].compose(all_frames_to_bb)

    return all_frames_to_global


def _frames_and_literature_positions_to_atom14_pos(
    t,
    f,
    default_frames,
    group_idx,
    atom_mask,
    lit_positions,
):

    # [*, N, 14, 4, 4] 
    default_4x4 = default_frames[f,...]

    # [*, N, 14]
    group_mask = group_idx[f,...]

    # [N, 14, 8]
    group_mask = nn.functional.one_hot(
        group_mask, num_classes=default_frames.shape[-3],
    )
 
    # [*, N, 14, 8]
    t_atoms_to_global = t[..., None, :] * group_mask

    # [*, N, 14]
    t_atoms_to_global = t_atoms_to_global.map_tensor_fn(
        lambda x: torch.sum(x, dim=-1)
    )

    # [N, 14, 1]
    atom_mask = atom_mask[f,...].unsqueeze(-1)

    # [N, 14, 3]
    lit_positions = lit_positions[f,...]
    pred_positions = t_atoms_to_global.apply(lit_positions)
    pred_positions *= atom_mask

    return pred_positions


class StructureModuleTransitionLayer(nn.Module):
    def __init__(self, c):
        super(StructureModuleTransitionLayer, self).__init__()

        self.c = c
        
        self.linear_1 = Linear(self.c, self.c, init="relu")
        self.linear_2 = Linear(self.c, self.c, init="relu")
        self.linear_3 = Linear(self.c, self.c, init="final")

        self.relu = nn.ReLU()

    def forward(self, s):
        s_initial = s
        s = self.linear_1(s)
        s = self.relu(s)
        s = self.linear_2(s)
        s = self.relu(s)
        s = self.linear_3(s)

        s = s + s_initial

        return s


class StructureModuleTransition(nn.Module):
    def __init__(self, c, num_layers, dropout_rate):
        super(StructureModuleTransition, self).__init__()

        self.c = c
        self.num_layers = num_layers
        self.dropout_rate = dropout_rate

        self.layers = nn.ModuleList()
        for _ in range(self.num_layers):
            l = StructureModuleTransitionLayer(self.c)
            self.layers.append(l)

        self.dropout = nn.Dropout(self.dropout_rate)
        self.layer_norm = nn.LayerNorm(self.c)

    def forward(self, s):
        for l in self.layers:
            s = l(s)

        s = self.dropout(s)
        s = self.layer_norm(s)

        return s


class StructureModule(nn.Module):
    def __init__(self, 
        c_s, 
        c_z,
        c_ipa,
        c_resnet,
        no_heads_ipa,
        no_qk_points,
        no_v_points,
        dropout_rate,
        no_blocks,
        no_transition_layers,
        no_resnet_blocks,
        no_angles,
        trans_scale_factor,
        epsilon,
        inf,
        **kwargs,
    ):
        """
            Args:
                c_s:
                    Single representation channel dimension
                c_z:
                    Pair representation channel dimension
                c_ipa:
                    IPA hidden channel dimension
                c_resnet:
                    Angle resnet (Alg. 23 lines 11-14) hidden channel dimension
                no_heads_ipa:
                    Number of IPA heads
                no_qk_points:
                    Number of query/key points to generate during IPA
                no_v_points:
                    Number of value points to generate during IPA
                dropout_rate:
                    Dropout rate used throughout the layer
                no_blocks:
                    Number of structure module blocks
                no_transition_layers:
                    Number of layers in the single representation transition
                    (Alg. 23 lines 8-9)
                no_resnet_blocks:
                    Number of blocks in the angle resnet
                no_angles:
                    Number of angles to generate in the angle resnet
                trans_scale_factor:
                    Scale of single representation transition hidden dimension
                epsilon:
                    Small number used in angle resnet normalization
                inf:
                    Large number used for attention masking
        """             
        super(StructureModule, self).__init__()

        self.c_s = c_s
        self.c_z = c_z
        self.c_ipa = c_ipa
        self.c_resnet = c_resnet
        self.no_heads_ipa = no_heads_ipa
        self.no_qk_points = no_qk_points
        self.no_v_points = no_v_points
        self.dropout_rate = dropout_rate
        self.no_blocks = no_blocks
        self.no_transition_layers = no_transition_layers
        self.no_resnet_blocks = no_resnet_blocks
        self.no_angles = no_angles
        self.trans_scale_factor = trans_scale_factor
        self.epsilon = epsilon
        self.inf = inf

        # To be lazily initialized later
        self.default_frames = None
        self.group_idx = None
        self.atom_mask = None
        self.lit_positions = None

        self.layer_norm_s = nn.LayerNorm(self.c_s)
        self.layer_norm_z = nn.LayerNorm(self.c_z)

        self.linear_in = Linear(self.c_s, self.c_s)

        self.ipa = InvariantPointAttention(
            self.c_s,
            self.c_z,
            self.c_ipa,
            self.no_heads_ipa,
            self.no_qk_points,
            self.no_v_points,
            inf=self.inf,
        )

        self.ipa_dropout = nn.Dropout(self.dropout_rate)
        self.layer_norm_ipa = nn.LayerNorm(self.c_s)

        self.transition = StructureModuleTransition(
            self.c_s,
            self.no_transition_layers,
            self.dropout_rate,
        )

        self.bb_update = BackboneUpdate(self.c_s)

        self.angle_resnet = AngleResnet(
            self.c_s, 
            self.c_resnet, 
            self.no_resnet_blocks,
            self.no_angles,
            self.epsilon,
        )

    def forward(self, 
        s,
        z,
        f,
        mask=None,
    ):
        """
            Args:
                s:
                    [*, N_res, C_s] single representation
                z:
                    [*, N_res, N_res, C_z] pair representation
                f:
                    [*, N_res] amino acid indices
                mask:
                    Optional [*, N_res] sequence mask
            Returns:
                A dictionary of outputs
        """
        if(mask is None):
            # [*, N]
            mask = s.new_ones(s.shape[:-1], requires_grad=False)

        # [*, N, C_s]
        s = self.layer_norm_s(s)

        # [*, N, N, C_z]
        z = self.layer_norm_z(z)

        # [*, N, C_s]
        s_initial = s
        s = self.linear_in(s)

        # [*, N]
        t = T.identity(s.shape[:-1], s.dtype, s.device, self.training)
 
        outputs = []
        for l in range(self.no_blocks):
            # [*, N, C_s]
            s = s + self.ipa(s, z, t, mask)
            s = self.ipa_dropout(s)
            s = self.layer_norm_ipa(s)
            s = self.transition(s)

            # [*, N]
            t = t.compose(self.bb_update(s))

            # [*, N, 7, 2]
            a = self.angle_resnet(s, s_initial)

            all_frames_to_global = self.torsion_angles_to_frames(
                t.scale_translation(self.trans_scale_factor), a, f,
            )

            pred_xyz = self.frames_and_literature_positions_to_atom14_pos(
                all_frames_to_global,
                f,
            )

            preds = {
                "transformations": 
                    t.scale_translation(self.trans_scale_factor).to_4x4(),
                "angles": a,
                "positions": pred_xyz,
            }

            outputs.append(preds)

            if(l < self.no_blocks - 1):
                t = t.stop_rot_gradient()

        outputs = stack_tensor_dicts(outputs)
        outputs["single_act"] = s

        return outputs

    def _init_residue_constants(self, device):
        if(self.default_frames is None):
            self.default_frames = torch.tensor(
                restype_rigid_group_default_frame, device=device,
            )
        if(self.group_idx is None):
            self.group_idx = torch.tensor(
                restype_atom14_to_rigid_group, device=device,
            )
        if(self.atom_mask is None):
            self.atom_mask = torch.tensor(
                restype_atom14_mask, device=device,
            )
        if(self.lit_positions is None):
            self.lit_positions = torch.tensor(
                restype_atom14_rigid_group_positions, device=device,
            )

    def torsion_angles_to_frames(self, t, alpha, f):
        # Lazily initialize the residue constants on the correct device
        self._init_residue_constants(f.device)
        # Separated purely to make testing less annoying
        return _torsion_angles_to_frames(t, alpha, f, self.default_frames)

    def frames_and_literature_positions_to_atom14_pos(self, 
        t,    # [*, N, 8]
        f     # [*, N]
    ):
        # Lazily initialize the residue constants on the correct device
        # TODO: Maybe this stuff should be done on CPU instead (so these
        # arrays 
        self._init_residue_constants(f.device)
        
        return _frames_and_literature_positions_to_atom14_pos(
            t, 
            f, 
            self.default_frames, 
            self.group_idx,
            self.atom_mask,
            self.lit_positions, 
        )


if __name__ == "__main__":
    c_m = 11
    c_z = 13
    c_hidden = 17
    no_heads = 3
    no_qp = 5
    no_vp = 7

    batch_size = 2

    s = torch.rand((batch_size, n_res, c_m))
    z = torch.rand((batch_size, n_res, n_res, c_z))

    rots = torch.rand((batch_size, n_res, 3, 3))
    trans = torch.rand((batch_size, n_res, 3))

    t = (rots, trans)

    ipa = InvariantPointAttention(c_m, c_z, c_hidden, no_heads, no_qp, no_vp)

    shape_before = s.shape
    s = ipa(s, z, t)

    assert(s.shape == shape_before)
