from openfold.model.structure_module import InvariantPointAttention, BackboneUpdate, AngleResnet
from openfold.utils.rigid_utils import Rigid, Rotation
from openfold.model.primitives import LayerNorm
from typing import Optional, Tuple, Sequence
from openfold.model.dropout import Dropout
import torch.nn as nn
import torch

from openfold.model.third_track_util_module import rbf, get_seqsep, init_lecun_normal, FeedForwardLayer, BiasedAxialAttention

from openfold.utils.feats import (
    frames_and_literature_positions_to_atom14_pos,
    torsion_angles_to_frames,
)

from openfold.np.residue_constants import (
    restype_rigid_group_default_frame,
    restype_atom14_to_rigid_group,
    restype_atom14_mask,
    restype_atom14_rigid_group_positions,
)

#ADD MODULE

class PairStr2Pair(nn.Module):
    def __init__(self, c_z=128, c_head=4, c_hidden=32, c_rbf=36, p_drop=0.15):
        super(PairStr2Pair, self).__init__()

        #Dropout(broadcast_dim=1, p_drop=p_drop)
        self.drop_row = Dropout(batch_dim=1, r=p_drop)
        self.drop_col = Dropout(batch_dim=2, r=p_drop)

        self.row_attn = BiasedAxialAttention(c_z, c_rbf, c_head, c_hidden, p_drop=p_drop, is_row=True)
        self.col_attn = BiasedAxialAttention(c_z, c_rbf, c_head, c_hidden, p_drop=p_drop, is_row=False)

        self.ff = FeedForwardLayer(c_z, 2)
        
    def forward(self, z, rbf_feat):
        z = z + self.drop_row(self.row_attn(z, rbf_feat))
        z = z + self.drop_col(self.col_attn(z, rbf_feat))
        z = z + self.ff(z)
        return z

class Str2Str(nn.Module):
    def __init__(self, c_m=64, c_z=128, c_s=16, 
            rbf_scale=1.0, p_drop=0.1
    ):
        super(Str2Str, self).__init__()
        
        # initial node & z feature process
        self.norm_m = nn.LayerNorm(c_m)
        self.norm_z = nn.LayerNorm(c_z)
        self.norm_s = nn.LayerNorm(c_s)
    
        self.embed_x = nn.Linear(c_m+c_s, c_m)
        self.embed_e = nn.Linear(c_z+36+1, c_z)
        
        self.norm_node = nn.LayerNorm(c_m)
        
        self.rbf_scale = rbf_scale

        self.ipa = InvariantPointAttention(
            c_s,
            c_z,
            c_hidden=16,
            no_heads=12,
            no_qk_points=4,
            no_v_points=8,
            inf=1e5,
            eps=1e-12,
        )

        self.ipa_dropout = nn.Dropout(p_drop)
        self.layer_norm_ipa = LayerNorm(c_s)

        self.bb_update = BackboneUpdate(c_s)

        self.angle_resnet = AngleResnet(
            c_s,
            c_hidden=128,
            no_blocks=2,
            no_angles=7,
            epsilon=1e-12,
        )

        self.reset_parameter()

    def reset_parameter(self):
        # initialize weights to normal distribution
        self.embed_x = init_lecun_normal(self.embed_x)
        self.embed_e = init_lecun_normal(self.embed_e)

        # initialize bias to zeros
        nn.init.zeros_(self.embed_x.bias)
        nn.init.zeros_(self.embed_e.bias)
    
    @torch.cuda.amp.autocast(enabled=False)
    def forward(self, 
                m, 
                z: torch.Tensor, 
                state: torch.Tensor, 
                rbf_feat: torch.Tensor, 
                aatype,
                rigids,
                mask: torch.Tensor,
                inplace_safe: bool = False,
                _offload_inference: bool = False,
                _z_reference_list: Optional[Sequence[torch.Tensor]] = None,
                ):
        # process msa & z features
        B, N, L = m.shape[:3]
        node = self.norm_m(m[..., 0, :, :])
        z = self.norm_z(z)
        state = self.norm_s(state)

        state_initial = state

        node = torch.cat((node, state), dim=-1)
        node = self.norm_node(self.embed_x(node))
        
        neighbor = get_seqsep(aatype)
        z = torch.cat((z, rbf_feat, neighbor), dim=-1)
        z = self.norm_edge(self.embed_e(z))
        
        if mask is None:
            # [*, N]
            mask = state.new_ones(state.shape[:-1])

        # [*, N, C_s]
        state = state + self.ipa(
            state, 
            z, 
            rigids, 
            mask, 
            inplace_safe=inplace_safe,
            _offload_inference=_offload_inference, 
            _z_reference_list=None,
        )
        state = self.ipa_dropout(state)
        state = self.layer_norm_ipa(state)
        state = self.transition(state)
        
        # [*, N]
        rigids = rigids.compose_q_update_vec(self.bb_update(state))

        # To hew as closely as possible to AlphaFold, we convert our
        # quaternion-based transformations to rotation-matrix ones
        # here
        backb_to_global = Rigid(
            Rotation(
                rot_mats=rigids.get_rots().get_rot_mats(), 
                quats=None
            ),
            rigids.get_trans(),
        )

        backb_to_global = backb_to_global.scale_translation(
            trans_scale_factor=10,
        )

        # [*, N, 7, 2]
        unnormalized_angles, angles = self.angle_resnet(state, state_initial)

        all_frames_to_global = self.torsion_angles_to_frames(
            backb_to_global,
            angles,
            aatype,
        )

        pred_xyz = self.frames_and_literature_positions_to_atom14_pos(
            all_frames_to_global,
            aatype,
        )

        rigids = rigids.stop_rot_gradient()

        return pred_xyz, rigids, state
    
    def _init_residue_constants(self, float_dtype, device):
        if not hasattr(self, "default_frames"):
            self.register_buffer(
                "default_frames",
                torch.tensor(
                    restype_rigid_group_default_frame,
                    dtype=float_dtype,
                    device=device,
                    requires_grad=False,
                ),
                persistent=False,
            )
        if not hasattr(self, "group_idx"):
            self.register_buffer(
                "group_idx",
                torch.tensor(
                    restype_atom14_to_rigid_group,
                    device=device,
                    requires_grad=False,
                ),
                persistent=False,
            )
        if not hasattr(self, "atom_mask"):
            self.register_buffer(
                "atom_mask",
                torch.tensor(
                    restype_atom14_mask,
                    dtype=float_dtype,
                    device=device,
                    requires_grad=False,
                ),
                persistent=False,
            )
        if not hasattr(self, "lit_positions"):
            self.register_buffer(
                "lit_positions",
                torch.tensor(
                    restype_atom14_rigid_group_positions,
                    dtype=float_dtype,
                    device=device,
                    requires_grad=False,
                ),
                persistent=False,
            )
    
    def torsion_angles_to_frames(self, r, alpha, f):
        # Lazily initialize the residue constants on the correct device
        self._init_residue_constants(alpha.dtype, alpha.device)
        # Separated purely to make testing less annoying
        return torsion_angles_to_frames(r, alpha, f, self.default_frames)

    def frames_and_literature_positions_to_atom14_pos(
        self, r, f  # [*, N, 8]  # [*, N]
    ):
        # Lazily initialize the residue constants on the correct device
        self._init_residue_constants(r.get_rots().dtype, r.get_rots().device)
        return frames_and_literature_positions_to_atom14_pos(
            r,
            f,
            self.default_frames,
            self.group_idx,
            self.atom_mask,
            self.lit_positions,
        )