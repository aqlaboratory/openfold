from openfold.model.structure_module import InvariantPointAttention, BackboneUpdate, AngleResnet
from openfold.utils.rigid_utils import Rigid, Rotation
from openfold.model.primitives import LayerNorm
from typing import Optional, Tuple, Sequence
from dropout import Dropout
import torch.nn as nn
import torch

from third_track_util_module import rbf, get_seqsep, init_lecun_normal, FeedForwardLayer, BiasedAxialAttention

#ADD MODULE

class PairStr2Pair(nn.Module):
    def __init__(self, c_z=128, c_head=4, c_hidden=32, c_rbf=36, p_drop=0.15):
        super(PairStr2Pair, self).__init__()

        self.drop_row = Dropout(broadcast_dim=1, p_drop=p_drop)
        self.drop_col = Dropout(broadcast_dim=2, p_drop=p_drop)

        self.row_attn = BiasedAxialAttention(c_z, c_rbf, c_head, c_hidden, p_drop=p_drop, is_row=True)
        self.col_attn = BiasedAxialAttention(c_z, c_rbf, c_head, c_hidden, p_drop=p_drop, is_row=False)

        self.ff = FeedForwardLayer(c_z, 2)
        
    def forward(self, z, rbf_feat):
        z = z + self.drop_row(self.row_attn(z, rbf_feat))
        z = z + self.drop_col(self.col_attn(z, rbf_feat))
        z = z + self.ff(z)
        return z

class Str2Str(nn.Module):
    def __init__(self, c_m=256, c_z=128, c_s=16, 
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
            self.c_s,
            self.c_z,
            self.c_ipa,
            self.no_heads_ipa,
            self.no_qk_points,
            self.no_v_points,
            inf=self.inf,
            eps=self.epsilon,
        )

        self.ipa_dropout = nn.Dropout(self.dropout_rate)
        self.layer_norm_ipa = LayerNorm(self.c_s)

        self.bb_update = BackboneUpdate(self.c_s)

        self.angle_resnet = AngleResnet(
            self.c_s,
            self.c_resnet,
            self.no_resnet_blocks,
            self.no_angles,
            self.epsilon,
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
                mask: torch.Tensor,
                inplace_safe: bool = False,
                _offload_inference: bool = False,
                _z_reference_list: Optional[Sequence[torch.Tensor]] = None,
                ):
        # process msa & z features
        B, N, L = m.shape[:3]
        node = self.norm_m(m[:,0])
        z = self.norm_z(z)
        state = self.norm_s(state)

        state_initial = state

        node = torch.cat((node, state), dim=-1)
        node = self.norm_node(self.embed_x(node))
        
        neighbor = get_seqsep(aatype)
        z = torch.cat((z, rbf_feat, neighbor), dim=-1)
        z = self.norm_edge(self.embed_e(z))

        rigids = Rigid.identity(
            state.shape[:-1], 
            state.dtype, 
            state.device, 
            self.training,
            fmt="quat",
        )

        for i in range(self.no_blocks):
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
                self.trans_scale_factor
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

        return pred_xyz, state