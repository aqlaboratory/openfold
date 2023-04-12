from structure_module import InvariantPointAttention, BackboneUpdate
from openfold.utils.rigid_utils import Rigid
from openfold.model.primitives import LayerNorm
from typing import Optional, Tuple, Sequence
from dropout import Dropout
import torch.nn as nn
import torch

from third_track_util_module import rbf, get_seqsep, init_lecun_normal, FeedForwardLayer, BiasedAxialAttention

#ADD MODULE

class Str2Str(nn.Module):
    def __init__(self, d_msa=256, d_z=128, d_state=16, 
            rbf_sigma=1.0, p_drop=0.1
    ):
        super(Str2Str, self).__init__()
        
        # initial node & z feature process
        self.norm_m = nn.LayerNorm(d_msa)
        self.norm_z = nn.LayerNorm(d_z)
        self.norm_s = nn.LayerNorm(d_state)
    
        self.embed_x = nn.Linear(d_msa+d_state, d_msa)
        self.embed_e = nn.Linear(d_z+36+1, d_z)
        
        self.norm_node = nn.LayerNorm(d_msa)
        
        self.rbf_sigma = rbf_sigma

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
                z: Optional[torch.Tensor], 
                r: Rigid,
                s: torch.Tensor, 
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
        s = self.norm_s(s)

        node = torch.cat((node, s), dim=-1)
        node = self.norm_node(self.embed_x(node))
        
        neighbor = get_seqsep(aatype)
        cas = rigids[:,:,1].contiguous()
        rbf_feat = rbf(torch.cdist(cas, cas), self.rbf_sigma)
        z = torch.cat((z, rbf_feat, neighbor), dim=-1)
        z = self.norm_edge(self.embed_e(z))

        rigids = Rigid.identity(
            s.shape[:-1], 
            s.dtype, 
            s.device, 
            self.training,
            fmt="quat",
        )

        for i in range(self.no_blocks):
            # [*, N, C_s]
            s = s + self.ipa(
                s, 
                z, 
                rigids, 
                mask, 
                inplace_safe=inplace_safe,
                _offload_inference=_offload_inference, 
                _z_reference_list=None,
            )
            s = self.ipa_dropout(s)
            s = self.layer_norm_ipa(s)
            s = self.transition(s)
           
            # [*, N]
            rigids = rigids.compose_q_update_vec(self.bb_update(s))

        return rigids, s

class PairStr2Pair(nn.Module):
    def __init__(self, d_pair=128, n_head=4, d_hidden=32, d_rbf=36, p_drop=0.15):
        super(PairStr2Pair, self).__init__()

        self.drop_row = Dropout(broadcast_dim=1, p_drop=p_drop)
        self.drop_col = Dropout(broadcast_dim=2, p_drop=p_drop)

        self.row_attn = BiasedAxialAttention(d_pair, d_rbf, n_head, d_hidden, p_drop=p_drop, is_row=True)
        self.col_attn = BiasedAxialAttention(d_pair, d_rbf, n_head, d_hidden, p_drop=p_drop, is_row=False)

        self.ff = FeedForwardLayer(d_pair, 2)
        
    def forward(self, z, rbf_feat):
        z = z + self.drop_row(self.row_attn(z, rbf_feat))
        z = z + self.drop_col(self.col_attn(z, rbf_feat))
        z = z + self.ff(z)
        return z