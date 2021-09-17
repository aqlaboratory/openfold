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

import torch
import torch.nn as nn
from typing import Tuple, Optional
from functools import partial

from alphafold.model.primitives import Linear
from alphafold.utils.deepspeed import checkpoint_blocks 
from alphafold.model.dropout import DropoutRowwise, DropoutColumnwise
from alphafold.model.msa import (
    MSARowAttentionWithPairBias,
    MSAColumnAttention,
    MSAColumnGlobalAttention,
)
from alphafold.model.outer_product_mean import OuterProductMean
from alphafold.model.pair_transition import PairTransition
from alphafold.model.triangular_attention import (
    TriangleAttentionStartingNode,
    TriangleAttentionEndingNode,
)
from alphafold.model.triangular_multiplicative_update import (
    TriangleMultiplicationOutgoing,
    TriangleMultiplicationIncoming,
)
from alphafold.utils.tensor_utils import chunk_layer


class MSATransition(nn.Module):
    """
        Feed-forward network applied to MSA activations after attention.

        Implements Algorithm 9
    """
    def __init__(self, c_m, n, chunk_size):
        """
            Args:
                c_m:
                    MSA channel dimension
                n:
                    Factor multiplied to c_m to obtain the hidden channel 
                    dimension
        """
        super(MSATransition, self).__init__()

        self.c_m = c_m
        self.n = n
        self.chunk_size = chunk_size

        self.layer_norm = nn.LayerNorm(self.c_m)
        self.linear_1 = Linear(self.c_m, self.n * self.c_m, init="relu")
        self.relu = nn.ReLU()
        self.linear_2 = Linear(self.n * self.c_m, self.c_m, init="final")
        
    def _transition(self, m, mask):
        m = self.linear_1(m)
        m = self.relu(m)
        m = self.linear_2(m) * mask
        return m

    def forward(self, 
        m: torch.Tensor,
        mask: torch.Tensor = None,
    ) -> torch.Tensor:
        """
            Args:
                m:
                    [*, N_seq, N_res, C_m] MSA activation
                mask:
                    [*, N_seq, N_res, C_m] MSA mask
            Returns:
                m:
                    [*, N_seq, N_res, C_m] MSA activation update
        """
        # DISCREPANCY: DeepMind forgets to apply the MSA mask here.
        if(mask is None):
            mask = m.new_ones(m.shape[:-1])

        mask = mask.unsqueeze(-1)

        m = self.layer_norm(m)

        inp = {"m": m, "mask": mask}
        if(not self.training and self.chunk_size is not None):
            m = chunk_layer(
                self._transition,
                inp,
                chunk_size=self.chunk_size,
                no_batch_dims=len(m.shape[:-2]),
            )
        else:
            m = self._transition(**inp)

        return m


class EvoformerBlock(nn.Module):
    def __init__(self,
        c_m: int,
        c_z: int,
        c_hidden_msa_att: int,
        c_hidden_opm: int,
        c_hidden_mul: int,
        c_hidden_pair_att: int,
        no_heads_msa: int,
        no_heads_pair: int,
        transition_n: int,
        msa_dropout: float,
        pair_dropout: float,
        chunk_size: int,
        inf: float,
        eps: float,
        _is_extra_msa_stack: bool = False,
    ):
        super(EvoformerBlock, self).__init__()
       
        self.msa_att_row = MSARowAttentionWithPairBias(
            c_m=c_m,
            c_z=c_z,
            c_hidden=c_hidden_msa_att,
            no_heads=no_heads_msa,
            inf=inf,
        )

        if(_is_extra_msa_stack):
            self.msa_att_col = MSAColumnGlobalAttention(
                c_in=c_m,
                c_hidden=c_hidden_msa_att,
                no_heads=no_heads_msa,
                chunk_size=chunk_size,
                inf=inf,
                eps=eps,
            )
        else:
            self.msa_att_col = MSAColumnAttention(
                c_m,
                c_hidden_msa_att,
                no_heads_msa,
                chunk_size=chunk_size,
                inf=inf,
            )

        self.msa_transition = MSATransition(
            c_m=c_m,
            n=transition_n,
            chunk_size=chunk_size,
        )

        self.outer_product_mean = OuterProductMean(
            c_m,
            c_z,
            c_hidden_opm,
            chunk_size=chunk_size,
        )

        self.tri_mul_out = TriangleMultiplicationOutgoing(
            c_z,
            c_hidden_mul,
        )
        self.tri_mul_in = TriangleMultiplicationIncoming(
            c_z,
            c_hidden_mul,
        )

        self.tri_att_start = TriangleAttentionStartingNode(
            c_z,
            c_hidden_pair_att,
            no_heads_pair,
            chunk_size=chunk_size,
            inf=inf,
        )
        self.tri_att_end = TriangleAttentionEndingNode(
            c_z,
            c_hidden_pair_att,
            no_heads_pair,
            chunk_size=chunk_size,
            inf=inf,
        )

        self.pair_transition = PairTransition(
            c_z,
            transition_n,
            chunk_size=chunk_size,
        )
        
        self.msa_dropout_layer = DropoutRowwise(msa_dropout)
        self.ps_dropout_row_layer = DropoutRowwise(pair_dropout)
        self.ps_dropout_col_layer = DropoutColumnwise(pair_dropout)

    def forward(self, 
        m: torch.Tensor, 
        z: torch.Tensor, 
        msa_mask: torch.Tensor, 
        pair_mask: torch.Tensor, 
        _mask_trans: bool = True,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        # DeepMind doesn't mask these transitions in the source, so _mask_trans
        # should be disabled to better approximate the exact activations of
        # the original.
        msa_trans_mask = msa_mask if _mask_trans else None
        pair_trans_mask = pair_mask if _mask_trans else None

        m = m + self.msa_dropout_layer(self.msa_att_row(m, z, mask=msa_mask))
        m = m + self.msa_att_col(m, mask=msa_mask)
        m = m + self.msa_transition(m, mask=msa_trans_mask)
        z = z + self.outer_product_mean(m, mask=msa_mask)
        z = z + self.ps_dropout_row_layer(self.tri_mul_out(z, mask=pair_mask))
        z = z + self.ps_dropout_row_layer(self.tri_mul_in(z, mask=pair_mask))
        z = z + self.ps_dropout_row_layer(self.tri_att_start(z, mask=pair_mask))
        z = z + self.ps_dropout_col_layer(self.tri_att_end(z, mask=pair_mask))
        z = z + self.pair_transition(z, mask=pair_trans_mask)

        return m, z


class EvoformerStack(nn.Module):
    """
        Main Evoformer trunk.

        Implements Algorithm 6.
    """
    def __init__(self,
        c_m: int,
        c_z: int,
        c_hidden_msa_att: int,
        c_hidden_opm: int,
        c_hidden_mul: int,
        c_hidden_pair_att: int,
        c_s: int,
        no_heads_msa: int,
        no_heads_pair: int,
        no_blocks: int,
        transition_n: int,
        msa_dropout: float,
        pair_dropout: float,
        blocks_per_ckpt: int,
        chunk_size: int, 
        inf: float,
        eps: float,
        _is_extra_msa_stack: bool = False,
        **kwargs,
    ):
        """
            Args:
                c_m:
                    MSA channel dimension
                c_z:
                    Pair channel dimension
                c_hidden_msa_att:
                    Hidden dimension in MSA attention
                c_hidden_opm:
                    Hidden dimension in outer product mean module
                c_hidden_mul:
                    Hidden dimension in multiplicative updates
                c_hidden_pair_att:
                    Hidden dimension in triangular attention
                c_s:
                    Channel dimension of the output "single" embedding
                no_heads_msa:
                    Number of heads used for MSA attention
                no_heads_pair:
                    Number of heads used for pair attention
                no_blocks:
                    Number of Evoformer blocks in the stack
                transition_n:
                    Factor by which to multiply c_m to obtain the MSATransition 
                    hidden dimension
                msa_dropout:
                    Dropout rate for MSA activations
                pair_dropout:
                    Dropout used for pair activations
                blocks_per_ckpt:
                    Number of Evoformer blocks in each activation checkpoint
        """
        super(EvoformerStack, self).__init__()

        self.blocks_per_ckpt = blocks_per_ckpt
        self._is_extra_msa_stack = _is_extra_msa_stack

        self.blocks = nn.ModuleList()

        for _ in range(no_blocks):
            block = EvoformerBlock(
                c_m=c_m,
                c_z=c_z,
                c_hidden_msa_att=c_hidden_msa_att,
                c_hidden_opm=c_hidden_opm,
                c_hidden_mul=c_hidden_mul,
                c_hidden_pair_att=c_hidden_pair_att,
                no_heads_msa=no_heads_msa,
                no_heads_pair=no_heads_pair,
                transition_n=transition_n,
                msa_dropout=msa_dropout,
                pair_dropout=pair_dropout,
                chunk_size=chunk_size,
                inf=inf,
                eps=eps,
                _is_extra_msa_stack=_is_extra_msa_stack,
            )
            self.blocks.append(block)

        if(not self._is_extra_msa_stack):
            self.linear = Linear(c_m, c_s)

    def forward(self, 
        m: torch.Tensor, 
        z: torch.Tensor,
        msa_mask: torch.Tensor,
        pair_mask: torch.Tensor,
        _mask_trans: bool = True,
    ) -> Tuple[torch.Tensor, torch.Tensor, Optional[torch.Tensor]]:
        """
            Args:
                m:
                    [*, N_seq, N_res, C_m] MSA embedding
                z:
                    [*, N_res, N_res, C_z] pair embedding
                msa_mask:
                    [*, N_seq, N_res] MSA mask
                pair_mask:
                    [*, N_res, N_res] pair mask
            Returns:
                m:
                    [*, N_seq, N_res, C_m] MSA embedding
                z:
                    [*, N_res, N_res, C_z] pair embedding
                s:
                    [*, N_res, C_s] single embedding
        """
        m, z = checkpoint_blocks(
            blocks=[
                partial(
                    b, 
                    msa_mask=msa_mask, 
                    pair_mask=pair_mask,
                    _mask_trans=_mask_trans,
                ) for b in self.blocks
            ], 
            args=(m, z),
            blocks_per_ckpt=self.blocks_per_ckpt,
        )

        s = None
        if(not self._is_extra_msa_stack):
            seq_dim = -3
            index = torch.tensor([0], device=m.device)
            s = self.linear(torch.index_select(m, dim=seq_dim, index=index))
            s = s.squeeze(seq_dim)

        return m, z, s


class ExtraMSAStack(nn.Module):
    """ 
        Implements Algorithm 18.
    """
    def __init__(self,
        c_m: int,
        c_z: int,
        c_hidden_msa_att: int,
        c_hidden_opm: int,
        c_hidden_mul: int,
        c_hidden_pair_att: int,
        no_heads_msa: int,
        no_heads_pair: int,
        no_blocks: int,
        transition_n: int,
        msa_dropout: float,
        pair_dropout: float,
        blocks_per_ckpt: int,
        chunk_size: int,
        inf: float,
        eps: float,
        **kwargs,
    ):
        super(ExtraMSAStack, self).__init__()

        c_s = None
        self.stack = EvoformerStack(
            c_m=c_m,
            c_z=c_z,
            c_hidden_msa_att=c_hidden_msa_att,
            c_hidden_opm=c_hidden_opm,
            c_hidden_mul=c_hidden_mul,
            c_hidden_pair_att=c_hidden_pair_att,
            c_s=c_s,
            no_heads_msa=no_heads_msa,
            no_heads_pair=no_heads_pair,
            no_blocks=no_blocks,
            transition_n=transition_n,
            msa_dropout=msa_dropout,
            pair_dropout=pair_dropout,
            blocks_per_ckpt=blocks_per_ckpt,
            chunk_size=chunk_size,
            inf=inf,
            eps=eps,
            _is_extra_msa_stack=True, 
        )

    def forward(self, 
        m: torch.Tensor, 
        z: torch.Tensor, 
        msa_mask: Optional[torch.Tensor] = None, 
        pair_mask: Optional[torch.Tensor] = None, 
        _mask_trans: bool = True
    ) -> torch.Tensor:
        """
            Args:
                m:
                    [*, N_extra, N_res, C_m] extra MSA embedding
                z:
                    [*, N_res, N_res, C_z] pair embedding
                msa_mask:
                    Optional [*, N_extra, N_res] MSA mask
                pair_mask:
                    Optional [*, N_res, N_res] pair mask
            Returns:
                [*, N_res, N_res, C_z] pair update
        """ 
        _, z, _ = self.stack(
            m, 
            z,
            msa_mask=msa_mask, 
            pair_mask=pair_mask, 
            _mask_trans=_mask_trans
        )
        return z


if __name__ == "__main__":
    batch_size = 2
    s_t = 3
    n_res = 100
    c_m = 128
    c_z = 64
    c_hidden_att = 32
    c_hidden_opm = 31
    c_hidden_mul = 30
    c_s = 29
    no_heads_msa = 4
    no_heads_pair = 8
    no_blocks = 2
    transition_n = 5
    msa_dropout = 0.15
    pair_dropout = 0.25

    es = EvoformerStack(
        c_m,
        c_z,
        c_hidden_att,
        c_hidden_opm,
        c_hidden_mul,
        c_s,
        no_heads_msa,
        no_heads_pair,
        no_blocks,
        transition_n,
        msa_dropout,
        pair_dropout,
    )

    m = torch.rand((batch_size, s_t, n_res, c_m))
    z = torch.rand((batch_size, n_res, n_res, c_z))

    shape_m_before = m.shape
    shape_z_before = z.shape

    m, z, s = es(m, z)

    assert(m.shape == shape_m_before)
    assert(z.shape == shape_z_before)
    assert(s.shape == (batch_size, n_res, c_s))

    batch_size = 2
    s = 5
    n_res = 100
    c_m = 256
    c = 32
    c_z = 128

    opm = OuterProductMean(c_m, c_z, c)

    m = torch.rand((batch_size, s, n_res, c_m))
    m = opm(m)

    assert(m.shape == (batch_size, n_res, n_res, c_z))
