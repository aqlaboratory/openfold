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

from typing import Optional, Sequence, Tuple

import torch
import torch.nn as nn

from openfold.model.dropout import (
    DropoutRowwise,
    DropoutColumnwise,
)
from openfold.model.evoformer import (
    EvoformerBlock,
    EvoformerStack,
)
from openfold.model.outer_product_mean import OuterProductMean
from openfold.model.msa import (
    MSARowAttentionWithPairBias, 
    MSAColumnAttention,
    MSAColumnGlobalAttention,
)
from openfold.model.pair_transition import PairTransition
from openfold.model.primitives import Attention, GlobalAttention
from openfold.model.structure_module import (
    InvariantPointAttention,
    BackboneUpdate,
)
from openfold.model.template import TemplatePairStackBlock
from openfold.model.triangular_attention import (
    TriangleAttentionStartingNode,
    TriangleAttentionEndingNode,
)
from openfold.model.triangular_multiplicative_update import (
    TriangleMultiplicationOutgoing,
    TriangleMultiplicationIncoming,
)


def script_preset_(model: torch.nn.Module):
    """
    TorchScript a handful of low-level but frequently used submodule types 
    that are known to be scriptable.

    Args:
        model: 
            A torch.nn.Module. It should contain at least some modules from 
            this repository, or this function won't do anything.
    """
    script_submodules_(
        model, 
        [
            nn.Dropout,
            Attention,
            GlobalAttention,
            EvoformerBlock,
            #TemplatePairStackBlock,
        ], 
        attempt_trace=False,
        batch_dims=None,
    ) 

    
def _get_module_device(module: torch.nn.Module) -> torch.device:
    """
    Fetches the device of a module, assuming that all of the module's
    parameters reside on a single device

    Args:
        module: A torch.nn.Module
    Returns:
        The module's device
    """
    return next(module.parameters()).device


def _trace_module(module, batch_dims=None):
    if(batch_dims is None):
        batch_dims = ()

    # Stand-in values
    n_seq = 10
    n_res = 10

    device = _get_module_device(module)

    def msa(channel_dim):
        return torch.rand(
            (*batch_dims, n_seq, n_res, channel_dim),
            device=device,
        )

    def pair(channel_dim):
        return torch.rand(
            (*batch_dims, n_res, n_res, channel_dim),
            device=device,
        )

    if(isinstance(module, MSARowAttentionWithPairBias)):
        inputs = {
            "forward": (
                msa(module.c_in), # m
                pair(module.c_z), # z
                torch.randint(
                    0, 2, 
                    (*batch_dims, n_seq, n_res)
                ), # mask
            ),
        }
    elif(isinstance(module, MSAColumnAttention)):
        inputs = {
            "forward": (
                msa(module.c_in), # m
                torch.randint(
                    0, 2, 
                    (*batch_dims, n_seq, n_res)
                ), # mask
            ),
        }
    elif(isinstance(module, OuterProductMean)):
        inputs = {
            "forward": (
                msa(module.c_m),
                torch.randint(
                    0, 2,
                    (*batch_dims, n_seq, n_res)
                )
            )
        }
    else:
        raise TypeError(
            f"tracing is not supported for modules of type {type(module)}"
        )

    return torch.jit.trace_module(module, inputs)


def _script_submodules_helper_(
    model,
    types,
    attempt_trace,
    to_trace,
):
    for name, child in model.named_children():
        if(types is None or any(isinstance(child, t) for t in types)):
            try:
                scripted = torch.jit.script(child)
                setattr(model, name, scripted)
                continue
            except (RuntimeError, torch.jit.frontend.NotSupportedError) as e:
                if(attempt_trace):
                    to_trace.add(type(child))
                else:
                    raise e
        
        _script_submodules_helper_(child, types, attempt_trace, to_trace)


def _trace_submodules_(
    model,
    types,
    batch_dims=None,
):
    for name, child in model.named_children():
        if(any(isinstance(child, t) for t in types)):
            traced = _trace_module(child, batch_dims=batch_dims)
            setattr(model, name, traced)
        else:
            _trace_submodules_(child, types, batch_dims=batch_dims)


def script_submodules_(
    model: nn.Module,
    types: Optional[Sequence[type]] = None,
    attempt_trace: Optional[bool] = True,
    batch_dims: Optional[Tuple[int]] = None,
):
    """
    Convert all submodules whose types match one of those in the input 
    list to recursively scripted equivalents in place. To script the entire
    model, just call torch.jit.script on it directly.

    When types is None, all submodules are scripted.

    Args:
        model: 
            A torch.nn.Module
        types: 
            A list of types of submodules to script
        attempt_trace: 
            Whether to attempt to trace specified modules if scripting 
            fails. Recall that tracing eliminates all conditional 
            logic---with great tracing comes the mild responsibility of 
            having to remember to ensure that the modules in question 
            perform the same computations no matter what.
    """
    to_trace = set()

    # Aggressively script as much as possible first...
    _script_submodules_helper_(model, types, attempt_trace, to_trace)
  
    # ... and then trace stragglers.
    if(attempt_trace and len(to_trace) > 0):
        _trace_submodules_(model, to_trace, batch_dims=batch_dims)
