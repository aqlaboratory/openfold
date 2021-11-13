from typing import Optional, Sequence

import torch
import torch.nn as nn


def script_submodules_(
    model: nn.Module,
    types: Optional[Sequence[type]] = None,
):
    """
        Convert all submodules whose types match one of those in the input 
        list to recursively scripted equivalents in place. To script the entire
        model, just call torch.jit.script on it directly.

        Args:
            model: A torch.nn.Module
            types: A list of types of submodules to script
    """
    for name, child in model.named_children():
        if(any(isinstance(child, t) for t in types)):
            setattr(model, name, torch.jit.script(child))
        else:
            script_submodules_(child, types)
