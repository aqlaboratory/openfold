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

import deepspeed
import torch
from typing import Any, Tuple, List, Callable

BLOCK_ARG = Any
BLOCK_ARGS = Tuple[BLOCK_ARG, ...]

def checkpoint_blocks(
    blocks: List[Callable[BLOCK_ARGS, BLOCK_ARGS]], 
    args: BLOCK_ARGS, 
    blocks_per_ckpt: int,
) -> BLOCK_ARGS:
    """
        Chunk a list of blocks and run each chunk with activation 
        checkpointing. We define a "block" as a callable whose only inputs are 
        the outputs of the previous block.

        This function assumes that deepspeed has already been initialized.

        Implements Subsection 1.11.8

        Args:
            blocks:
                List of blocks
            args:
                Tuple of arguments for the first block.
            blocks_per_ckpt:
                Size of each chunk. A higher value corresponds to higher memory
                consumption but fewer checkpoints. If None, no checkpointing is
                performed.
        Returns:
            The output of the final block
    """

    def wrap(a):
        return (a,) if type(a) is not tuple else a

    def exec(b, a):
        for block in b:
            a = wrap(block(*a))
        return a

    def chunker(s, e):
        def exec_sliced(a):
            return exec(blocks[s:e], a)
        return exec_sliced

    # Avoids mishaps when the blocks take just one argument
    args = wrap(args)

    if(blocks_per_ckpt is None):
        return exec(blocks, args)
    elif(blocks_per_ckpt < 1 or blocks_per_ckpt > len(blocks)):
        raise ValueError("blocks_per_ckpt must be between 1 and len(blocks)")

    for s in range(0, len(blocks), blocks_per_ckpt):
        e = s + blocks_per_ckpt
        args = deepspeed.checkpointing.checkpoint(chunker(s, e), args)

    return args
