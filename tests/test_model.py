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

import torch
import numpy as np
import unittest
from config import *
from alphafold.model.model import *
from alphafold.utils.utils import my_tree_map
from tests.alphafold.utils.utils import (
    random_template_feats,
    random_extra_msa_feats,
)


class TestModel(unittest.TestCase):
    def test_dry_run(self):
        batch_size = 2
        n_seq = 5
        n_templ = 7
        n_res = 11
        n_extra_seq = 13

        c = model_config("model_1").model
        c.no_cycles = 2  
        c.evoformer_stack.no_blocks = 4 # no need to go overboard here
        c.evoformer_stack.blocks_per_ckpt = None # don't want to set up 
                                                 # deepspeed for this test

        model = AlphaFold(c)

        batch = {}
        tf = torch.randint(
            c.input_embedder.tf_dim - 1, size=(batch_size, n_res)
        )
        batch["target_feat"] = nn.functional.one_hot(
            tf, c.input_embedder.tf_dim).float()
        batch["aatype"] = torch.argmax(batch["target_feat"], dim=-1)
        batch["residue_index"] = torch.arange(n_res)
        batch["msa_feat"] = torch.rand(
            (batch_size, n_seq, n_res, c.input_embedder.msa_dim)
        )
        t_feats = random_template_feats(n_templ, n_res, batch_size=batch_size)
        batch.update({k:torch.tensor(v) for k, v in t_feats.items()})
        extra_feats = random_extra_msa_feats(
            n_extra_seq, n_res, batch_size=batch_size
        )
        batch.update({k:torch.tensor(v) for k, v in extra_feats.items()}) 
        batch["msa_mask"] = torch.randint(
            low=0, high=2, size=(batch_size, n_seq, n_res)
        )
        batch["seq_mask"] = torch.randint(
            low=0, high=2, size=(batch_size, n_res)
        )

        add_recycling_dims = lambda t: (
            t.unsqueeze(-1).expand(*t.shape, c.no_cycles)
        )
        batch = my_tree_map(add_recycling_dims, batch, torch.Tensor)

        with torch.no_grad(): 
            out = model(batch)


if __name__ == "__main__":
    unittest.main()

