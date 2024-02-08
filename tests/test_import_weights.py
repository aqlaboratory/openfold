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

import os
import torch
import numpy as np
import unittest
from pathlib import Path

from tests.config import consts
from openfold.config import model_config
from openfold.model.model import AlphaFold
from openfold.utils.import_weights import import_jax_weights_, import_openfold_weights_


class TestImportWeights(unittest.TestCase):
    def test_import_jax_weights_(self):
        npz_path = Path(__file__).parent.resolve() / f"../openfold/resources/params/params_{consts.model}.npz"

        c = model_config(consts.model)
        c.globals.blocks_per_ckpt = None
        model = AlphaFold(c)
        model.eval()

        import_jax_weights_(
            model,
            npz_path,
            version=consts.model
        )

        data = np.load(npz_path)
        prefix = "alphafold/alphafold_iteration/"

        test_pairs = [
            # Normal linear weight
            (
                torch.as_tensor(
                    data[
                        prefix + "structure_module/initial_projection//weights"
                    ]
                ).transpose(-1, -2),
                model.structure_module.linear_in.weight,
            ),
            # Normal layer norm param
            (
                torch.as_tensor(
                    data[prefix + "evoformer/prev_pair_norm//offset"],
                ),
                model.recycling_embedder.layer_norm_z.bias,
            ),
            # From a stack
            (
                torch.as_tensor(
                    data[
                        prefix
                        + (
                            "evoformer/evoformer_iteration/outer_product_mean/"
                            "left_projection//weights"
                        )
                    ][1].transpose(-1, -2)
                ),
                model.evoformer.blocks[1].outer_product_mean.linear_1.weight,
            ),
        ]

        for w_alpha, w_repro in test_pairs:
            self.assertTrue(torch.all(w_alpha == w_repro))

    def test_import_openfold_weights_(self):
        model_name = 'initial_training'
        pt_path = Path(__file__).parent.resolve() / f"../openfold/resources/openfold_params/{model_name}.pt"

        if os.path.exists(pt_path):
            c = model_config(model_name)
            c.globals.blocks_per_ckpt = None
            model = AlphaFold(c)
            model.eval()

            d = torch.load(pt_path)

            import_openfold_weights_(
                model=model,
                state_dict=d,
            )
