# Copyright 2026 AlQuraishi Laboratory
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

import copy
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

from openfold.data import data_pipeline


class TestDataPipelineMultimer(unittest.TestCase):
    def test_duplicate_sequence_mmcif_chains_get_distinct_structural_features(self):
        pipeline = data_pipeline.DataPipelineMultimer(
            monomer_data_pipeline=object()
        )
        mmcif = SimpleNamespace(
            file_id="fake",
            chain_to_seqres={"A": "AC", "B": "AC"},
        )
        processed_chain_ids = []

        def fake_process_single_chain(chain_id, sequence, description, **_):
            processed_chain_ids.append(chain_id)
            num_res = len(sequence)
            aatype = np.zeros((num_res, 21), dtype=np.float32)
            aatype[np.arange(num_res), np.arange(num_res) % 21] = 1.0
            return {
                "aatype": aatype,
                "sequence": np.array([sequence.encode("utf-8")], dtype=object),
                "domain_name": np.array(
                    [description.encode("utf-8")], dtype=object
                ),
                "num_alignments": np.array(
                    [1] * num_res, dtype=np.int32
                ),
                "seq_length": np.array([num_res] * num_res, dtype=np.int32),
            }

        atom_positions = {
            "A": np.full((2, 37, 3), 1.0, dtype=np.float32),
            "B": np.full((2, 37, 3), 2.0, dtype=np.float32),
        }
        atom_masks = {
            "A": np.ones((2, 37), dtype=np.float32),
            "B": np.zeros((2, 37), dtype=np.float32),
        }
        requested_chain_ids = []

        def fake_get_mmcif_features(mmcif_object, chain_id):
            requested_chain_ids.append(chain_id)
            return {
                "all_atom_positions": atom_positions[chain_id],
                "all_atom_mask": atom_masks[chain_id],
                "resolution": np.array(1.0, dtype=np.float32),
                "release_date": np.array([b"2026-05-10"], dtype=object),
                "is_distillation": np.array(0.0, dtype=np.float32),
            }

        captured = {}

        def fake_pair_and_merge(all_chain_features):
            captured["all_chain_features"] = copy.deepcopy(all_chain_features)
            return {"msa": np.zeros((1, 1), dtype=np.int32)}

        pipeline._process_single_chain = fake_process_single_chain
        pipeline.get_mmcif_features = fake_get_mmcif_features

        with mock.patch.object(
            data_pipeline.feature_processing_multimer,
            "pair_and_merge",
            side_effect=fake_pair_and_merge,
        ), mock.patch.object(
            data_pipeline,
            "pad_msa",
            side_effect=lambda np_example, _: np_example,
        ):
            pipeline.process_mmcif(mmcif=mmcif, alignment_dir="/unused")

        self.assertEqual(processed_chain_ids, ["fake_A"])
        self.assertEqual(requested_chain_ids, ["A", "B"])

        all_chain_features = captured["all_chain_features"]
        np.testing.assert_array_equal(
            all_chain_features["A_1"]["all_atom_positions"],
            atom_positions["A"],
        )
        np.testing.assert_array_equal(
            all_chain_features["A_2"]["all_atom_positions"],
            atom_positions["B"],
        )
        np.testing.assert_array_equal(
            all_chain_features["A_1"]["all_atom_mask"],
            atom_masks["A"],
        )
        np.testing.assert_array_equal(
            all_chain_features["A_2"]["all_atom_mask"],
            atom_masks["B"],
        )
        self.assertFalse(
            np.array_equal(
                all_chain_features["A_1"]["all_atom_positions"],
                all_chain_features["A_2"]["all_atom_positions"],
            )
        )
        self.assertEqual(
            all_chain_features["A_2"]["auth_chain_id"].item(),
            "fake_B",
        )


if __name__ == "__main__":
    unittest.main()
