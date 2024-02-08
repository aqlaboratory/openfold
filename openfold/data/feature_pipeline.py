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

import copy
from typing import Mapping, Tuple, List, Optional, Dict, Sequence

import ml_collections
import numpy as np
import torch

from openfold.data import input_pipeline, input_pipeline_multimer


FeatureDict = Mapping[str, np.ndarray]
TensorDict = Dict[str, torch.Tensor]


def np_to_tensor_dict(
    np_example: Mapping[str, np.ndarray],
    features: Sequence[str],
) -> TensorDict:
    """Creates dict of tensors from a dict of NumPy arrays.

    Args:
        np_example: A dict of NumPy feature arrays.
        features: A list of strings of feature names to be returned in the dataset.

    Returns:
        A dictionary of features mapping feature names to features. Only the given
        features are returned, all other ones are filtered out.
    """
    # torch generates warnings if feature is already a torch Tensor
    to_tensor = lambda t: torch.tensor(t) if type(t) != torch.Tensor else t.clone().detach()
    tensor_dict = {
        k: to_tensor(v) for k, v in np_example.items() if k in features
    }

    return tensor_dict


def make_data_config(
    config: ml_collections.ConfigDict,
    mode: str,
    num_res: int,
) -> Tuple[ml_collections.ConfigDict, List[str]]:
    cfg = copy.deepcopy(config)
    mode_cfg = cfg[mode]
    with cfg.unlocked():
        if mode_cfg.crop_size is None:
            mode_cfg.crop_size = num_res

    feature_names = cfg.common.unsupervised_features

    # Add seqemb related features if using seqemb mode.
    if cfg.seqemb_mode.enabled:
        feature_names += cfg.common.seqemb_features

    if cfg.common.use_templates:
        feature_names += cfg.common.template_features

    if cfg[mode].supervised:
        feature_names += cfg.supervised.supervised_features

    return cfg, feature_names


def np_example_to_features(
    np_example: FeatureDict,
    config: ml_collections.ConfigDict,
    mode: str,
    is_multimer: bool = False
):
    np_example = dict(np_example)

    seq_length = np_example["seq_length"]
    num_res = int(seq_length[0]) if seq_length.ndim != 0 else int(seq_length)
    cfg, feature_names = make_data_config(config, mode=mode, num_res=num_res)
 
    if "deletion_matrix_int" in np_example:
        np_example["deletion_matrix"] = np_example.pop(
            "deletion_matrix_int"
        ).astype(np.float32)

    tensor_dict = np_to_tensor_dict(
        np_example=np_example, features=feature_names
    )

    with torch.no_grad():
        if is_multimer:
            features = input_pipeline_multimer.process_tensors_from_config(
                tensor_dict,
                cfg.common,
                cfg[mode],
            )
        else:
            features = input_pipeline.process_tensors_from_config(
                tensor_dict,
                cfg.common,
                cfg[mode],
            )

    if mode == "train":
        p = torch.rand(1).item()
        use_clamped_fape_value = float(p < cfg.supervised.clamp_prob)
        features["use_clamped_fape"] = torch.full(
            size=[cfg.common.max_recycling_iters + 1],
            fill_value=use_clamped_fape_value,
            dtype=torch.float32,
        )
    else:
        features["use_clamped_fape"] = torch.full(
            size=[cfg.common.max_recycling_iters + 1],
            fill_value=0.0,
            dtype=torch.float32,
        )

    return {k: v for k, v in features.items()}


class FeaturePipeline:
    def __init__(
        self,
        config: ml_collections.ConfigDict,
    ):
        self.config = config

    def process_features(
        self,
        raw_features: FeatureDict,
        mode: str = "train",
        is_multimer: bool = False,
    ) -> FeatureDict:
        # if(is_multimer and mode != "predict"):
        #     raise ValueError("Multimer mode is not currently trainable")
        
        return np_example_to_features(
            np_example=raw_features,
            config=self.config,
            mode=mode,
            is_multimer=is_multimer,
        )
