import copy

import ml_collections
import torch
from typing import Mapping, Tuple, List, Optional, Dict, Sequence
import numpy as np

from features import input_pipeline

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
    tensor_dict = {k: torch.tensor(v) for k, v in np_example.items() if k in features}
    return tensor_dict

def make_data_config(
        config: ml_collections.ConfigDict,
        num_res: int,
        ) -> Tuple[ml_collections.ConfigDict, List[str]]:
    cfg = copy.deepcopy(config.data)

    feature_names = cfg.common.unsupervised_features
    if cfg.common.use_templates:
        feature_names += cfg.common.template_features

    with cfg.unlocked():
        cfg.eval.crop_size = num_res

    return cfg, feature_names

def np_example_to_features(np_example: FeatureDict,
                           config: ml_collections.ConfigDict,
                           random_seed: int = 0):
    np_example = dict(np_example)
    num_res = int(np_example['seq_length'][0])
    cfg, feature_names = make_data_config(config, num_res=num_res)

    if 'deletion_matrix_int' in np_example:
        np_example['deletion_matrix'] = (
            np_example.pop('deletion_matrix_int').astype(np.float32))

    torch.manual_seed(random_seed)
    tensor_dict = np_to_tensor_dict(
        np_example=np_example, features=feature_names)
    features = input_pipeline.process_tensors_from_config(tensor_dict, cfg)

    return {k: v for k, v in features.items()}


class FeaturePipeline:
    def __init__(self,
                 config: ml_collections.ConfigDict,
                 params: Optional[Mapping[str, Mapping[str, np.ndarray]]] = None):
        self.config = config
        self.params = params

    def process_features(self,
                         raw_features: FeatureDict,
                         random_seed: int) -> FeatureDict:
        return np_example_to_features(
            np_example=raw_features,
            config=self.config,
            random_seed=random_seed
        )