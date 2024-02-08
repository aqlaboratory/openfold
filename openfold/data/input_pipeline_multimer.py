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

import random
import torch

from openfold.data import (
    data_transforms,
    data_transforms_multimer,
)


def groundtruth_transforms_fns():
    transforms = [data_transforms.make_atom14_masks,
                  data_transforms.make_atom14_positions,
                  data_transforms.atom37_to_frames,
                  data_transforms.atom37_to_torsion_angles(""),
                  data_transforms.make_pseudo_beta(""),
                  data_transforms.get_backbone_frames,
                  data_transforms.get_chi_angles]
    return transforms


def nonensembled_transform_fns():
    """Input pipeline data transformers that are not ensembled."""
    transforms = [
        data_transforms.cast_to_64bit_ints,
        data_transforms_multimer.make_msa_profile,
        data_transforms_multimer.create_target_feat,
        data_transforms.make_atom14_masks
    ]

    return transforms


def ensembled_transform_fns(common_cfg, mode_cfg, ensemble_seed):
    """Input pipeline data transformers that can be ensembled and averaged."""
    transforms = []

    pad_msa_clusters = mode_cfg.max_msa_clusters
    max_msa_clusters = pad_msa_clusters
    max_extra_msa = mode_cfg.max_extra_msa

    msa_seed = None
    if(not common_cfg.resample_msa_in_recycling):
        msa_seed = ensemble_seed
    
    transforms.append(
        data_transforms_multimer.sample_msa(
            max_msa_clusters, 
            max_extra_msa,
            seed=msa_seed,
        )
    )

    if "masked_msa" in common_cfg:
        # Masked MSA should come *before* MSA clustering so that
        # the clustering and full MSA profile do not leak information about
        # the masked locations and secret corrupted locations.
        transforms.append(
            data_transforms_multimer.make_masked_msa(
                common_cfg.masked_msa, 
                mode_cfg.masked_msa_replace_fraction,
                seed=(msa_seed + 1) if msa_seed else None,
            )
        )

    transforms.append(data_transforms_multimer.nearest_neighbor_clusters())
    transforms.append(data_transforms_multimer.create_msa_feat)

    crop_feats = dict(common_cfg.feat)

    if mode_cfg.fixed_size:
        transforms.append(data_transforms.select_feat(list(crop_feats)))

        if mode_cfg.crop:
            transforms.append(
                data_transforms_multimer.random_crop_to_size(
                    crop_size=mode_cfg.crop_size,
                    max_templates=mode_cfg.max_templates,
                    shape_schema=crop_feats,
                    spatial_crop_prob=mode_cfg.spatial_crop_prob,
                    interface_threshold=mode_cfg.interface_threshold,
                    subsample_templates=mode_cfg.subsample_templates,
                    seed=ensemble_seed + 1,
                )
            )
        transforms.append(
            data_transforms.make_fixed_size(
                shape_schema=crop_feats,
                msa_cluster_size=pad_msa_clusters,
                extra_msa_size=mode_cfg.max_extra_msa,
                num_res=mode_cfg.crop_size,
                num_templates=mode_cfg.max_templates,
            )
        )
    else:
        transforms.append(
            data_transforms.crop_templates(mode_cfg.max_templates)
        )

    return transforms


def prepare_ground_truth_features(tensors):
    """Prepare ground truth features that are only needed for loss calculation during training"""

    gt_features = ['all_atom_mask', 'all_atom_positions', 'asym_id', 'sym_id', 'entity_id']
    gt_tensors = {k: v for k, v in tensors.items() if k in gt_features}
    gt_tensors['aatype'] = tensors['aatype'].to(torch.long)
    gt_tensors = compose(groundtruth_transforms_fns())(gt_tensors)
    return gt_tensors


def process_tensors_from_config(tensors, common_cfg, mode_cfg):
    """Based on the config, apply filters and transformations to the data."""

    process_gt_feats = mode_cfg.supervised
    gt_tensors = {}
    if process_gt_feats:
        gt_tensors = prepare_ground_truth_features(tensors)

    ensemble_seed = random.randint(0, torch.iinfo(torch.int32).max)
    tensors['aatype'] = tensors['aatype'].to(torch.long)
    nonensembled = nonensembled_transform_fns()
    tensors = compose(nonensembled)(tensors)
    if("no_recycling_iters" in tensors):
        num_recycling = int(tensors["no_recycling_iters"])
    else:
        num_recycling = common_cfg.max_recycling_iters

    def wrap_ensemble_fn(data, i):
        """Function to be mapped over the ensemble dimension."""
        d = data.copy()
        fns = ensembled_transform_fns(
            common_cfg, 
            mode_cfg, 
            ensemble_seed,
        )
        fn = compose(fns)
        d["ensemble_index"] = i
        return fn(d)

    tensors = map_fn(
        lambda x: wrap_ensemble_fn(tensors, x), torch.arange(num_recycling + 1)
    )

    if process_gt_feats:
        tensors['gt_features'] = gt_tensors

    return tensors

@data_transforms.curry1
def compose(x, fs):
    for f in fs:
        x = f(x)
    return x


def map_fn(fun, x):
    ensembles = [fun(elem) for elem in x]
    features = ensembles[0].keys()
    ensembled_dict = {}
    for feat in features:
        ensembled_dict[feat] = torch.stack(
            [dict_i[feat] for dict_i in ensembles], dim=-1
        )
    return ensembled_dict
