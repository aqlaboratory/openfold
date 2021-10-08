import torch

from features import data_transforms


def nonensembled_transform_fns(data_config):
    """Input pipeline data transformers that are not ensembled."""
    common_cfg = data_config.common

    transforms = [
        data_transforms.cast_to_64bit_ints,
        data_transforms.correct_msa_restypes,
        data_transforms.add_distillation_flag(False),
        data_transforms.squeeze_features,
        data_transforms.randomly_replace_msa_with_unknown(0.0),
        data_transforms.make_seq_mask,
        data_transforms.make_msa_mask,
        data_transforms.make_hhblits_profile,
    ]
    if common_cfg.use_templates:
        transforms.extend([
            data_transforms.fix_templates_aatype,
            data_transforms.make_template_mask,
            data_transforms.make_pseudo_beta('template_')
        ])
    transforms.extend([
        data_transforms.make_atom14_masks,
    ])

    return transforms


def ensembled_transform_fns(data_config):
    """Input pipeline data transformers that can be ensembled and averaged."""
    common_cfg = data_config.common
    eval_cfg = data_config.eval
    transforms = []

    if common_cfg.reduce_msa_clusters_by_max_templates:
        pad_msa_clusters = eval_cfg.max_msa_clusters - eval_cfg.max_templates
    else:
        pad_msa_clusters = eval_cfg.max_msa_clusters

    max_msa_clusters = pad_msa_clusters
    max_extra_msa = common_cfg.max_extra_msa

    transforms.append(
        data_transforms.sample_msa(max_msa_clusters, keep_extra=True)
    )

    if 'masked_msa' in common_cfg:
        # Masked MSA should come *before* MSA clustering so that
        # the clustering and full MSA profile do not leak information about
        # the masked locations and secret corrupted locations.
        transforms.append(
            data_transforms.make_masked_msa(common_cfg.masked_msa,
                                            eval_cfg.masked_msa_replace_fraction)
        )

    if common_cfg.msa_cluster_features:
        transforms.append(data_transforms.nearest_neighbor_clusters())
        transforms.append(data_transforms.summarize_clusters())

    # Crop after creating the cluster profiles.
    if max_extra_msa:
        transforms.append(data_transforms.crop_extra_msa(max_extra_msa))
    else:
        transforms.append(data_transforms.delete_extra_msa)

    transforms.append(data_transforms.make_msa_feat())

    crop_feats = dict(eval_cfg.feat)

    if eval_cfg.fixed_size:
        transforms.append(data_transforms.select_feat(list(crop_feats)))

        transforms.append(data_transforms.make_fixed_size(
            crop_feats,
            pad_msa_clusters,
            common_cfg.max_extra_msa,
            eval_cfg.crop_size,
            eval_cfg.max_templates
        ))
    else:
        transforms.append(data_transforms.crop_templates(eval_cfg.max_templates))

    return transforms

def process_tensors_from_config(tensors, data_config):
    """Based on the config, apply filters and transformations to the data."""

    def wrap_ensemble_fn(data, i):
        """Function to be mapped over the ensemble dimension."""
        d = data.copy()
        fns = ensembled_transform_fns(data_config)
        fn = compose(fns)
        d['ensemble_index'] = i
        return fn(d)

    eval_cfg = data_config.eval
    tensors = compose(
        nonensembled_transform_fns(data_config)
    )(tensors)

    tensors_0 = wrap_ensemble_fn(tensors, 0)
    num_ensemble = eval_cfg.num_ensemble
    if data_config.common.resample_msa_in_recycling:
        # Separate batch per ensembling & recycling step.
        num_ensemble *= data_config.common.num_recycle + 1

    if isinstance(num_ensemble, torch.Tensor) or num_ensemble > 1:
        tensors = map_fn(lambda x: wrap_ensemble_fn(tensors, x),
                         torch.arange(num_ensemble))
    else:
        tensors = tree.map_structure(lambda x: x[None], tensors_0)

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
        ensembled_dict[feat] = torch.stack([dict_i[feat] for dict_i in ensembles])
    return ensembled_dict
