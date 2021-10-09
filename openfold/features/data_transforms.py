import itertools
from functools import reduce

import numpy as np
import torch
from operator import add

from openfold.config import NUM_RES, NUM_EXTRA_SEQ, NUM_TEMPLATES, NUM_MSA_SEQ
from openfold.np import residue_constants

MSA_FEATURE_NAMES = [
    'msa', 'deletion_matrix', 'msa_mask', 'msa_row_mask', 'bert_mask', 'true_msa'
]

def cast_to_64bit_ints(protein):
    # We keep all ints as int64
    for k, v in protein.items():
        if v.dtype == torch.int32:
            protein[k] = v.type(torch.int64)
    return protein

def make_one_hot(x, num_classes):
    x_one_hot = torch.zeros(*x.shape, num_classes)
    x_one_hot.scatter_(-1, x.unsqueeze(-1), 1)
    return x_one_hot

def make_seq_mask(protein):
    protein['seq_mask'] = torch.ones(protein['aatype'].shape, dtype=torch.float32)
    return protein

def make_template_mask(protein):
    protein['template_mask'] = torch.ones(
        protein['template_aatype'].shape[0], dtype=torch.float32
    )
    return protein

def curry1(f):
  """Supply all arguments but the first."""

  def fc(*args, **kwargs):
    return lambda x: f(x, *args, **kwargs)

  return fc

@curry1
def add_distillation_flag(protein, distillation):
    protein['is_distillation'] = torch.tensor(
        float(distillation), dtype=torch.float32
    )
    return protein

def make_all_atom_aatype(protein):
    protein['all_atom_aatype'] = protein['aatype']
    return protein

def fix_templates_aatype(protein):
    # Map one-hot to indices
    num_templates = protein['template_aatype'].shape[0]
    protein['template_aatype'] = torch.argmax(protein['template_aatype'], dim=-1)
    # Map hhsearch-aatype to our aatype.
    new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
    new_order = torch.tensor(
        new_order_list, dtype=torch.int32
    ).expand(num_templates, -1)
    protein['template_aatype'] = torch.gather(
        new_order, 1, index=protein['template_aatype']
    )
    return protein

def correct_msa_restypes(protein):
    """Correct MSA restype to have the same order as residue_constants."""
    new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
    new_order = torch.tensor(
        [new_order_list]*protein['msa'].shape[1], dtype=protein['msa'].dtype
    ).transpose(0,1)
    protein['msa'] = torch.gather(new_order, 0, protein['msa'])

    perm_matrix = np.zeros((22, 22), dtype=np.float32)
    perm_matrix[range(len(new_order_list)), new_order_list] = 1.

    for k in protein:
        if 'profile' in k:
            num_dim = protein[k].shape.as_list()[-1]
            assert num_dim in [20,21,22], (
                'num_dim for %s out of expected range: %s' % (k, num_dim))
            protein[k] = torch.dot(protein[k], perm_matrix[:num_dim, :num_dim])
    return protein

def squeeze_features(protein):
    """Remove singleton and repeated dimensions in protein features."""
    protein['aatype'] = torch.argmax(protein['aatype'], dim=-1)
    for k in [
            'domain_name', 'msa', 'num_alignments', 'seq_length', 'sequence',
            'superfamily', 'deletion_matrix', 'resolution',
            'between_segment_residues', 'residue_index', 'template_all_atom_masks']:
        if k in protein:
            final_dim = protein[k].shape[-1]
            if isinstance(final_dim, int) and final_dim == 1:
                protein[k] = torch.squeeze(protein[k], dim=-1)

    for k in ['seq_length', 'num_alignments']:
        if k in protein:
            protein[k] = protein[k][0]
    return protein

def make_protein_crop_to_size_seed(protein):
    protein['random_crop_to_size_seed'] = torch.distributions.Uniform(
        low=torch.int32, high=torch.int32).sample((2)
    )
    return protein

@curry1
def randomly_replace_msa_with_unknown(protein, replace_proportion):
    """Replace a portion of the MSA with 'X'."""
    msa_mask = (torch.rand(protein['msa'].shape) < replace_proportion)
    x_idx = 20
    gap_idx = 21
    msa_mask = torch.logical_and(msa_mask, protein['msa'] != gap_idx)
    protein['msa'] = torch.where(msa_mask, torch.ones_like(protein['msa'])*x_idx,
                                 protein['msa'])
    aatype_mask = (
        torch.rand(protein['aatype'].shape) < replace_proportion
    )

    protein['aatype'] = torch.where(
        aatype_mask, torch.ones_like(protein['aatype']) * x_idx,
        protein['aatype']
    )
    return protein

@curry1
def sample_msa(protein, max_seq, keep_extra):
    """Sample MSA randomly, remaining sequences are stored are stored as `extra_*`.
    """
    num_seq = protein['msa'].shape[0]
    shuffled = torch.randperm(num_seq-1)+1
    index_order = torch.cat((torch.tensor([0]), shuffled), dim=0)
    num_sel = min(max_seq, num_seq)
    sel_seq, not_sel_seq = torch.split(index_order, [num_sel, num_seq-num_sel])

    for k in MSA_FEATURE_NAMES:
        if k in protein:
            if keep_extra:
                protein['extra_'+k] = torch.index_select(protein[k], 0, not_sel_seq)
            protein[k] = torch.index_select(protein[k], 0, sel_seq)
    return protein

@curry1
def crop_extra_msa(protein, max_extra_msa):
    num_seq = protein['extra_msa'].shape[0]
    num_sel = min(max_extra_msa, num_seq)
    select_indices = torch.randperm(num_seq)[:num_sel]
    for k in MSA_FEATURE_NAMES:
        if 'extra_' + k in protein:
            protein['extra_'+k] = torch.index_select(protein['extra_'+k], 0, select_indices)
    return protein

def delete_extra_msa(protein):
    for k in MSA_FEATURE_NAMES:
        if 'extra_' + k in protein:
            del protein['extra_' + k]
    return protein

# Not used in inference
@curry1
def block_delete_msa(protein, config):
    num_seq = protein['msa'].shape[0]
    block_num_seq = torch.floor(
        torch.tensor(
            num_seq, dtype=torch.float32
        ) * config.msa_fraction_per_block
    ).to(torch.int32)

    if config.randomize_num_blocks:
        nb = torch.distributions.uniform.Uniform(0, config.num_blocks+1).sample()
    else:
        nb = config.num_blocks

    del_block_starts = torch.distributions.Uniform(0, num_seq).sample(nb)
    del_blocks = del_block_starts[:, None] + torch.range(block_num_seq)
    del_blocks = torch.clip(del_blocks, 0, num_seq-1)
    del_indices = torch.unique(torch.sort(torch.reshape(del_blocks, [-1])))[0]

    # Make sure we keep the original sequence
    combined = torch.cat((torch.range(1, num_seq)[None], del_indices[None]))
    uniques, counts = combined.unique(return_counts=True)
    difference = uniques[counts == 1]
    intersection = uniques[counts > 1]
    keep_indices = torch.squeeze(difference, 0)

    for k in MSA_FEATURE_NAMES:
        if k in protein:
            protein[k] = torch.gather(protein[k], keep_indices)

    return protein

@curry1
def nearest_neighbor_clusters(protein, gap_agreement_weight=0.):
    weights = torch.cat([
        torch.ones(21),
        gap_agreement_weight * torch.ones(1),
        torch.zeros(1)
    ], 0)

    # Make agreement score as weighted Hamming distance
    msa_one_hot = make_one_hot(protein['msa'], 23)
    sample_one_hot = (protein['msa_mask'][:,:,None] * msa_one_hot)
    extra_msa_one_hot = make_one_hot(protein['extra_msa'], 23)
    extra_one_hot = (protein['extra_msa_mask'][:,:,None] * extra_msa_one_hot)

    num_seq, num_res, _ = sample_one_hot.shape
    extra_num_seq, _, _ = extra_one_hot.shape

    # Compute tf.einsum('mrc,nrc,c->mn', sample_one_hot, extra_one_hot, weights)
    # in an optimized fashion to avoid possible memory or computation blowup.
    agreement = torch.matmul(
        torch.reshape(extra_one_hot, [extra_num_seq, num_res*23]),
        torch.reshape(
            sample_one_hot * weights, [num_seq, num_res * 23]
        ).transpose(0, 1),
    )

    # Assign each sequence in the extra sequences to the closest MSA sample
    protein['extra_cluster_assignment'] = torch.argmax(agreement, dim=1).to(torch.int64)

    return protein

def unsorted_segment_sum(data, segment_ids, num_segments):
    """
    Computes the sum along segments of a tensor. Analogous to tf.unsorted_segment_sum.

    :param data: A tensor whose segments are to be summed.
    :param segment_ids: The segment indices tensor.
    :param num_segments: The number of segments.
    :return: A tensor of same data type as the data argument.
    """
    # segment_ids.shape should be a prefix of data.shape
    assert all([i in data.shape for i in segment_ids.shape])

    # segment_ids is a 1-D tensor repeat it to have the same shape as data
    if len(segment_ids.shape) == 1:
        s = torch.prod(torch.tensor(data.shape[1:])).long()
        segment_ids = segment_ids.repeat_interleave(s).view(
            segment_ids.shape[0], *data.shape[1:]
        )

    # data.shape and segment_ids.shape should be equal
    assert data.shape == segment_ids.shape

    shape = [num_segments] + list(data.shape[1:])
    tensor = torch.zeros(*shape).scatter_add(0, segment_ids, data.float())
    tensor = tensor.type(data.dtype)
    return tensor

@curry1
def summarize_clusters(protein):
    """Produce profile and deletion_matrix_mean within each cluster."""
    num_seq = protein['msa'].shape[0]
    def csum(x):
        return unsorted_segment_sum(
            x, protein['extra_cluster_assignment'], num_seq
        )

    mask = protein['extra_msa_mask']
    mask_counts = 1e-6 + protein['msa_mask'] + csum(mask)   # Include center

    msa_sum = csum(mask[:, :, None] * make_one_hot(protein['extra_msa'], 23))
    msa_sum += make_one_hot(protein['msa'], 23)  # Original sequence
    protein['cluster_profile'] = msa_sum / mask_counts[:, :, None]

    del msa_sum

    del_sum = csum(mask * protein['extra_deletion_matrix'])
    del_sum += protein['deletion_matrix'] # Original sequence
    protein['cluster_deletion_mean'] = del_sum / mask_counts
    del del_sum

    return protein

def make_msa_mask(protein):
    """Mask features are all ones, but will later be zero-padded."""
    protein['msa_mask'] = torch.ones(protein['msa'].shape, dtype=torch.float32)
    protein['msa_row_mask'] = torch.ones(protein['msa'].shape[0], dtype=torch.float32)
    return protein

def pseudo_beta_fn(aatype, all_atom_positions, all_atom_masks):
    """Create pseudo beta features."""
    is_gly = torch.eq(aatype, residue_constants.restype_order['G'])
    ca_idx = residue_constants.atom_order['CA']
    cb_idx = residue_constants.atom_order['CB']
    pseudo_beta = torch.where(
        torch.tile(is_gly[..., None], [1] * len(is_gly.shape) + [3]),
        all_atom_positions[..., ca_idx, :],
        all_atom_positions[..., cb_idx, :])

    if all_atom_masks is not None:
        pseudo_beta_mask = torch.where(
            is_gly, all_atom_masks[..., ca_idx], all_atom_masks[..., cb_idx])
        return pseudo_beta, pseudo_beta_mask
    else:
        return pseudo_beta

@curry1
def make_pseudo_beta(protein, prefix=''):
    """Create pseudo-beta (alpha for glycine) position and mask."""
    assert prefix in ['', 'template_']
    protein[prefix + 'pseudo_beta'], protein[prefix + 'pseudo_beta_mask'] = (
        pseudo_beta_fn(
            protein['template_aatype' if prefix else 'all_atom_aatype'],
            protein[prefix + 'all_atom_positions'],
            protein['template_all_atom_masks' if prefix else 'all_atom_mask']))
    return protein

@curry1
def add_constant_field(protein, key, value):
    protein[key] = torch.tensor(value)
    return protein

def shaped_categorical(probs, epsilon=1e-10):
    ds = probs.shape
    num_classes = ds[-1]
    distribution = torch.distributions.categorical.Categorical(
        torch.reshape(probs+epsilon,[-1, num_classes])
    )
    counts = distribution.sample()
    return torch.reshape(counts, ds[:-1])

def make_hhblits_profile(protein):
    """Compute the HHblits MSA profile if not already present."""
    if 'hhblits_profile' in protein:
        return protein

    # Compute the profile for every residue (over all MSA sequences).
    msa_one_hot = make_one_hot(protein['msa'], 22)

    protein['hhblits_profile'] = torch.mean(msa_one_hot, dim=0)
    return protein

@curry1
def make_masked_msa(protein, config, replace_fraction):
    """Create data for BERT on raw MSA."""
    # Add a random amino acid uniformly.
    random_aa = torch.tensor([0.05] * 20 + [0., 0.], dtype=torch.float32)

    categorical_probs = (
        config.uniform_prob * random_aa +
        config.profile_prob * protein['hhblits_profile'] +
        config.same_prob * make_one_hot(protein['msa'], 22))

    # Put all remaining probability on [MASK] which is a new column
    pad_shapes = list(reduce(add, [(0, 0) for _ in range(len(categorical_probs.shape))]))
    pad_shapes[1] = 1
    mask_prob = 1. - config.profile_prob - config.same_prob - config.uniform_prob
    assert mask_prob >= 0.
    categorical_probs = torch.nn.functional.pad(
        categorical_probs, pad_shapes, value=mask_prob
    )

    sh = protein['msa'].shape
    mask_position = torch.rand(sh) < replace_fraction

    bert_msa = shaped_categorical(categorical_probs)
    bert_msa = torch.where(mask_position, bert_msa, protein['msa'])

    # Mix real and masked MSA
    protein['bert_mask'] = mask_position.to(torch.float32)
    protein['true_msa'] = protein['msa']
    protein['msa'] = bert_msa

    return protein

@curry1
def make_fixed_size(
    protein, 
    shape_schema, 
    msa_cluster_size, 
    extra_msa_size, 
    num_res=0, 
    num_templates=0
):
    """Guess at the MSA and sequence dimension to make fixed size."""

    pad_size_map = {
        NUM_RES: num_res,
        NUM_MSA_SEQ: msa_cluster_size,
        NUM_EXTRA_SEQ: extra_msa_size,
        NUM_TEMPLATES: num_templates,
    }

    for k, v in protein.items():
        # Don't transfer this to the accelerator.
        if k == 'extra_cluster_assignment':
            continue
        shape = list(v.shape)
        schema = shape_schema[k]
        msd = "Rank mismatch between shape and shape schema for"
        assert len(shape) == len(schema), (
            f'{msg} {k}: {shape} vs {schema}'
        )
        pad_size = [
            pad_size_map.get(s2, None) or s1 for (s1, s2) in zip(shape, schema)
        ]

        padding = [(0, p - v.shape[i]) for i, p in enumerate(pad_size)]
        padding.reverse()
        padding = list(itertools.chain(*padding))
        if padding:
            protein[k] = torch.nn.functional.pad(v, padding)
            protein[k] = torch.reshape(protein[k], pad_size)

    return protein

@curry1
def make_msa_feat(protein):
    """Create and concatenate MSA features."""
    # Whether there is a domain break. Always zero for chains, but keeping for 
    # compatibility with domain datasets.
    has_break = torch.clip(
        protein['between_segment_residues'].to(torch.float32), 0, 1
    )
    aatype_1hot = make_one_hot(protein['aatype'], 21)

    target_feat = [
        torch.unsqueeze(has_break, dim=-1),
        aatype_1hot, # Everyone gets the original sequence.
    ]

    msa_1hot = make_one_hot(protein['msa'], 23)
    has_deletion = torch.clip(protein['deletion_matrix'], 0., 1.)
    deletion_value = torch.atan(protein['deletion_matrix'] / 3.) * (2. / np.pi)

    msa_feat = [
        msa_1hot,
        torch.unsqueeze(has_deletion, dim=-1),
        torch.unsqueeze(deletion_value, dim=-1),
    ]

    if 'cluster_profile' in protein:
        deletion_mean_value = (
            torch.atan(protein['cluster_deletion_mean'] / 3.) * (2. / np.pi)
        )
        msa_feat.extend([protein['cluster_profile'],
                         torch.unsqueeze(deletion_mean_value, dim=-1),
        ])

    if 'extra_deletion_matrix' in protein:
        protein['extra_has_deletion'] = torch.clip(
            protein['extra_deletion_matrix'], 0., 1.
        )
        protein['extra_deletion_value'] = torch.atan(
            protein['extra_deletion_matrix'] / 3.
        ) * (2. / np.pi)

    protein['msa_feat'] = torch.cat(msa_feat, dim=-1)
    protein['target_feat'] = torch.cat(target_feat, dim=-1)
    return protein

@curry1
def select_feat(protein, feature_list):
    return {k: v for k, v in protein.items() if k in feature_list}

@curry1
def crop_templates(protein, max_templates):
    for k, v in protein.items():
        if k.startswith('template_'):
            protein[k] = v[:max_templates]
    return protein

def make_atom14_masks(protein):
    """Construct denser atom positions (14 dimensions instead of 37)."""
    restype_atom14_to_atom37 = []
    restype_atom37_to_atom14 = []
    restype_atom14_mask = []

    for rt in residue_constants.restypes:
        atom_names = residue_constants.restype_name_to_atom14_names[
            residue_constants.restype_1to3[rt]
        ]
        restype_atom14_to_atom37.append([
            (residue_constants.atom_order[name] if name else 0) 
            for name in atom_names
        ])
        atom_name_to_idx14 = {name: i for i, name in enumerate(atom_names)}
        restype_atom37_to_atom14.append([
            (atom_name_to_idx14[name] if name in atom_name_to_idx14 else 0)
            for name in residue_constants.atom_types
        ])
        # Since all 14 atoms are not present in every residue, use this mask to 
        # tell which atom is there in this residue
        restype_atom14_mask.append([(1. if name else 0.) for name in atom_names])

    # Add dummy mapping for restype 'UNK'
    restype_atom14_to_atom37.append([0] * 14)
    restype_atom37_to_atom14.append([0] * 37)
    restype_atom14_to_atom37 = torch.tensor(
        restype_atom14_to_atom37, dtype=torch.int32
    )
    restype_atom37_to_atom14 = torch.tensor(
        restype_atom37_to_atom14, dtype=torch.int32
    )
    restype_atom14_mask = torch.tensor(
        restype_atom14_mask, dtype=torch.float32
    )

    # create the mapping for (residx, atom14) --> atom37, i.e. an array
    # with shape (num_res, 14) containing the atom37 indices for this protein
    residx_atom14_to_atom37 = torch.index_select(
        restype_atom14_to_atom37, 0, protein['aatype']
    )
    residx_atom14_mask = torch.index_select(
        restype_atom14_mask, 0, protein['aatype']
    )

    protein['atom14_atom_exists'] = residx_atom14_mask
    protein['residx_atom14_to_atom37'] = residx_atom14_to_atom37

    # create the gather indices for mapping back
    residx_atom37_to_atom14 = torch.index_select(
        restype_atom37_to_atom14, 0, protein['aatype']
    )
    protein['residx_atom37_to_atom14'] = residx_atom37_to_atom14

    # create the corresponding mask
    restype_atom37_mask = torch.zeros([21, 37], dtype=torch.float32)
    for restype, restype_letter in enumerate(residue_constants.restypes):
        restype_name = residue_constants.restype_1to3[restype_letter]
        atom_names = residue_constants.residue_atoms[restype_name]
        for atom_name in atom_names:
            atom_type = residue_constants.atom_order[atom_name]
            restype_atom37_mask[restype, atom_type] = 1

    residx_atom37_mask = torch.index_select(
        restype_atom37_mask, 0, protein['aatype']
    )
    protein['atom37_atom_exists'] = residx_atom37_mask

    return protein
