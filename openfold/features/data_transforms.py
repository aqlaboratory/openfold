import numpy as np
import torch

from np import residue_constants

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
    protein['template_mask'] = torch.ones(protein['template_aatype'].shape[0], dtype=torch.float32)
    return protein

def curry1(f):
  """Supply all arguments but the first."""

  def fc(*args, **kwargs):
    return lambda x: f(x, *args, **kwargs)

  return fc

@curry1
def add_distillation_flag(protein, distillation):
    protein['is_distillation'] = torch.tensor(float(distillation), dtype=torch.float32)
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
    new_order = torch.tensor(new_order_list, dtype=torch.int32).expand(num_templates, -1)
    protein['template_aatype'] = torch.gather(new_order, 1, index=protein['template_aatype'])
    return protein

def correct_msa_restypes(protein):
    """Correct MSA restype to have the same order as residue_constants."""
    new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
    new_order = torch.tensor([new_order_list]*protein['msa'].shape[1], dtype=protein['msa'].dtype).transpose(0,1)
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
    protein['random_crop_to_size_seed'] = torch.distributions.Uniform(low=torch.int32, high=torch.int32).sample((2))
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

    protein['aatype'] = torch.where(aatype_mask, torch.ones_like(protein['aatype']) * x_idx,
                                    protein['aatype'])
    return protein

@curry1
def sample_msa(protein, max_seq, keep_extra):
    """Sample MSA randomly, remaining sequences are stored are stored as `extra_*`.
    """
    num_seq = protein['msa'].shape[0]
    shuffled = torch.randperm(num_seq-1)+1
    index_order = torch.cat((torch.tensor([0]), shuffled), dim=0)
    num_sel = min(max_seq, num_seq)
    print('sample_msa num_sel', num_sel, ' num_seq', num_seq)
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
    print('select_indices', select_indices)
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
    block_num_seq = torch.floor(torch.tensor(num_seq, dtype=torch.float32) * config.msa_fraction_per_block).to(torch.int32)

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
    print('msa_one_hot shape', msa_one_hot.shape)
    sample_one_hot = (protein['msa_mask'][:,:,None] * msa_one_hot)
    extra_msa_one_hot = make_one_hot(protein['extra_msa'], 23)
    print('extra_msa_one_hot shape', extra_msa_one_hot.shape)
    extra_one_hot = (protein['extra_msa_mask'][:,:,None] * extra_msa_one_hot)

    num_seq, num_res, _ = sample_one_hot.shape
    extra_num_seq, _, _ = extra_one_hot.shape

    # Compute tf.einsum('mrc,nrc,c->mn', sample_one_hot, extra_one_hot, weights)
    # in an optimized fashion to avoid possible memory or computation blowup.
    agreement = torch.matmul(torch.reshape(extra_one_hot, [extra_num_seq, num_res*23]),
                             torch.reshape(sample_one_hot * weights, [num_seq, num_res * 23]).transpose(0, 1),
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
    assert all([i in data.shape for i in segment_ids.shape]), "segment_ids.shape should be a prefix of data.shape"

    # segment_ids is a 1-D tensor repeat it to have the same shape as data
    if len(segment_ids.shape) == 1:
        s = torch.prod(torch.tensor(data.shape[1:])).long()
        segment_ids = segment_ids.repeat_interleave(s).view(segment_ids.shape[0], *data.shape[1:])

    assert data.shape == segment_ids.shape, "data.shape and segment_ids.shape should be equal"

    shape = [num_segments] + list(data.shape[1:])
    tensor = torch.zeros(*shape).scatter_add(0, segment_ids, data.float())
    tensor = tensor.type(data.dtype)
    return tensor

@curry1
def summarize_clusters(protein):
    """Produce profile and deletion_matrix_mean within each cluster."""
    num_seq = protein['msa'].shape[0]
    def csum(x):
        return unsorted_segment_sum(x, protein['extra_cluster_assignment'], num_seq)

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
