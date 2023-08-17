from typing import Sequence

import torch

from openfold.data.data_transforms import curry1
from openfold.utils.tensor_utils import masked_mean


def gumbel_noise(
    shape: Sequence[int], 
    device: torch.device, 
    eps=1e-6,
    generator=None,
) -> torch.Tensor:
    """Generate Gumbel Noise of given Shape.

    This generates samples from Gumbel(0, 1).

    Args:
        shape: Shape of noise to return.

    Returns:
        Gumbel noise of given shape.
    """
    uniform_noise = torch.rand(
        shape, dtype=torch.float32, device=device, generator=generator
    )
    gumbel = -torch.log(-torch.log(uniform_noise + eps) + eps)
    return gumbel


def gumbel_max_sample(logits: torch.Tensor, generator=None) -> torch.Tensor:
    """Samples from a probability distribution given by 'logits'.

    This uses Gumbel-max trick to implement the sampling in an efficient manner.

    Args:
        logits: Logarithm of probabilities to sample from, probabilities can be
            unnormalized.

    Returns:
        Sample from logprobs in one-hot form.
    """
    z = gumbel_noise(logits.shape, device=logits.device, generator=generator)
    return torch.nn.functional.one_hot(
        torch.argmax(logits + z, dim=-1),
        logits.shape[-1],
    )


def gumbel_argsort_sample_idx(
    logits: torch.Tensor, 
    generator=None
) -> torch.Tensor:
    """Samples with replacement from a distribution given by 'logits'.

    This uses Gumbel trick to implement the sampling an efficient manner. For a
    distribution over k items this samples k times without replacement, so this
    is effectively sampling a random permutation with probabilities over the
    permutations derived from the logprobs.

    Args:
        logits: Logarithm of probabilities to sample from, probabilities can be
            unnormalized.

    Returns:
        Sample from logprobs in one-hot form.
    """
    z = gumbel_noise(logits.shape, device=logits.device, generator=generator)
    return torch.argsort(logits + z, dim=-1, descending=True)


@curry1
def make_masked_msa(batch, config, replace_fraction, seed, eps=1e-6):
    """Create data for BERT on raw MSA."""
    # Add a random amino acid uniformly.
    random_aa = torch.Tensor(
        [0.05] * 20 + [0., 0.], 
        device=batch['msa'].device
    )

    categorical_probs = (
        config.uniform_prob * random_aa +
        config.profile_prob * batch['msa_profile'] +
        config.same_prob * torch.nn.functional.one_hot(batch['msa'], 22)
    )

    # Put all remaining probability on [MASK] which is a new column.
    mask_prob = 1. - config.profile_prob - config.same_prob - config.uniform_prob
   
    categorical_probs = torch.nn.functional.pad(
        categorical_probs, [0,1], value=mask_prob
    )

    sh = batch['msa'].shape
    mask_position = torch.rand(sh, device=batch['msa'].device) < replace_fraction
    mask_position *= batch['msa_mask'].to(mask_position.dtype)
    
    logits = torch.log(categorical_probs + eps)
    
    g = torch.Generator(device=batch["msa"].device)
    if seed is not None:
        g.manual_seed(seed)

    bert_msa = gumbel_max_sample(logits, generator=g)
    
    bert_msa = torch.where(
        mask_position,
        torch.argmax(bert_msa, dim=-1), 
        batch['msa']
    )
    bert_msa *= batch['msa_mask'].to(bert_msa.dtype)

    # Mix real and masked MSA.
    if 'bert_mask' in batch:
        batch['bert_mask'] *= mask_position.to(torch.float32)
    else:
        batch['bert_mask'] = mask_position.to(torch.float32)
    batch['true_msa'] = batch['msa']
    batch['msa'] = bert_msa

    return batch


@curry1
def nearest_neighbor_clusters(batch, gap_agreement_weight=0.):
    """Assign each extra MSA sequence to its nearest neighbor in sampled MSA."""
    device = batch["msa_mask"].device

    # Determine how much weight we assign to each agreement.    In theory, we could
    # use a full blosum matrix here, but right now let's just down-weight gap
    # agreement because it could be spurious.
    # Never put weight on agreeing on BERT mask.

    weights = torch.Tensor(
        [1.] * 21 + [gap_agreement_weight] + [0.], 
        device=device,
    )

    msa_mask = batch['msa_mask']
    msa_one_hot = torch.nn.functional.one_hot(batch['msa'], 23)

    extra_mask = batch['extra_msa_mask']
    extra_one_hot = torch.nn.functional.one_hot(batch['extra_msa'], 23)

    msa_one_hot_masked = msa_mask[:, :, None] * msa_one_hot
    extra_one_hot_masked = extra_mask[:, :, None] * extra_one_hot

    agreement = torch.einsum(
        'mrc, nrc->nm', 
        extra_one_hot_masked,
        weights * msa_one_hot_masked
    )

    cluster_assignment = torch.nn.functional.softmax(1e3 * agreement, dim=0)
    cluster_assignment *= torch.einsum('mr, nr->mn', msa_mask, extra_mask)

    cluster_count = torch.sum(cluster_assignment, dim=-1)
    cluster_count += 1.    # We always include the sequence itself.

    msa_sum = torch.einsum('nm, mrc->nrc', cluster_assignment, extra_one_hot_masked)
    msa_sum += msa_one_hot_masked

    cluster_profile = msa_sum / cluster_count[:, None, None]

    extra_deletion_matrix = batch['extra_deletion_matrix']
    deletion_matrix = batch['deletion_matrix']

    del_sum = torch.einsum(
        'nm, mc->nc', 
        cluster_assignment,
        extra_mask * extra_deletion_matrix
    )
    del_sum += deletion_matrix    # Original sequence.
    cluster_deletion_mean = del_sum / cluster_count[:, None]

    batch['cluster_profile'] = cluster_profile
    batch['cluster_deletion_mean'] = cluster_deletion_mean

    return batch


def create_target_feat(batch):
    """Create the target features"""
    batch["target_feat"] = torch.nn.functional.one_hot(
        batch["aatype"], 21
    ).to(torch.float32)
    return batch


def create_msa_feat(batch):
    """Create and concatenate MSA features."""
    device = batch["msa"]
    msa_1hot = torch.nn.functional.one_hot(batch['msa'], 23)
    deletion_matrix = batch['deletion_matrix']
    has_deletion = torch.clamp(deletion_matrix, min=0., max=1.)[..., None]
    pi = torch.acos(torch.zeros(1, device=deletion_matrix.device)) * 2
    deletion_value = (torch.atan(deletion_matrix / 3.) * (2. / pi))[..., None]

    deletion_mean_value = (
        torch.atan(
            batch['cluster_deletion_mean'] / 3.) *
            (2. / pi)
    )[..., None]

    msa_feat = torch.cat(
        [
            msa_1hot,
            has_deletion,
            deletion_value,
            batch['cluster_profile'],
            deletion_mean_value
        ], 
        dim=-1,
    )

    batch["msa_feat"] = msa_feat

    return batch 


def build_extra_msa_feat(batch):
    """Expand extra_msa into 1hot and concat with other extra msa features.

    We do this as late as possible as the one_hot extra msa can be very large.

    Args:
        batch: a dictionary with the following keys:
         * 'extra_msa': [num_seq, num_res] MSA that wasn't selected as a cluster
             centre. Note - This isn't one-hotted.
         * 'extra_deletion_matrix': [num_seq, num_res] Number of deletions at given
                position.
        num_extra_msa: Number of extra msa to use.

    Returns:
        Concatenated tensor of extra MSA features.
    """
    # 23 = 20 amino acids + 'X' for unknown + gap + bert mask
    extra_msa = batch['extra_msa']
    deletion_matrix = batch['extra_deletion_matrix']
    msa_1hot = torch.nn.functional.one_hot(extra_msa, 23)
    has_deletion = torch.clamp(deletion_matrix, min=0., max=1.)[..., None]
    pi = torch.acos(torch.zeros(1, device=deletion_matrix.device)) * 2
    deletion_value = (
        (torch.atan(deletion_matrix / 3.) * (2. / pi))[..., None]
    )
    extra_msa_mask = batch['extra_msa_mask']
    catted = torch.cat([msa_1hot, has_deletion, deletion_value], dim=-1)
  
    return catted


@curry1
def sample_msa(batch, max_seq, max_extra_msa_seq, seed, inf=1e6):
    """Sample MSA randomly, remaining sequences are stored as `extra_*`.

    Args:
        batch: batch to sample msa from.
        max_seq: number of sequences to sample.
    Returns:
        Protein with sampled msa.
    """
    g = torch.Generator(device=batch["msa"].device)
    if seed is not None:
        g.manual_seed(seed)

    # Sample uniformly among sequences with at least one non-masked position.
    logits = (torch.clamp(torch.sum(batch['msa_mask'], dim=-1), 0., 1.) - 1.) * inf
    # The cluster_bias_mask can be used to preserve the first row (target
    # sequence) for each chain, for example.
    if 'cluster_bias_mask' not in batch:
        cluster_bias_mask = torch.nn.functional.pad(
            batch['msa'].new_zeros(batch['msa'].shape[0] - 1), 
            (1, 0), 
            value=1.
        )
    else:
        cluster_bias_mask = batch['cluster_bias_mask']

    logits += cluster_bias_mask * inf
    index_order = gumbel_argsort_sample_idx(logits, generator=g)
    sel_idx = index_order[:max_seq]
    extra_idx = index_order[max_seq:][:max_extra_msa_seq]

    for k in ['msa', 'deletion_matrix', 'msa_mask', 'bert_mask']:
        if k in batch:
            batch['extra_' + k] = batch[k][extra_idx]
            batch[k] = batch[k][sel_idx]

    return batch


def make_msa_profile(batch):
    """Compute the MSA profile."""

    # Compute the profile for every residue (over all MSA sequences).
    batch["msa_profile"] = masked_mean(
        batch['msa_mask'][..., None], 
        torch.nn.functional.one_hot(batch['msa'], 22), 
        dim=-3,
    )

    return batch
