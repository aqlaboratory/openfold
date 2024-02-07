import logging
import random
import torch

from openfold.np import residue_constants as rc

logger = logging.getLogger(__name__)


def compute_rmsd(
    true_atom_pos: torch.Tensor,
    pred_atom_pos: torch.Tensor,
    atom_mask: torch.Tensor = None,
    eps: float = 1e-6,
) -> torch.Tensor:
    sq_diff = torch.square(true_atom_pos - pred_atom_pos).sum(dim=-1, keepdim=False)
    if atom_mask is not None:
        sq_diff = torch.masked_select(sq_diff, atom_mask.to(sq_diff.device))
    msd = torch.mean(sq_diff)
    msd = torch.nan_to_num(msd, nan=1e8)
    return torch.sqrt(msd + eps)  # prevent sqrt 0


def kabsch_rotation(P, Q):
    """
    Calculate the best rotation that minimises the RMSD between P and Q.

    The optimal rotation matrix was calculated using Kabsch algorithm:
    https://en.wikipedia.org/wiki/Kabsch_algorithm

    Args:
    P: [N * 3] Nres is the number of atoms and each row corresponds to the atom's x,y,z coordinates
    Q: [N * 3] the same dimension as P

    return:
    A 3*3 rotation matrix
    """
    assert P.shape == torch.Size([Q.shape[0], Q.shape[1]])

    # Firstly, compute SVD of P.T * Q
    u, _, vt = torch.linalg.svd(torch.matmul(P.to(torch.float32).T,
                                             Q.to(torch.float32)))
    # Then construct s matrix
    s = torch.eye(P.shape[1], device=P.device)
    # correct the rotation matrix to ensure a right-handed coordinate
    s[-1, -1] = torch.sign(torch.linalg.det(torch.matmul(u, vt)))
    # finally compute the rotation matrix
    r_opt = torch.matmul(torch.matmul(u, s), vt)
    assert r_opt.shape == torch.Size([3,3])
    return r_opt.to(device=P.device, dtype=P.dtype)


def get_optimal_transform(
    src_atoms: torch.Tensor,
    tgt_atoms: torch.Tensor,
    mask: torch.Tensor = None,
):
    """
    src_atoms: predicted CA positions, shape:[num_res,3]
    tgt_atoms: ground-truth CA positions, shape:[num_res,3]
    mask: a vector of boolean values, shape:[num_res]
    """
    assert src_atoms.shape == tgt_atoms.shape, (src_atoms.shape, tgt_atoms.shape)
    assert src_atoms.shape[-1] == 3
    if mask is not None:
        assert len(mask.shape) == 1, "mask should have the shape of [num_res]"
    if torch.isnan(src_atoms).any() or torch.isinf(src_atoms).any():
        #
        # sometimes using fake test inputs generates NaN in the predicted atom positions
        # #
        logging.warning(f"src_atom has nan or inf")
        src_atoms = torch.nan_to_num(src_atoms, nan=0.0, posinf=1.0, neginf=1.0)

    if mask is not None:
        assert mask.dtype == torch.bool
        assert mask.shape[-1] == src_atoms.shape[-2]
        if mask.sum() == 0:
            src_atoms = torch.zeros((1, 3), device=src_atoms.device, dtype=src_atoms.dtype)
            tgt_atoms = src_atoms
        else:
            src_atoms = src_atoms[mask, :]
            tgt_atoms = tgt_atoms[mask, :]
    src_center = src_atoms.mean(-2, keepdim=True, dtype=src_atoms.dtype)
    tgt_center = tgt_atoms.mean(-2, keepdim=True, dtype=src_atoms.dtype)
    r = kabsch_rotation(src_atoms, tgt_atoms)

    x = tgt_center - src_center @ r
    return r, x


def get_least_asym_entity_or_longest_length(batch, input_asym_id):
    """
    First check how many subunit(s) one sequence has. Select the subunit that is less
    common, e.g. if the protein was AABBB then select one of the A as anchor

    If there is a tie, e.g. AABB, first check which sequence is the longer/longest,
    then choose one of the corresponding subunits as anchor

    Args:
    batch: in this function batch is the full ground truth features
    input_asym_id: A list of asym_ids that are in the cropped input features

    Return:
    anchor_gt_asym_id: Tensor(int) selected ground truth asym_id
    anchor_pred_asym_ids: list(Tensor(int)) a list of all possible pred anchor candidates
    """
    entity_2_asym_list = get_entity_2_asym_list(batch)
    unique_entity_ids = torch.unique(batch["entity_id"])
    entity_asym_count = {}
    entity_length = {}

    for entity_id in unique_entity_ids:
        asym_ids = torch.unique(batch["asym_id"][batch["entity_id"] == entity_id])

        # Make sure some asym IDs associated with ground truth entity ID exist in cropped prediction
        asym_ids_in_pred = [a for a in asym_ids if a in input_asym_id]
        if not asym_ids_in_pred:
            continue

        entity_asym_count[int(entity_id)] = len(asym_ids)

        # Calculate entity length
        entity_mask = (batch["entity_id"] == entity_id)
        entity_length[int(entity_id)] = entity_mask.sum().item()

    min_asym_count = min(entity_asym_count.values())
    least_asym_entities = [entity for entity, count in entity_asym_count.items() if count == min_asym_count]

    # If multiple entities have the least asym_id count, return those with the longest length
    if len(least_asym_entities) > 1:
        max_length = max([entity_length[entity] for entity in least_asym_entities])
        least_asym_entities = [entity for entity in least_asym_entities if entity_length[entity] == max_length]

    # If still multiple entities, return a random one
    if len(least_asym_entities) > 1:
        least_asym_entities = [random.choice(least_asym_entities)]

    assert len(least_asym_entities) == 1
    least_asym_entities = least_asym_entities[0]

    anchor_gt_asym_id = random.choice(entity_2_asym_list[least_asym_entities])
    anchor_pred_asym_ids = [asym_id for asym_id in entity_2_asym_list[least_asym_entities] if asym_id in input_asym_id]

    return anchor_gt_asym_id, anchor_pred_asym_ids


def greedy_align(
    batch,
    per_asym_residue_index,
    entity_2_asym_list,
    pred_ca_pos,
    pred_ca_mask,
    true_ca_poses,
    true_ca_masks,
):
    """
    Implement Algorithm 4 in the Supplementary Information of AlphaFold-Multimer paper:
    Evans,R et al., 2022 Protein complex prediction with AlphaFold-Multimer, bioRxiv 2021.10.04.463034; doi: https://doi.org/10.1101/2021.10.04.463034
    """
    used = [False for _ in range(len(true_ca_poses))]
    align = []
    unique_asym_ids = [i for i in torch.unique(batch["asym_id"]) if i != 0]
    for cur_asym_id in unique_asym_ids:
        i = int(cur_asym_id - 1)
        asym_mask = batch["asym_id"] == cur_asym_id
        cur_entity_ids = batch["entity_id"][asym_mask][0]
        best_rmsd = torch.inf
        best_idx = None
        cur_asym_list = entity_2_asym_list[int(cur_entity_ids)]
        cur_residue_index = per_asym_residue_index[int(cur_asym_id)]
        cur_pred_pos = pred_ca_pos[asym_mask]
        cur_pred_mask = pred_ca_mask[asym_mask]
        for next_asym_id in cur_asym_list:
            j = int(next_asym_id - 1)
            if not used[j]:  # possible candidate
                cropped_pos = torch.index_select(true_ca_poses[j], 1, cur_residue_index)
                mask = torch.index_select(true_ca_masks[j], 1, cur_residue_index)
                rmsd = compute_rmsd(
                    torch.squeeze(cropped_pos, 0), torch.squeeze(cur_pred_pos, 0),
                    (cur_pred_mask * mask).bool()
                )
                if rmsd is not None and rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_idx = j

        assert best_idx is not None
        used[best_idx] = True
        align.append((i, best_idx))
    return align


def pad_features(feature_tensor, nres_pad, pad_dim):
    """Pad input feature tensor"""
    pad_shape = list(feature_tensor.shape)
    pad_shape[pad_dim] = nres_pad
    padding_tensor = feature_tensor.new_zeros(pad_shape, device=feature_tensor.device)
    return torch.concat((feature_tensor, padding_tensor), dim=pad_dim)


def merge_labels(per_asym_residue_index, labels, align, original_nres):
    """
    Merge ground truth labels according to the permutation results

    labels: list of original ground truth feats
    align: list of tuples, each entry specify the corresponding label of the asym.

    modified based on UniFold:
    https://github.com/dptech-corp/Uni-Fold/blob/b1c89a2cebd4e4ee4c47b4e443f92beeb9138fbb/unifold/losses/chain_align.py#L176C1-L176C1
    """
    outs = {}
    for k, v in labels[0].items():
        cur_out = {}
        for i, j in align:
            label = labels[j][k]
            # to 1-based
            cur_residue_index = per_asym_residue_index[i + 1]
            if len(v.shape) <= 1 or "template" in k or "row_mask" in k:
                continue
            else:
                dimension_to_merge = 1
                cur_out[i] = label.index_select(dimension_to_merge, cur_residue_index)
        cur_out = [x[1] for x in sorted(cur_out.items())]
        if len(cur_out) > 0:
            new_v = torch.concat(cur_out, dim=dimension_to_merge)
            # below check whether padding is needed
            if new_v.shape[dimension_to_merge] != original_nres:
                nres_pad = original_nres - new_v.shape[dimension_to_merge]
                new_v = pad_features(new_v, nres_pad, pad_dim=dimension_to_merge)
            outs[k] = new_v
    return outs


def split_ground_truth_labels(gt_features):
    """
    Splits ground truth features according to chains

    Returns:
    a list of feature dictionaries with only necessary ground truth features
    required to finish multi-chain permutation
    """
    unique_asym_ids, asym_id_counts = torch.unique(gt_features["asym_id"], sorted=True, return_counts=True)
    n_res = gt_features["asym_id"].shape[-1]

    def split_dim(shape):
        return next(iter(i for i, size in enumerate(shape) if size == n_res), None)

    labels = list(map(dict, zip(*[[(k, v) for v in torch.split(v_all, asym_id_counts.tolist(),
                                                               dim=split_dim(v_all.shape))]
                                  for k, v_all in gt_features.items()
                                  if n_res in v_all.shape])))
    return labels


def get_per_asym_residue_index(features):
    unique_asym_ids = [i for i in torch.unique(features["asym_id"]) if i != 0]
    per_asym_residue_index = {}
    for cur_asym_id in unique_asym_ids:
        asym_mask = (features["asym_id"] == cur_asym_id).bool()
        per_asym_residue_index[int(cur_asym_id)] = torch.masked_select(features["residue_index"], asym_mask)

    return per_asym_residue_index


def get_entity_2_asym_list(batch):
    """
    Generates a dictionary mapping unique entity IDs to lists of unique asymmetry IDs (asym_id) for each entity.

    Args:
        batch (dict): A dictionary containing data batches, including "entity_id" and "asym_id" tensors.

    Returns:
        entity_2_asym_list (dict): A dictionary where keys are unique entity IDs, and values are lists of unique asymmetry IDs
        associated with each entity.
    """
    entity_2_asym_list = {}
    unique_entity_ids = torch.unique(batch["entity_id"])
    for cur_ent_id in unique_entity_ids:
        ent_mask = batch["entity_id"] == cur_ent_id
        cur_asym_id = torch.unique(batch["asym_id"][ent_mask])
        entity_2_asym_list[int(cur_ent_id)] = cur_asym_id
    return entity_2_asym_list


def calculate_input_mask(true_ca_masks, anchor_gt_idx, anchor_gt_residue,
                         asym_mask, pred_ca_mask):
    """
    Calculate an input mask for downstream optimal transformation computation

    Args:
        true_ca_masks (Tensor): ca mask from ground truth.
        anchor_gt_idx (Tensor): The index of selected ground truth anchor.
        asym_mask (Tensor): Boolean tensor indicating which regions are selected predicted anchor.
        pred_ca_mask (Tensor): ca mask from predicted structure.

    Returns:
        input_mask (Tensor): A boolean mask
    """
    pred_ca_mask = torch.squeeze(pred_ca_mask, 0)
    asym_mask = torch.squeeze(asym_mask, 0)
    anchor_pred_mask = pred_ca_mask[asym_mask]
    anchor_true_mask = torch.index_select(true_ca_masks[anchor_gt_idx], 1, anchor_gt_residue)
    input_mask = (anchor_true_mask * anchor_pred_mask).bool()
    return input_mask


def calculate_optimal_transform(true_ca_poses,
                                anchor_gt_idx, anchor_gt_residue,
                                true_ca_masks, pred_ca_mask,
                                asym_mask,
                                pred_ca_pos):
    input_mask = calculate_input_mask(true_ca_masks,
                                      anchor_gt_idx,
                                      anchor_gt_residue,
                                      asym_mask,
                                      pred_ca_mask)
    input_mask = torch.squeeze(input_mask, 0)
    pred_ca_pos = torch.squeeze(pred_ca_pos, 0)
    asym_mask = torch.squeeze(asym_mask, 0)
    anchor_true_pos = torch.index_select(true_ca_poses[anchor_gt_idx], 1, anchor_gt_residue)
    anchor_pred_pos = pred_ca_pos[asym_mask]
    r, x = get_optimal_transform(
        anchor_pred_pos, torch.squeeze(anchor_true_pos, 0),
        mask=input_mask
    )

    return r, x


def compute_permutation_alignment(out, features, ground_truth):
    """
    A class method that first permutate chains in ground truth first
    before calculating the loss.

    Details are described in Section 7.3 in the Supplementary of AlphaFold-Multimer paper:
    https://www.biorxiv.org/content/10.1101/2021.10.04.463034v2
    """
    unique_asym_ids = set(torch.unique(features['asym_id']).tolist())
    unique_asym_ids.discard(0)  # Remove padding asym_id
    is_monomer = len(unique_asym_ids) == 1

    per_asym_residue_index = get_per_asym_residue_index(features)

    if is_monomer:
        best_align = list(enumerate(range(len(per_asym_residue_index))))
        return best_align, per_asym_residue_index

    best_rmsd = float('inf')
    best_align = None
    # First select anchors from predicted structures and ground truths
    anchor_gt_asym, anchor_pred_asym_ids = get_least_asym_entity_or_longest_length(ground_truth,
                                                                                   features['asym_id'])
    entity_2_asym_list = get_entity_2_asym_list(ground_truth)
    labels = split_ground_truth_labels(ground_truth)
    assert isinstance(labels, list)
    anchor_gt_idx = int(anchor_gt_asym) - 1
    # Then calculate optimal transform by aligning anchors
    ca_idx = rc.atom_order["CA"]
    pred_ca_pos = out["final_atom_positions"][..., ca_idx, :]  # [bsz, nres, 3]
    pred_ca_mask = out["final_atom_mask"][..., ca_idx].to(dtype=pred_ca_pos.dtype)  # [bsz, nres]

    true_ca_poses = [
        l["all_atom_positions"][..., ca_idx, :] for l in labels
    ]  # list([nres, 3])
    true_ca_masks = [
        l["all_atom_mask"][..., ca_idx].long() for l in labels
    ]  # list([nres,])
    for candidate_pred_anchor in anchor_pred_asym_ids:
        asym_mask = (features["asym_id"] == candidate_pred_anchor).bool()
        anchor_gt_residue = per_asym_residue_index[candidate_pred_anchor.item()]
        r, x = calculate_optimal_transform(true_ca_poses,
                                           anchor_gt_idx,
                                           anchor_gt_residue,
                                           true_ca_masks,
                                           pred_ca_mask,
                                           asym_mask,
                                           pred_ca_pos)
        aligned_true_ca_poses = [ca.to(r.dtype) @ r + x for ca in true_ca_poses]  # apply transforms
        align = greedy_align(
            features,
            per_asym_residue_index,
            entity_2_asym_list,
            pred_ca_pos,
            pred_ca_mask,
            aligned_true_ca_poses,
            true_ca_masks,
        )
        merged_labels = merge_labels(per_asym_residue_index, labels, align,
                                     original_nres=features['aatype'].shape[-1])
        rmsd = compute_rmsd(
            true_atom_pos=merged_labels['all_atom_positions'][..., ca_idx, :].to(r.dtype) @ r + x,
            pred_atom_pos=pred_ca_pos,
            atom_mask=(pred_ca_mask * merged_labels['all_atom_mask'][..., ca_idx].long()).bool())
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_align = align

    return best_align, per_asym_residue_index


def multi_chain_permutation_align(out, features, ground_truth):
    """Compute multi-chain permutation alignment.

    Args:
        out: The output of model.forward()
        features: Input features
        ground_truth: Ground truth features
    """

    labels = split_ground_truth_labels(ground_truth)

    # Then permute ground truth chains before calculating the loss
    align, per_asym_residue_index = compute_permutation_alignment(out=out,
                                                                  features=features,
                                                                  ground_truth=ground_truth)

    # reorder ground truth labels according to permutation results
    labels = merge_labels(per_asym_residue_index, labels, align,
                          original_nres=features['aatype'].shape[-1])

    features.update(labels)
    return features
