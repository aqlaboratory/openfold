import torch
from openfold.np import residue_constants as rc
import logging
logger = logging.getLogger(__name__)
import sys

def kabsch_rotation(P, Q):
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """

    # Computation of the covariance matrix
    C = P.transpose(-1, -2) @ Q

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, _, W = torch.linalg.svd(C)
    d = (torch.linalg.det(V) * torch.linalg.det(W)) < 0.0

    if d:
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = V @ W
    return U

def get_optimal_transform(
    src_atoms: torch.Tensor,
    tgt_atoms: torch.Tensor,
    mask: torch.Tensor = None,
):
    assert src_atoms.shape == tgt_atoms.shape, (src_atoms.shape, tgt_atoms.shape)
    assert src_atoms.shape[-1] == 3
    if mask is not None:
        assert mask.dtype == torch.bool
        assert mask.shape[-1] == src_atoms.shape[-2]
        if mask.sum() == 0:
            src_atoms = torch.zeros((1, 3), device=src_atoms.device).float()
            tgt_atoms = src_atoms
        else:
            src_atoms = src_atoms[mask, :]
            tgt_atoms = tgt_atoms[mask, :]
    src_center = src_atoms.mean(-2, keepdim=True)
    tgt_center = tgt_atoms.mean(-2, keepdim=True)
    r = kabsch_rotation(src_atoms - src_center, tgt_atoms - tgt_center)
    x = tgt_center - src_center @ r
    return r, x


def compute_rmsd(
    true_atom_pos: torch.Tensor,
    pred_atom_pos: torch.Tensor,
    atom_mask: torch.Tensor = None,
    eps: float = 1e-6,
) -> torch.Tensor:
    # shape check
    sq_diff = torch.square(true_atom_pos - pred_atom_pos).sum(dim=-1, keepdim=False)
    if atom_mask is not None:
        sq_diff = sq_diff[atom_mask]
    msd = torch.mean(sq_diff)
    msd = torch.nan_to_num(msd, nan=1e8)
    return torch.sqrt(msd + eps)

def kabsch_rmsd(
    true_atom_pos: torch.Tensor,
    pred_atom_pos: torch.Tensor,
    atom_mask: torch.Tensor,
):
    r, x = get_optimal_transform(
        true_atom_pos,
        pred_atom_pos,
        atom_mask,
    )
    aligned_true_atom_pos = true_atom_pos @ r + x
    return compute_rmsd(aligned_true_atom_pos, pred_atom_pos, atom_mask)



def multi_chain_perm_align(out, batch, labels, shuffle_times=2):
    assert isinstance(labels, list)
    ca_idx = rc.atom_order["CA"]
    pred_ca_pos = out["final_atom_positions"][..., ca_idx, :].float()  # [bsz, nres, 3]
    pred_ca_mask = out["final_atom_mask"][..., ca_idx].float()  # [bsz, nres]
    true_ca_poses = [
        l["all_atom_positions"][..., ca_idx, :].float() for l in labels
    ]  # list([nres, 3])
    true_ca_masks = [
        l["all_atom_mask"][..., ca_idx].float() for l in labels
    ]  # list([nres,])

    unique_asym_ids = torch.unique(batch["asym_id"])

    per_asym_residue_index = {}
    for cur_asym_id in unique_asym_ids:
        asym_mask = (batch["asym_id"] == cur_asym_id).bool()
        asym_mask = asym_mask[:,0] # somehow need to adjust the asym_mask shape
        per_asym_residue_index[int(cur_asym_id)] = batch["residue_index"][asym_mask]
    anchor_gt_asym, anchor_pred_asym = get_anchor_candidates(
        batch, per_asym_residue_index, true_ca_masks
    )
    print(f"anchor_gt_asym is {anchor_gt_asym}, anchor_pred_asym is {anchor_pred_asym}")
    anchor_gt_idx = int(anchor_gt_asym) - 1

    best_rmsd = 1e9
    best_labels = None

    unique_entity_ids = torch.unique(batch["entity_id"])
    entity_2_asym_list = {}
    for cur_ent_id in unique_entity_ids:
        ent_mask = batch["entity_id"] == cur_ent_id
        cur_asym_id = torch.unique(batch["asym_id"][ent_mask])
        entity_2_asym_list[int(cur_ent_id)] = cur_asym_id

    for cur_asym_id in anchor_pred_asym:
        asym_mask = (batch["asym_id"] == cur_asym_id).bool()
        anchor_residue_idx = per_asym_residue_index[int(cur_asym_id)]
        print(f"anchor_residue_idx:{anchor_residue_idx},anchor_gt_idx:{anchor_gt_idx}\n")

        anchor_true_pos = true_ca_poses[anchor_gt_idx][anchor_residue_idx]
        asym_mask = asym_mask[:,0] # somehow need to adjust the asym_mask shape
        anchor_pred_pos = pred_ca_pos[asym_mask]
        print(f"true_ca_masks:\n")
        print(true_ca_masks[anchor_gt_idx].bool())
        print(f"pred_ca_mask\n")
        print(pred_ca_mask.bool())
        anchor_true_mask = true_ca_masks[anchor_gt_idx][anchor_residue_idx[:,0]]
        anchor_pred_mask = pred_ca_mask[asym_mask]
        print(f"anchor_true_mask:\n")
        print(anchor_true_mask.shape)
        print(f"anchor_pred_mask:\n")
        print(anchor_pred_mask.shape)
        r, x = get_optimal_transform(
            anchor_true_pos,
            anchor_pred_pos,
            (anchor_true_mask.to('cpu') * anchor_pred_mask.to('cpu')).bool(),
        )

        aligned_true_ca_poses = [ca @ r + x for ca in true_ca_poses]  # apply transforms
        for _ in range(shuffle_times):
            shuffle_idx = torch.randperm(
                unique_asym_ids.shape[0], device=unique_asym_ids.device
            )
            shuffled_asym_ids = unique_asym_ids[shuffle_idx]
            align = greedy_align(
                batch,
                per_asym_residue_index,
                shuffled_asym_ids,
                entity_2_asym_list,
                pred_ca_pos,
                pred_ca_mask,
                aligned_true_ca_poses,
                true_ca_masks,
            )

            merged_labels = merge_labels(
                batch,
                per_asym_residue_index,
                labels,
                align,
            )

            rmsd = kabsch_rmsd(
                merged_labels["all_atom_positions"][..., ca_idx, :] @ r + x,
                pred_ca_pos,
                (pred_ca_mask * merged_labels["all_atom_mask"][..., ca_idx]).bool(),
            )

            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_labels = merged_labels
    return best_labels


def get_anchor_candidates(batch, per_asym_residue_index, true_masks):
    def find_by_num_sym(min_num_sym):
        best_len = -1
        best_gt_asym = None
        asym_ids = torch.unique(batch["asym_id"][batch["num_sym"] == min_num_sym])

        for cur_asym_id in asym_ids:
            assert cur_asym_id > 0
            cur_residue_index = per_asym_residue_index[int(cur_asym_id)]
            j = int(cur_asym_id - 2) # somehow have to change from -1 to -2
            cur_true_mask = true_masks[j][cur_residue_index]
            cur_len = cur_true_mask.sum()
            if cur_len > best_len:
                best_len = cur_len
                best_gt_asym = cur_asym_id
        print(f"finished selected the best anchor\n best_gt_asym is {best_gt_asym} and best_len is {best_len}")
        return best_gt_asym, best_len

    sorted_num_sym = batch["num_sym"][batch["num_sym"] > 0].sort()[0]
    best_gt_asym = None
    best_len = -1
    for cur_num_sym in sorted_num_sym:
        if cur_num_sym <= 0:
            continue
        cur_gt_sym, cur_len = find_by_num_sym(cur_num_sym)
        if cur_len > best_len:
            best_len = cur_len
            best_gt_asym = cur_gt_sym
        if best_len >= 3:
            break
    best_entity = batch["entity_id"][batch["asym_id"] == best_gt_asym][0]
    print(f"best_entity is {best_entity}\n")
    best_pred_asym = torch.unique(batch["asym_id"][batch["entity_id"] == best_entity])
    return best_gt_asym, best_pred_asym


def greedy_align(
    batch,
    per_asym_residue_index,
    unique_asym_ids,
    entity_2_asym_list,
    pred_ca_pos,
    pred_ca_mask,
    true_ca_poses,
    true_ca_masks,
):
    used = [False for _ in range(len(true_ca_poses))]
    align = []
    for cur_asym_id in unique_asym_ids:
        # skip padding
        if cur_asym_id == 0:
            continue
        i = int(cur_asym_id - 1)
        asym_mask = batch["asym_id"] == cur_asym_id
        num_sym = batch["num_sym"][asym_mask][0]
        # don't need to align
        if (num_sym) == 1:
            align.append((i, i))
            assert used[i] == False
            used[i] = True
            continue
        cur_entity_ids = batch["entity_id"][asym_mask][0]
        best_rmsd = 1e10
        best_idx = None
        cur_asym_list = entity_2_asym_list[int(cur_entity_ids)]
        cur_residue_index = per_asym_residue_index[int(cur_asym_id)]
        cur_pred_pos = pred_ca_pos[asym_mask]
        cur_pred_mask = pred_ca_mask[asym_mask]
        for next_asym_id in cur_asym_list:
            if next_asym_id == 0:
                continue
            j = int(next_asym_id - 1)
            if not used[j]:  # posesible candidate
                cropped_pos = true_ca_poses[j][cur_residue_index]
                mask = true_ca_masks[j][cur_residue_index]
                rmsd = compute_rmsd(
                    cropped_pos, cur_pred_pos, (cur_pred_mask * mask).bool()
                )
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_idx = j

        assert best_idx is not None
        used[best_idx] = True
        align.append((i, best_idx))

    return align


def merge_labels(batch, per_asym_residue_index, labels, align):
    """
    batch:
    labels: list of label dicts, each with shape [nk, *]
    align: list of int, such as [2, None, 0, 1], each entry specify the corresponding label of the asym.
    """
    num_res = batch["msa_mask"].shape[-1]
    outs = {}
    for k, v in labels[0].items():
        if k in [
            "resolution",
        ]:
            continue
        cur_out = {}
        for i, j in align:
            label = labels[j][k]
            # to 1-based
            cur_residue_index = per_asym_residue_index[i + 1]
            cur_out[i] = label[cur_residue_index]
        cur_out = [x[1] for x in sorted(cur_out.items())]
        new_v = torch.concat(cur_out, dim=0)
        merged_nres = new_v.shape[0]
        assert (
            merged_nres <= num_res
        ), f"bad merged num res: {merged_nres} > {num_res}. something is wrong."
        if merged_nres < num_res:  # must pad
            pad_dim = new_v.shape[1:]
            pad_v = new_v.new_zeros((num_res - merged_nres, *pad_dim))
            new_v = torch.concat((new_v, pad_v), dim=0)
        outs[k] = new_v
    return outs