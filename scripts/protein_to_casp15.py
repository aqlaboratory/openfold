from np import residue_constants
from np.protein import Protein


def protein_to_casp15(prot: Protein, target_name: str, template_names: str) -> str:
    """Convert a `Protein` object into the CASP15 submission format.

    Args:
        prot: The protein object

    Returns:
        CASP15 submission string
    """
    restypes = residue_constants.restypes + ["X"]
    res_1to3 = lambda r: residue_constants.restype_1to3.get(restypes[r], "UNK")
    atom_types = residue_constants.atom_types

    atom_mask = prot.atom_mask
    aatype = prot.aatype
    atom_positions = prot.atom_positions
    residue_index = prot.residue_index.astype(np.int32)
    b_factors = prot.b_factors

    lines = []
    lines.append("PFRMAT TS")
    lines.append(f"TARGET {target_name}")
    lines.append("AUTHOR {Openfold group registration code}")
    lines.append("REMARK ")
    lines.append("METHOD ")

    lines.append("MODEL 1")
    lines.append(f"PARENT \t {template_names}")
    atom_index = 1
    chain_id = 1
    for i in range(aatype.shape[0]):
        res_name_3 = res_1to3(aatype[i])
        for atom_name, pos, mask, b_factor in zip(
                atom_types, atom_positions[i], atom_mask[i], b_factors[i]):
            if mask < 0.5:
                continue

            record_type = "ATOM"
            name = atom_name if len(atom_name) == 4 else f" {atom_name}"
            alt_loc = ""
            insertion_code = ""
            occupancy = 1.00
            element = atom_name[0]
            charge = ""
            atom_line = (
                f"{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}"
                f"{res_name_3:>3} {chain_id:>1}"
                f"{residue_index[i]:>4}{insertion_code:>1}   "
                f"{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}"
                f"{occupancy:>6.2f}{b_factor:>6.2f}          "
                f"{element:>2}{charge:>2}"
            )
            lines.append(atom_line)
            atom_index += 1

    lines.append("TER")
    lines.append("END")
    lines.append("")
    return "\n".join(lines)
