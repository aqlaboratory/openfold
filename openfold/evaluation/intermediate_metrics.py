from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict

import numpy as np

# Optional dependencies: these will be used if installed, otherwise we fall back to None.
try:
    import tmscoring  # TM-score implementation
except ImportError:
    tmscoring = None

try:
    import biotite.structure as struc
    import biotite.structure.io as strucio
except ImportError:
    struc = None
    strucio = None


@dataclass
class StructureMetrics:
    """
    Container for structure comparison metrics between a predicted and reference structure.
    """
    tm_score: Optional[float]
    lddt_mean: Optional[float]
    lddt_per_residue: Optional[np.ndarray]
    n_residues: int


def _tm_score_from_pdb(pred_pdb: str, ref_pdb: str) -> Optional[float]:
    """
    Compute TM-score using tmscoring if available.
    Returns None if tmscoring is not installed.
    """
    if tmscoring is None:
        return None

    alignment = tmscoring.TMscoring(pred_pdb, ref_pdb)
    alignment.optimise()
    # Use the built-in scoring API
    return float(alignment.tmscore(**alignment.get_current_values()))


def _lddt_from_pdb(
    pred_pdb: str,
    ref_pdb: str,
    atom_name: str = "CA",
) -> Dict[str, Optional[np.ndarray]]:
    """
    Compute per-residue and mean lDDT using Biotite if available.

    Returns:
        {
            "per_residue": np.ndarray of shape (L,) or None,
            "mean": float or None,
        }
    """
    if struc is None or strucio is None:
        return {"per_residue": None, "mean": None}

    ref = strucio.load_structure(ref_pdb)
    pred = strucio.load_structure(pred_pdb)

    # For now, assume residue mapping is aligned and just filter by atom name.
    ref_ca = ref[ref.atom_name == atom_name]
    pred_ca = pred[pred.atom_name == atom_name]

    # Safety: length mismatch → bail out
    if len(ref_ca) == 0 or len(pred_ca) == 0 or len(ref_ca) != len(pred_ca):
        return {"per_residue": None, "mean": None}

    lddt_per_res = struc.lddt(
        ref_ca,
        pred_ca,
        aggregation="residue",
    )
    lddt_mean = float(lddt_per_res.mean())

    return {"per_residue": lddt_per_res, "mean": lddt_mean}


def evaluate_structures(pred_pdb: str, ref_pdb: str) -> StructureMetrics:
    """
    High-level helper that computes TM-score and lDDT between
    a predicted and reference structure.

    Args:
        pred_pdb: Path to predicted structure PDB file.
        ref_pdb: Path to reference structure PDB file.

    Returns:
        StructureMetrics with TM-score, mean lDDT, per-residue lDDT and residue count.
        If a metric cannot be computed (missing optional deps), its value is None.
    """
    tm = _tm_score_from_pdb(pred_pdb, ref_pdb)
    lddt = _lddt_from_pdb(pred_pdb, ref_pdb)

    per_res = lddt["per_residue"]
    mean_lddt = lddt["mean"]
    n_res = int(per_res.shape[0]) if per_res is not None else 0

    return StructureMetrics(
        tm_score=tm,
        lddt_mean=mean_lddt,
        lddt_per_residue=per_res,
        n_residues=n_res,
    )

