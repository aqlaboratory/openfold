#!/usr/bin/env python

import argparse
from pathlib import Path

from openfold.evaluation.intermediate_metrics import evaluate_structures


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate predicted vs reference structures with TM-score and lDDT."
    )
    parser.add_argument("--pred-pdb", required=True, help="Path to predicted PDB")
    parser.add_argument("--ref-pdb", required=True, help="Path to reference PDB")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("metrics.txt"),
        help="Where to save a simple text summary",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    metrics = evaluate_structures(args.pred_pdb, args.ref_pdb)

    lines = [
        f"Predicted: {args.pred_pdb}",
        f"Reference: {args.ref_pdb}",
        "",
        f"TM-score: {metrics.tm_score if metrics.tm_score is not None else 'N/A'}",
        f"mean lDDT: {metrics.lddt_mean if metrics.lddt_mean is not None else 'N/A'}",
        f"n_residues: {metrics.n_residues}",
    ]

    args.out.write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

