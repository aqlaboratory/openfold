#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import json
import os
import re
import requests

from openfold.data import mmcif_parsing


VALID_PERIODS = [
    "1-year",
    "6-months",
    "3-months",
    "1-month",
    "1-week",
]


def generate_url(period, end_date):
    return '/'.join([
        "https://www.cameo3d.org/",
        "modeling",
        "targets",
        period,
        "ajax",
        f"?to_date={end_date}",
    ])


def main(args):
    data_dir_path = os.path.join(args.output_dir, "data_dir")
    fasta_dir_path = os.path.join(args.output_dir, "fasta_dir")
    
    os.makedirs(data_dir_path, exist_ok=True)
    os.makedirs(fasta_dir_path, exist_ok=True)

    url = generate_url(args.period, args.end_date)
    raw_data = requests.get(url).text
    parsed_data = json.loads(raw_data)

    chain_data = parsed_data["aaData"]
    for chain in chain_data:
        pdb_id = chain["pdbid"]
        chain_id = chain["pdbid_chain"]

        pdb_url = f"https://files.rcsb.org/view/{pdb_id.upper()}.cif"
        pdb_file = requests.get(pdb_url).text

        parsed_cif = mmcif_parsing.parse(
            file_id=pdb_id, mmcif_string=pdb_file
        )
        mmcif_object = parsed_cif.mmcif_object
        if(mmcif_object is None):
            raise list(parsed_cif.errors.values())[0]

        seq = mmcif_object.chain_to_seqres[chain_id]

        if(args.max_seqlen > 0):
            if(len(seq) > len(seq)):
                continue

        fasta_file = '\n'.join([
            f">{pdb_id}_{chain_id}",
            seq,
        ])

        fasta_filename = f"{pdb_id}_{chain_id}.fasta"
        with open(os.path.join(fasta_dir_path, fasta_filename), "w") as fp:
            fp.write(fasta_file)

        cif_filename = f"{pdb_id}.cif"
        with open(os.path.join(data_dir_path, cif_filename), "w") as fp:
            fp.write(pdb_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "period", type=str,
        help=f"""The length of the period from which to draw CAMEO proteins. 
             Choose from {VALID_PERIODS}"""
    )
    parser.add_argument(
        "end_date", type=str,
        help="The date marking the end of the period (YYYY-MM-DD)"
    )
    parser.add_argument("output_dir")
    parser.add_argument(
        "--max_seqlen", type=int, default=700,
        help="The maximum length in residues of downloaded proteins (or -1)"
    )

    args = parser.parse_args()

    if(args.period not in VALID_PERIODS):
        raise ValueError(f"Invalid period. Choose from {VALID_PERIODS}")

    date_regex = re.compile("^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
    if(not date_regex.match(args.end_date)):
        raise ValueError(f"Invalid end_date: {args.end_date}. Use YYYY-MM-DD format")

    main(args)
