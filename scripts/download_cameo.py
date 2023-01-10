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
        "https://www.cameo3d.org",
        "modeling",
        period,
        "ajax",
        f"?to_date={end_date}",
    ])


def main(args):
    data_dir_path = os.path.join(args.output_dir, "data_dir")
    fasta_dir_path = os.path.join(args.output_dir, "fasta_dir")
    pdb_dir_path = os.path.join(args.output_dir, "pdb_dir")
    os.makedirs(data_dir_path, exist_ok=True)
    os.makedirs(fasta_dir_path, exist_ok=True)
    os.makedirs(pdb_dir_path, exist_ok=True)

    url = generate_url(args.period, args.end_date)
    # raw_data = requests.get(url).text
    # parsed_data = json.loads(raw_data)
    # chain_data = parsed_data["aaData"]

    raw_data = ["7RPY:A", "7PB4:I", "7FJS:L", "7U2R:A", "7AAL:A", "7KO9:A", "7RKC:A", "7EQH:A", "7U4H:A", "7ETS:B", "7E4S:A", "7F3A:A", "7LXK:A", "7EAD:A", "7DKI:A", "7ESH:A", "7VNO:A", "7Z79:B", "7RW4:A", "7W5M:A", "7WNW:B", "7OPB:D", "7EQE:A", "7N0E:A", "7T4Z:A", "7ED6:A", "7NUV:A", "7TV9:C", "7ZCL:B", "7VWT:A", "7PC6:A", "7NQD:B", "7TXP:A", "7PC4:A", "7QRY:B", "7FEV:A", "7FIW:B", "7RPR:A", "7OA7:A", "7EBQ:A", "7YWG:B", "7UGH:A", "7F0A:A", "7U5Y:A", "7NDE:A", "7QIL:A", "7X4E:A", "7OD9:C", "7TBU:A", "7W26:A", "7X4O:B", "7NMI:B", "7WRK:A", "7QSS:A", "7LI0:A", "7RGV:A", "7VSP:C", "7X8J:A", "7QDV:A", "7E52:A", "7RCW:A", "7TNI:C", "7PC3:A", "7N29:B", "7F2Y:A", "7ZGF:A", "7T03:A", "7MYV:B", "7BLL:A", "7MQ4:A", "7X9E:A", "7F6J:C", "7EJG:C", "7V4S:A", "7QAP:A", "7ACY:B", "7MLA:B", "7QAO:A", "7WWR:A", "7QSU:A", "7PZJ:A", "7V1K:A", "7SGN:C", "7Z5P:A", "7N3T:C", "7EGT:B", "7O4O:A", "7CTX:B", "7VNX:A", "7YXG:A", "7QS5:A", "7X8V:A", "7MKU:A", "7RPS:A", "7MS2:A", "7QBP:A", "7QS2:A", "7EHG:E", "7PRQ:B", "7S2R:B", "7R74:B", "7W5U:A", "7O0B:A", "7TA5:A", "7WWX:A", "7Q51:A", "7SXB:A", "7WCJ:A", "7LXS:A", "7OSW:A", "7WRP:A", "7PNO:D", "7WJ9:A", "7RCZ:A", "7U5F:D", "7WME:A", "7RXE:A", "7B0K:A", "7ERN:C", "7R63:C", "7SJL:A", "7BI4:A", "7W1F:B", "7ED1:A", "7RAW:A", "7SO5:H", "7VOH:A", "7Q05:E", "7QBZ:A", "7PSG:C", "7P0H:A", "7MHW:A", "7P3I:B", "7ULH:A", "7R09:A", "7F0O:B", "7EQB:A", "7EFS:D", "7TZE:C", "7W74:A", "7PXY:A", "7PW1:A", "7E5J:A", "7V8E:B", "7ERP:B", "7R5Y:E"]
    chain_data = []
    for r in raw_data:
        chain_data.append({"pdbid": r[:4], "pdbid_chain": r[5:]})
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
        # with open(os.path.join(pdb_dir_path, f"{pdb_id}.pdb"), "w") as fp:
        #     fp.write(pdb_file)


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
