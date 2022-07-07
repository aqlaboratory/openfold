#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import urllib.request

""" Downloads CAMEO proteins from PDB. The "cameo_table_path" should be a file
    containing a CAMEO target table (exported using the "Copy to clipboard"
    option). Useful for constructing validation sets.
    E.g. https://www.cameo3d.org/modeling/targets/3-months/?to_date=2022-07-02
"""


def main(args):
    with open(args.cameo_table_path, "r") as fp:
        lines = [l.strip() for l in fp.readlines()]

    splits = [l.split() for l in lines]
    prots, chain_ids = zip(*[s[5:7] for s in splits])
    chain_ids = [chain_id.strip('[').strip(']') for chain_id in chain_ids]

    for prot in prots:
        url = f"https://files.rcsb.org/view/{prot.upper()}.cif"
        out_path = os.path.join(args.output_dir, f"{prot}.cif")
        if(not os.path.exists(out_path)):
            urllib.request.urlretrieve(url, out_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("cameo_table_path", type=str)
    parser.add_argument("output_dir", type=str)

    args = parser.parse_args()
    
    main(args)

