#!/usr/bin/env python3

import argparse
import os
import pandas as pd

from foycc import pmap
from openfold.data.parsers import parse_fasta, parse_a3m, parse_hhr, parse_stockholm

USE_OK_FILES = False

STATUS_OK         = "OK"
STATUS_ERROR      = "ERR"
STATUS_INCOMPLETE = "INC"
STATUS_NOT_FOUND  = "NA"


def check_file(path, parse_count=None):
    if not os.path.exists(path):
        return STATUS_NOT_FOUND, None
    elif USE_OK_FILES and not os.path.exists(path+"_ok"):
        return STATUS_INCOMPLETE, None

    if parser is not None:
        try:
            with open(path, "r") as f:
                return STATUS_OK, parse_count(f.read())

        except:
            return STATUS_ERROR, None

    return STATUS_OK, None


def check_dir(path):
    ret = {"name": os.path.basename(path)}

    exists = os.path.exists(path)
    ret["exists"] = exists

    parse_count_a3m = lambda x: len(parse_a3m(x)[0])
    parse_count_sto = lambda x: len(parse_stockholm(x)[0])
    parse_count_hhr = lambda x: len(parse_hhr(x)) # REVIEW

    types = [
        ("uniref90",  parse_count_a3m, "uniref90_hits.a3m"    ),
        ("mgnify",    parse_count_a3m, "mgnify_hits.a3m"      ),
        ("bfd",       parse_count_a3m, "bfd_uniclust_hits.a3m"),
        ("pdb70",     parse_count_hhr, "pdb70_hits.hhr"       ),
        ("small_bfd", parse_count_sto, "small_bfd_hits.sto"   ),
    ]

    for db, parse_count, filename in types:
        status, count = check_file(
            os.path.join(path, filename),
            parse_count=parse_count,
        )
        if status == STATUS_OK:
            assert count is not None
            ret[db] = count
        else:
            assert count is None
            ret[db] = status

    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="The input .fasta file")
    parser.add_argument("output_dir", help="The output directory")
    args = parser.parse_args()

    input_file = args.input_file
    with open(input_file, 'r') as fp:
        fasta_str = fp.read()
        input_seqs, input_chains = parse_fasta(fasta_str)
        input_chains_set = input_chains
    if len(input_chains) != len(input_chains_set):
        sys.stderr.write("The input chain names must be unique.\n")
        exit()

    input_chains = sorted(input_chains_set)

    print(f"Checking {len(input_chains)} directories...")
    results = pmap(check_dir, [os.path.join(args.output_dir, x) for x in input_chains])
    df = pd.DataFrame(results)
    print(df)
    df.to_csv("results.csv", index=False)
