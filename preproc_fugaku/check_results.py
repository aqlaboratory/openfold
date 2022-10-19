#!/usr/bin/env python3

import argparse
import os
import pandas as pd

from foycc import pmap
from openfold.data.parsers import parse_a3m, parse_hhr, parse_stockholm


STATUS_OK         = "OK"
STATUS_ERROR      = "ERR"
STATUS_INCOMPLETE = "INC"
STATUS_NOT_FOUND  = "NA"


def check_file(path, parser=None):
    if not os.path.exists(path):
        return STATUS_NOT_FOUND
    elif not os.path.exists(path+"_ok"):
        return STATUS_INCOMPLETE

    if parser is not None:
        try:
            with open(path, "r") as f:
                _ = parser(f.read().strip().plit("\n"))

            return STATUS_OK

        except:
            return STATUS_ERROR

    return None


def check_dir(path):
    ret = {"name": os.path.basename(path)}

    ret["uniref90"] = check_file(
        os.path.join(path, "uniref90_hits.a3m"),
        parser=parse_a3m)
    ret["mgnify"] = check_file(
        os.path.join(path, "mgnify_hits.a3m"),
        parser=parse_a3m)
    ret["bfd"] = check_file(
        os.path.join(path, "bfd_uniclust_hits.a3m"),
        parser=parse_a3m)
    ret["pdb70"] = check_file(
        os.path.join(path, "pdb70_hits.hhr"),
        parser=parse_hhr)
    ret["small_bfd"] = check_file(
        os.path.join(path, "small_bfd_hits.sto"),
        parser=parse_stockholm)

    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="The directory path")
    args = parser.parse_args()

    root = args.path

    dirs = os.listdir(root)

    print(f"Checking {len(dirs)} directories...")
    results = pmap(check_dir, [os.path.join(root, x) for x in dirs])
    df = pd.DataFrame(results)
    print(df)
