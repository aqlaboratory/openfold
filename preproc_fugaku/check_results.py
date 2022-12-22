#!/usr/bin/env python3

import argparse
import os
import re
import glob

import pandas as pd

from foycc import pmap

USE_OK_FILES = True

STATUS_OK         = "OK"
STATUS_ERROR      = "ERR"
STATUS_INCOMPLETE = "INC"
STATUS_NOT_FOUND  = "NOF"


def parse_output(lines):
    def grep_lines(pattern, multiple=False):
        matcher = re.compile(pattern)
        hits = []
        for l in lines:
            m = matcher.match(l)
            if m is not None:
                hits.append(m.groups())

        if not multiple:
            if len(hits) == 0:
                raise Exception(f"Pattern not found: {matcher}")
            elif len(hits) > 1:
                raise Exception(f"Multiple hits found: {hits}")

            return hits[0]

        return hits

    def grep_env(env):
        return grep_lines(f"{env}=(.+)")[0]

    num_nodes   = int(grep_env("NumNodes"))
    num_procs   = int(grep_env("NumProcs"))
    num_threads = int(grep_env("NumThreads"))
    mode        = grep_env("Mode")
    times        = grep_lines(
        ".*:Finished .+ query \(([^)]+)\) in ([^ ]+) seconds",
        multiple=True)

    ret = []
    for chain, time in times:
        ret.append({
            "mode":        mode,
            "chain":       chain,
            "time":        float(time),
            "num_nodes":   num_nodes,
            "num_procs":   num_procs,
            "num_threads": num_threads,
        })

    return ret


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


def check_dir(input_chain):
    path, seq = input_chain

    ret = {"name": os.path.basename(path)}

    exists = os.path.exists(path)
    ret["exists"] = exists

    parse_count_a3m = lambda x: len(parse_a3m(x)[0])
    parse_count_sto = lambda x: len(parse_stockholm(x)[0])
    parse_count_hhr = lambda x: len(parse_hhr(x))

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
    os.environ["OPENFOLD_IGNORE_IMPORT"] = "1"
    from openfold.data.parsers import parse_fasta, parse_a3m, parse_hhr, parse_stockholm

    parser = argparse.ArgumentParser()
    parser.add_argument("mode", help="output or log")
    parser.add_argument("directory", help="The output/log directory")
    parser.add_argument(
        "--fasta",
        default=None,
        help="The input .fasta file")
    parser.add_argument("--weak_scale",
                    dest="weak_scale", action="store_const",
                    const=True, default=False,
                    help="Show the maximum time for each chain")
    args = parser.parse_args()

    if args.mode == "output":
        input_file = args.fasta
        if input_file is not None:
            with open(input_file, 'r') as fp:
                fasta_str = fp.read()

            input_seqs, input_chains = parse_fasta(fasta_str)
            print(f"Loaded {len(input_chains)} from the fast file ({input_file})")

        else:
            input_chains = os.listdir(args.directory)
            input_seqs = [None]*len(input_chains)

        if len(input_chains) != len(set(input_chains)):
            sys.stderr.write("The input chain names must be unique.\n")
            exit()

        input_chains = zip(input_chains, input_seqs)
        input_chains = sorted(input_chains, key=lambda x: x[0])
        input_chains = [(os.path.join(args.directory, x[0]), x[1]) for x in input_chains]

        print(f"Checking {len(input_chains)} directories...")
        results = pmap(check_dir, input_chains)

        df = pd.DataFrame(results)
        print(df)
        df.to_csv("results.csv", index=False)

    elif args.mode == "log":
        # Fugaku
        outputs = glob.glob(os.path.join(args.directory, "script.*.out"))
        if len(outputs) > 0:
            for i, o in enumerate(outputs):
                jobid = int(re.compile(".*/script\.(\d+)\.out").match(o).group(1))
                assert jobid > 0
                stderrs = glob.glob(os.path.join(
                    args.directory, f"output.{jobid}", "*", "*", "stderr.*"))
                if len(stderrs) > 0:
                    outputs[i] = [o]+stderrs
                else:
                    outputs[i] = None

            outputs = [x for x in outputs if x is not None]

        else:
            outputs = [[x] for x in glob.glob(os.path.join(args.directory, "*", "output"))]

        times = []
        for output in outputs:
            lines = []
            for o in output:
                with open(o, "r") as f:
                    lines.extend(f.read().strip().split("\n"))

            times.extend(parse_output(lines))

        if args.weak_scale:
            def get_merge_key(time):
                time = dict(time.items())
                del time["time"]
                chain = re.compile("(.+)_\d+").match(time["chain"]).group(1)
                time["chain"] = chain
                time = sorted(time.items(), key=lambda x: x[0])
                return str(time), chain

            uniq_keys = set()
            for t in times:
                key, chain = get_merge_key(t)
                t["key"] = key
                t["chain"] = chain
                uniq_keys.add(key)

            merged_times = []
            for key in uniq_keys:
                ts = [x for x in times if x["key"] == key]
                time = max([x["time"] for x in ts])
                merged = dict(ts[0].items())
                merged["time"] = time
                merged["count"] = len(ts)
                del merged["key"]
                merged_times.append(merged)

            times = merged_times

        times = pd.DataFrame(times)
        print(times)
        times.to_csv("times.csv", index=False)

    else:
        raise Exception(f"Invalid mode: {args.mode}")
