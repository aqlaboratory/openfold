#!/usr/bin/env python3

import os
import re
import argparse
from datetime import datetime
from typing import List, Dict, Any

try:
    import pandas as pd
except:
    print("Error: Importing the pandas package failed. Make sure that pandas is alredy installed, or just run `pip install pandas`.")
    exit(1)


MODE_ORDER = [
    "before_complete",
    "before_incomplete",
    "success",
    "failure",
    "after_complete",
    "after_incomplete",
]

def get_log_info(path: str) -> pd.DataFrame:
    # database
    match = re.compile(r".+\.([^.]+)").match(os.path.basename(path))
    if match is not None:
        database = match.group(1)
    else:
        database = None

    # job_id
    environ_path = os.path.join(path, "environ")
    if os.path.exists(environ_path):
        with open(environ_path) as f:
            lines = f.read().strip().split("\n")

        matcher = re.compile(r"PJM_JOBID=(\d+)")
        job_id = [matcher.match(x) for x in lines]
        job_id = [x.group(1) for x in job_id if x is not None]
        job_id = job_id[0] if len(job_id) > 0 else None
    else:
        job_id = None

    # chains
    matcher = re.compile(r"chains_(.+)_(\d+)_(.+)\.csv")
    chain_files = [(x, matcher.match(x)) for x in os.listdir(path)]
    chain_files = [(x[0], x[1].groups()) for x in chain_files if x[1] is not None]
    chain_stats = []
    for chain_file, (_, step, mode) in chain_files:
        with open(os.path.join(path, chain_file)) as f:
            chains = f.read().strip().split("\n")

        if "" in chains:
            chains.remove("")

        chain_stats.append(
            {
                "step": int(step),
                "mode": mode,
                "count": len(chains),
                "time": datetime.fromtimestamp(os.path.getmtime(os.path.join(path, chain_file))),
            }
        )

    # Per-step # of processes and threads
    output_path = os.path.join(path, "output")
    if os.path.exists(output_path):
        with open(output_path) as f:
            lines = f.read().strip().split("\n")

        matcher = re.compile(r"Starting the script:\s*NumProcs=(\d+),\s*NumTaskProcs=(\d+),\s*NumThreads=(\d+),.*")
        proc_configs = [matcher.match(x) for x in lines]
        proc_configs = [x.groups() for x in proc_configs if x is not None]
    else:
        proc_configs = []

    # Returns if job has not started
    if len(chain_stats) > 0:
        df = []
        for step in sorted(set([x["step"] for x in chain_stats])):
            row = {}
            row["job_id"] = job_id
            row["database"] = database
            row["step"] = step
            row["time"] = None
            row["proc"]   = proc_configs[step][0] if step < len(proc_configs) else None
            row["thread"] = proc_configs[step][2] if step < len(proc_configs) else None
            for stat in chain_stats:
                if stat["step"] == step:
                    row[stat["mode"]] = stat["count"]
                    if row["time"] is None or stat["time"] > row["time"]:
                        row["time"] = stat["time"]

            if "after_complete" in row.keys() and \
               "after_incomplete" in row.keys():
               comp = row["after_complete"]
               incomp = row["after_incomplete"]
               row["progress"] = float(comp)/(comp+incomp)*100
            elif "before_complete" in row.keys() and \
               "before_incomplete" in row.keys() and \
               "success" in row.keys():
               comp = row["before_complete"] + row["success"]
               incomp = row["before_incomplete"] - row["success"]
               row["progress"] = float(comp)/(comp+incomp)*100

            df.append(row)

    else:
        if database is None:
            return None

        df = [{
            "job_id": job_id,
            "database": database,
        }]

    df = pd.DataFrame(df)
    column_order = ["job_id", "database", "step", "proc", "thread", "time"] + MODE_ORDER + ["progress"]
    df = df.reindex(columns=column_order)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root-dir",
        type=str,
        default="log",
        help="Path to the root log directory",
    )
    parser.add_argument(
        "--sort-by-time", "-t",
        dest="sort_time",
        action="store_const",
        const=True,
        default=False,
        help="Sort by last update time"
    )
    parser.add_argument(
        "--sort-by-database", "-d",
        dest="sort_database",
        action="store_const",
        const=True,
        default=False,
        help="Sort by database name"
    )

    args = parser.parse_args()

    log_dirs = [x for x in os.listdir(args.root_dir) if os.path.isdir(os.path.join(args.root_dir, x))]
    df = None
    for log_dir in log_dirs:
        ret = get_log_info(os.path.join(args.root_dir, log_dir))
        if ret is None:
            continue

        if df is None:
            df = ret
        else:
            df = pd.concat([df, ret])

    if df is None or len(df.index) == 0:
        print(f"No active log directories found. Please check the \"{args.root_dir}\" directory contains at least one directory.")
        exit()

    if args.sort_time:
        sort_by = "time"
        na_position = "first"
    elif args.sort_database:
        sort_by = ["database", "job_id", "step"]
        na_position = "last"
    else:
        sort_by = ["job_id", "step"]
        na_position = "last"

    df = df.sort_values(sort_by, na_position=na_position)

    df = df.rename(
        columns=
        {
            "job_id": "Job ID",
            "database": "DB",
            "step": "Step",
            "proc": "#Procs.",
            "thread": "#Threads",
            "time": "Last update",
            "before_complete": "#Compl.(b)",
            "before_incomplete": "#Incompl.(b)",
            "success": "#Success",
            "failure": "#Failure",
            "after_complete": "#Compl.(a)",
            "after_incomplete": "#Incompl.(a)",
            "progress": "Progress [%]",
        }
    )
    df = df.fillna("-")
    print(df.to_string(index=False))
