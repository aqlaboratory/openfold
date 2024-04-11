#!/usr/bin/env python3

import os
import re
import argparse
from datetime import datetime
from typing import List, Dict, Any

pat = re.compile('job_([0-9]+).csv')

try:
    import pandas as pd
except:
    print("Error: Importing the pandas package failed. Make sure that pandas is alredy installed, or just run `pip install pandas`.")
    exit(1)


def get_chains(csv_file):
    if not os.path.isfile(csv_file):
        return None

    with open(csv_file) as f:
        chains = f.read().strip().split("\n")

    return set(filter(None, chains))

def get_log_info(path: str) -> pd.DataFrame:

    job_id = os.path.basename(path)

    # set directory's update time
    last_update =  datetime.fromtimestamp(os.path.getmtime(path))

    complete_path   = os.path.join(path, 'before_complete.csv')
    incomplete_path = os.path.join(path, 'before_incomplete.csv')
    noalign_path    = os.path.join(path, 'before_noalign.csv')
    skip_path       = os.path.join(path, 'before_skip.csv')
    processed_path  = os.path.join(path, 'processed.csv')

    complete_chains   = get_chains(complete_path)
    incomplete_chains = get_chains(incomplete_path)
    noalign_chains    = get_chains(noalign_path)
    skip_chains       = get_chains(skip_path)

    if (complete_chains is None) or \
       (incomplete_chains is None) or \
       (noalign_chains is None) or \
       (skip_chains is None):
        data = {'Job ID'      : job_id,
                'Last update' : last_update,
                '#Compl.(b)'  : len(complete_chains) if complete_chains else None,
                '#Incompl.(b)': len(incomplete_chains) if incomplete_chains else None,
                '#NoAlign.'   : len(noalign_chains) if noalign_chains else None,
                '#Skip'       : len(skip_chains) if noalign_chains else None,
                '#Success'    : None,
                '#Failure'    : None}
        return pd.DataFrame([data])

    statuses = ['OK', 'NG_timeout', 'NG_unknown', 'NG_noalignment']
    if os.path.isfile(processed_path):
        last_update =  datetime.fromtimestamp(os.path.getmtime(processed_path))

        df = pd.read_csv(processed_path,
                         names=['chain', 'seq_len', 'status', 'time_all', 'time_infer', 'time_relax'],
                         usecols=['chain', 'seq_len', 'status'],
                         dtype = {'chain':'str', 'seq_len':'int32', 'status':'str'})

        status_count = { st: (df['status'] == st).sum() for st in statuses}
    else:
        status_count = { st: 0 for st in statuses}

    n_compl = len(complete_chains)
    n_incompl = len(incomplete_chains)
    data = {'Job ID'      : job_id,
            'Last update' : last_update,
            '#Compl.(b)'  : n_compl,
            '#Incompl.(b)': n_incompl,
            '#NoAlign.'   : len(noalign_chains) + status_count['NG_noalignment'],
            '#Skip'       : len(skip_chains),
            '#Success'    : status_count['OK'],
            '#Failure'    : status_count['NG_timeout'] + status_count['NG_unknown'],
    }

    return pd.DataFrame([data])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root-dir",
        type=str,
        default="log",
        help="Path to the root log directory",
    )

    args = parser.parse_args()

    log_dir = os.path.join(args.root_dir, 'result')

    if not os.path.isdir(log_dir):
        raise Exception(f'There is no result directory in specified root directory. Please check --root-dir option.')

    job_dirs = [x for x in os.listdir(log_dir) if os.path.isdir(os.path.join(log_dir, x))]

    df = None
    for job_dir in job_dirs:
        ret = get_log_info(os.path.join(log_dir, job_dir))
        if ret is None:
            continue

        if df is None:
            df = ret
        else:
            df = pd.concat([df, ret])

    if df is not None:
        df = df.sort_values('Job ID')
        df['Progress[%]'] = ((df['#Compl.(b)'] + df['#Success']) * 100.0 / (df['#Compl.(b)'] + df['#Incompl.(b)']))
        df['Progress[%]'] = df['Progress[%]'].astype('float').round(1)
        df = df.fillna('-')
        print(df.to_string(index=False))
    else:
        print("No data!")
