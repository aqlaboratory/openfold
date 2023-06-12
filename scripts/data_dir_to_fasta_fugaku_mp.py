# Copyright 2023 RIKEN & Fujitsu Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import logging
import os
#import threading
from multiprocessing import Pool

from openfold.data import mmcif_parsing
from openfold.np import protein, residue_constants

from functools import partial

def parse_all(data_dir, fnames, raise_errors = False):
    fasta = []
    num_files = len(fnames)
    for i, fname in enumerate(fnames):
        if i % 10 == 0:
            print("{}/{}".format(i, num_files))
        basename, ext = os.path.splitext(fname)
        basename = basename.upper()
        fpath = os.path.join(data_dir, fname)
        if(ext == ".cif"):
            with open(fpath, 'r') as fp:
                mmcif_str = fp.read()

            mmcif = mmcif_parsing.parse(
                file_id=basename, mmcif_string=mmcif_str
            )
            if(mmcif.mmcif_object is None):
                logging.warning(f'Failed to parse {fname}...')
                if(raise_errors):
                    raise list(mmcif.errors.values())[0]
                else:
                    continue

            mmcif = mmcif.mmcif_object
            for chain, seq in mmcif.chain_to_seqres.items():
                chain_id = '_'.join([basename, chain])
                fasta.append(f">{chain_id}")
                fasta.append(seq)
        elif(ext == ".core"):
            with open(fpath, 'r') as fp:
                core_str = fp.read()

            core_protein = protein.from_proteinnet_string(core_str)
            aatype = core_protein.aatype
            seq = ''.join([
                residue_constants.restypes_with_x[aatype[i]]
                for i in range(len(aatype))
            ])
            fasta.append(f">{basename}")
            fasta.append(seq)

    return fasta


def main(args):
    fnames = os.listdir(args.data_dir)
    # fasta = parse_all(args.data_dir, fnames, args.raise_errors)

    pool = Pool(args.no_tasks)

    f = partial(parse_all,
                args.data_dir,
                raise_errors = args.raise_errors,
                )
    
    ret_fasta = pool.map(f, [fnames[i::args.no_tasks] for i in range(args.no_tasks)])
    
    #procs = []
    #ret_fasta = [None]*args.no_tasks
    #for i in range(args.no_tasks):
    #    ret_fasta[i] = []
    #    print(f"Started process {i}...")
    #    targs = [
    #        args.data_dir,
    #        fnames[i::args.no_tasks],
    #        args.raise_errors,
    #        ret_fasta[i],
    #    ]
    #    p = Process(target=parse_all, args=targs)
    #    procs.append(p)
    #    p.start()
    #
    #for p in procs:
    #    p.join()
    #

    fasta = []
    for ret in ret_fasta:
        fasta.extend(ret)
    

    with open(args.output_path, "w") as fp:
        fp.write('\n'.join(fasta))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_dir", type=str,
        help="Path to a directory containing mmCIF or .core files"
    )
    parser.add_argument(
        "output_path", type=str,
        help="Path to output FASTA file"
    )
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )
    parser.add_argument(
        "--no_tasks", type=int, default=1,
    )

    args = parser.parse_args()
    if args.no_tasks < 1:
        raise ValueError("Use data_dir_to_fasta.py if you don't use multiple processes.")

    main(args)
