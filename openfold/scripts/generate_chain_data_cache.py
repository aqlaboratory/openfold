import argparse
from functools import partial
import json
import logging
from multiprocessing import Pool
import os

import sys
sys.path.append(".") # an innocent hack to get this to run from the top level

from tqdm import tqdm

from openfold.data.mmcif_parsing import parse 
from openfold.np import protein, residue_constants


def parse_file(
    f, 
    args,
    chain_cluster_size_dict
):
    file_id, ext = os.path.splitext(f)
    if(ext == ".cif"):
        with open(os.path.join(args.data_dir, f), "r") as fp:
            mmcif_string = fp.read()
        mmcif = parse(file_id=file_id, mmcif_string=mmcif_string)
        if mmcif.mmcif_object is None:
            logging.info(f"Could not parse {f}. Skipping...")
            return {}
        else:
            mmcif = mmcif.mmcif_object

        out = {}
        for chain_id, seq in mmcif.chain_to_seqres.items():
            full_name = "_".join([file_id, chain_id])
            out[full_name] = {}
            local_data = out[full_name]
            local_data["release_date"] = mmcif.header["release_date"]
            local_data["seq"] = seq
            local_data["resolution"] = mmcif.header["resolution"]
           
            if(chain_cluster_size_dict is not None):
                cluster_size = chain_cluster_size_dict.get(
                    full_name.upper(), -1
                )
                local_data["cluster_size"] = cluster_size
    elif(ext == ".pdb"):
        with open(os.path.join(args.data_dir, f), "r") as fp:
            pdb_string = fp.read()
          
        protein_object = protein.from_pdb_string(pdb_string, None)

        chain_dict = {} 
        chain_dict["seq"] = residue_constants.aatype_to_str_sequence(
            protein_object.aatype,
        )
        chain_dict["resolution"] = 0.
        
        if(chain_cluster_size_dict is not None):
            cluster_size = chain_cluster_size_dict.get(
                full_name.upper(), -1
            )
            chain_dict["cluster_size"] = cluster_size

        out = {file_id: chain_dict}

    return out


def main(args):
    chain_cluster_size_dict = None
    if(args.cluster_file is not None):
        chain_cluster_size_dict = {}
        with open(args.cluster_file, "r") as fp:
            clusters = [l.strip() for l in fp.readlines()]

        for cluster in clusters:
            chain_ids = cluster.split()
            cluster_len = len(chain_ids)
            for chain_id in chain_ids:
                chain_id = chain_id.upper()
                chain_cluster_size_dict[chain_id] = cluster_len
   
    accepted_exts = [".cif", ".pdb"]
    files = list(os.listdir(args.data_dir))
    files = [f for f in files if os.path.splitext(f)[-1] in accepted_exts]
    fn = partial(
        parse_file, 
        args=args,
        chain_cluster_size_dict=chain_cluster_size_dict,
    )
    data = {}
    with Pool(processes=args.no_workers) as p:
        with tqdm(total=len(files)) as pbar:
            for d in p.imap_unordered(fn, files, chunksize=args.chunksize):
                data.update(d)
                pbar.update()

    with open(args.output_path, "w") as fp:
        fp.write(json.dumps(data, indent=4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_dir", type=str, help="Directory containing mmCIF or PDB files"
    )
    parser.add_argument(
        "output_path", type=str, help="Path for .json output"
    )
    parser.add_argument(
        "--cluster_file", type=str, default=None,
        help=(
            "Path to a cluster file (e.g. PDB40), one cluster "
            "({PROT1_ID}_{CHAIN_ID} {PROT2_ID}_{CHAIN_ID} ...) per line. "
            "Chains not in this cluster file will NOT be filtered by cluster "
            "size."
        )
    )
    parser.add_argument(
        "--no_workers", type=int, default=4,
        help="Number of workers to use for parsing"
    )
    parser.add_argument(
        "--chunksize", type=int, default=10,
        help="How many files should be distributed to each worker at a time"
    )

    args = parser.parse_args()

    main(args)
