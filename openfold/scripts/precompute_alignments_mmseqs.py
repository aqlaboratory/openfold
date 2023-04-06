import argparse
import logging
import os
from pathlib import Path
import subprocess

from openfold.data.tools import hhsearch


def _split_a3ms(output_dir):
    for fname in os.listdir(output_dir):
        if(not os.path.splitext(fname)[-1] == ".a3m"):
            continue

        fpath = os.path.join(output_dir, fname)

        with open(fpath, "r") as fp:
            a3ms = fp.read()

        # Split by the null byte, excluding the terminating null byte
        a3ms = a3ms.split('\x00')[:-1]

        for a3m in a3ms:
            name = a3m.split('\n', 1)[0][1:]
            prot_dir = os.path.join(output_dir, name)
            Path(prot_dir).mkdir(parents=True, exist_ok=True)
            with open(os.path.join(prot_dir, fname), "w") as fp:
                fp.write(a3m)

        os.remove(fpath)
        os.remove(fpath + ".dbtype")
        os.remove(fpath + ".index")


def main(args):
    with open(args.input_fasta, "r") as f:
        lines = [l.strip() for l in f.readlines()]

    names = lines[::2]
    seqs =  lines[1::2]

    if(args.fasta_chunk_size is None):
        chunk_size = len(seqs)
    else:
        chunk_size = args.fasta_chunk_size

    # Make the output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)


    s = 0
    while(s < len(seqs)):
        e = s + chunk_size
        chunk_fasta = [el for tup in zip(names[s:e], seqs[s:e]) for el in tup] 
        s = e
        
        prot_dir = os.path.join(args.output_dir, chunk_fasta[0][1:].upper())
        if(os.path.exists(prot_dir)):
            # We've already computed this chunk
            continue

        chunk_fasta_path = os.path.join(args.output_dir, "tmp.fasta")
        with open(chunk_fasta_path, "w") as f:
            f.write('\n'.join(chunk_fasta) + '\n')

        cmd = [
            "scripts/colabfold_search.sh",
            args.mmseqs_binary_path,
            chunk_fasta_path,
            args.mmseqs_db_dir,
            args.output_dir,
            args.uniref_db,
            '""',
            '""' if args.env_db is None else args.env_db,
            "0" if args.env_db is None else "1",
            "0", # compute templates
            "1", # filter
            "1", # use precomputed index 
            "0", # db-load-mode
        ]

        logging.info('Launching subprocess "%s"', " ".join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        stdout, stderr = process.communicate()
        retcode = process.wait()
        
        if retcode:
            raise RuntimeError(
                "MMseqs failed\nstdout:\n%s\n\nstderr:\n%s\n"
                % (stdout.decode("utf-8"), stderr.decode("utf-8"))
            )

        _split_a3ms(args.output_dir)

        # Clean up temporary files
        os.remove(chunk_fasta_path)


    hhsearch_pdb70_runner = hhsearch.HHSearch(
        binary_path=args.hhsearch_binary_path, databases=[args.pdb70]
    )


    for d in os.listdir(args.output_dir):
        dpath = os.path.join(args.output_dir, d)
        if(not os.path.isdir(dpath)):
            continue
        for fname in os.listdir(dpath):
            fpath = os.path.join(dpath, fname)
            if(not "uniref" in fname or 
                not os.path.splitext(fname)[-1] == ".a3m"):
                continue

            with open(fpath, "r") as fp:
                a3m = fp.read()

            hhsearch_result = hhsearch_pdb70_runner.query(a3m)
            pdb70_out_path = os.path.join(dpath, "pdb70_hits.hhr")
            with open(pdb70_out_path, "w") as f:
                f.write(hhsearch_result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_fasta", type=str,
        help="Path to input FASTA file. Can contain one or more sequences."
    )
    parser.add_argument(
        "mmseqs_db_dir", type=str,
        help="""Path to directory containing pre-processed MMSeqs2 DBs 
                (see README)"""
    )
    parser.add_argument(
        "uniref_db", type=str,
        help="Basename of uniref database"
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Output directory"
    )
    parser.add_argument(
        "mmseqs_binary_path", type=str,
        help="Path to mmseqs binary"
    )
    parser.add_argument(
        "--hhsearch_binary_path", type=str, default=None,
        help="""Path to hhsearch binary (for template search). In future 
                versions, we'll also use mmseqs for this"""
    )
    parser.add_argument(
        "--pdb70", type=str, default=None,
        help="Basename of the pdb70 database"
    )
    parser.add_argument(
        "--env_db", type=str, default=None,
        help="Basename of environmental database"
    )
    parser.add_argument(
        "--fasta_chunk_size", type=int, default=None,
        help="""How many sequences should be processed at once. All sequences 
                processed at once by default."""
    )

    args = parser.parse_args()

    if(args.hhsearch_binary_path is not None and args.pdb70 is None):
        raise ValueError(
            "pdb70 must be specified along with hhsearch_binary_path"
        )

    main(args)
