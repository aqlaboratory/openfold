# Some functions borrowed from [ESM](https://www.github.com/facebookresearch/esm)
import argparse
import logging
import os

import torch

from openfold.data import parsers

logging.basicConfig(level=logging.INFO)

class SequenceDataset(object):
    def __init__(self, labels, sequences) -> None:
        self.labels = labels
        self.sequences = sequences
    
    @classmethod
    def from_file(cls, fasta_file):
        labels, sequences = [], []

        with open(fasta_file, "r") as infile:
            fasta_str = infile.read()
            sequences, labels = parsers.parse_fasta(fasta_str)
        
        assert len(set(labels)) == len(labels),\
            "Sequence labels need to be unique. Duplicates found!"
        
        return cls(labels, sequences)
    
    def __len__(self):
        return len(self.labels)
    
    def __getitem__(self, idx):
        return self.labels[idx], self.sequences[idx]
    
    def get_batch_indices(self, toks_per_batch, extra_toks_per_seq):
        sizes = [(len(s), i) for i, s in enumerate(self.sequences)]
        sizes.sort()
        batches = []
        buf = []
        max_len = 0

        def _flush_current_buf():
            nonlocal max_len, buf
            if len(buf) == 0:
                return
            batches.append(buf)
            buf = []
            max_len = 0
        
        for sz, i in sizes:
            sz += extra_toks_per_seq
            if max(sz, max_len) * (len(buf)+1) > toks_per_batch:
                _flush_current_buf()
            max_len = max(max_len, sz)
            buf.append(i)
        
        _flush_current_buf()
        return batches

def main(args):
    labels = []
    seqs = []

    # Generate a single bulk file
    for f in os.listdir(args.fasta_dir):
        f_name, ext = os.path.splitext(f)
        if ext != '.fasta' and ext != '.fa':
            logging.warning(f"Ignoring non-FASTA file: {f}")
            continue
        with open(os.path.join(args.fasta_dir, f), 'r') as infile:
            seq = infile.readlines()[1].strip()
        labels.append(f_name)
        seqs.append(seq)
    
    lines = []
    for label, seq in zip(labels, seqs):
        lines += f'>{label}\n'
        lines += f'{seq}\n'
    os.makedirs(args.output_dir, exist_ok=True)
    temp_fasta_file = os.path.join(args.output_dir, 'temp.fasta')
    with open(temp_fasta_file, 'w') as outfile:
        outfile.writelines(lines)

    # Generate embeddings in bulk
    if args.use_local_esm:
        model, alphabet = torch.hub.load(args.use_local_esm, "esm1b_t33_650M_UR50S", source='local')
    else:
        model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
    if torch.cuda.is_available() and not args.nogpu:
        model = model.to(device="cuda")
    dataset = SequenceDataset.from_file(temp_fasta_file)
    batches = dataset.get_batch_indices(args.toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(), batch_sampler=batches
    )
    logging.info("Loaded all sequences")
    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in args.repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in args.repr_layers]

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            logging.info(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
            if torch.cuda.is_available() and not args.nogpu:
                toks = toks.to(device="cuda", non_blocking=True)
            
            if args.truncate:
                toks = toks[:1022]
            
            out = model(toks, repr_layers=repr_layers, return_contacts=False)

            logits = out["logits"].to(device="cpu")
            representations = {
                layer: t.to(device="cpu") for layer, t in out["representations"].items()
            }

            for i, label in enumerate(labels):
                os.makedirs(os.path.join(args.output_dir, label), exist_ok=True)
                result = {"label": label}

                if "per_tok" in args.include:
                    result["representations"] = {
                        layer: t[i, 1: len(strs[i]) + 1].clone()
                        for layer, t in representations.items()
                    }
                torch.save(
                    result,
                    os.path.join(args.output_dir, label, label+".pt")
                )
    
    os.remove(temp_fasta_file)
    logging.info("Completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_dir", type=str,
        help="""Path to directory containing FASTA files."""
    )
    parser.add_argument(
        "output_dir", type=str,
        help="Directory in which to output embeddings"
    )
    parser.add_argument(
        "--toks_per_batch", type=int, default=4096, 
        help="maximum tokens in a batch"
    )
    parser.add_argument(
        "--repr_layers", type=int, default=[-1], nargs="+",
        help="Layer indices from which to extract representations (0 to num_layers, inclusive)"
    )
    parser.add_argument(
        "--include", type=str, default=["per_tok"], nargs="+",
        choices=["mean", "per_tok", "bos", "contacts"],
        help="Specify which representations to return"
    )
    parser.add_argument(
        "--truncate", action="store_true", default=True,
        help="Truncate sequences longer than 1022 (ESM restriction). Default: True"
    )
    parser.add_argument(
        "--use_local_esm", type=str, default=None,
        help="Use a local ESM repository instead of cloning from Github"
    )
    parser.add_argument(
        "--nogpu", action="store_true",
        help="Do not use GPU"
    )

    args = parser.parse_args()

    main(args)
