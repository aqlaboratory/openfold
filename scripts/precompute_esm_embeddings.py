import argparse
import pathlib

import torch

from openfold.data import parsers


class SequenceDataset(object):
    def __init__(self, labels, sequences):
        self.labels = labels
        self.sequences = sequences

    @classmethod
    def from_file(cls, fasta_file):
        labels, sequences = [], []

        with open(fasta_file, "r") as infile:
            fasta_str = infile.read()
            sequences, labels = parsers.parse_fasta(fasta_str)

        assert len(set(labels)) == len(labels),\
            "Sequence lables need to be unique. Duplicates found!"

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
    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
    model.eval()
    if torch.cuda.is_available() and not args.nogpu:
        print("Using device=GPU")
        model = model.cuda()

    fasta_file_path = args.fasta_path
    dataset = SequenceDataset.from_file(fasta_file_path)
    batches = dataset.get_batch_indices(args.toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(), batch_sampler=batches
    )
    print(f"Processed {fasta_file_path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in range(args.repr_layers))
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in args.repr_layers]

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
        if torch.cuda.is_available() and not args.nogpu:
            toks = toks.to(device="cuda", non_blocking=True)

        # ESM cannot handle >1022 residue sequences
        if args.truncate:
            toks = toks[:1022]

        out = model(toks, repr_layers=repr_layers)


        logits = out["logits"].to(device="cpu")
        representations = {
            layer: t.to(device="cpu") for layer, t in out["representations"].items()
        }

        for i, label in enumerate(labels):
            args.output_file = args.output_dir / label / f"{label}.pt"
            args.output_file.parent.mkdir(parents=True, exist_ok=True)
            result = {"label": label}

            if "per_tok" in args.include:
                result["representations"] = {
                    layer: t[i, 1 : len(strs[i] + 1).clone()]
                    for layer, t in representations.items()
                }
            torch.save(
                result,
                args.output_file,
            )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract per-residue ESM embeddings for sequences in a FASTA file, "
                    "or a dir of FASTA files."
    )
    parser.add_argument("model",
        type=str,
        help="Either the name of a pretrained model to download via PyTorch Hub, "
             "or path to a PyTorch model file.",
        default="esm1b_t33_650M_UR50S"
    )
    parser.add_argument("fasta_path",
        type=pathlib.Path,
        help="Path to a FASTA file, or to a dir containing FASTA files."
    )
    parser.add_argument("output_path",
        type=pathlib.Path,
        help="Output directory where ESM representations will be stored."
    )
    parser.add_argument(
        "--toks_per_batch",
        type=int,
        default=4096,
        help="Maximum batch size, i.e. number of tokens in batch."
    )
    parser.add_argument(
        "--repr_layers",
        type=int,
        default=[-1],
        nargs="+",
        help="Layers indices from which to extract representations (0 to num_layers, inclusive)."
    )
    parser.add_argument(
        "--truncate",
        action="store_true",
        default=True,
        help="Truncate sequences longer than 1022 (ESM restriction). Default: True"
    )
    parser.add_argument("--nogpu",
        action="store_true",
        help="Do not use GPU"
    )

    args = parser.parse_args()
    main(args)