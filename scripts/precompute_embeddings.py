# Some functions borrowed from [ESM](https://www.github.com/facebookresearch/esm)
import argparse
import logging
import os

import torch

from openfold.data import parsers

logging.basicConfig(level=logging.INFO)

ESM_MAX_RESIDUES = 1022

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


class EmbeddingGenerator:
    """Generates the ESM-1b embeddings for the single sequence model"""
    def __init__(self,
        toks_per_batch: int = 4096,
        truncate: bool = True,
        use_local_esm: str = None,
        nogpu: bool = False,
    ):
        self.toks_per_batch = toks_per_batch
        self.truncate = truncate
        self.use_local_esm = use_local_esm
        self.nogpu = nogpu
        self.max_residues = ESM_MAX_RESIDUES
        
        # Generate embeddings in bulk
        if self.use_local_esm:
            self.model, self.alphabet = torch.hub.load(self.use_local_esm, "esm1b_t33_650M_UR50S", source='local')
        else:
            self.model, self.alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
        if torch.cuda.is_available() and not self.nogpu:
            self.model = self.model.to(device="cuda")

    def parse_sequences(self, fasta_dir, output_dir):
        labels = []
        seqs = []

        # Generate a single bulk file
        for f in os.listdir(fasta_dir):
            f_name, ext = os.path.splitext(f)
            if ext != '.fasta' and ext != '.fa':
                logging.warning(f"Ignoring non-FASTA file: {f}")
                continue
            with open(os.path.join(fasta_dir, f), 'r') as infile:
                seq = infile.readlines()[1].strip()
            labels.append(f_name)
            seqs.append(seq)
        
        lines = []
        for label, seq in zip(labels, seqs):
            lines += f'>{label}\n'
            lines += f'{seq}\n'
        os.makedirs(output_dir, exist_ok=True)
        temp_fasta_file = os.path.join(output_dir, 'temp.fasta')
        with open(temp_fasta_file, 'w') as outfile:
            outfile.writelines(lines)
        return temp_fasta_file
    
    def run(
        self,
        fasta_file,
        output_dir,
    ):

        dataset = SequenceDataset.from_file(fasta_file)
        batches = dataset.get_batch_indices(self.toks_per_batch, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(), batch_sampler=batches
        )
        logging.info("Loaded all sequences")
        repr_layers = [33]

        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                logging.info(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
                effective_lengths = []
                truncated_labels = []
                for label, seq in zip(labels, strs):
                    seq_len = len(seq)
                    if seq_len > self.max_residues:
                        if not self.truncate:
                            raise ValueError(
                                f"Sequence '{label}' has length {seq_len}, exceeding "
                                f"the ESM-1b limit of {self.max_residues} residues. "
                                "Re-run with truncation enabled."
                            )
                        truncated_labels.append((label, seq_len))
                    effective_lengths.append(min(seq_len, self.max_residues))

                if truncated_labels:
                    for label, seq_len in truncated_labels:
                        logging.warning(
                            "Truncating sequence '%s' from %d to %d residues for ESM-1b.",
                            label,
                            seq_len,
                            self.max_residues,
                        )

                max_effective_length = max(effective_lengths)
                toks = toks[:, : max_effective_length + 2]
                if torch.cuda.is_available() and not self.nogpu:
                    toks = toks.to(device="cuda", non_blocking=True)

                out = self.model(toks, repr_layers=repr_layers, return_contacts=False)

                representations = {
                    33: out["representations"][33].to(device="cpu")
                }

                for i, (label, effective_length) in enumerate(zip(labels, effective_lengths)):
                    os.makedirs(os.path.join(output_dir, label), exist_ok=True)
                    result = {"label": label}

                    result["representations"] = {
                        33: representations[33][i, 1: effective_length + 1].clone()
                    }
                    torch.save(
                        result,
                        os.path.join(output_dir, label, label+".pt")
                    )
        
        

def main(args):
    logging.info("Loading the model...")
    embedding_generator = EmbeddingGenerator(
        args.toks_per_batch,
        args.truncate,
        args.use_local_esm,
        args.nogpu)
    logging.info("Loading the sequences and running the inference...")
    temp_fasta_file = embedding_generator.parse_sequences(
        args.fasta_dir,
        args.output_dir
    )
    embedding_generator.run(
        temp_fasta_file,
        args.output_dir
    )
    os.remove(temp_fasta_file)
    logging.info("Completed.")


def create_parser():
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
        "--truncate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            f"Truncate sequences longer than {ESM_MAX_RESIDUES} residues "
            "to satisfy the ESM-1b limit. Default: True"
        ),
    )
    parser.add_argument(
        "--use_local_esm", type=str, default=None,
        help="Use a local ESM repository instead of cloning from Github"
    )
    parser.add_argument(
        "--nogpu", action="store_true",
        help="Do not use GPU"
    )

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
