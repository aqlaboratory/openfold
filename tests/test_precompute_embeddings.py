import os
import tempfile
import unittest
from unittest import mock

import torch

from scripts import precompute_embeddings


class _FakeModel:
    def __init__(self):
        self.last_toks_shape = None

    def to(self, device="cpu"):
        return self

    def __call__(self, toks, repr_layers=None, return_contacts=False):
        self.last_toks_shape = tuple(toks.shape)
        batch, tok_len = toks.shape
        reps = torch.arange(batch * tok_len * 3, dtype=torch.float32)
        reps = reps.reshape(batch, tok_len, 3)
        return {"representations": {33: reps}}


class _FakeAlphabet:
    @staticmethod
    def get_batch_converter():
        def _convert(batch):
            labels = [label for label, _ in batch]
            sequences = [seq for _, seq in batch]
            max_len = max(len(seq) for seq in sequences)
            toks = torch.zeros((len(batch), max_len + 2), dtype=torch.int64)
            return labels, sequences, toks

        return _convert


class TestPrecomputeEmbeddings(unittest.TestCase):
    def test_parser_defaults_and_no_truncate_flag(self):
        parser = precompute_embeddings.create_parser()

        args = parser.parse_args(["input", "output"])
        self.assertTrue(args.truncate)

        args = parser.parse_args(["input", "output", "--no-truncate"])
        self.assertFalse(args.truncate)

    @mock.patch("scripts.precompute_embeddings.torch.hub.load")
    def test_long_sequences_are_truncated_on_token_axis(self, mock_load):
        model = _FakeModel()
        mock_load.return_value = (model, _FakeAlphabet())

        generator = precompute_embeddings.EmbeddingGenerator(
            truncate=True,
            nogpu=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            output_dir = os.path.join(tmpdir, "embeddings")
            long_seq = "A" * 1030
            short_seq = "G" * 12
            with open(fasta_path, "w") as fp:
                fp.write(f">long\n{long_seq}\n>short\n{short_seq}\n")

            generator.run(fasta_path, output_dir)

            long_result = torch.load(os.path.join(output_dir, "long", "long.pt"))
            short_result = torch.load(os.path.join(output_dir, "short", "short.pt"))

            self.assertEqual(long_result["representations"][33].shape[0], 1022)
            self.assertEqual(short_result["representations"][33].shape[0], len(short_seq))
            self.assertEqual(model.last_toks_shape, (2, 1024))

    @mock.patch("scripts.precompute_embeddings.torch.hub.load")
    def test_long_sequences_fail_when_truncation_disabled(self, mock_load):
        mock_load.return_value = (_FakeModel(), _FakeAlphabet())

        generator = precompute_embeddings.EmbeddingGenerator(
            truncate=False,
            nogpu=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "input.fasta")
            with open(fasta_path, "w") as fp:
                fp.write(f">long\n{'A' * 1030}\n")

            with self.assertRaisesRegex(ValueError, "exceeding the ESM-1b limit"):
                generator.run(fasta_path, tmpdir)


if __name__ == "__main__":
    unittest.main()
