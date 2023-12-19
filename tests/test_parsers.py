import os
import tempfile
import unittest

from openfold.data.parsers import \
    convert_stockholm_to_a3m, \
    convert_stockholm_to_a3m_stream


class TestParsers(unittest.TestCase):
    def test_sto_to_a3m_stream_equal(self):
        self._test_sto_to_a3m_stream_equal()

    def test_sto_to_a3m_stream_equal_max_sequences_1(self):
        self._test_sto_to_a3m_stream_equal(1)

    def test_sto_to_a3m_stream_equal_max_sequences_2(self):
        self._test_sto_to_a3m_stream_equal(2)

    def _test_sto_to_a3m_stream_equal(self, max_sequences: int=None):
        """
        Test if the output of convert_stockholm_to_a3m_stream is equals to
        that of convert_stockholm_to_a3m.
        """

        kwargs = {}
        if max_sequences:
            kwargs["max_sequences"] = max_sequences

        sto_paths = [
            "tests/test_data/alignments/mgnify_hits.sto",
            "tests/test_data/alignments/uniref90_hits.sto",
            "tests/test_data/sto_to_a3m_alignments/hits.sto",
        ]
        for sto_path in sto_paths:
            with open(sto_path, "r") as f:
                sto = f.read()
            a3m = convert_stockholm_to_a3m(sto, **kwargs)
            self.assertTrue(len(a3m) > 0)

            a3m_stream_path = tempfile.mktemp()
            convert_stockholm_to_a3m_stream(sto_path, a3m_stream_path, **kwargs)
            self.assertTrue(os.path.exists(a3m_stream_path))
            with open(a3m_stream_path, "r") as f:
                a3m_stream = f.read()

            self.assertEqual(a3m, a3m_stream)
            os.remove(a3m_stream_path)

if __name__ == '__main__':
    unittest.main()
