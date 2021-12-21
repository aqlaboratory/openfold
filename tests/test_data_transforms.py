import copy
import gzip

import os

import pickle

import numpy
import torch
import unittest

from data.data_transforms import make_seq_mask, add_distillation_flag, make_all_atom_aatype
from openfold.config import model_config


class TestDataTransforms(unittest.TestCase):
    def test_make_seq_mask(self):
        seq = torch.tensor([range(20)], dtype=torch.int64).transpose(0,1)
        seq_one_hot = torch.FloatTensor(seq.shape[0], 20).zero_()
        seq_one_hot.scatter_(1, seq, 1)
        protein_aatype = torch.tensor(seq_one_hot)
        protein = {'aatype': protein_aatype}
        protein = make_seq_mask(protein)
        print(protein)
        assert 'seq_mask' in protein
        assert protein['seq_mask'].shape == torch.Size((seq.shape[0], 20))

    def test_add_distillation_flag(self):
        protein = {}
        protein = add_distillation_flag.__wrapped__(protein, True)
        print(protein)
        assert 'is_distillation' in protein
        assert protein['is_distillation'] is True

    def test_make_all_atom_aatype(self):
        seq = torch.tensor([range(20)], dtype=torch.int64).transpose(0, 1)
        seq_one_hot = torch.FloatTensor(seq.shape[0], 20).zero_()
        seq_one_hot.scatter_(1, seq, 1)
        protein_aatype = torch.tensor(seq_one_hot)
        protein = {'aatype': protein_aatype}
        protein = make_all_atom_aatype(protein)
        print(protein)
        assert 'all_atom_aatype' in protein
        assert protein['all_atom_aatype'].shape == protein['aatype'].shape


if __name__ == '__main__':
    unittest.main()

