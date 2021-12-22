import copy
import gzip

import os

import pickle

import numpy as np
import torch
import unittest

from data.data_transforms import make_seq_mask, add_distillation_flag, make_all_atom_aatype, fix_templates_aatype, \
    correct_msa_restypes, squeeze_features, randomly_replace_msa_with_unknown, MSA_FEATURE_NAMES, sample_msa, \
    crop_extra_msa, delete_extra_msa, nearest_neighbor_clusters, make_msa_mask, make_hhblits_profile
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

    def test_fix_templates_aatype(self):
        template_seq = torch.tensor(list(range(20))*2, dtype=torch.int64)
        template_seq = template_seq.unsqueeze(0).transpose(0, 1)
        template_seq_one_hot = torch.FloatTensor(template_seq.shape[0], 20).zero_()
        template_seq_one_hot.scatter_(1, template_seq, 1)
        template_aatype = torch.tensor(template_seq_one_hot).unsqueeze(0)
        protein = {'template_aatype': template_aatype}
        protein = fix_templates_aatype(protein)
        print(protein)
        template_seq_ours = torch.tensor([[0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18]*2])
        assert torch.all(torch.eq(protein['template_aatype'], template_seq_ours))

    def test_correct_msa_restypes(self):
        with open("../test_data/features.pkl", 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa'], dtype=torch.int64)}
        protein = correct_msa_restypes(protein)
        print(protein)
        assert torch.all(torch.eq(torch.tensor(features['msa'].shape), torch.tensor(protein['msa'].shape)))

    def test_squeeze_features(self):
        with open("../test_data/features.pkl", "rb") as file:
            features = pickle.load(file)
            print(os.path.realpath(file.name), 'Keys: ', features.keys())

        features_list = [
            'domain_name', 'msa', 'num_alignments', 'seq_length', 'sequence',
            'superfamily', 'deletion_matrix', 'resolution',
            'between_segment_residues', 'residue_index', 'template_all_atom_mask']

        protein = {'aatype': torch.tensor(features['aatype'])}
        for k in features_list:
            if k in features:
                print(k, features[k].dtype)
                if k in ['domain_name', 'sequence']:
                    protein[k] = np.expand_dims(features[k], -1)
                else:
                    protein[k] = torch.tensor(features[k]).unsqueeze(-1)

        for k in ['seq_length', 'num_alignments']:
            if k in protein:
                protein[k] = torch.tensor(protein[k]).unsqueeze(0)

        protein_squeezed = squeeze_features(protein)
        print(protein)
        for k in features_list:
            if k in protein:
                print(k, protein_squeezed[k].shape, features[k].shape)
                assert protein_squeezed[k].shape == features[k].shape

    def test_randomly_replace_msa_with_unknown(self):
        with open('../test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa']),
                   'aatype': torch.argmax(torch.tensor(features['aatype']), dim=1)}
        replace_proportion = 0.15
        x_idx = 20
        protein = randomly_replace_msa_with_unknown.__wrapped__(protein, replace_proportion)
        unknown_proportion_in_msa = torch.bincount(protein['msa'].flatten()) / torch.numel(protein['msa'])
        unknown_proportion_in_seq = torch.bincount(protein['aatype'].flatten()) / torch.numel(protein['aatype'])
        print(protein)
        print('Proportion of X in MSA: ', unknown_proportion_in_msa[x_idx])
        print('Proportion of X in sequence: ', unknown_proportion_in_seq[x_idx])

    def test_sample_msa(self):
        with open('../test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        max_seq = 1000
        keep_extra = True
        protein = {}
        for k in MSA_FEATURE_NAMES:
            if k in features:
                protein[k] = torch.tensor(features[k])

        protein_processed = sample_msa.__wrapped__(protein.copy(), max_seq, keep_extra)
        print(protein)
        for k in MSA_FEATURE_NAMES:
            if k in protein and keep_extra:
                assert protein_processed[k].shape[0] == min(protein[k].shape[0], max_seq)
                assert 'extra_'+k in protein_processed
                print('extra_'+str(k), protein_processed['extra_'+k].shape)
                print('msa', protein[k].shape[0] - min(protein[k].shape[0], max_seq))
                assert protein_processed['extra_'+k].shape[0] == protein[k].shape[0] - min(protein[k].shape[0], max_seq)

    def test_crop_extra_msa(self):
        with open('../test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        max_extra_msa = 10
        protein = {'extra_msa': torch.tensor(features['msa'])}
        num_seq = protein["extra_msa"].shape[0]

        protein = crop_extra_msa.__wrapped__(protein, max_extra_msa)
        print(protein)
        for k in MSA_FEATURE_NAMES:
            if "extra_" + k in protein:
                assert protein["extra_" + k].shape[0] == min(max_extra_msa, num_seq)

    def test_delete_extra_msa(self):
        protein = {'extra_msa': torch.rand((512, 100, 23))}
        extra_msa_has_deletion_shape = list(protein['extra_msa'].shape)
        extra_msa_has_deletion_shape[2] = 1
        protein['extra_deletion_matrix'] = torch.rand(extra_msa_has_deletion_shape)
        protein = delete_extra_msa(protein)
        print(protein)
        for k in MSA_FEATURE_NAMES:
            assert 'extra_' + k not in protein
        assert 'extra_msa' not in protein

    def test_nearest_neighbor_clusters(self):
        with gzip.open('../test_data/sample_feats.pickle.gz', 'rb') as f:
            features = pickle.load(f)

        protein = {'msa': torch.tensor(features['true_msa'][0], dtype=torch.int64),
                   'msa_mask': torch.tensor(features['msa_mask'][0], dtype=torch.int64),
                   'extra_msa': torch.tensor(features['extra_msa'][0], dtype=torch.int64),
                   'extra_msa_mask': torch.tensor(features['extra_msa_mask'][0], dtype=torch.int64)}
        protein = nearest_neighbor_clusters.__wrapped__(protein, 0)
        print(protein)
        assert 'extra_cluster_assignment' in protein

    def test_make_msa_mask(self):
        with open('../test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        msa_mat = torch.tensor(features['msa'])
        protein = {'msa': msa_mat}
        protein = make_msa_mask(protein)
        print(protein)
        assert 'msa_row_mask' in protein
        assert protein['msa_row_mask'].shape[0] == msa_mat.shape[0]

    def test_make_hhblits_profile(self):
        with open('../test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa'], dtype=torch.int64)}
        protein = make_hhblits_profile(protein)
        assert 'hhblits_profile' in protein
        assert protein['hhblits_profile'].shape == torch.Size((protein['msa'].shape[1], 22))


if __name__ == '__main__':
    unittest.main()

