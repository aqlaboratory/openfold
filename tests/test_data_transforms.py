import copy
import gzip

import os

import pickle

import numpy as np
import torch
import unittest

from openfold.data.data_transforms import make_seq_mask, add_distillation_flag, make_all_atom_aatype, fix_templates_aatype, \
    correct_msa_restypes, squeeze_features, randomly_replace_msa_with_unknown, MSA_FEATURE_NAMES, sample_msa, \
    crop_extra_msa, delete_extra_msa, nearest_neighbor_clusters, make_msa_mask, make_hhblits_profile, make_masked_msa, \
    make_msa_feat, crop_templates, make_atom14_masks
from tests.config import config


class TestDataTransforms(unittest.TestCase):
    def test_make_seq_mask(self):
        seq = torch.tensor([range(20)], dtype=torch.int64).transpose(0,1)
        seq_one_hot = torch.FloatTensor(seq.shape[0], 20).zero_()
        seq_one_hot.scatter_(1, seq, 1)
        protein_aatype = seq_one_hot.clone().detach()
        protein = {'aatype': protein_aatype}
        protein = make_seq_mask(protein)
        assert 'seq_mask' in protein
        assert protein['seq_mask'].shape == torch.Size((seq.shape[0], 20))

    def test_add_distillation_flag(self):
        protein = {}
        protein = add_distillation_flag.__wrapped__(protein, True)
        assert 'is_distillation' in protein
        assert protein['is_distillation'] is True

    def test_make_all_atom_aatype(self):
        seq = torch.tensor([range(20)], dtype=torch.int64).transpose(0, 1)
        seq_one_hot = torch.FloatTensor(seq.shape[0], 20).zero_()
        seq_one_hot.scatter_(1, seq, 1)
        protein_aatype = seq_one_hot.clone().detach()
        protein = {'aatype': protein_aatype}
        protein = make_all_atom_aatype(protein)
        assert 'all_atom_aatype' in protein
        assert protein['all_atom_aatype'].shape == protein['aatype'].shape

    def test_fix_templates_aatype(self):
        template_seq = torch.tensor(list(range(20))*2, dtype=torch.int64)
        template_seq = template_seq.unsqueeze(0).transpose(0, 1)
        template_seq_one_hot = torch.FloatTensor(template_seq.shape[0], 20).zero_()
        template_seq_one_hot.scatter_(1, template_seq, 1)
        template_aatype = template_seq_one_hot.clone().detach().unsqueeze(0)
        protein = {'template_aatype': template_aatype}
        protein = fix_templates_aatype(protein)
        template_seq_ours = torch.tensor([[0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18]*2])
        assert torch.all(torch.eq(protein['template_aatype'], template_seq_ours))

    def test_correct_msa_restypes(self):
        with open("tests/test_data/features.pkl", 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa'], dtype=torch.int64)}
        protein = correct_msa_restypes(protein)
        assert torch.all(torch.eq(torch.tensor(features['msa'].shape), torch.tensor(protein['msa'].shape)))

    def test_squeeze_features(self):
        with open("tests/test_data/features.pkl", "rb") as file:
            features = pickle.load(file)

        features_list = [
            'domain_name', 'msa', 'num_alignments', 'seq_length', 'sequence',
            'superfamily', 'deletion_matrix', 'resolution',
            'between_segment_residues', 'residue_index', 'template_all_atom_mask']

        protein = {'aatype': torch.tensor(features['aatype'])}
        for k in features_list:
            if k in features:
                if k in ['domain_name', 'sequence']:
                    protein[k] = np.expand_dims(features[k], -1)
                else:
                    protein[k] = torch.tensor(features[k]).unsqueeze(-1)

        for k in ['seq_length', 'num_alignments']:
            if k in protein:
                protein[k] = protein[k].clone().detach().unsqueeze(0)

        protein_squeezed = squeeze_features(protein)
        for k in features_list:
            if k in protein:
                assert protein_squeezed[k].shape == features[k].shape

    def test_randomly_replace_msa_with_unknown(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa']),
                   'aatype': torch.argmax(torch.tensor(features['aatype']), dim=1)}
        replace_proportion = 0.15
        x_idx = 20
        protein = randomly_replace_msa_with_unknown.__wrapped__(protein, replace_proportion)
        unknown_proportion_in_msa = torch.bincount(protein['msa'].flatten()) / torch.numel(protein['msa'])
        unknown_proportion_in_seq = torch.bincount(protein['aatype'].flatten()) / torch.numel(protein['aatype'])

    def test_sample_msa(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        max_seq = 1000
        keep_extra = True
        protein = {}
        for k in MSA_FEATURE_NAMES:
            if k in features:
                protein[k] = torch.tensor(features[k])

        protein_processed = sample_msa.__wrapped__(protein.copy(), max_seq, keep_extra)
        for k in MSA_FEATURE_NAMES:
            if k in protein and keep_extra:
                assert protein_processed[k].shape[0] == min(protein[k].shape[0], max_seq)
                assert 'extra_'+k in protein_processed
                assert protein_processed['extra_'+k].shape[0] == protein[k].shape[0] - min(protein[k].shape[0], max_seq)

    def test_crop_extra_msa(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        max_extra_msa = 10
        protein = {'extra_msa': torch.tensor(features['msa'])}
        num_seq = protein["extra_msa"].shape[0]

        protein = crop_extra_msa.__wrapped__(protein, max_extra_msa)
        for k in MSA_FEATURE_NAMES:
            if "extra_" + k in protein:
                assert protein["extra_" + k].shape[0] == min(max_extra_msa, num_seq)

    def test_delete_extra_msa(self):
        protein = {'extra_msa': torch.rand((512, 100, 23))}
        extra_msa_has_deletion_shape = list(protein['extra_msa'].shape)
        extra_msa_has_deletion_shape[2] = 1
        protein['extra_deletion_matrix'] = torch.rand(extra_msa_has_deletion_shape)
        protein = delete_extra_msa(protein)
        for k in MSA_FEATURE_NAMES:
            assert 'extra_' + k not in protein
        assert 'extra_msa' not in protein

    def test_nearest_neighbor_clusters(self):
        with gzip.open('tests/test_data/sample_feats.pickle.gz', 'rb') as f:
            features = pickle.load(f)

        protein = {'msa': torch.tensor(features['true_msa'][0], dtype=torch.int64),
                   'msa_mask': torch.tensor(features['msa_mask'][0], dtype=torch.int64),
                   'extra_msa': torch.tensor(features['extra_msa'][0], dtype=torch.int64),
                   'extra_msa_mask': torch.tensor(features['extra_msa_mask'][0], dtype=torch.int64)}
        protein = nearest_neighbor_clusters.__wrapped__(protein, 0)
        assert 'extra_cluster_assignment' in protein

    def test_make_msa_mask(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        msa_mat = torch.tensor(features['msa'])
        protein = {'msa': msa_mat}
        protein = make_msa_mask(protein)
        assert 'msa_row_mask' in protein
        assert protein['msa_row_mask'].shape[0] == msa_mat.shape[0]

    def test_make_hhblits_profile(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa'], dtype=torch.int64)}
        protein = make_hhblits_profile(protein)
        assert 'hhblits_profile' in protein
        assert protein['hhblits_profile'].shape == torch.Size((protein['msa'].shape[1], 22))

    def test_make_masked_msa(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        protein = {'msa': torch.tensor(features['msa'], dtype=torch.int64)}
        protein = make_hhblits_profile(protein)
        masked_msa_config = config.data.common.masked_msa
        protein = make_masked_msa.__wrapped__(protein, masked_msa_config, replace_fraction=0.15)
        assert 'bert_mask' in protein
        assert 'true_msa' in protein
        assert 'msa' in protein
        assert protein['bert_mask'].sum() >= 0
        assert torch.all(torch.eq(
            protein['true_msa'] * (1-protein['bert_mask']), protein['msa'] * (1-protein['bert_mask'])))

    def test_make_msa_feat(self):
        with open('tests/test_data/features.pkl', 'rb') as file:
            features = pickle.load(file)

        protein = {'between_segment_residues': torch.tensor(features['between_segment_residues']),
                   'msa': torch.tensor(features['msa'], dtype=torch.int64),
                   'deletion_matrix': torch.tensor(features['deletion_matrix_int']),
                   'aatype': torch.argmax(torch.tensor(features['aatype']), dim=1)}
        protein = make_msa_feat.__wrapped__(protein)
        assert 'msa_feat' in protein
        assert 'target_feat' in protein
        assert protein['target_feat'].shape == torch.Size((protein['msa'].shape[1], 22))
        assert protein['msa_feat'].shape == torch.Size((*protein['msa'].shape, 25))

    def test_crop_templates(self):
        with gzip.open('tests/test_data/sample_feats.pickle.gz', 'rb') as f:
            features = pickle.load(f)

        protein = {'template_aatype': torch.tensor(features['true_msa'][0]),
                   'template_all_atom_masks': torch.tensor(features['msa_mask'][0])}
        max_templates = 2
        protein = crop_templates.__wrapped__(protein, max_templates)
        assert protein['template_aatype'].shape[0] == max_templates
        assert protein['template_all_atom_masks'].shape[0] == max_templates

    def test_make_atom14_masks(self):
        with gzip.open('tests/test_data/sample_feats.pickle.gz', 'rb') as file:
            features = pickle.load(file)

        protein = {'aatype': torch.tensor(features['aatype'][0])}
        protein = make_atom14_masks(protein)
        assert 'atom14_atom_exists' in protein
        assert 'residx_atom14_to_atom37' in protein
        assert 'residx_atom37_to_atom14' in protein
        assert 'atom37_atom_exists' in protein


if __name__ == '__main__':
    unittest.main()
