import numpy as np
import torch

from np import residue_constants


def cast_to_64bit_ints(protein):
    # We keep all ints as int64
    for k, v in protein.items():
        if v.dtype == torch.int32:
            protein[k] = v.type(torch.int64)
    return protein

def make_seq_mask(protein):
    protein['seq_mask'] = torch.ones(protein['aatype'].shape, dtype=torch.float32)
    return protein

def make_template_mask(protein):
    protein['template_mask'] = torch.ones(protein['template_aatype'].shape[0], dtype=torch.float32)
    return protein

def curry1(f):
  """Supply all arguments but the first."""

  def fc(*args, **kwargs):
    return lambda x: f(x, *args, **kwargs)

  return fc

@curry1
def add_distillation_flag(protein, distillation):
    protein['is_distillation'] = torch.tensor(float(distillation), dtype=torch.float32)
    return protein

def make_all_atom_aatype(protein):
    protein['all_atom_aatype'] = protein['aatype']
    return protein

def fix_templates_aatype(protein):
    # Map one-hot to indices
    num_templates = protein['template_aatype'].shape[0]
    protein['template_aatype'] = torch.argmax(protein['template_aatype'], dim=-1)
    # Map hhsearch-aatype to our aatype.
    new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
    new_order = torch.tensor(new_order_list, dtype=torch.int32).expand(num_templates, -1)
    protein['template_aatype'] = torch.gather(new_order, 1, index=protein['template_aatype'])
    return protein

def correct_msa_restypes(protein):
    """Correct MSA restype to have the same order as residue_constants."""
    new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
    new_order = torch.tensor([new_order_list]*protein['msa'].shape[1], dtype=protein['msa'].dtype).transpose(0,1)
    protein['msa'] = torch.gather(new_order, 0, protein['msa'])

    perm_matrix = np.zeros((22, 22), dtype=np.float32)
    perm_matrix[range(len(new_order_list)), new_order_list] = 1.

    for k in protein:
        if 'profile' in k:
            num_dim = protein[k].shape.as_list()[-1]
            assert num_dim in [20,21,22], (
                'num_dim for %s out of expected range: %s' % (k, num_dim))
            protein[k] = torch.dot(protein[k], perm_matrix[:num_dim, :num_dim])
    return protein

def squeeze_features(protein):
    """Remove singleton and repeated dimensions in protein features."""
    protein['aatype'] = torch.argmax(protein['aatype'], dim=-1)
    for k in [
            'domain_name', 'msa', 'num_alignments', 'seq_length', 'sequence',
            'superfamily', 'deletion_matrix', 'resolution',
            'between_segment_residues', 'residue_index', 'template_all_atom_masks']:
        if k in protein:
            final_dim = protein[k].shape[-1]
            if isinstance(final_dim, int) and final_dim == 1:
                protein[k] = torch.squeeze(protein[k], dim=-1)

    for k in ['seq_length', 'num_alignments']:
        if k in protein:
            protein[k] = protein[k][0]
    return protein

def make_protein_crop_to_size_seed(protein):
    protein['random_crop_to_size_seed'] = torch.distributions.Uniform(low=torch.int32, high=torch.int32).sample((2))
    return protein
