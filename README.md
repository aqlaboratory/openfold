# Permutation code README

## Overview:

NB: before running the test codes,please download the procrustes package first:
from https://github.com/theochem/procrustes


The only file that has been changed is:
[```openfold/utils/loss.py```](https://github.com/dingquanyu/openfold/blob/permutation/openfold/utils/loss.py), in which the forward function is modified in 
original ```AlphaFoldLoss``` class; create a child class called ```AlphaFoldMultimerLoss``` that not only inherited all the loss calculations but also 
has multi-chain permutation codes.

These files are newly added:
[```tests/test_permutation.py```](https://github.com/dingquanyu/openfold/blob/permutation/tests/test_permutation.py): A unittest script 
that tests permutation functions.

[```tests/test_data/label_1.pkl```](https://github.com/dingquanyu/openfold/blob/permutation/tests/test_data/label_1.pkl) 
and [```tests/test_data/label_2.pkl```](https://github.com/dingquanyu/openfold/blob/permutation/tests/test_data/label_2.pkl) are 2 fake ground truth structures.
```label_1.pkl``` has 9 residues and ```label_2.pkl``` has 13 residues
29/06/23 Fill NaN in the lddt scores with the matrix mean for now because the test data are randomly generated and it gives NaN in the lddt score somehow.
**Delete** this step before merging to Multimer branch
