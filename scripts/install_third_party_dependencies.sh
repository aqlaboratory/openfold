#!/bin/bash

# Download folding resources
wget -N --no-check-certificate -P openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Certain tests need access to this file
mkdir -p tests/test_data/alphafold/common
ln -rs openfold/resources/stereo_chemical_props.txt tests/test_data/alphafold/common

# echo "Downloading OpenFold parameters..."
# bash scripts/download_openfold_params.sh openfold/resources
# 
# echo "Downloading AlphaFold parameters..."
# bash scripts/download_alphafold_params.sh openfold/resources

# Decompress test data
gunzip -c tests/test_data/sample_feats.pickle.gz > tests/test_data/sample_feats.pickle

python setup.py install
