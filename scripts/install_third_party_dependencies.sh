#!/bin/bash
CONDA_INSTALL_URL=${CONDA_INSTALL_URL:-"https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"}

source scripts/vars.sh

# Install Miniconda locally
rm -rf lib/conda
rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh
wget -P /tmp \
    "${CONDA_INSTALL_URL}" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p lib/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

# Grab conda-only packages
export PATH=lib/conda/bin:$PATH
lib/conda/bin/python3 -m pip install nvidia-pyindex
conda env create --name=${ENV_NAME} -f environment.yml
source scripts/activate_conda_env.sh

echo "Attempting to install FlashAttention"
pip install git+https://github.com/HazyResearch/flash-attention.git@5b838a8bef78186196244a4156ec35bbb58c337d && echo "Installation successful"

# Install DeepMind's OpenMM patch
OPENFOLD_DIR=$PWD
pushd lib/conda/envs/$ENV_NAME/lib/python3.7/site-packages/ \
    && patch -p0 < $OPENFOLD_DIR/lib/openmm.patch \
    && popd

# Download folding resources
wget --no-check-certificate -P openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Certain tests need access to this file
mkdir -p tests/test_data/alphafold/common
ln -rs openfold/resources/stereo_chemical_props.txt tests/test_data/alphafold/common

echo "Downloading OpenFold parameters..."
bash scripts/download_openfold_params.sh openfold/resources

echo "Downloading AlphaFold parameters..."
bash scripts/download_alphafold_params.sh openfold/resources

# Decompress test data
gunzip tests/test_data/sample_feats.pickle.gz
