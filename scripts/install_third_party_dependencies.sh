#!/bin/bash

source scripts/vars.sh

apt update
apt install -y aria2
conda env update -n openfold_venv --file environment.yml && conda clean --all -y
conda install -y -c bioconda aria2
patch -p0 -d /opt/conda/lib/python3.8/site-packages/ < lib/openmm.patch

# install mmseqs2
pushd ~/
# git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2
#mkdir build
cd build
#cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
#make
make install
export PATH=$(pwd)/bin/:$PATH
echo PATH=$(pwd)/bin/:'$PATH' >> ~/.bashrc
popd


# lib/conda/bin/python3 -m pip install nvidia-pyindex

echo "Attempting to install FlashAttention"
git clone https://github.com/HazyResearch/flash-attention
CUR_DIR=$PWD
cd flash-attention
git checkout 5b838a8bef
python3 setup.py install
cd $CUR_DIR

# Download folding resources
# wget --no-check-certificate -P openfold/resources \
#     https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# # Certain tests need access to this file
# mkdir -p tests/test_data/alphafold/common
# ln -rs openfold/resources/stereo_chemical_props.txt tests/test_data/alphafold/common

# echo "Downloading OpenFold parameters..."
# bash scripts/download_openfold_params.sh openfold/resources

# echo "Downloading AlphaFold parameters..."
# bash scripts/download_alphafold_params.sh openfold/resources

# # Decompress test data
# gunzip tests/test_data/sample_feats.pickle.gz
