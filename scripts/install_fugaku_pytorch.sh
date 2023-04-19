#!/bin/bash -eux
#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "elapse=4:0:00"
#PJM --rsc-list "retention_state=0"
#PJM -j

. scripts/setenv

set -eu

cd $WORKDIR

rm -rf pytorch
git clone -b fujitsu_v1.10.1_for_a64fx https://github.com/fujitsu/pytorch.git

# apply gather&scatter improve patch
pushd pytorch
git show 0c49800 --unified=0 > gather_scatter.patch
sed -i "s|@@ -134|@@ -133|g" gather_scatter.patch
patch -p1 < gather_scatter.patch
popd

pushd pytorch/scripts/fujitsu

# set config
sed -i "s|TCSDS_PATH=/opt/FJSVstclanga/cp-1.0.21.01|TCSDS_PATH=/opt/FJSVxtclanga/tcsds-1.2.35|g" env.src
sed -i "s|PREFIX=~/prefix|PREFIX=$PREFIX|g" env.src
sed -i "s|VENV_PATH=~/venv||g" env.src
# force flushing denormal numbers to zero
sed -i "s|CFLAGS=-O3 CXXFLAGS=-O3|CFLAGS=-Kfast|g" 5_pytorch.sh

export fjenv_use_venv=false
bash 1_python.sh
bash 3_venv.sh
bash 4_numpy_scipy.sh
bash 5_pytorch.sh
bash 6_vision.sh
bash 8_libtcmalloc.sh

# make symbolic link (python -> python3)
pushd $PREFIX/bin
ln -s python3 python
popd

echo "Finished"
