#!/bin/bash -eux
#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "elapse=2:0:00"
#PJM --rsc-list "retention_state=0"
#PJM -j
#
# Copyright 2023 RIKEN & Fujitsu Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

. scripts/setenv

set -eu

mkdir -p $WORKDIR
cd $WORKDIR

# Bazel
wget -O $PREFIX/bin/bazel https://github.com/bazelbuild/bazelisk/releases/download/v1.11.0/bazelisk-linux-arm64
chmod +x $PREFIX/bin/bazel
export USE_BAZEL_VERSION=5.4.0

# Doxygen
rm -rf doxygen
git clone https://github.com/doxygen/doxygen.git -b Release_1_9_4
pushd doxygen
mkdir build
cd build
cmake -D CMAKE_INSTALL_PREFIX=$PREFIX -G "Unix Makefiles" ..
make -j48
make install
popd

# OpenMM
rm -rf openmm
git clone https://github.com/openmm/openmm -b 7.7.0
pushd openmm
mkdir build
cd build
(
    export CC="fcc -Nclang -Kfast -Knolargepage -lpthread -Kopenmp"
    export CXX="FCC -Nclang -Kfast -Knolargepage -lpthread -Kopenmp"
    cmake -D CMAKE_INSTALL_PREFIX=$PREFIX -G "Unix Makefiles" ..
    make -j48
    make install
    make PythonInstall
)
popd

# HMMER
rm -rf hmmer
git clone https://github.com/EddyRivasLab/hmmer.git -b h3-arm
pushd hmmer
git clone https://github.com/EddyRivasLab/easel -b develop
pushd easel
git checkout 367f817 # 2023/09/22
popd
sed -i "s|AC_PREREQ(\[2.71\])|AC_PREREQ(\[2.69\])|g" configure.ac
autoconf
./configure --prefix=$PREFIX
make -j48
make check
make install
popd

# HH suite
rm -rf hh-suite
git clone https://github.com/soedinglab/hh-suite.git -b v3.3.0
mkdir -p hh-suite/build
pushd hh-suite/build
cmake -D CMAKE_INSTALL_PREFIX=$PREFIX -G "Unix Makefiles" ..
make -j48
make install
popd

# kalign2
rm -rf kalign2
mkdir kalign2
pushd kalign2
wget http://msa.sbc.su.se/downloads/kalign/current.tar.gz
tar xf current.tar.gz
sed -i 's/=\ gcc/=\ @CC@/g' Makefile.in
sed -i 's/-O9\ \ -Wall/@CFLAGS@/g' Makefile.in
./configure
make -j48
mkdir -p $PREFIX/bin
cp kalign $PREFIX/bin/
popd

# pip packages
pip3 install -r $OPENFOLDDIR/scripts/install_fugaku_requirements.txt

# patch pytorch-lightning
pushd $(pip3 show pytorch-lightning | sed -n -e "s/Location: \(.*\)/\1/p")/pytorch_lightning
patch -p0 < $OPENFOLDDIR/lib/pytorch_lightning.patch
popd

# DeepSpeed
rm -rf DeepSpeed
git clone https://github.com/microsoft/DeepSpeed.git
pushd DeepSpeed
git checkout b4e5826a
patch -p1 < $OPENFOLDDIR/lib/deepspeed.patch
DS_BUILD_UTILS=1 pip3 install .
popd

# mpi4py
env MPICC=`which mpifcc` pip3 install mpi4py

# OpenFold
pushd $OPENFOLDDIR
pip3 install -e .
popd

echo "Finished"
