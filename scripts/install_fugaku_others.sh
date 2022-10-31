#!/bin/bash -eux
#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "elapse=2:0:00"
#PJM --rsc-list "retention_state=0"
#PJM -j

. scripts/setenv

cd $WORKDIR

# Bazel
wget -O $PREFIX/bin/bazel https://github.com/bazelbuild/bazelisk/releases/download/v1.11.0/bazelisk-linux-arm64
chmod +x $PREFIX/bin/bazel

# Doxygen
git clone https://github.com/doxygen/doxygen.git -b Release_1_9_4
pushd doxygen
mkdir build
cd build
cmake -D CMAKE_INSTALL_PREFIX=$PREFIX -G "Unix Makefiles" ..
make -j48
make install
popd

# OpenMM
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
git clone https://github.com/EddyRivasLab/hmmer.git -b h3-arm
pushd hmmer
git clone https://github.com/EddyRivasLab/easel -b develop
sed -i "s|AC_PREREQ(\[2.71\])|AC_PREREQ(\[2.69\])|g" configure.ac
autoconf
./configure --prefix=$PREFIX
make -j48
make check
make install
popd

# HH suite
git clone https://github.com/soedinglab/hh-suite.git -b v3.3.0
mkdir -p hh-suite/build
pushd hh-suite/build
cmake -D CMAKE_INSTALL_PREFIX=$PREFIX -G "Unix Makefiles" ..
make -j48
make install
popd

# kalign2
mkdir kalign2
pushd kalign2
wget http://msa.sbc.su.se/downloads/kalign/current.tar.gz
tar xf current.tar.gz
patch -p1 < $OPENFOLDDIR/lib/kalign2.patch
./configure --prefix=$PREFIX
make -j48
make install
popd

# pip packages
pip3 install -r $OPENFOLDDIR/scripts/install_fugaku_requirements.txt

# patch pytorch-lightning
pushd $(pip show pytorch-lightning | sed -n -e "s/Location: \(.*\)/\1/p")/pytorch_lightning
patch -p0 < $OPENFOLDDIR/lib/pytorch_lightning.patch
popd

# DeepSpeed
pushd DeepSpeedCPU
DS_BUILD_UTILS=1 pip3 install .
popd

# OpenFold
pushd $OPENFOLDDIR
python3 setup.py install
popd

echo "Finished"
