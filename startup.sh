#!/usr/bin/env bash

# Export mmseqs path
export PATH=/MMseqs2/build/bin/:$PATH

# Install DeepMind's OpenMM patch
patch -p0 -d /tools/miniconda/envs/ottf-openfold/lib/python3.7/site-packages/ < /openfold/lib/openmm.patch

