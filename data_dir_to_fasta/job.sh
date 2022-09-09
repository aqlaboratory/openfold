#!/bin/bash

source ../scripts/setenv
python3 ../scripts/data_dir_to_fasta.py \
    /vol0003/share/hp210263/OpenFold/data/pdb_mmcif/mmcif_files \
    pdb_mmcif.fasta
