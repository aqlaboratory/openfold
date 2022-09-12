#!/bin/bash

source ../scripts/setenv

set -eu

PDB_MMCIF=/vol0003/share/hp210263/OpenFold/data/pdb_mmcif/mmcif_files
# PDB_MMCIF=/vol0003/hp210263/data/oyama/repo/openfold/data_dir_to_fasta_fugaku/test_dir
DirTransfer=/home/system/tool/dir_transfer

echo "Executing LLIO transfer"
time $DirTransfer $PDB_MMCIF

echo "Starting the script"
python3 ../scripts/data_dir_to_fasta_fugaku.py \
	--no_tasks 48 \
	$PDB_MMCIF \
	pdb_mmcif.fasta
