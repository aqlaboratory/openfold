#!/bin/bash

. ../scripts/setenv

InputDir=$1
OutputDir=$2
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

python3 ../scripts/precompute_alignments_smallbfd.py $InputDir $OutputDir \
	--uniref90_database_path $DataDir/uniref90/uniref90.fasta \
	--mgnify_database_path $DataDir/mgnify/mgy_clusters_2018_12.fa \
	--pdb70_database_path $DataDir/pdb70/pdb70 \
	--uniclust30_database_path $DataDir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
	--cpus_per_task 48 \
	--jackhmmer_binary_path $BinDir/jackhmmer \
	--hhblits_binary_path $BinDir/hhblits \
	--hhsearch_binary_path $BinDir/hhsearch \
	--kalign_binary_path $BinDir/kalign \
	--bfd_database_path $DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta 
