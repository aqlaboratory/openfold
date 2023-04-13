#!/bin/bash

. ../scripts/setenv

hostname

InputDir=$1
OutputDir=$2
InputFileList=$3
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

export SLURM_JOB_NUM_NODES=$GLOBAL_SIZE  #$OMPI_MCA_orte_ess_num_procs
export SLURM_NODEID=$((GLOBAL_RANK_OFFSET + PMIX_RANK))

if [ "$InputFileList" == "" ]; then
    FileListArg=""
else
    FileListArg="--mmcif_file_list $InputFileList"
fi

python3 ../scripts/precompute_alignments_smallbfd.py $InputDir $OutputDir \
	--uniref90_database_path $DataDir/uniref90/uniref90.fasta \
	--mgnify_database_path $DataDir/mgnify/mgy_clusters_2018_12.fa \
	--pdb70_database_path $DataDir/pdb70/pdb70 \
	--uniclust30_database_path $DataDir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
	--cpus_per_task 48 \
	--no_tasks 1 \
	--jackhmmer_binary_path $BinDir/jackhmmer \
	--hhblits_binary_path $BinDir/hhblits \
	--hhsearch_binary_path $BinDir/hhsearch \
	--kalign_binary_path $BinDir/kalign \
	--bfd_database_path $DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta \
	$FileListArg
