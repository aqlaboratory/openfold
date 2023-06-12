#!/bin/bash
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
