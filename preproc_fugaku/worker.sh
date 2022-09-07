#!/bin/bash

. ../scripts/setenv

set -eu

InputDir=$1
OutputDir=$2
Mode=$3
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

Mgnify=$DataDir/mgnify/mgy_clusters_2018_12.fa

DirTransfer=/home/system/tool/dir_transfer

DatabaseArgs=""
if [[ ${Mode} == "mgnify" ]]; then
    # Move DB from 2ndfs to cache
    $DirTransfer $Mgnify

    # Move cache to local temp
    # FIXME: This works only when one node is used
    MgnifyLocal=$PJM_LOCALTMP/mgnify
    cp -r $Mgnify=$MgnifyLocal
    MgnifyLocal=$Mgnify

    DatabaseArgs="--mgnify_database_path $Mgnify"
else
    echo "Invalid Mode: ${Mode}" >&2
    exit
fi

mpiexec \
    -n 1 \
    python3 ../scripts/precompute_alignments_fugaku.py $InputDir $OutputDir \
    --cpus_per_task 48 \
    --jackhmmer_binary_path $BinDir/jackhmmer \
    --hhblits_binary_path $BinDir/hhblits \
    --hhsearch_binary_path $BinDir/hhsearch \
    --kalign_binary_path $BinDir/kalign \
    --use_small_bfd \
    $DatabaseArgs

# --uniref90_database_path $DataDir/uniref90/uniref90.fasta \
    # --mgnify_database_path $Mgnify \
    # --pdb70_database_path $DataDir/pdb70/pdb70 \
    # --uniclust30_database_path $DataDir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    # --bfd_database_path $DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta \
