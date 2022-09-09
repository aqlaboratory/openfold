#!/bin/bash

. ../scripts/setenv

set -eu

# FIXME
mpiexec hostname

InputDir=$1
OutputDir=$2
Mode=$3
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

Uniref90=$DataDir/uniref90/uniref90.fasta
Pdb70=$DataDir/pdb70
Mgnify=$DataDir/mgnify/mgy_clusters_2018_12.fa

LLIOTransfer=llio_transfer
DirTransfer=/home/system/tool/dir_transfer

if [[ $Mode = "uniref90" ]]; then
    Database=$Uniref90
    DatabaseIsFile=1
elif [[ $Mode = "pdb70" ]]; then
    Database=$Pdb70
    DatabaseIsFile=0
elif [[ $Mode = "mgnify" ]]; then
    Database=$Mgnify
    DatabaseIsFile=1
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi

echo "Executing LLIO transfer for OpenFold"
time $DirTransfer $PREFIX
time $DirTransfer $OPENFOLDDIR/openfold
time $DirTransfer $OPENFOLDDIR/openfold.egg-info

DatabaseArgs=""
if [[ $DatabaseIsFile = 1 ]]; then
    # Move DB from 2ndfs to cache
    echo "Executing LLIO transfer for the database"
    time $LLIOTransfer $Database

    # Move cache to local temp
    # FIXME: This works only when one node is used
    echo "Copying DB to local temp"
    DatabaseLocal=$PJM_LOCALTMP/`basename $Database`
    time cp $Database $DatabaseLocal
    Database=$DatabaseLocal

else
    echo "Executing LLIO transfer for the database"
    time $DirTransfer $Database

    # Move cache to local temp
    # FIXME: This works only when one node is used
    echo "Copying DB to local temp"
    DatabaseLocal=$PJM_LOCALTMP/`basename $Database`
    time cp -r $Database $DatabaseLocal
    Database=$DatabaseLocal
fi

# FIXME: remove
echo Database=$Database

if [[ $Mode = "uniref90" ]]; then
    DatabaseArgs="--uniref90_database_path $Database"
elif [[ $Mode = "pdb70" ]]; then
    DatabaseArgs="--pdb70_database_path $Database"
elif [[ $Mode = "mgnify" ]]; then
    DatabaseArgs="--mgnify_database_path $Database"
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi


echo "Starting the script"
mpiexec \
    -n 1 \
    python3 ../scripts/precompute_alignments_fugaku.py $InputDir $OutputDir \
    --cpus_per_task 48 \
    --jackhmmer_binary_path $BinDir/jackhmmer \
    --hhblits_binary_path $BinDir/hhblits \
    --hhsearch_binary_path $BinDir/hhsearch \
    --kalign_binary_path $BinDir/kalign \
    --use_small_bfd True \
    $DatabaseArgs

# --uniref90_database_path $DataDir/uniref90/uniref90.fasta \
    # --mgnify_database_path $Mgnify \
    # --pdb70_database_path $DataDir/pdb70/pdb70 \
    # --uniclust30_database_path $DataDir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    # --bfd_database_path $DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta \
