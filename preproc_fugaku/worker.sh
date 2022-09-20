#!/bin/bash

. ../scripts/setenv

set -eu

InputDir=$1
OutputDir=$2
Mode=$3
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

# DBs
Uniref90=$DataDir/uniref90/uniref90.fasta
Pdb70=$DataDir/pdb70
Mgnify=$DataDir/mgnify/mgy_clusters_2018_12.fa
Uniclust30=$DataDir/uniclust30/uniclust30_2018_08

SmallBfd=$DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta
Bfd=$DataDir/bfd
BfdCs219Data=$Bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt_cs219.ffdata

# commands
Script=../scripts/precompute_alignments_fugaku.py
LLIOTransfer=llio_transfer
DirTransfer=/home/system/tool/dir_transfer
MpiExecEachNode=mpiexec

if [[ $Mode = "uniref90" ]]; then
    Database=$Uniref90
    DatabaseIsFile=1
elif [[ $Mode = "pdb70" ]]; then
    Database=$Pdb70
    DatabaseIsFile=0
elif [[ $Mode = "mgnify" ]]; then
    Database=$Mgnify
    DatabaseIsFile=1
elif [[ $Mode = "small_bfd" ]]; then
    Database=$SmallBfd
    DatabaseIsFile=1
elif [[ $Mode = "bfd" ]]; then
    Database=$Uniclust30
    DatabaseIsFile=0
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi

echo "Executing LLIO transfer for OpenFold"
if [[ -v IS_FUGAKU ]]; then
    time $DirTransfer $PREFIX/lib
    time $DirTransfer $PREFIX/bin
    time $DirTransfer $OPENFOLDDIR/scripts
    time $DirTransfer $OPENFOLDDIR/openfold
    time $DirTransfer $OPENFOLDDIR/openfold.egg-info
fi

DatabaseArgs=""
if [[ $Mode = "bfd" ]]; then
    echo "Executing LLIO transfer for a BFD file"
    time $DirTransfer $Bfd

    $MpiExecEachNode \
	./prepare_bfd.sh \
	$Bfd \
	$BfdLocal \
	$BfdCs219Data

    Bfd=$BfdLocal
fi

if [[ -v IS_FUGAKU ]]; then
    if [[ $DatabaseIsFile = 1 ]]; then
	echo "Executing LLIO transfer for the DB file"
	time $LLIOTransfer $Database

	echo "Copying the DB file to local temp"
	DatabaseLocal=$PJM_LOCALTMP/`basename $Database`
	time $MpiExecEachNode cp $Database $DatabaseLocal

	Database=$DatabaseLocal

    else
	echo "Executing LLIO transfer for the DB directory"
	time $DirTransfer $Database

	echo "Copying the DB directory to local temp"
	DatabaseLocal=$PJM_LOCALTMP/`basename $Database`
	time $MpiExecEachNode cp -r $Database $DatabaseLocal

	Database=$DatabaseLocal
    fi
fi

# FIXME: remove
echo Database=$Database

if [[ $Mode = "uniref90" ]]; then
    DatabaseArgs="--uniref90_database_path $Database"
elif [[ $Mode = "pdb70" ]]; then
    DatabaseArgs="--pdb70_database_path $Database"
elif [[ $Mode = "mgnify" ]]; then
    DatabaseArgs="--mgnify_database_path $Database"
elif [[ $Mode = "small_bfd" ]]; then
    DatabaseArgs="--bfd_database_path $Database --use_small_bfd True"
elif [[ $Mode = "bfd" ]]; then
    DatabaseArgs="--uniclust30_database_path $Database --bfd_database_path $Bfd"
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi


echo "Starting the script"
mpiexec \
    -n 1 \
    python3 $Script \
    $InputDir $OutputDir \
    --cpus_per_task 48 \
    --jackhmmer_binary_path $BinDir/jackhmmer \
    --hhblits_binary_path $BinDir/hhblits \
    --hhsearch_binary_path $BinDir/hhsearch \
    --kalign_binary_path $BinDir/kalign \
    $DatabaseArgs

# --uniref90_database_path $DataDir/uniref90/uniref90.fasta \
    # --mgnify_database_path $Mgnify \
    # --pdb70_database_path $DataDir/pdb70/pdb70 \
    # --uniclust30_database_path $DataDir/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    # --bfd_database_path $DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta \
