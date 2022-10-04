#!/bin/bash

. ../scripts/setenv

set -eu

InputDir=$1
OutputDir=$2
NumNodes=$3
NumThreads=$4
Mode=$5
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

echo InputDir=$InputDir
echo OutputDir=$OutputDir
echo NumNodes=$NumNodes
echo NumThreads=$NumThreads
echo Mode=$Mode

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
MpiExecEachNode="mpiexec -np $NumNodes -npernode 1"
NumProcsPerNode=`expr $NUM_CORES / $NumThreads`
# NumProcsPerNode=1
NumTaskProcs=`expr $NumProcsPerNode \* $NumNodes`
MpiExecTask="mpiexec -np $NumTaskProcs -npernode $NumProcsPerNode"

echo NumProcsPerNode=$NumProcsPerNode
echo NumTaskProcs=$NumTaskProcs
echo MpiExecEachNode=$MpiExecEachNode
echo MpiExecTask=$MpiExecTask

# Temp locations
if [[ ${IS_FUGAKU} == 1 ]]; then
    LOCALTMP=$PJM_LOCALTMP # Shared SSD, 87 GiB/node
elif [[ ${IS_ABCI} == 1 ]]; then
    rm -rf /dev/shm/* 2>/dev/null || true
    LOCALTMP=/dev/shm # Private memory, 384 GiB/node (V)
    # LOCALTMP=$SGE_LOCALDIR # Private SSD, 1.6 TB
fi

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

if [[ ${IS_FUGAKU} == 1 ]]; then
    echo "Executing LLIO transfer for OpenFold"
    time $DirTransfer $PREFIX/lib
    time $DirTransfer $PREFIX/bin
    time $DirTransfer $OPENFOLDDIR/scripts
    time $DirTransfer $OPENFOLDDIR/openfold
    time $DirTransfer $OPENFOLDDIR/openfold.egg-info
fi

DatabaseArgs=""
if [[ $Mode = "bfd" ]]; then
    if [[ ${IS_FUGAKU} == 1 ]]; then
	echo "Executing LLIO transfer for a BFD file"
	time $DirTransfer $Bfd
    fi

    BfdLocal=$LOCALTMP/`basename $Bfd`

    # TODO: remove
    echo Bfd=$Bfd
    echo BfdLocal=$BfdLocal
    echo BfdCs219Data=$BfdCs219Data

    time $MpiExecEachNode \
	 ./prepare_bfd.sh \
	 $Bfd \
	 $BfdLocal \
	 $BfdCs219Data
    Bfd=$BfdLocal
fi

if [[ $DatabaseIsFile = 1 ]]; then
    if [[ ${IS_FUGAKU} == 1 ]]; then
	echo "Executing LLIO transfer for the DB file"
	time $LLIOTransfer $Database
    fi

    echo "Copying the DB file to local temp"
    DatabaseLocal=$LOCALTMP/`basename $Database`
    time $MpiExecEachNode cp $Database $DatabaseLocal
    Database=$DatabaseLocal

else
    if [[ ${IS_FUGAKU} == 1 ]]; then
	echo "Executing LLIO transfer for the DB directory"
	time $DirTransfer $Database
    fi

    echo "Copying the DB directory to local temp"
    DatabaseLocal=$LOCALTMP/`basename $Database`
    time $MpiExecEachNode cp -r $Database $DatabaseLocal
    Database=$DatabaseLocal
fi

if [[ $Mode = "uniref90" ]]; then
    DatabaseArgs="--uniref90_database_path $Database"
elif [[ $Mode = "pdb70" ]]; then
    DatabaseArgs="--pdb70_database_path $Database/pdb70"
elif [[ $Mode = "mgnify" ]]; then
    DatabaseArgs="--mgnify_database_path $Database"
elif [[ $Mode = "small_bfd" ]]; then
    DatabaseArgs="--bfd_database_path $Database --use_small_bfd True"
elif [[ $Mode = "bfd" ]]; then
    DatabaseArgs="--uniclust30_database_path $Database/uniclust30_2018_08 --bfd_database_path $Bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi

mkdir -p $OutputDir

echo "Starting the script"
export OMP_NUM_THREADS=$NumThreads # Define just in case it is used
$MpiExecTask \
    --bind-to none \
    -x OMP_NUM_THREADS \
    python3 $Script \
    $InputDir $OutputDir \
    --cpus_per_task $NumThreads \
    --jackhmmer_binary_path $BinDir/jackhmmer \
    --hhblits_binary_path $BinDir/hhblits \
    --hhsearch_binary_path $BinDir/hhsearch \
    --kalign_binary_path $BinDir/kalign \
    $DatabaseArgs \
    --ignore_unique # TODO: remove in real runs
