#!/bin/bash

. ../scripts/setenv

set -eu

InputDir=$1
OutputDir=$2
NumNodes=$3
NumProcs=$4
NumThreads=$5
Mode=$6
ScriptArgs=$7
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

echo InputDir=$InputDir
echo OutputDir=$OutputDir
echo NumNodes=$NumNodes
echo NumProcs=$NumProcs
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
NumTaskProcs=`expr $NumProcs \* $NumNodes`

if [[ ${IS_FUGAKU} == 1 ]]; then
    MpiExecEachNode="STUB"
    MpiExecTask="mpiexec -n $NumTaskProcs"
    MpiArgs=""
else
    MpiExecEachNode="mpiexec -np $NumNodes -npernode 1"
    MpiExecTask="mpiexec -np $NumTaskProcs -npernode $NumProcs"
    MpiArgs="--bind-to none"
fi

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

echo LOCALTMP=$LOCALTMP

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

    # FIXME: The following dirs fail because some files include a space
    # time $DirTransfer $PREFIX/lib
    # time $DirTransfer $PREFIX/bin

    set +eu
    time $DirTransfer $OPENFOLDDIR/scripts
    time $DirTransfer $OPENFOLDDIR/openfold
    time $DirTransfer $OPENFOLDDIR/openfold.egg-info
    set -eu
fi

DatabaseArgs=""
if [[ $Mode = "bfd" ]]; then
    if [[ ${IS_FUGAKU} == 1 ]]; then
	echo "Executing LLIO transfer for BFD files"
	for ffindex in `find $Bfd -name *ffindex`; do
	    time $LLIOTransfer $ffindex
	done
	time $LLIOTransfer $BfdCs219Data
	time $DirTransfer $Bfd

    else
	BfdLocal=$LOCALTMP/`basename $Bfd`

	echo Bfd=$Bfd
	echo BfdLocal=$BfdLocal
	echo BfdCs219Data=$BfdCs219Data

	echo "Executing prepare_bfd"
	time $MpiExecEachNode \
	     ./prepare_bfd.sh \
	     $Bfd \
	     $BfdLocal \
	     $BfdCs219Data
	Bfd=$BfdLocal
    fi
fi

if [[ $DatabaseIsFile = 1 ]]; then
    if [[ ${IS_FUGAKU} == 1 ]]; then
	echo "Executing LLIO transfer for the DB file"
	time $LLIOTransfer $Database

    else
	echo "Copying the DB file to local temp"
	DatabaseLocal=$LOCALTMP/`basename $Database`
	time $MpiExecEachNode cp $Database $DatabaseLocal
	Database=$DatabaseLocal
    fi

else
    if [[ ${IS_FUGAKU} == 1 ]]; then
	echo "Executing LLIO transfer for the DB directory"
	time $DirTransfer $Database

    else
	echo "Copying the DB directory to local temp"
	DatabaseLocal=$LOCALTMP/`basename $Database`
	time $MpiExecEachNode cp -r $Database $DatabaseLocal
	Database=$DatabaseLocal
    fi
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
export PARALLEL=$OMP_NUM_THREADS
$MpiExecTask \
    $MpiArgs \
    -x PARALLEL \
    -x OMP_NUM_THREADS \
    python3 $Script \
    $InputDir $OutputDir \
    --cpus_per_task $NumThreads \
    --jackhmmer_binary_path $BinDir/jackhmmer \
    --hhblits_binary_path $BinDir/hhblits \
    --hhsearch_binary_path $BinDir/hhsearch \
    --kalign_binary_path $BinDir/kalign \
    $DatabaseArgs \
    $ScriptArgs
