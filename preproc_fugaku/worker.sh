#!/bin/bash

. ../scripts/setenv

set -eu

InputDir=$1
OutputDir=$2
ToolTimeLimit=$3
NumNodes=$4
NumProcs=$5
NumThreads=$6
Mode=$7
ScriptArgs=$8
DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

echo InputDir=$InputDir
echo OutputDir=$OutputDir
echo ToolTimeLimit=$ToolTimeLimit
echo NumNodes=$NumNodes
echo NumProcs=$NumProcs
echo NumThreads=$NumThreads
echo Mode=$Mode
echo ScriptArgs=$ScriptArgs

# DBs
Uniref90=$DataDir/uniref90/uniref90.fasta
Pdb70=$DataDir/pdb70
Mgnify=$DataDir/mgnify/mgy_clusters_2018_12.fa

SmallBfd=$DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta

if [[ $Mode = "bfd_deflate" ]]; then
    if [[ ${IS_FUGAKU} == 1 ]]; then
	# TODO: Move to share
	BfdA3M=/vol0006/ra000018/data/u10130/hhblits_deflate/bfd
    else
	# TODO: Do staging for a3m.ffindex
	BfdA3M=$DataDir/hhblits_deflate/bfd
    fi
else
    BfdA3M=""
fi
Bfd=$DataDir/bfd
Uniclust30=$DataDir/uniclust30/uniclust30_2018_08
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
elif [[ $Mode = "bfd" || $Mode = "bfd_deflate" ]]; then
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
if [[ $Mode = "bfd" || $Mode = "bfd_deflate" ]]; then
    if [[ ${IS_FUGAKU} == 1 ]]; then

	if [[ $Mode = "bfd" ]]; then
	    echo "Mode must be bfd_deflate on Fugaku" >&2
	    exit
	fi

	echo "Executing LLIO transfer for BFD files"
	echo "Bfd=${Bfd}"
	set +e
	time $DirTransfer $Bfd
	set -e

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
	Exts=(ffdata ffindex)
	for Ext in ${Exts[@]}; do
	    set +e
	    time $LLIOTransfer $Database/*.$Exts
	    set -e
	done

    else
	echo "Copying the DB directory to local temp"
	DatabaseLocal=$LOCALTMP/`basename $Database`
	time $MpiExecEachNode mkdir -p $DatabaseLocal
	time $MpiExecEachNode cp -r $Database/*.ffdata $DatabaseLocal/
	time $MpiExecEachNode cp -r $Database/*.ffindex $DatabaseLocal/
	Database=$DatabaseLocal

	# FIXME: remove
	echo "DatabaseLocal:"
	ls $DatabaseLocal/

    fi
fi

echo "Setting DatabaseArgs"
if [[ $Mode = "uniref90" ]]; then
    DatabaseArgs="--uniref90_database_path $Database"
elif [[ $Mode = "pdb70" ]]; then
    DatabaseArgs="--pdb70_database_path $Database/pdb70"
elif [[ $Mode = "mgnify" ]]; then
    DatabaseArgs="--mgnify_database_path $Database"
elif [[ $Mode = "small_bfd" ]]; then
    DatabaseArgs="--bfd_database_path $Database --use_small_bfd True"
elif [[ $Mode = "bfd" || $Mode = "bfd_deflate" ]]; then
    DatabaseArgs="--uniclust30_database_path $Database/uniclust30_2018_08 --bfd_database_path $Bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    if [[ -n $BfdA3M ]]; then
	DatabaseArgs="$DatabaseArgs --bfd_deflate_a3m_path $BfdA3M/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    fi
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi

echo "DatabaseArgs=$DatabaseArgs"

echo "Creating output directory"
mkdir -p $OutputDir


if [[ $Mode = "bfd_deflate" ]]; then
    # TODO: Do not use hhblits in the home directory
    if [[ ${IS_FUGAKU} == 1 ]]; then
	export PATH=/vol0006/ra000018/data/u10130/repo/hh-suite/build/src:${PATH}
    else
	export PATH=/home/acd13359ld/repo/hh-suite/build/src:${PATH}
    fi
    echo hhblits=`which hhblits`
fi

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
    --hhblits_binary_path `which hhblits` \
    --hhsearch_binary_path $BinDir/hhsearch \
    --kalign_binary_path $BinDir/kalign \
    --timeout $ToolTimeLimit \
    $DatabaseArgs \
    $ScriptArgs
