#!/bin/bash
#
# Copyright 2023 Fujitsu Limited
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

set -eu

DataDir=$DATADIR
BinDir=$PREFIX/bin
Log=log_`date +%Y%m%d_%H%M%S`

NumTotalThreads=$(($NumProcs * $NumThreads))

echo "---- worker.sh arguments -----"
echo InputDir=$InputDir
echo OutputDir=$OutputDir
echo ToolTimeLimit=$ToolTimeLimit
echo NumNodes=$NumNodes
echo NumProcs=$NumProcs
echo NumThreads=$NumThreads
echo Mode=$Mode
echo ScriptArgs=$ScriptArgs
echo NumTotalThreads=$NumTotalThreads
echo DoStaging=$DoStaging
echo LimitMaxMem=$LimitMaxMem
echo "--- worker.sh arguments end ---"

# Database path in $DataDir
Uniref90=$DataDir/uniref90/uniref90.fasta
Pdb70=$DataDir/pdb70
Mgnify=$DataDir/mgnify/mgy_clusters_2018_12.fa
SmallBfd=$DataDir/small_bfd/bfd-first_non_consensus_sequences.fasta
Bfd=$DataDir/bfd
Uniclust30=$DataDir/uniclust30/uniclust30_2018_08
BfdCs219Data=$Bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt_cs219.ffdata

if [[ $Mode = "bfd_deflate" ]]; then
    if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	# TODO: Move to share
	BfdA3M=/vol0006/ra000018/data/u10130/hhblits_deflate/bfd
	Uniclust30A3M=/vol0006/ra000018/data/u10130/hhblits_deflate/uniclust30
    else
	# TODO: Do staging for .ffindex
	BfdA3M=$DataDir/hhblits_deflate/bfd
	Uniclust30A3M=$DataDir/hhblits_deflate/uniclust30
    fi
else
    BfdA3M=""
    Uniclust30A3M=""
fi

# commands
Script=scripts/precompute_alignments_fugaku.py
LLIOTransfer=llio_transfer
DirTransfer=./scripts/my_dir_transfer

if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
    MpiExecEachNode="STUB"
    MpiArgs=""
    export LD_PRELOAD=/usr/lib/FJSVtcs/ple/lib64/libpmix.so # Required to use mpi4py
else
    MpiExecEachNode="mpiexec -np $NumNodes -npernode 1"
    MpiArgs="--bind-to none"
fi

# Temp locations
if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
    LOCALTMP=$PJM_LOCALTMP # Shared SSD, 87 GiB/node
elif [[ ${OPENFOLD_MACHINE} == "pbs" ]]; then
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
elif [[ $Mode = "bfd" || $Mode = "bfd_deflate" ]]; then
    Database=$Uniclust30
    DatabaseIsFile=0
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi

if [[ $DoStaging = 1 ]]; then
    if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	echo "Executing LLIO transfer for OpenFold"

	set -x
	time $DirTransfer \
	     $PREFIX/bin \
	     $PREFIX/lib
	time $LLIOTransfer \
	     $OPENFOLDDIR/*.py
	set +e # Allow errors as some files (such as scripts/setenv) might be already cached
	time $DirTransfer \
	     $OPENFOLDDIR/scripts \
	     $OPENFOLDDIR/openfold \
	     $OPENFOLDDIR/openfold.egg-info
	set -e
	set +x
    fi

    if [[ $Mode = "bfd" || $Mode = "bfd_deflate" ]]; then
	if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then

	    if [[ $Mode = "bfd" ]]; then
		echo "Mode must be bfd_deflate on Fugaku" >&2
		exit
	    fi

	    echo "Executing LLIO transfer for BFD files"
	    echo "Bfd=${Bfd}"
	    set +e
	    set -x
	    time $LLIOTransfer \
		 $Bfd/*_cs219*
	    set +x
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
	if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	    echo "Executing LLIO transfer for the DB file"
	    set +e
	    set -x
	    time $LLIOTransfer \
		 $Database
	    set +x
	    set -e

	else
	    echo "Copying the DB file to local temp"
	    DatabaseLocal=$LOCALTMP/`basename $Database`
	    time $MpiExecEachNode cp $Database $DatabaseLocal
	    Database=$DatabaseLocal
	fi

    else
	if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	    echo "Executing LLIO transfer for the DB directory"
	    if [[ $Mode = "bfd_deflate" ]]; then
		set +e
		set -x
		time $LLIOTransfer \
		     $Database/*_cs219.ffdata \
		     $Database/*_cs219.ffindex
		set +x
		set -e

	    else
		Exts=(ffdata ffindex)
		for Ext in ${Exts[@]}; do
		    set +e
		    set -x
		    time $LLIOTransfer \
			 $Database/*.$Ext
		    set +x
		    set -e
		done
	    fi

	else
	    echo "Copying the DB directory to local temp"
	    DatabaseLocal=$LOCALTMP/`basename $Database`
	    time $MpiExecEachNode mkdir -p $DatabaseLocal
	    if [[ $Mode = "bfd_deflate" ]]; then
		time $MpiExecEachNode cp -r $Database/*_cs219.ffdata $DatabaseLocal/
		time $MpiExecEachNode cp -r $Database/*_cs219.ffindex $DatabaseLocal/
	    else
		time $MpiExecEachNode cp -r $Database/*.ffdata $DatabaseLocal/
		time $MpiExecEachNode cp -r $Database/*.ffindex $DatabaseLocal/
	    fi
	    Database=$DatabaseLocal
	fi
    fi

    if [[ -n $BfdA3M ]]; then
	if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	    echo "Executing LLIO transfer for deflate BFD"
	    set +e
	    set -x
	    time $LLIOTransfer \
		 $BfdA3M/*.ffdata \
		 $BfdA3M/*.ffindex
	    set +x
	    set -e
	fi
    fi

    if [[ -n $Uniclust30A3M ]]; then
	if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	    echo "Executing LLIO transfer for deflate Uniclust30"
	    set +e
	    set -x
	    time $LLIOTransfer \
		 $Uniclust30A3M/*.ffdata \
		 $Uniclust30A3M/*.ffindex
	    set +x
	    set -e
	fi
    fi
fi

DatabaseArgs=""
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
	DatabaseArgs="$DatabaseArgs --bfd_deflate_a3m_hhm_path $BfdA3M/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    fi
    if [[ -n $Uniclust30A3M ]]; then
	DatabaseArgs="$DatabaseArgs --uniclust30_deflate_a3m_hhm_path $Uniclust30A3M/uniclust30_2018_08"
    fi
else
    echo "Invalid Mode: $Mode" >&2
    exit
fi

mkdir -p $OutputDir

if [[ $Mode = "bfd_deflate" ]]; then
    # TODO: Do not use hhblits in the home directory
    if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	export PATH=/vol0006/ra000018/data/u10130/repo/hh-suite/build/src:${PATH}
    else
	export PATH=/home/acd13359ld/repo/hh-suite/build/src:${PATH}
    fi
    echo hhblits=`which hhblits`
fi


while (( $NumProcs > 0 )); do

    MaxMem=""
    MaxMemArg=""
    if [[ $LimitMaxMem = 1 ]]; then
	MaxMem=$(($OPENFOLD_MAX_MEM / $NumProcs))
	MaxMemArg="--max_memory ${MaxMem}"
    fi

    export OMP_NUM_THREADS=$NumThreads # Define just in case it is used
    export PARALLEL=$OMP_NUM_THREADS

    NumTaskProcs=$(($NumProcs * $NumNodes))
    if [[ ${OPENFOLD_MACHINE} == "fugaku" ]]; then
	MpiExecTask="mpiexec -n $NumTaskProcs"
    else
	MpiExecTask="mpiexec -np $NumTaskProcs -npernode $NumProcs"
    fi

    echo ""
    echo "Starting the script: NumProcs=${NumProcs}, NumTaskProcs=${NumTaskProcs}, NumThreads=${NumThreads}, MaxMem=${MaxMem}"

    ReportOutPath=`mktemp`

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
	--report_out_path $ReportOutPath \
	$DatabaseArgs \
	$MaxMemArg \
	$ScriptArgs

    RemainingCount=`cat $ReportOutPath`
    rm ${ReportOutPath}

    echo "Remaining number of sequences: $RemainingCount"

    NumProcs=$(($NumProcs / 2))

    NumTaskProcs=$(($NumProcs * $NumNodes))
    if (( $NumTaskProcs > $RemainingCount )); then
	NumProcs=$(( ($RemainingCount + $NumNodes - 1) / $NumNodes ))
    fi

    if (( $NumProcs > 0 )); then
	NumThreads=$(($NumTotalThreads / $NumProcs))
    else
	NumThreads=0
    fi

done
