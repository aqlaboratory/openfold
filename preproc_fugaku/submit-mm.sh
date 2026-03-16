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

set -eu

#---------- I/O configurations ----------#

# The path to the directory containing the fasta files
InputDir=input_example_for_mm

# The path to the output directory
OutputDir=output

# The root path of temporary directories.
# Each job creates one or more sub-directories in $TempDir.
# It is safe to remove all temporary directories once all jobs are completed.
# /tmp is NOT recommended on Fugaku as OpenFold consumes GBs of disk space.
TempDir=temp

# Whether or not LLIO transfer is performed for Python modules, executables, and databases
DoStaging=1

# Convert small BFD MSAs from the sto format to a3m format to improve I/O performance in training/inference
ConvertSmallBFDToA3M=0

# Create output directories just before alignment file is created
CreateDirOnDemand=1

# If SubDirectorySize > 0, create subdirectories for each number of sequences specified by this
SubDirectorySize=0

#---------- Job configurations ----------#

# Whether or not each job is submitted
Do_small_bfd=1
Do_mgnify=1
Do_uniref90=1
Do_pdb_seqres=1
Do_uniprot=1

# The number of nodes of each search job
NumNodes=10

# The initial number of processes/node and threads/process for each search tool
NumProcsJackhmmer=1
NumThreadsJackhmmer=48
NumProcsHMMsearch=1
NumThreadsHMMsearch=48

# Limit max memory
LimitMaxMemJackhmmer=1
LimitMaxMemHMMsearch=0

# Time limit for each job in the pjsub time format
JobTime_small_bfd=2:00:00
JobTime_mgnify=2:00:00
JobTime_uniref90=1:00:00
JobTime_pdb_seqres=1:00:00
JobTime_uniprot=1:00:00

# Time limit for each tool in seconds
Timeout_small_bfd=3600
Timeout_mgnify=7200
Timeout_uniref90=7200
Timeout_pdb_seqres=3600
Timeout_uniprot=3600

# Use the stream version of sto-to-a3m conversion if sto file size is
# larger than this size, which keeps memory usage small
StreamSTOSize=1073741824 # 1 GiB

#------- Output size configurations -------#

# The maximum number of MSA hits ("": unlimited)
MaxHits_small_bfd=""
MaxHits_mgnify=5000
MaxHits_uniref90=10000
MaxHits_uniprot=50000

#----------- Configurations end -----------#

if [[ -n $MaxHits_small_bfd && $ConvertSmallBFDToA3M != 1 ]]; then
    echo "Warning: MaxHits_small_bfd is not effective because this is applied only if the output format is a3m (ConvertSmallBFDToA3M=1)." >&2
fi

SubmitScript=scripts/Submit_preproc_fugaku

# Submit Jackhmmer jobs
StepNameBase=openfold_preproc_`date +%Y%m%d%H%M%S`
JackhmmerDatabases=(small_bfd mgnify uniprot uniref90)
for Database in ${JackhmmerDatabases[@]}; do
    JobTimeVar="JobTime_${Database}"
    TimeoutVar="Timeout_${Database}"
    MaxHitsVar="MaxHits_${Database}"
    DoVar="Do_${Database}"
    if (( $DoVar == 1 )); then
        TimeLimit=${!JobTimeVar} \
                 ToolTimeLimit=${!TimeoutVar} \
                 NumNodes=$NumNodes \
                 NumProcs=$NumProcsJackhmmer \
                 NumThreads=$NumThreadsJackhmmer \
                 InputFile=$InputDir \
                 OutputDir=$OutputDir \
                 TempDir=$TempDir \
                 Mode=$Database \
                 StepName="${StepNameBase}_${Database}" \
                 DoStaging=$DoStaging \
                 LimitMaxMem=$LimitMaxMemJackhmmer \
                 StreamSTOSize=$StreamSTOSize \
                 ConvertSmallBFDToA3M=$ConvertSmallBFDToA3M \
                 CreateDirOnDemand=$CreateDirOnDemand \
                 SubDirectorySize=$SubDirectorySize \
                 MaxHits=${!MaxHitsVar} \
                 ScriptArgs="" \
                 $SubmitScript
    fi
done

# Submit a HMMsearch job
if (( $Do_pdb_seqres  == 1 )); then
    if (( $Do_uniref90  == 1 )); then
        # Follow the corresponding uniref90 job
        StepName="${StepNameBase}_uniref90"
    else
        # Run independently because the corresponding job is not submitted
        StepName="${StepNameBase}_pdb_seqres"
    fi

    TimeLimit=$JobTime_pdb_seqres \
             ToolTimeLimit=$Timeout_pdb_seqres \
             NumNodes=$NumNodes \
             NumProcs=$NumProcsHMMsearch \
             NumThreads=$NumThreadsHMMsearch \
                   InputFile=$InputDir \
             OutputDir=$OutputDir \
             TempDir=$TempDir \
             Mode=pdb_seqres \
             StepName=$StepName \
             DoStaging=$DoStaging \
             LimitMaxMem=$LimitMaxMemHMMsearch \
             StreamSTOSize=$StreamSTOSize \
             ConvertSmallBFDToA3M=$ConvertSmallBFDToA3M \
             CreateDirOnDemand=$CreateDirOnDemand \
             SubDirectorySize=$SubDirectorySize \
             MaxHits="" \
             ScriptArgs="" \
             $SubmitScript
fi
