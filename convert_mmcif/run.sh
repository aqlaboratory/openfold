#!/bin/bash -eux
#PJM -L "node=16"
#PJM -L "rscgrp=small"
#PJM -L "elapse=2:00:00"
#PJM --rsc-list "retention_state=0"
#PJM -j
#PJM --mpi max-proc-per-node=24
#PJM -S
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


source ../scripts/setenv

## set output directory where converted .pkl files stored
OutputDir=

## set input directory which contains .cif files
InputDir=$DATADIR/pdb_mmcif/mmcif_files

[ -d "$InputDir" ] && echo "InputDir($InputDir) is not exist" && exit 1
[ -d "$OutputDir" ] && echo "OutputDir($OutputDir) is not exist" && exit 1

MyName=run
Time=`date "+%y%m%d%H%M%S%3N"`
JobName="$MyName.$Time"

LogDir="log"/"$JobName"
mkdir -p $LogDir

ListFile=$LogDir/file_list
find $InputDir -name *.cif > $ListFile

NumProcs=$(( $PJM_PROC_BY_NODE * $PJM_NODE))

mkdir -p $OutputDir

$FJSVXTCLANGA/bin/mpiexec -np $NumProcs \
                          --mca orte_abort_print_stack 1 \
                          --of-proc $LogDir/output/%/1000r/out \
                          -x PYTHONPATH=$OPENFOLDDIR \
                          python3 convert.py $ListFile $OutputDir $NumProcs
