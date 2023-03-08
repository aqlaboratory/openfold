#!/bin/bash -eux
#PJM -L "node=16"
#PJM -L "rscgrp=small"
#PJM -L "elapse=2:00:00"
#PJM --rsc-list "retention_state=0"
#PJM -j
#PJM --mpi max-proc-per-node=24
#PJM -S


## set input directory which contains .cif files
InputDir=
## set output directory where converted .pkl files stored
OutputDir=



[ -z "$InputDir" ] && echo "InputDir is not defined" && exit 1
[ -z "$OutputDir" ] && echo "OutputDir is not defined" && exit 1

source ../scripts/setenv

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
