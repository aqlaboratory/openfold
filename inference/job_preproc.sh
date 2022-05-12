#!/bin/bash
#PJM -L "rscunit=rscunit_ft01,rscgrp=small"
#PJM -L elapse=2:0:0
#PJM -L proc-core=48
#PJM -L "node=1"
##:strict,freq=2200
#PJM --rsc-list "retention_state=0"
#PJM -j

source /home/ra010011/a04350/data/OpenFold/setenv

## without bfd
#./preproc.sh input alignment

## with bfd
./preproc_bfd.sh input alignment_bfd
