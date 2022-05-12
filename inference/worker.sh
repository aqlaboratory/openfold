#!/bin/bash -ex

echo `date +%s` `date` -- worker.sh start

SIZE=${NumNodes}
RANK=${PMIX_RANK}
echo `date` "#start hostname: " `hostname` "JOBID: "${PJM_SUBJOBID}

### Copy opt.tgz file and cpdata via /var/crash/MLP
export SCRIPT_DIR=${ShareDir}

export CPSCRIPT="cpdata_decomp8K.sh"

#-- OPT file(SRC_FILE) defined in BatchBase
#export OPT_PATH=${ShareDir}/TF220
#export OPT_PATH=${ShareDir}/TF220-33

#-- setenv
module list

#time -p tar -I pigz -xf ${SRC_FILE} -C ${ShareDir}
#cp ./${CPSCRIPT} ${ShareDir}/cpdata.sh


LogDir=${LOGDIR}

. "$ParameterFile"

ulimit -s 16384

if [ $RANK -eq "0" ]; then
    export DNNL_VERBOSE=$dnnlverbose
    env > ${LogDir}/rank0_env
    ulimit -a > ${LogDir}/rank0_ulimit
    #(vmstat -t 1 > ${LogDir}/vmstat.log) &
    #PID_VMSTAT=$!
fi

time -p numactl --cpunodebind 4-7 --membind 4-7 \
    "${PARAMS[@]}"

if [ $RANK -eq "0" ]; then
    #kill -9 $PID_VMSTAT
    :
fi

unset LD_PRELOAD

echo `date +%s` `date` -- ALL DONE!
