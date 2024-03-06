#!/bin/bash -ex
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

echo `date +%s` `date` -- worker.sh start

SIZE=${NumNodes}
RANK=${PMIX_RANK}
echo `date` "#start hostname: " `hostname` "JOBID: "${PJM_SUBJOBID}

### Copy opt.tgz file and cpdata via /var/crash/MLP
export SCRIPT_DIR=${ShareDir}

#-- OPT file(SRC_FILE) defined in BatchBase
#export OPT_PATH=${ShareDir}/TF220
#export OPT_PATH=${ShareDir}/TF220-33

#-- setenv
module list

#time -p tar -I pigz -xf ${SRC_FILE} -C ${ShareDir}
#cp ./${CPSCRIPT} ${ShareDir}/cpdata.sh


LogDir=${LOGDIR}

. "$ParameterFile"

export LD_PRELOAD=/usr/lib/FJSVtcs/ple/lib64/libpmix.so:$LD_PRELOAD

ulimit -s 16384
ulimit -c 0

if [ $RANK -eq "0" ]; then
    export DNNL_VERBOSE=$dnnlverbose
    env > ${LogDir}/rank0_env
    ulimit -a > ${LogDir}/rank0_ulimit
    #(vmstat -t 1 > ${LogDir}/vmstat.log) &
    #PID_VMSTAT=$!
fi

#         strace -ff -e trace=open,openat -o ${LogDir}/strace.${PMIX_RANK} 
time -p numactl --cpunodebind 4-7 --membind 4-7 \
    "${PARAMS[@]}"

if [ $RANK -eq "0" ]; then
    #kill -9 $PID_VMSTAT
    :
fi

unset LD_PRELOAD

echo `date +%s` `date` -- ALL DONE!
