#!/bin/bash

LogDir=$1

if [ ! -d $LogDir ]; then
    echo "LogDir ($LogDir) does not exist"
    exit 1
fi

grep "inference_stat" ${LogDir}/output/0/* | grep "NG" | sed -e 's/.*inference_stat \([^ ]*\).*/\1/'

