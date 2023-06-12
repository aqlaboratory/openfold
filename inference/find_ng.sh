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

LogDir=$1

if [ ! -d $LogDir ]; then
    echo "LogDir ($LogDir) does not exist"
    exit 1
fi

grep "inference_stat" ${LogDir}/output/0/* | grep "NG" | sed -e 's/.*inference_stat \([^ ]*\).*/\1/'

