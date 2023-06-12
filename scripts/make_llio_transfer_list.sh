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

if [ ! -d $1 ]; then
    echo "Please specify a directory to 1st argument"
    exit 1
fi

if [ "$PREFIX" = "" -o "$OPENFOLDDIR" = "" ]; then
    echo "Please specify PREFIX and OPENFOLDDIR"
    exit 1
fi

egrep -v '= \-1 ENOENT|O_DIRECTORY' $1/strace.* | egrep O_RDONLY | cut -d\" -f 2  |egrep ^/vol.... | sort | uniq | grep -v -E "\.(fasta|cif|a3m|hhr|sto|npz)$" | sed -e 's|'"$PREFIX|"'$PREFIX|' | sed -e 's|'"$OPENFOLDDIR|"'$OPENFOLDDIR|' > $2

