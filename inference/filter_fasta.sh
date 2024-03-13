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

#set -x

InputFasta=$1
OutputFasta=$2
AlignmentDir=$3
ProcessedDir=$4
NGList=$5

if [ ! -f $InputFasta ]; then
    echo "InputFasta ($InputFasta) does not exist"
    exit 1
fi

if [ -e $OutputFasta ]; then
    echo "OutputFasta ($OutputFasta) already exist"
    exit 2
fi

if [ ! -d $AlignmentDir ]; then
    echo "AlignmenDir ($AlignmentDir) does not exist"
    exit 3
fi

PredictionDir=${ProcessedDir}/predictions
if [ "$ProcessedDir" != "" ] && [ ! -d $PredictionDir ]; then
    echo "ProcessedDir/predictions ($PredictionDir) does not exist"
    exit 4
fi

if [ "$NGList" != "" ] && [ ! -f $NGList ]; then
    echo "NGList ($NGList) does not exist"
    exit 5
fi

function num_exist_output() {
    ListFile=$1
    Name=$2
    Result=$(grep -c "^${Name}_.*\.pdb$" $ListFile)
    echo $Result
}

NumNoAlignment=0
NumContainX=0
NumProcessed=0
NumNG=0
NumOutput=0

Names=()
Seqs=()

while true ; do
    read NameLine || break
    if [[ ${NameLine} =~ ^\>.*$ ]]; then
        Name=${NameLine:1}
        read Seq || break

        #echo "Name: $Name, Seq: $Seq"
        Names+=("$Name")
        Seqs+=("$Seq")
    fi
done < $InputFasta

if [ "${#Names[@]}" != "${#Seqs[@]}" ]; then
    echo "The number of Names ${#Names[@]} is not match with that of Seqs ${#Seqs[@]}"
    exit 1
fi

NumInput=${#Names[@]}

OutputFastaNG="${OutputFasta}.ng"

if [ "$ProcessedDir" != "" ] ; then
    ProcessedPDBFiles=$(mktemp)
    ls -1 $PredictionDir >> $ProcessedPDBFiles
fi

for ((i=0; i<$NumInput; i++)); do
    Name=${Names[$i]}
    Seq=${Seqs[$i]}

    NameFirst=${Name%%_*}
    NameLast=${Name##*_}
    NameLower=${NameFirst,,}_${NameLast}

    if [[ "$Seq" = "" ]] ; then
        echo -e ">${Name}\n${Seq}" >> $OutputFastaNG
        ((NumEmpty++))
        continue
    fi
    if [[ "$Seq" =~ "X" ]] ; then
        echo -e ">${Name}\n${Seq}" >> $OutputFastaNG
        ((NumContainX++))
        continue
    fi

    if [ "$NGList" != "" ]; then
        grepNG=`grep -E "^${NameLower}$" $NGList`
        if [ "$grepNG" != "" ]; then
            echo -e ">${Name}\n${Seq}" >> $OutputFastaNG
            ((NumNG++))
            continue
        fi
    fi

    if [ ! -d $AlignmentDir/$NameLower ] ||\
       [ ! -f $AlignmentDir/$NameLower/mgnify_hits.a3m ] ||\
       [ ! -f $AlignmentDir/$NameLower/pdb70_hits.hhr ] ||\
       [ ! -f $AlignmentDir/$NameLower/small_bfd_hits.sto ] ||\
       [ ! -f $AlignmentDir/$NameLower/uniref90_hits.a3m ]; then
        #echo "no alignment:" $NameLower
        echo -e ">${Name}\n${Seq}" >> $OutputFastaNG
        ((NumNoAlignment++))
        continue
    fi
    if [ "$ProcessedDir" != "" ]; then
        if [ `num_exist_output $ProcessedPDBFiles $NameLower` -eq 2 ]; then
            #echo "already exist output:" $NameLower
            echo -e ">${Name}\n${Seq}" >> $OutputFastaNG
            ((NumProcessed++))
            continue
        fi
    fi
    echo -e ">${Name}\n${Seq}" >> $OutputFasta
    ((NumOutput++))
done

if [ "$ProcessedDir" != "" ] ; then
    rm -f $ProcessedPDBFiles
fi

echo "NumInput:" $NumInput
echo "NumEmpty:" $NumEmpty
echo "NumContainX:" $NumContainX
echo "NumNoAlignment:" $NumNoAlignment
echo "NumProcessed:" $NumProcessed
echo "NumNG:" $NumNG
echo "NumOutput:" $NumOutput
