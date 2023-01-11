#!/bin/bash

InputFasta=$1
OutputFasta=$2
AlignmentDir=$3
PredictionDir=$4
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

if [ "$PredictionDir" != "" ] && [ ! -d $PredictionDir ]; then
    echo "PredictionDir ($PredictionDir) does not exist"
    exit 4
fi

if [ "$NGList" != "" ] && [ ! -f $NGList ]; then
    echo "NGList ($NGList) does not exist"
    exit 5
fi

function num_exist_output() {
    Dir=$1
    Name=$2
    #1awi_A_model_1_unrelaxed.pdb
    Result=(`find $Dir -maxdepth 1 -name $Name_*.pdb`)
    echo ${#Result[@]}
}

NumNoAlignment=0
NumContainX=0
NumProcessed=0
NumNG=0
NumOutput=0

NameLines=(`grep -E "^>" $InputFasta`)
Seqs=(`grep -E "^[^>]" $InputFasta`)

NumInput=${#NameLines[@]}

for ((i=0; i<$NumInput; i++)); do
    NameLine=${NameLines[$i]}
    Name=${NameLine:1}
    Seq=${Seqs[$i]}

    NameFirst=${Name%%_*}
    NameLast=${Name##*_}
    NameLower=${NameFirst,,}_${NameLast}

    if [[ "$Seq" =~ "X" ]] ; then
	((NumContainX++))
	continue
    fi

    grepNG=`grep -E "^${NameLower}$" $NGList`
    if [ "$grepNG" != "" ]; then
	((NumNG++))
	continue
    fi

    if [ ! -d $AlignmentDir/$NameLower ] ||\
       [ ! -f $AlignmentDir/$NameLower/mgnify_hits.a3m_ok ] ||\
       [ ! -f $AlignmentDir/$NameLower/pdb70_hits.hhr_ok ] ||\
       [ ! -f $AlignmentDir/$NameLower/small_bfd_hits.sto_ok ] ||\
       [ ! -f $AlignmentDir/$NameLower/uniref90_hits.a3m_ok ]; then
	#echo "no alignment:" $NameLower
	((NumNoAlignment++))
	continue
    fi
    if [ "$PredictionDir" != "" ]; then
	if [ `num_exist_output $PredictionDir $NameLower` -eq 2 ]; then
	    #echo "already exist output:" $NameLower
	    ((NumProcessed++))
	    continue
	fi
    fi
    grep -E "^>$Name$" -A 1 $InputFasta >> $OutputFasta
    ((NumOutput++))
done

echo "NumInput:" $NumInput
echo "NumContainX:" $NumContainX
echo "NumNoAlignment:" $NumNoAlignment
echo "NumProcessed:" $NumProcessed
echo "NumNG:" $NumNG
echo "NumOutput:" $NumOutput
