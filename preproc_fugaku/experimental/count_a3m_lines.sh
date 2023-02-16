#!/bin/bash

DB=uniref90
DataDir=output_202212_all

oks=`cd $DataDir && find . -type f | grep $DB | grep a3m$ | sort | uniq`
for x in ${oks[@]}; do
    a3m=`echo $x | sed "s/_ok$//g"`
    chain=`echo $a3m | sed "s/^\.\///g" | sed "s/\/.*//g"`
    count=`wc -l $DataDir/$a3m | sed "s/ .*//g"`
    count=$(( $count / 2 )) # because it's FASTA
    echo $chain $count
done
