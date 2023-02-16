#!/bin/bash

DB=pdb70
DataDir=output_202212_all

oks=`cd $DataDir && find . -type f | grep $DB | grep hhr$ | sort | uniq`
for x in ${oks[@]}; do
    hhr=`echo $x | sed "s/_ok$//g"`
    chain=`echo $hhr | sed "s/^\.\///g" | sed "s/\/.*//g"`
    count=`cat $DataDir/$hhr | grep "^No " | wc -l | sed "s/ .*//g"`
    echo $chain $count
done
