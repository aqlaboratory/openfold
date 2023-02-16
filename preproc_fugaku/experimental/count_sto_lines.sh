#!/bin/bash

DB=small_bfd
DataDir=output_202212_all

oks=`cd $DataDir && find . -type f | grep $DB | grep sto$ | sort | uniq`
for x in ${oks[@]}; do
    sto=`echo $x | sed "s/_ok$//g"`
    chain=`echo $sto | sed "s/^\.\///g" | sed "s/\/.*//g"`
    count=`cat $DataDir/$sto | grep "^#=GS " | wc -l | sed "s/ .*//g"`
    echo $chain $count
done
