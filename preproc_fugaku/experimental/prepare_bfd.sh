#!/bin/bash

set -eu

Bfd=$1
BfdLocal=$2
BfdCs219Data=$3

echo "prepare_bfd: Copying the DB file to local temp"
mkdir $BfdLocal
time cp $BfdCs219Data $BfdLocal
time cp $Bfd/*.ffindex $BfdLocal

echo "prepare_bfd: Adding links of the DB files to local temp"
for x in `ls $Bfd`; do
    if [[ ! -e $BfdLocal/$x ]]; then
	ln -s `realpath $Bfd/$x` $BfdLocal/$x
    fi
done
