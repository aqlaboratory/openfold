#!/bin/bash

if [ ! -d $1 ]; then
    echo "Please specify a directory to 1st argument"
    exit 1
fi

egrep -v '= \-1 ENOENT|O_DIRECTORY' $1/strace.* | egrep O_RDONLY | cut -d\" -f 2  |egrep ^/vol.... | sort | uniq | grep -v -E "\.(fasta|cif|a3m|hhr|sto|npz)$" > $2

