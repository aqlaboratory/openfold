#!/bin/bash
#
# Copyright 2021 AlQuraishi Laboratories
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
#
# Downloads .cif files matching the RODA alignments. Outputs a list of 
# RODA alignments for which .cif files could not be found..
if [[ $# != 2 ]]; then
    echo "usage: ./download_roda_pdbs.sh <out_dir> <roda_pdb_alignment_dir>"
    exit 1
fi

OUT_DIR=$1
RODA_ALIGNMENT_DIR=$2

if [[ -d $OUT_DIR ]]; then
    echo "${OUT_DIR} already exists. Download failed..."
    exit 1
fi

SERVER=snapshotrsync.rcsb.org                       # RCSB server name
PORT=873                                           # port RCSB server is using

rsync -rlpt -v -z --delete --port=$PORT $SERVER::20220103/pub/pdb/data/structures/divided/mmCIF/ $OUT_DIR 2>&1 > /dev/null

for f in $(find $OUT_DIR -mindepth 2 -type f); do
    mv $f $OUT_DIR
    BASENAME=$(basename $f)
    gunzip "${OUT_DIR}/${BASENAME}"
done

find $OUT_DIR -mindepth 1 -type d,l -delete

for d in $(find $RODA_ALIGNMENT_DIR -mindepth 1 -maxdepth 1 -type d); do
    BASENAME=$(basename $d)
    PDB_ID=$(echo $BASENAME | cut -d '_' -f 1)
    CIF_PATH="${OUT_DIR}/${PDB_ID}.cif"
    if [[ ! -f $CIF_PATH ]]; then
        echo $d
    fi
done
