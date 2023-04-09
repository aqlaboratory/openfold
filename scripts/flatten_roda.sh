#!/usr/bin/env sh
#
# Flattens a downloaded RODA database into the format expected by OpenFold
# Args:
#     roda_dir: 
#           The path to the database you want to flatten. E.g. "roda/pdb" 
#           or "roda/uniclust30". Note that, to save space, this script
#           will empty this directory.
#     output_dir:
#           The directory in which to construct the reformatted data

if [[ $# != 2 ]]; then
    echo "usage: ./flatten_roda.sh <roda_dir> <output_dir>"
    exit 1
fi

RODA_DIR=$1
OUTPUT_DIR=$2

DATA_DIR="${OUTPUT_DIR}/data"
ALIGNMENT_DIR="${OUTPUT_DIR}/alignments"

mkdir -p "${DATA_DIR}"
mkdir -p "${ALIGNMENT_DIR}"

for chain_dir in $(ls "${RODA_DIR}"); do
    CHAIN_DIR_PATH="${RODA_DIR}/${chain_dir}"
    for subdir in $(ls "${CHAIN_DIR_PATH}"); do
        if [[ $subdir = "pdb" ]] || [[ $subdir = "cif" ]]; then
            mv "${CHAIN_DIR_PATH}/${subdir}"/* "${DATA_DIR}"
        else
            CHAIN_ALIGNMENT_DIR="${ALIGNMENT_DIR}/${chain_dir}"
            mkdir -p "${CHAIN_ALIGNMENT_DIR}"
            mv "${CHAIN_DIR_PATH}/${subdir}"/* "${CHAIN_ALIGNMENT_DIR}"
        fi
    done
done

NO_DATA_FILES=$(find "${DATA_DIR}" -type f | wc -l)
if [[ $NO_DATA_FILES = 0 ]]; then
    rm -rf ${DATA_DIR}
fi
