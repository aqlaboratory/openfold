#!/bin/bash

# Generates uniclust30 all-against-all alignments on a SLURM cluster.
# Thanks to Milot Mirdita for help & feedback on this script.

set -e

if [[ $# != 3 ]]; then
    echo "usage: ./run_uniclust30_search.sh <uniclust30_path> <scratch_dir> <out_dir>"
    exit
fi

UNICLUST_PATH=$1
SCRATCH_DIR_BN=$2
OUT_DIR=$3

CPUS_PER_TASK=4
MAX_SIZE=10000000000 # 10GB

SCRATCH_DIR="${SCRATCH_DIR_BN}_${SLURM_NODEID}"

mkdir -p ${SCRATCH_DIR}
mkdir -p ${OUT_DIR}

# copy database to local ssd
DB_BN=$(basename $UNICLUST_PATH)
DB_DIR="/dev/shm/uniclust30"
mkdir -p $DB_DIR
cp ${UNICLUST_PATH}*.ff* $DB_DIR
DB="${DB_DIR}/${DB_BN}"

for f in $(ls $OUT_DIR/*.zip)
do 
    zipinfo -1 $f '*/' | awk -F/ '{print $(NF-1)}' >> ${DB_DIR}/already_searched.txt
done

python3 filter_ffindex.py ${DB}_a3m.ffindex ${DB_DIR}/already_searched.txt ${DB_DIR}/filtered_a3m.ffindex 

TARGET="${DB}_a3m_${SLURM_NODEID}.ffindex"
split -n "l/$((SLURM_NODEID + 1))/${SLURM_JOB_NUM_NODES}" "${DB_DIR}/filtered_a3m.ffindex" > $TARGET

open_sem() {
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for ((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock() {
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
        ( "$@"; )
        # push the return code of the command to the semaphore
        printf '%.3d' $? >&3
    )&
}

task() {
    dd if="${DB}_a3m.ffdata" ibs=1 skip="${OFF}" count="${LEN}" status=none | \
	hhblits -i stdin \
            -oa3m "${SCRATCH_DIR}/${KEY}/uniclust30.a3m" \
            -v 0 \
            -o /dev/null \
            -cpu $CPUS_PER_TASK \
            -d $DB \
            -n 3 \
            -e 0.001
}

zip_or_not() {
    SIZE=$(du -hbs $SCRATCH_DIR | sed 's/|/ /' | awk '{print $1}')
    #if [[ "$SIZE" -gt "$MAX_SIZE" ]]
    if [[ "2" -gt "1" ]]
    then
        wait
        RANDOM_NAME=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32)
        zip -r "${OUT_DIR}/${RANDOM_NAME}.zip" $SCRATCH_DIR
        find $SCRATCH_DIR -mindepth 1 -type d -exec rm -rf {} +
    fi
}

N=$(($(nproc) / ${CPUS_PER_TASK}))
open_sem $N
while read -r KEY OFF LEN; do
    PROT_DIR="${SCRATCH_DIR}/${KEY}"
    
    if [[ -d $PROT_DIR ]]
    then
        continue
    fi
    
    mkdir -p $PROT_DIR
    run_with_lock task "${KEY}" "${OFF}" "${LEN}"
    zip_or_not
done < $TARGET

wait

zip_or_not

wait
