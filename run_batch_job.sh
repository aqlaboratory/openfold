#!/bin/bash

################
# Example CMD
./run_batch_job.sh \
    s3://sagemaker-us-east-2-032243382548/openfold_testing/T1084.fasta \
    /data/fasta_dir \
    "python3 /opt/openfold/run_pretrained_openfold.py \
    /data/fasta_dir \
    /database/pdb_mmcif/mmcif_files/ \
    --uniref90_database_path /database/uniref90/uniref90.fasta \
    --mgnify_database_path /database/mgnify/mgy_clusters_2018_12.fa \
    --pdb70_database_path /database/pdb70/pdb70 \
    --uniclust30_database_path /database/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --output_dir /data \
    --model_device cuda:0 \
    --jackhmmer_binary_path /opt/conda/bin/jackhmmer \
    --hhblits_binary_path /opt/conda/bin/hhblits \
    --hhsearch_binary_path /opt/conda/bin/hhsearch \
    --kalign_binary_path /opt/conda/bin/kalign \
    --jax_param_path /database/params/params_model_1.npz" \
    /data \
    s3://sagemaker-us-east-2-032243382548/openfold_testing/

input_source=$1
input_destination=$2
script=$3
output_source=$4
output_destination=$5

aws s3 cp $input_source $input_destination
$script
a2s s3 cp $output_source $output_destination