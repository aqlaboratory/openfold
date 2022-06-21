#!/bin/bash

# the alphafold script has a few frustrating behaviors in how it parses input flags, this script reduces the need for
# code duplication in the argo workflow
set -x

database_path=$1
fasta_path=$2

base_arguments=(
    ${fasta_path}
    "${database_path}/pdb_mmcif/mmcif_files/"
    "--uniref90_database_path" "${database_path}/uniref90/uniref90.fasta"
    "--mgnify_database_path" "${database_path}/mgnify/mgy_clusters_2018_12.fa"
    "--pdb70_database_path" "${database_path}/pdb70/pdb70"
    "--uniclust30_database_path" "${database_path}/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
    "--output_dir" "out"
    "--bfd_database_path" "${database_path}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    "--model_device" "cuda:0"
    "--jackhmmer_binary_path" "/home/samdeluca/miniconda3/envs/openfold_venv/bin/jackhmmer"
    "--hhblits_binary_path" "/home/samdeluca/miniconda3/envs/openfold_venv/bin/hhblits"
    "--hhsearch_binary_path" "/home/samdeluca/miniconda3/envs/openfold_venv/bin/hhsearch"
    "--kalign_binary_path" "/home/samdeluca/miniconda3/envs/openfold_venv/bin/kalign"
    "--openfold_checkpoint_path" "/mnt/openfold_params/101-80999.ckpt,/mnt/openfold_params/116-84749.ckpt,/mnt/openfold_params/94-79249.ckpt"
)

python3 run_pretrained_openfold.py  ${base_arguments[*]}
