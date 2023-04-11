#!/bin/bash
# cd openfold
# python3 setup.py install
python3 train_openfold.py /bgfs02/PSIVANT/External/Openfold/Common/Datasets/pdb_mmcif/mmcif_files /bgfs02/PSIVANT/External/Openfold/Common/Datasets/openfold_retrain/OpenProteinSet/roda_pdb_flatten/alignments /bgfs02/PSIVANT/External/Openfold/Common/Datasets/pdb_mmcif/mmcif_files /bgfs02/PSIVANT/External/Openfold/Common/Datasets/train_output \
	2021-10-01 \
	--template_release_dates_cache_path /bgfs02/PSIVANT/External/Openfold/Users/milen.ferev/mmcif_cache.json \
        --seed 6102022 \
        --gpus 4 \
        --replace_sampler_ddp=True \
        --precision bf16 \
        --train_epoch_len 33000 \
        --train_chain_data_cache_path /bgfs02/PSIVANT/External/Openfold/Users/milen.ferev/chain_data_cache.json \
        --deepspeed_config /bgfs02/PSIVANT/External/Openfold/Users/milen.ferev/openfold/deepspeed_config.json \
        --obsolete_pdbs_file_path /bgfs02/PSIVANT/External/Openfold/Common/Datasets/pdb_mmcif/obsolete.dat \
        --checkpoint_every_epoch
