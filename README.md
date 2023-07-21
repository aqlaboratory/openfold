# Permutation code README

## Overview: before running training script

NB: before running the test codes,please download the procrustes package first:
from https://github.com/theochem/procrustes

Make sure that the product of running ```scripts/generate_mmcif_cache.py``` is ready and available in ```tests/test_data```
I have uploaded the json file to owncloud [here](https://oc.embl.de/index.php/s/wVUwc1IHiJUt9sP)

To test the train multimer codes:
```bash
python3 train_openfold.py /g/alphafold/AlphaFold_DBs/pdb_mmcif/mmcif_files/ \
    tests/test_data/alignments/ \
    /g/alphafold/AlphaFold_DBs/pdb_mmcif/mmcif_files/ \
    /scratch/gyu/train_openfold_output \
    2500-01-01 \
    --train_mmcif_data_cache_path=/tests/test_data/train_mmcifs_cache.json \
    --template_release_dates_cache_path=tests/test_data/mmcif_cache.json \
    --config_preset=model_1_multimer_v3 --seed=42 --gpus=1
```

Unlike training monomer, chain_cache_data is not required but the train_mmcifs_cache is required. In this case, I selected these 9 mmcifs that are already in the previous test_data folder as a training set. ```./tests/test_data/train_mmcifs_cache.json``` in the command above record the information of these 9 structures and is needed to run the training code. 

## Issues
Testing the codes on cpu works fine but when running it on a gpu, it causes ```RuntimeError: CUDA error: device-side assert triggered``` at unexpected steps.

For example, this error was raised while calculating the best rotation matrix that aligns selected anchors during multi-chain permutation steps, I have to use 
```torch.masked_select``` and ```torch.index_select``` in https://github.com/dingquanyu/openfold/blob/a1ef4c8fa99da5cff9501051de71be440ca3cedf/openfold/utils/loss.py#L2043 and https://github.com/dingquanyu/openfold/blob/a1ef4c8fa99da5cff9501051de71be440ca3cedf/openfold/utils/loss.py#L2060 instead of simply slicing the matrix like ```matrix[index]```.

Later on the same ```CUDA error: device-side assert triggered``` error was raised while adding dimensions to the ```atom_pred_positions``` in https://github.com/dingquanyu/openfold/blob/a1ef4c8fa99da5cff9501051de71be440ca3cedf/openfold/utils/loss.py#L989

I've dumped the matrices in a pickle and load them individually outside the programme to a GPU then the indexing steps worked without the CUDA error.
