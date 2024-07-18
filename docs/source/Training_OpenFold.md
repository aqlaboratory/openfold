# Training OpenFold
## Background

This guide covers how to train an OpenFold model for monomers. Some additional instructions are provided at the end for fine-tuning your model.

### Pre-requisites: 

This guide requires the following:
- [Installation of OpenFold and dependencies](Installation.md) (Including jackhmmer and hhblits depedencies)
- A preprocessed dataset:
	- For this guide, we will use the original OpenFold dataset which is available on RODA, processed with [these instructions](OpenFold_Training_Setup.md).
- GPUs configured with CUDA. Training OpenFold with CPUs only is not supported. 

## Training a new OpenFold model 

#### Basic command

For a dataset that has the default alignment file structure, e.g. 

```
-$DATA_DIR
  └── pdb_data 
 	   ├──  mmcifs 
 	   	    ├── 3lrm.cif
 	   	    └── 6kwc.cif
			...
	   ├── obsolete.dat 	
	   ├── duplicate_pdb_chains.txt 
	   └── data_caches 
	   		├── duplicate_pdb_chains.txt 
	   		└── data_caches 
  └── alignment_data
	  └── alignments
   	  		├── 3lrm_A/ 
   	  		├── 3lrm_B/ 
	  	    └── 6kwc_A/
	  		...
```

The basic command to train a new OpenFold model is: 

```
python3 train_openfold.py $DATA_DIR/pdb/mmcifs $DATA_DIR/alignment_data/alignments $TEMPLATE_MMCIF_DIR $OUTPUT_DIR \
    --max_template_date 2021-10-10 \ 
    --train_chain_data_cache_path $DATA_DIR/pdb_data/data_caches/chain_data_cache.json \
    --template_release_dates_cache_path $DATA_DIR/pdb_data/data_caches/mmcif_cache.json \ 
	--config_preset initial_training \
    --seed 42 \
    --obsolete_pdbs_file_path $DATA_DIR/pdb_data/obsolete.dat \
    --num_nodes 1 \
    --gpus 4 \
    --num_workers 4
```

The required arguments are:
- `mmcif_dir` : Mmcif files for the training set.
- `alignments_dir`: Alignments for the sequences in `mmcif_dir`, see expected directory structure 
- `template_mmcif_dir`:  Template mmcif files with structures, which can be the same directory as mmcif_dir. The `max_template_date` and `template_release_dates_cache_path` will specify which templates will be allowed based on a date cutoff
- `output_dir` : Where model checkpoint files and other outputs will be saved. 

Commonly used flags include:
- `config_preset`: Specifies which selection of hyperparameters should be used for initial model training. Commonly used configs are defined in [`openfold/config.py`](https://github.com/aqlaboratory/openfold)
- `num_nodes` and `gpus`:  Specifies number of nodes and GPUs available to train OpenFold.
- `seed` - Specifies random seed
- `num_workers`: Number of CPU workers to assign for creating dataset examples
- `obsolete_pdbs_file_path`: Specifies obsolete pdb IDs that should be excluded from training.
- `val_data_dir` and `val_alignment_dir`: Specifies data directory and alignments for validation dataset. 

```{note}
Note that `--seed` must be specified to correctly configure training examples on multi-GPU training runs
```

#### Train with OpenFold Dataset Configuration

If the [OpenFold alignment database](OpenFold_Training_Setup.md#2-creating-alignment-dbs-optional) setup is used, resulting in a data directory such as:
```
- $DATA_DIR 
  ├── duplicate_pdb_chains.txt
  ├── pdb_data
	  └── mmcifs 
		  ├── 3lrm.cif
		  └── 6kwc.cif
  └── alignment_data 
	  └── alignment_db
		  ├── alignment_db_0.db 
          ├── alignment_db_1.db
          ...
          ├── alignment_db_9.db
		  └── alignment_db.index 
```

The training command will use the `alignment_index_path` argument to specify `db.index` files, e.g.: 

```
python3 train_openfold.py $DATA_DIR/pdb_data/mmcifs $DATA_DIR/alignment_data/alignment_db $TEMPLATE_MMCIF_DIR $OUTPUT_DIR \
    --max_template_date 2021-10-10 \ 
    --train_chain_data_cache_path $DATA_DIR/pdb_data/data_caches/chain_data_cache.json \
    --template_release_dates_cache_path $DATA_DIR/pdb_data/data_caches/mmcif_cache.json \ 
	--alignment_index_path $DATA_DIR/pdb/alignment_db.index 
	--config_preset initial_training \
    --seed 42 \
    --obsolete_pdbs_file_path $DATA_DIR/pdb/obsolete.dat \
    --num_nodes 1 \
    --gpus 4 \
    --num_workers 4
```

#### Additional command line flag options:

Here we provide brief descriptions for customizing your training run of OpenFold. A full description of all flags can be accessed by using the `--help` option in the script 

- **Use Deepspeed acceleration strategy:** `--deepspeed_config` This option configures OpenFold to use custom Deepspeed kernels. This option requires a deepspeed_config.json, you can create your own, or use the one in the OpenFold directory 

- **Use a validation dataset:** Specify validation database paths with `--val_data_dir` + `--val_alignment_dir`. Validation metrics will be evaluated on these datasets.

- **Use a self-distillation dataset:**  Specify paths with `--distillation_data_dir` and `--distillation_alignment_dir` flags

- **Change specific parameters in the model or data setup:**  `--experiment_config_json`. These parameters must be defined in the [`openfold/config.py`](https://github.com/aqlaboratory/openfold/blob/main/openfold/config.py). For example to change the crop size for training a model, you can write the following json:
	```cropsize.json
	{
			"data.train.crop_size": 128
	}
	```

- **Configure training settings with PyTorch Lightning** 
	
	Some flags e.g. `--precision`, `--max_epochs` configure training behavior. See the Pytorch Lightning Trainer args section in the `--help`  menu for more information and consult [Pytorch lightning documentation](https://lightning.ai/docs/pytorch/stable/)
	
	- Precision: On A100s, OpenFold training works best with bfloat 16 precision (e.g. `--precision bf16-mixed`) 
	
- **Restart training from an existing checkpoint:** Use the `--resume_from_ckpt` to restart training from an existing checkpoint.

## Advanced Training Configurations 
:::

### Fine tuning from existing model weights 

If you have existing model weights, you can fine tune the model by specifying a checkpoint path with `--resume_from_ckpt` and `--resume_model_weights_only` arguments, e.g. 

```
python3 train_openfold.py $DATA_DIR/mmcifs $DATA_DIR/alignment.db $TEMPLATE_MMCIF_DIR $OUTPUT_DIR \
    --max_template_date 2021-10-10 \ 
    --train_chain_data_cache_path chain_data_cache.json \
    --template_release_dates_cache_path mmcif_cache.json \ 
	--config_preset finetuning \
	--alignment_index_path $DATA_DIR/pdb/alignment_db.index \ 
    --seed 4242022 \
    --obsolete_pdbs_file_path obsolete.dat \
    --num_nodes 1 \
    --gpus 4 \
    --num_workers 4 \
	--resume_from_ckpt $CHECKPOINT_PATH \
	--resume_model_weights_only
```

If you have model parameters from OpenFold v1.x, you may need to convert your checkpoint file or parameter. See [Converting OpenFold v1 Weights](convert_of_v1_weights.md) for more details. 

### Using MPI

If MPI is configured on your system, and you would like to use MPI to train OpenFold models, you may do so with the following step:

 1. Add the `mpi4py` package, which are available through pip and conda. Please see [mpi4py documentation](https://pypi.org/project/mpi4py/) for more instructions on installation.
2. Add the `--mpi_plugin` flag to your training command.


### Training Multimer  models

```{note}
Coming soon.
```