# Training OpenFold
## Background

This guide covers how to train an OpenFold model. These instructions focus on training a model for predicting monomers, but additional instructions are provided for training a monomer / multimer model. 

### Pre-requisites: 

This guide requires the following:
- [Installation of OpenFold and dependencies](Installation.md) (Including jackhmmer and hhblits depedencies)
- A preprocessed dataset:
	- For this guide, we will use the original OpenFold dataset which is available on RODA (TODO: add link to processed dataset).
	- If you wish to construct your own dataset, [these instructions](OpenFold_Training_Setup.md) provide guidance for preprocessing alignments into an OpenFold format. 
- GPUs configured with CUDA. Training OpenFold with CPUs only is not supported. 

Expected data directory structure:
```
- OpenProteinSet 
  └── alignments 
	  └── 2x7l_M
		  └── mgnify_hits.a3m
		  └── bfd_uniclust_hits.a3m
		  └── uniref90_hits.a3m
		  └── pdb70_hits.hhr 
		...
  └── mmcifs 
	  └── 3u8d.cif
	  └── 3lrm.cif
	  ... 
  └── mmcif_cache.json 
  └── chain_data_cache.json 
```

The `mmcif_cache.json` and the `chain_data_cache.json` provide metadata for the mmcif and the protein chains in the dataset.

## Training a new OpenFold model 

#### Basic command
The basic command to train a new OpenFold model is 
```
python3 train_openfold.py $DATA_DIR/mmcifs/ $DATA_DIR/alignments/ template_mmcif_dir/ $OUTPUT_DIR \
    --max_template_date 2021-10-10 \ 
    --train_chain_data_cache_path chain_data_cache.json \
    --template_release_dates_cache_path mmcif_cache.json \ 
	--config_preset initial_training \
    --seed 42 \
    --obsolete_pdbs_file_path obsolete.dat \
    --num_nodes 1 \
    --gpus 4 \
    --num_workers 4 \
```

The required arguments are:
- `mmcif_dir` : Mmcif files for the training set.
- `alignment_dir`: Alignments for the sequences in `mmcif_dir`, see expected directory structure 
- `template_mmcif_dir`:  Template mmcif files with structures, which can be the same directory as mmcif_dir. The `max_template_date` and `template_release_dates_cache_path` will specify which templates will be allowed based on a date cutoff
- `$OUTPUT_DIR` : Where model checkpoint files and other outputs will be saved. 

Commonly used flags include:
- `config_preset`: Specifies which selection of hyperparameters should be used for initial model training. Commonly used configs are defined in `openfold/config.py` 
- `num_nodes` and `gpus`:  Specifies number of nodes and GPUs available to train OpenFold.
- `seed` - Specifies random seed
- `num_workers`: Number of CPU workers to assign for creating dataset examples
- `obsolete_pdbs_file_path`: Specifies obsolete pdb IDs that should be excluded from training.
- `val_data_dir` and `val_alignment_dir`: Specifies data directory and alignments for validation dataset. 

```{note}
Note that `--seed` must be specified to correctly configure training examples on multi-GPU training runs
```



#### Train OpenFold with Different Dataset Configurations

If the [OpenFold alignment database](OpenFold_Training_Setup.md#2-creating-alignment-dbs-optional) setup is used, the training command will instead look like this: 






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

### Training OpenFold Multimer 

At this time, we do not have a multimer training set available. To prepare your own multimer training set, please see the instructions at [Data Processing - multimer] 

The basic command for training a multimer model is then:

```
multimer training command here
```

The key differences are:
- Dataset configuration / preparation

### Fine tuning from existing model weights 

If you have existing model weights, you can fine tune the model using the following command:

```
python3 train_openfold.py mmcif_dir/ alignment_dir/ template_mmcif_dir/ $OUTPUT_DIR \
    --max_template_date 2021-10-10 \ 
    --train_chain_data_cache_path chain_data_cache.json \
    --template_release_dates_cache_path mmcif_cache.json \ 
	--config_preset finetuning \
    --seed 4242022 \
    --obsolete_pdbs_file_path obsolete.dat \
    --num_nodes 1 \
    --gpus 4 \
    --num_workers 4 \
	--resume_from_ckpt $CHECKPOINT_PATH
	--resume_model_weights_only
```

If you have model parameters from OpenFold v1.x, you may need to convert your checkpoint file or parameter. See [[Converting OpenFold v1 Weights]]  

### Using MPI

If MPI is configured on your system, and you would like to use MPI to train OpenFold models, you may do so with the following step:

 1. Add the `mpi4py` package, which are available through pip and conda. Please see [mpi4py documentation](https://pypi.org/project/mpi4py/) for more instructions on installation.
2. Add the `--mpi_plugin` flag to your training command.
