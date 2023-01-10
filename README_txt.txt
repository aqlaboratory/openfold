# OpenFold

A faithful PyTorch reproduction of DeepMind's 
[AlphaFold 2](https://github.com/deepmind/alphafold).

## Features

OpenFold carefully reproduces (almost) all of the features of the original open
source inference code. The sole exception is model ensembling, which fared
poorly in DeepMind's own ablation testing and is being phased out in future
DeepMind experiments. It is omitted here for the sake of reducing clutter. In 
cases where the Nature paper differs from the source, we always defer to the 
latter.

Unlike DeepMind's public code, OpenFold is also trainable. It can be trained 
with or without [DeepSpeed](https://github.com/microsoft/deepspeed) and with 
mixed precision. bfloat16 training is not currently supported, but will be 
soon. 

## Installation

Python dependencies available through `pip` are provided in `requirements.txt`. 
OpenFold also depends on `openmm==7.5.1` and `pdbfixer`, which are only
available via `conda`. 

For convenience, we provide a script that installs Miniconda locally, creates a 
`conda` virtual environment, installs all Python dependencies, and downloads
useful resources (including DeepMind's pretrained parameters). Run:

```bash
scripts/install_third_party_dependencies.sh
```

To activate the environment, run:

```bash
source scripts/activate_conda_venv.sh
```

To deactivate it, run:

```bash
source scripts/deactivate_conda_venv.sh
```

## Download Databases
There are multiple genetic databases that OpenFold uses to generate features.
We have provided scripts for downloading these databases. `aria2c` is a 
dependency for running them. 

There are two options for these databases.

* Full databases
```bash
scripts/download_all_data.sh <DOWNLOAD_DIRECTORY>
```

* Reduced databases
```bash
scripts/download_all_data.sh <DOWNLOAD_DIRECTORY> reduced_dbs
```

To process these databases, OpenFold also requires some external tools as a
dependency. See further in this document for how to set them up.

## Usage

To run inference on a sequence using a set of DeepMind's pretrained parameters, 
run e.g.

```bash
python3 run_pretrained_alphafold.py 
    --model model_1_ptm
    --preset reduced_dbs 
    --device cuda:1 
    --fasta_path target.fasta 
    --output_dir output_dir 
    --param_path openfold/resources/params/params_model_1.npz 
    --uniref90_database_path <DATABASE_DIR>/uniref90/uniref90.fasta 
    --mgnify_database_path <DATABASE_DIR>/mgnify/mgy_clusters_2018_12.fa 
    --pdb70_database_path <DATABASE_DIR>/pdb70/pdb70 
    --template_mmcif_dir <DATABASE_DIR>/pdb_mmcif/mmcif_files 
    --max_template_date 2100-01-01 
    --small_bfd_database_path <DATABASE_DIR>/small_bfd/bfd-first_non_consensus_sequences.fasta 
```

For running using the full databases:

```bash
python3 run_pretrained_alphafold.py
    --model model_1_ptm
    --preset full_dbs
    --device cuda:1
    --fasta_path target6_fulldbs.fasta 
    --output_dir output_dir 
    --param_path openfold/resources/params/params_model_1.npz
    --uniref90_database_path <DATABASE_DIR>/uniref90/uniref90.fasta 
    --mgnify_database_path <DATABASE_DIR>/mgnify/mgy_clusters_2018_12.fa 
    --pdb70_database_path <DATABASE_DIR>/pdb70/pdb70 
    --template_mmcif_dir <DATABASE_DIR>/pdb_mmcif/mmcif_files 
    --max_template_date 2100-01-01 
    --bfd_database_path <DATABASE_DIR>/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt 
    --uniclust30_database_path <DATABASE_DIR>/uniclust30/uniclust30_2018_08/uniclust30_2018_08 
```

## Copyright notice

While AlphaFold's and, by extension, OpenFold's source code is licensed under
the permissive Apache Licence, Version 2.0, DeepMind's pretrained parameters 
remain under the more restrictive CC BY-NC 4.0 license, a copy of which is 
downloaded to `openfold/resources/params` by the installation script. They are
thereby made unavailable for commercial use.
