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

## Usage

To run inference on a sequence using a set of DeepMind's pretrained parameters, 
run e.g.

```bash
python3 run_pretrained_openfold.py \
        test.fasta \
        data/uniref90/uniref90.fasta \
        data/mgnify/mgy_clusters_2018_12.fa \
        data/pdb70/pdb70 \
        data/pdb_mmcif/mmcif_files/ \
        --output_dir ./ \
        --bfd_database_path data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
        --uniclust30_database_path data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
        --device cuda:1
```

where `data` is a directory populated by `scripts/download_all_data.sh`. Run

```bash
python3 run_pretrained_openfold.py --help
```

for a full list of options.

## Copyright notice

While AlphaFold's and, by extension, OpenFold's source code is licensed under
the permissive Apache Licence, Version 2.0, DeepMind's pretrained parameters 
remain under the more restrictive CC BY-NC 4.0 license, a copy of which is 
downloaded to `openfold/resources/params` by the installation script. They are
thereby made unavailable for commercial use.
