# OpenFold

A faithful PyTorch reproduction of DeepMind's 
[AlphaFold 2](https://github.com/deepmind/alphafold).

## Features

OpenFold carefully reproduces (almost) all of the features of the original open
source inference code. The sole exception is model ensembling, which fared
poorly in DeepMind's own ablation testing and is being phased out in future
DeepMind experiments. It is omitted here for the sake of reducing clutter. In 
cases where the *Nature* paper differs from the source, we always defer to the 
latter. 

OpenFold is built to support inference with AlphaFold's original JAX weights.
Try it out with our [Colab notebook](https://colab.research.google.com/github/aqlaboratory/openfold/blob/main/notebooks/OpenFold.ipynb).

Unlike DeepMind's public code, OpenFold is also trainable. It can be trained 
with or without [DeepSpeed](https://github.com/microsoft/deepspeed) and with 
mixed precision. `bfloat16` training is not currently supported, but will be 
in the future.

## Installation (Linux)

Python dependencies available through `pip` are provided in `requirements.txt`. 
OpenFold depends on `openmm==7.5.1` and `pdbfixer`, which are only available 
via `conda`. For producing sequence alignments, you'll also need `jackhmmer`, 
`kalign`, and the [HH-suite](https://github.com/soedinglab/hh-suite) installed 
on your system.

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

To download the genetic databases used by AlphaFold/OpenFold, run:

```bash
scripts/download_all_data.sh data/
```

This script depends on `aria2c`.

### Inference

To run inference on a sequence using a set of DeepMind's pretrained parameters, 
run e.g.

```bash
python3 run_pretrained_openfold.py \
    target.fasta \
    data/uniref90/uniref90.fasta \
    data/mgnify/mgy_clusters_2018_12.fa \
    data/pdb70/pdb70 \
    data/pdb_mmcif/mmcif_files/ \
    data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --output_dir ./ \
    --bfd_database_path data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --device cuda:1
```

where `data` is the same directory as in the previous step.

### Training

To train the model, you will first need to precompute protein alignments. After
installing OpenFold using `setup.py`, do so with:

```bash
python3 scripts/precompute_alignments.py mmcif_dir/ alignment_dir/ \
    data/uniref90/uniref90.fasta \
    data/mgnify/mgy_clusters_2018_12.fa \
    data/pdb70/pdb70 \
    data/pdb_mmcif/mmcif_files/ \
    data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --bfd_database_path data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --cpus 16
```

Expect this step to take a very long time, even for small numbers of proteins.

Next, generate a cache of certain datapoints in the mmCIF files:

```bash
python3 scripts/generate_mmcif_cache.py \
    mmcif_dir/ \
    mmcif_cache.json \
    --no_workers 16
```

This cache is used to minimize the number of mmCIF parses performed during 
training-time data preprocessing. Finally, call the training script:

```bash
python3 train_openfold.py mmcif_dir/ alignment_dir/ template_mmcif_dir/ \
    2021-10-10 \ 
    --template_release_dates_cache_path mmcif_cache.json \ 
    --precision 16 \
    --gpus 8 --replace_sampler_ddp=True \
    --accelerator ddp \ 
    --seed 42 \ # in multi-gpu settings, the seed must be specified
    --deepspeed_config_path deepspeed_config.json
```

where `--template_release_dates_cache_path` is a path to the `.json` file
generated in the previous step. A suitable DeepSpeed configuration file can be 
generated with `scripts/build_deepspeed_config.py`. The training script is 
written with [PyTorch Lightning](https://github.com/PyTorchLightning/pytorch-lightning) 
and supports the full range of training options that entails, including 
multi-node distributed training. For more information, consult PyTorch 
Lightning documentation and the `--help` flag of the training script.

## Testing

To run unit tests, use

```bash
scripts/run_unit_tests.sh
```

The script is a thin wrapper around Python's `unittest` suite, and recognizes
`unittest` commands. E.g., to run a specific test verbosely:

```bash
scripts/run_unit_tests.sh -v tests.test_model
```

Certain tests require that AlphaFold be installed in the same Python
environment. These run components of AlphaFold and OpenFold side by side and
ensure that output activations are adequately similar. For most modules, we
target a maximum difference of `1e-4`.

## Copyright notice

While AlphaFold's and, by extension, OpenFold's source code is licensed under
the permissive Apache Licence, Version 2.0, DeepMind's pretrained parameters 
remain under the more restrictive CC BY-NC 4.0 license, a copy of which is 
downloaded to `openfold/resources/params` by the installation script. They are
thereby made unavailable for commercial use.

## Contributing

If you encounter problems using OpenFold, feel free to create an issue!
