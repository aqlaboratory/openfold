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
Try it out with our [Colab notebook](https://colab.research.google.com/github/aqlaboratory/openfold/blob/main/notebooks/OpenFold.ipynb)
(not yet visible from Colab because the repo is still private).

Unlike DeepMind's public code, OpenFold is also trainable. It can be trained 
with [DeepSpeed](https://github.com/microsoft/deepspeed) and with mixed 
precision. `bfloat16` training is not currently supported, but will be in the 
future.

## Installation (Linux)

Python dependencies available through `pip` are provided in `requirements.txt`. 
OpenFold depends on `openmm==7.5.1` and `pdbfixer`, which are only available 
via `conda`. For producing sequence alignments, you'll also need `jackhmmer`, 
`kalign`, and the [HH-suite](https://github.com/soedinglab/hh-suite) installed 
on your system. Finally, some download scripts require `aria2c`.

Note that the required version of PyTorch Lightning is 1.5.0, which has not
yet been released. Install that package from the nightly build.

For convenience, we provide a script that installs Miniconda locally, creates a 
`conda` virtual environment, installs all Python dependencies, and downloads
useful resources (including DeepMind's pretrained parameters). Run:

```bash
scripts/install_third_party_dependencies.sh
```

To activate the environment, run:

```bash
source scripts/activate_conda_env.sh
```

To deactivate it, run:

```bash
source scripts/deactivate_conda_env.sh
```

To install the HH-suite to `/usr/bin`, run

```bash
# scripts/install_hh_suite.sh
```

## Usage

To download the genetic databases used by AlphaFold/OpenFold, run:

```bash
scripts/download_all_data.sh data/
```

Alternatively, you can use raw MSAs from 
[ProteinNet](https://github.com/aqlaboratory/proteinnet). After downloading
the database, use `scripts/prepare_proteinnet_msas.py` to convert the data into
a format recognized by the OpenFold parser. The resulting directory becomes the
`alignment_dir` used in subsequent steps.

### Inference

To run inference on a sequence `target.fasta` (e.g., `wget https://www.rcsb.org/fasta/entry/4DSN`) using a set of DeepMind's pretrained parameters, run e.g.

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
    --device cuda:1 \
    --jackhmmer_binary_path lib/conda/envs/openfold_venv/bin/jackhmmer \
    --hhblits_binary_path lib/conda/envs/openfold_venv/bin/hhblits \
    --hhsearch_binary_path lib/conda/envs/openfold_venv/bin/hhsearch \
    --kalign_binary_path lib/conda/envs/openfold_venv/bin/kalign
```

where `data` is the same directory as in the previous step. If `jackhmmer`, `hhblits`, `hhsearch` and `kalign` are available at the default path of `/usr/bin`, their `binary_path` command-line arguments can be dropped. 

### Training

After activating the OpenFold environment with `source scripts/activate_conda_env.sh`, install OpenFold by running

```bash
python setup.py install
```

To train the model, you will first need to precompute protein alignments. Create `mmcif_dir/` and download `.cif` files from the PDB (e.g., `wget https://files.rcsb.org/download/4DSN.cif`). Then run:

```bash
python3 scripts/precompute_alignments.py mmcif_dir/ alignment_dir/ \
    data/uniref90/uniref90.fasta \
    data/mgnify/mgy_clusters_2018_12.fa \
    data/pdb70/pdb70 \
    data/pdb_mmcif/mmcif_files/ \
    data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --bfd_database_path data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --cpus 16 \
    --jackhmmer_binary_path lib/conda/envs/openfold_venv/bin/jackhmmer \
    --hhblits_binary_path lib/conda/envs/openfold_venv/bin/hhblits \
    --hhsearch_binary_path lib/conda/envs/openfold_venv/bin/hhsearch \
    --kalign_binary_path lib/conda/envs/openfold_venv/bin/kalign
```
As noted before, you can skip the `binary_path` arguments if these binaries are at `/usr/bin`.
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
    --seed 42 \ # in multi-gpu settings, the seed must be specified
    --deepspeed_config_path deepspeed_config.json \
    --resume_from_ckpt ckpt_dir/
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
