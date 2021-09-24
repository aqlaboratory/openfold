# OpenFold

A faithful PyTorch reproduction of DeepMind's 
[AlphaFold 2](https://github.com/deepmind/alphafold).

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
scripts/activate_conda_venv.sh
```

To deactivate it, run:

```bash
scripts/deactivate_conda_venv.sh
```

## Features

OpenFold supports all features of the original AlphaFold required for inference
with the sole exception of model ensembling, which fared poorly in DeepMind's 
ablation testing. It is even capable of importing AlphaFold's original 
pretrained model parameters. 

Future versions will support multi-GPU training with 
[DeepSpeed](https://github.com/microsoft/DeepSpeed).

## Copyright notice

While AlphaFold's and, by extension, OpenFold's source code is licensed under
the permissive Apache Licence, Version 2.0, DeepMind's pretrained parameters 
remain under the more restrictive CC BY-NC 4.0 license, a copy of which is 
downloaded to `openfold/resources/params` by the installation script. They are
thereby made unavailable for commercial use.
