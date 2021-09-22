# OpenFold

A faithful PyTorch reproduction of DeepMind's 
[AlphaFold 2](https://github.com/deepmind/alphafold).

## Installation

Python dependencies available through `pip` are provided in `requirements.txt`. 
OpenFold also depends on `openmm==7.5.1` and `pdbfixer`, which are only
available via `conda`. 

For convenience, we provide a script that installs Miniconda locally, creates a 
`conda` virtual environment, and installs all Python dependencies. Run:

```bash
scripts/install_third_party_dependencies.sh
```

To activate the environment, run:

```bash
scripts/activate_conda_venv.sh
```

To deactivate it, run

```bash
scripts/deactivate_conda_venv.sh
```

## Features

OpenFold reproduces AlphaFold's inference and data processing pipelines. With 
the exception of model ensembling, which fared poorly in DeepMind's ablation 
testing and is therefore omitted here, OpenFold supports all features of the 
original required for inference. It is even capable of importing AlphaFold's 
original pretrained weights. 

Future versions will support multi-GPU training with 
[DeepSpeed](https://github.com/microsoft/DeepSpeed).
