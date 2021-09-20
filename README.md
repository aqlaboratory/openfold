# OpenFold

A faithful PyTorch reproduction of DeepMind's 
[AlphaFold 2](https://github.com/deepmind/alphafold).

## Installation

1. Having installed Python 3.9, install dependencies.

```bash
pip3 install -r requirements.txt
```

2. Install third-party dependencies.

```bash
scripts/install_third_party_dependencies.py
```

## Features

OpenFold fully reproduces AlphaFold's inference and data processing pipelines, 
but is written entirely in PyTorch rather than JAX/TensorFlow. It is capable of 
importing AlphaFold's original pretrained weights.

Future versions will support multi-GPU training with 
[DeepSpeed](https://github.com/microsoft/DeepSpeed).
