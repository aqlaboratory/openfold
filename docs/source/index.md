# OpenFold

```{figure} ../imgs/of_banner.png
:width: 900px
:align: center
:alt: Comparison of OpenFold and AlphaFold2 predictions to the experimental structure of PDB 7KDX, chain B._
```
Welcome to the Documentation for OpenFold, the fully open source, trainable, PyTorch-based reproduction of DeepMind's 
[AlphaFold 2](https://github.com/deepmind/alphafold).

Here, you will find guides and documentation for:
- [Getting started with OpenFold](installation.md)!
- Learn how to [run inference with OpenFold](Inference.md)
- [Train your own OpenFold models](Training_OpenFold.md)
- Find guidance for setup and running OpenFold in the [FAQ](FAQ.md).

We also have a [Colab notebook](https://colab.research.google.com/github/aqlaboratory/openfold/blob/main/notebooks/OpenFold.ipynb) that can be used for single structure / multimer prediction.

Some portions of the documentation are still under migration from the original README, which can be found [here](original_readme.md).

# Features

OpenFold carefully reproduces (almost) all of the features of the original open
source monomer (v2.0.1) and multimer (v2.3.2) inference code. The sole exception is 
model ensembling, which fared poorly in DeepMind's own ablation testing and is being 
phased out in future DeepMind experiments. It is omitted here for the sake of reducing 
clutter. In cases where the *Nature* paper differs from the source, we always defer to the 
latter.

OpenFold is trainable in full precision, half precision, or `bfloat16` with or without DeepSpeed, 
and we've trained it from scratch, matching the performance of the original. 
We've publicly released model weights and our training data &mdash; some 400,000 
MSAs and PDB70 template hit files &mdash; under a permissive license. Model weights 
are available via scripts in this repository while the MSAs are hosted by the 
[Registry of Open Data on AWS (RODA)](https://registry.opendata.aws/openfold). 
Try out running inference for yourself with our [Colab notebook](https://colab.research.google.com/github/aqlaboratory/openfold/blob/main/notebooks/OpenFold.ipynb).

OpenFold also supports inference using AlphaFold's official parameters, and 
vice versa (see `scripts/convert_of_weights_to_jax.py`).

OpenFold has the following advantages over the reference implementation:

- **Faster inference** on GPU, sometimes by as much as 2x. The greatest speedups are achieved on Ampere or higher architecture GPUs.
- **Inference on extremely long chains**, made possible by our implementation of low-memory attention 
([Rabe & Staats 2021](https://arxiv.org/pdf/2112.05682.pdf)). OpenFold can predict the structures of
  sequences with more than 4000 residues on a single A100, and even longer ones with CPU offloading.
- **Custom CUDA attention kernels** modified from [FastFold](https://github.com/hpcaitech/FastFold)'s 
kernels support in-place attention during inference and training. They use 
4x and 5x less GPU memory than equivalent FastFold and stock PyTorch 
implementations, respectively.
- **Efficient alignment scripts** using the original AlphaFold HHblits/JackHMMER pipeline or [ColabFold](https://github.com/sokrypton/ColabFold)'s, which uses the faster MMseqs2 instead. We've used them to generate millions of alignments.
- **FlashAttention** support greatly speeds up MSA attention.
- **DeepSpeed DS4Sci_EvoformerAttention kernel** is a memory-efficient attention kernel developed as part of a collaboration between OpenFold and the DeepSpeed4Science initiative. The kernel provides substantial speedups for training and inference, and significantly reduces the model's peak device memory requirement by 13X. The model is 15% faster during the initial training and finetuning stages, and up to 4x faster during inference.

# Copyright Notice

While AlphaFold's and, by extension, OpenFold's source code is licensed under
the permissive Apache Licence, Version 2.0, DeepMind's pretrained parameters 
fall under the CC BY 4.0 license, a copy of which is downloaded to 
`openfold/resources/params` by the installation script. Note that the latter
replaces the original, more restrictive CC BY-NC 4.0 license as of January 2022.

## Contributing

If you encounter problems using OpenFold, feel free to create an issue! We also
welcome pull requests from the community.

## Citing this Work

Please cite our paper:

```bibtex
@article {Ahdritz2022.11.20.517210,
	author = {Ahdritz, Gustaf and Bouatta, Nazim and Floristean, Christina and Kadyan, Sachin and Xia, Qinghui and Gerecke, William and O{\textquoteright}Donnell, Timothy J and Berenberg, Daniel and Fisk, Ian and Zanichelli, Niccolò and Zhang, Bo and Nowaczynski, Arkadiusz and Wang, Bei and Stepniewska-Dziubinska, Marta M and Zhang, Shang and Ojewole, Adegoke and Guney, Murat Efe and Biderman, Stella and Watkins, Andrew M and Ra, Stephen and Lorenzo, Pablo Ribalta and Nivon, Lucas and Weitzner, Brian and Ban, Yih-En Andrew and Sorger, Peter K and Mostaque, Emad and Zhang, Zhao and Bonneau, Richard and AlQuraishi, Mohammed},
	title = {{O}pen{F}old: {R}etraining {A}lpha{F}old2 yields new insights into its learning mechanisms and capacity for generalization},
	elocation-id = {2022.11.20.517210},
	year = {2022},
	doi = {10.1101/2022.11.20.517210},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/10.1101/2022.11.20.517210},
	eprint = {https://www.biorxiv.org/content/early/2022/11/22/2022.11.20.517210.full.pdf},
	journal = {bioRxiv}
}
```
If you use OpenProteinSet, please also cite:

```bibtex
@misc{ahdritz2023openproteinset,
      title={{O}pen{P}rotein{S}et: {T}raining data for structural biology at scale}, 
      author={Gustaf Ahdritz and Nazim Bouatta and Sachin Kadyan and Lukas Jarosch and Daniel Berenberg and Ian Fisk and Andrew M. Watkins and Stephen Ra and Richard Bonneau and Mohammed AlQuraishi},
      year={2023},
      eprint={2308.05326},
      archivePrefix={arXiv},
      primaryClass={q-bio.BM}
}
```
Any work that cites OpenFold should also cite [AlphaFold](https://www.nature.com/articles/s41586-021-03819-2) and [AlphaFold-Multimer](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1) if applicable.


```{toctree}
:hidden: 
:caption: Guides
Installation.md
Inference.md
Single_Sequence_Inference.md
Multimer_Inference.md
OpenFold_Training_Setup.md
Training_OpenFold.md
```

```{toctree}
:hidden: 
:caption: Reference 
Aux_seq_files.md
OpenFold_Parameters.md
FAQ.md
original_readme.md
```