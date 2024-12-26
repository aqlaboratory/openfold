# Setting Up OpenFold 

In this guide, we will OpenFold and its dependencies.

**Pre-requisites**

This package is currently supported for CUDA 11 and Pytorch 1.12. All dependencies are listed in the [`environment.yml`](https://github.com/aqlaboratory/openfold/blob/main/environment.yml). To install OpenFold for CUDA 12, please refer to the [Environment specific modifications](#Environment-specific-modifications) section.

At this time, only Linux systems are supported.

## Instructions
:::

### Installation:
1. Clone the repository, e.g. `git clone https://github.com/aqlaboratory/openfold.git`
1. From the `openfold` repo:
    - Create a [Mamba]("https://github.com/conda-forge/miniforge/releases/latest/download/) environment, e.g.
        `mamba env create -n openfold_env -f environment.yml`
      Mamba is recommended as the dependencies required by OpenFold are quite large and mamba can speed up the process.
    - Activate the environment, e.g `conda activate openfold_env`
1. Run the setup script to configure kernels and folding resources.
	> scripts/install_third_party_dependencies.sh
1. Prepend the conda environment to the `$LD_LIBRARY_PATH` and `$LIBRARY_PATH`., e.g. 

	```
	export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
	export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
	```

	You may optionally set this as a conda environment variable according to the [conda docs](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables) to activate each time the environment is used.

1. Download parameters. We recommend using a destination as `openfold/resources` as our unittests will look for the weights there.
	-  For AlphaFold2 weights, use 
		> ./scripts/download_alphafold_params.sh <dest>
	 - For OpenFold weights, use : 
		>  ./scripts/download_openfold_params.sh <dest>
	 - For OpenFold SoloSeq weights, use: 
		> ./scripts/download_openfold_soloseq_params.sh <dest>

### Checking your build with unit tests: 

To test your installation, you can run OpenFold unit tests. Make sure that the OpenFold and AlphaFold parameters have been downloaded, and that they are located (or symlinked) in the directory `openfold/resources` 

Run with the following script:
> scripts/run_unit_tests.sh

The script is a thin wrapper around Python's `unittest` suite, and recognizes `unittest` arguments. E.g., to run a specific test verbosely:

> scripts/run_unit_tests.sh -v tests.test_model

**Alphafold Comparison tests:**
Certain tests perform equivalence comparisons with the AlphaFold implementation. Instructions to run this level of tests requires an environment with both AlphaFold 2.0.1 and OpenFold installed, and is not covered in this guide. These tests are skipped by default if no installation of AlphaFold is found. 

## Environment specific modifications 

### CUDA 12
To use OpenFold on CUDA 12 environment rather than a CUDA 11 environment.
	In step 1, use the branch [`pl_upgrades`](https://github.com/aqlaboratory/openfold/tree/pl_upgrades) rather than the main branch, i.e. replace the command in step 1 with `git clone -b pl_upgrades https://github.com/aqlaboratory/openfold.git`
	and follow the rest of the steps of [Installation Guide](#Installation)


### MPI
To use OpenFold with MPI support, you will need to add the package [`mpi4py`](https://pypi.org/project/mpi4py/). This can be done with pip in your OpenFold environment, e.g. `$ pip install mpi4py`. 


### Install OpenFold parameters without aws
If you don't have access to `aws` on your system, you can use a different download source:

- HuggingFace (requires `git-lts`):	`scripts/download_openfold_params_huggingface.sh`
- Google Drive: `scripts/download_openfold_params_gdrive.sh`

### Docker setup

A [`Dockerfile`] is provided to build an OpenFold Docker image. Additional notes for setting up a docker container for OpenFold and running inference can be found [here](original_readme.md#building-and-using-the-docker-container).
