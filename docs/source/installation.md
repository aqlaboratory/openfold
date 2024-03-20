# Setup Guide

In this guide, we will OpenFold and its dependencies.

**Pre-requisites**

This package is currently supported for CUDA 12 and Pytorch 2.1. All dependencies are listed in the `environment.yml`

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
	> scripts/install_third_party_dependencies.sh`
3. Prepend the conda environment to the $LD_LIBRARY_PATH., e.g. 
		`export $LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH``. You may optionally set this as a conda environment variable according to the [conda docs](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables) to activate each time the environment is used.
4. Download parameters. We recommend using a destination as `openfold/resources` as our unittests will look for the weights there.
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

## Modifications

### CUDA 11 environment
To use OpenFold on CUDA 11 environment rather than a CUDA 12 environment.
	In step 1, replace the github repository link with [OpenFold version v2](https://github.com/aqlaboratory/openfold/tree/v2.0.0)
	Follow the rest of the steps of [Installation Guide](#installation)

```{note}
Replace with link to last stable pre-update version
```

### Install OpenFold parameters without aws
If you don't have access to `aws` on your system, you can use a different download source:

- HuggingFace (requires `git-lts`):	`scripts/download_openfold_params_huggingface.sh`
- Google Drive: `scripts/download_openfold_params_gdrive.sh`

### Docker setup

```{note}
Add / check docker installation instructions 
```

## Troubleshooting FAQ

- In the unit tests, I see an error such as  
	```
	ImportError: version GLIBCXX_3.4.30 not found
	```

	> Solution: Make sure that the `$LD_LIBRARY_PATH` environment has been set to include the conda path, e.g. `export $LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`

- I see a CUDA mismatch error, eg. 
```
The detected CUDA version (11.8) mismatches the version that was used to compile
PyTorch (12.1). Please make sure to use the same CUDA versions.
```
 
 > 	Solution: Ensure that your system's CUDA driver and toolkit are is 12.x.  You can check the CUDA driver version with a command such as `nvidia-smi`

- I get some error involving `fatal error: cuda_runtime.h: No such file or directory` and or `ninja: build stopped: subcommand failed.`. 

> Solution: Something went wrong with setting up some of the custom kernels. Try running `install_third_party_dependencies.sh` again or try `python3 setup.py install` from inside the OpenFold folder. Make sure to prepend the conda environment as described above before running this.