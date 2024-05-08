# FAQ

Frequently asked questions or encountered issues when running OpenFold. 

## Setup

- When running unit tests (e.g. [`./scripts/run_unit_tests.sh`](https://github.com/aqlaboratory/openfold/blob/main/scripts/run_unit_tests.sh)), I see an error such as  
	```
	ImportError: version GLIBCXX_3.4.30 not found
	```

	> Solution: Make sure that the `$LD_LIBRARY_PATH` environment has been set to include the conda path, e.g. `export $LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`

- I see a CUDA mismatch error, eg. 
```
The detected CUDA version (11.8) mismatches the version that was used to compile
PyTorch (12.1). Please make sure to use the same CUDA versions.
```
 
 > 	Solution: Ensure that your system's CUDA driver and toolkit match your intended OpenFold installation (CUDA 11 by default).  You can check the CUDA driver version with a command such as `nvidia-smi`

- I get some error involving `fatal error: cuda_runtime.h: No such file or directory` and or `ninja: build stopped: subcommand failed.`. 

> Solution: Something went wrong with setting up some of the custom kernels. Try running `install_third_party_dependencies.sh` again or try `python3 setup.py install` from inside the OpenFold folder. Make sure to prepend the conda environment as described above before running this.

## Training

- My model training is hanging on the data loading step:
	 > Solution: While each system is different, a few general suggestions:
		 - Check your `$KMP_AFFINITY` environment setting and see if it is suitable for your system.
  		 - Adjust the number of data workers used to prepare data with the `--num_workers` setting. Increasing the number could help with dataset processing speed. However, to many workers could cause an OOM issue. 

- When I reload my pretrained model weights or checkpoints, I get `RuntimeError: Error(s) in loading state_dict for OpenFoldWrapper: Unexpected key(s) in state_dict:`
	> Solution: This suggests that your checkpoint / model weights are in OpenFold v1 format with outdated model layer names. Convert your weights/checkpoints following [this guide](convert_of_v1_weights.md).