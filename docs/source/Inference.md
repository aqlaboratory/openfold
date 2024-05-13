# OpenFold Inference 

In this guide, we will cover how to use OpenFold to make structure predictions.

## Background

We currently offer three modes of inference prediction:

- Monomer
- Multimer
- Single Sequence (Soloseq) 

This guide will focus on monomer prediction, the next sections will describe [Multimer](Multimer_Inference.md) and [Single Sequence](Single_Sequence_Inference.md) prediction. 
`
### Pre-requisites: 

- OpenFold Conda Environment. See [OpenFold Installation](Installation.md) for instructions on how to build this environment. 
- Downloading sequence databases for performing multiple sequence alignments. We provide a script to download the AlphaFold databases [here](https://github.com/aqlaboratory/openfold/blob/main/scripts/download_alphafold_dbs.sh).
   

## Running AlphaFold Model Inference 

The script [`run_pretrained_openfold.py`](https://github.com/aqlaboratory/openfold/blob/main/run_pretrained_openfold.py) performs model inference. We will go through the steps of how to use this script.

An example directory for performing infernce on [PDB:6KWC](https://www.rcsb.org/structure/6KWC) is provided [here](https://github.com/aqlaboratory/openfold/tree/main/examples/monomer). We refer to this example directory for the below examples.

### Download Model Parameters 

For monomer inference, you may either use the model parameters provided by Deepmind, or you may use the OpenFold trained parameters. Both models should give similar performance, please see [our main paper](https://www.biorxiv.org/content/10.1101/2022.11.20.517210v3) for further reference.

The model parameters provided by Deepmind can be downloaded with the following script located in this repository's `scripts/` directory:

```
$ bash scripts/download_alphafold_params.sh $PARAMS_DIR
```

To use the OpenFold trained parameters, you can use the following script

```
$ bash scripts/download_openfold_params.sh $PARAMS_DIR
```

We recommend selecting `openfold/resources` as the params directory as this is the default directory used by the `run_pretrained_openfold.py` to locate parameters. 

If you choose to use a different directory, you may make a symlink to the `openfold/resources` directory, or specify an alternate parameter path with the command line argument `--jax_param_path` for AlphaFold parameters or `--openfold_checkpoint_path` for OpenFold parameters. 


### Model Inference 

The input to [`run_pretrained_openfold.py`](https://github.com/aqlaboratory/openfold/blob/main/run_pretrained_openfold.py) is a directory of FASTA files. AlphaFold-style models also require a sequence alignment to perform inference.

If you do not have sequence alignments for your input sequences, you can compute them using the inference script directly by following the instructions for the following section [inference without pre-computed alignments](#model-inference-without-pre-computed-alignments).

Otherwise, if you already have alignments for your input FASTA sequences, skip ahead to the [inference with pre-computed alignments](#model-inference-with-pre-computed-alignments) section. 

#### Model inference without pre-computed alignments 
The following command performs a sequence alignment against the OpenProteinSet databases and performs model inference. 

```
python3 run_pretrained_openfold.py \
    $INPUT_FASTA_DIR \
    $TEMPLATE_MMCIF_DIR 
    --output_dir $OUTPUT_DIR \
    --config_preset model_1_ptm \
    --uniref90_database_path $BASE_DATA_DIR/uniref90 \
    --mgnify_database_path $BASE_DATA_DIR/mgnify/mgy_clusters_2018_12.fa \
    --pdb70_database_path $BASE_DATA_DIR/pdb70 \
    --uniclust30_database_path $BASE_DATA_DIR/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --bfd_database_path $BASE_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --model_device "cuda:0" 
```

**Required arguments:**
- `--output_dir`: specify the output directory
- `$INPUT_FASTA_DIR`: Directory of query fasta files, one sequence per file,e.g. `examples/monomer/fasta_dir`
- `$TEMPLATE_MMCIF_DIR`: MMCIF files to use for template matching. This directory is required even if using template free inference. 
- `*_database_path`: Paths to sequence databases for sequence alignment.
- `--model_device`: Specify to use a GPU is one is available.

#### Model inference with pre-computed alignments 
To perform model inference with pre-computed alignments, use the following command

```
python3 run_pretrained_openfold.py ${INPUT_FASTA_DIR} \
  $TEMPLATE_MMCIF_DIR \
  --output_dir $OUTPUT_DIR \
  --use_precomputed_alignments $PRECOMPUTED_ALIGNMENTS \
  --config_preset model_1_ptm \
  --model_device "cuda:0" \
```

where `$PRECOMPUTED_ALIGNMENTS` is a directory that contains alignments. A sample alignments directory structure for a single query is:

```
alignments
└── 6KWC_1 
    ├── bfd_uniclust_hits.a3m
    ├── hhsearch_output.hhr
    ├── mgnify_hits.sto
    └── uniref90_hits.sto
```

`bfd_uniclust_hits.a3m`, `mgnify_hits.sto`, and `uniref90_hits.sto` are all alignments of the query structure against the BFD, Mgnify, and Uniref90 datasets respsectively. `hhsearch_output.hhr` contains hits against the PDB70 database used for template matching. The example directory `examples/monomer/alignments` shows examples of expected directories.


#### Configuration settings for template modeling / pTM scoring 
There are a few configuration settings available for template based and template-free modeling, and for the option to estimate a predicted template modeling score (pTM). 

This table provides guidance on which setting to use for each set of predictions, as well as the parameters to select for each preset.  

|                    Setting |                           `config_preset` | AlphaFold params (match config name)                                              | OpenFold params (any are allowed)  |
| -------------------------: | ----------------------------------------: | :-------------------------------------------------------------------------------- | :--------------------------------- |
|      With template, no ptm |                        model_1<br>model_2 | `parms_model_1.npz`<br>`parms_model_2.npz`                                        | `finetuning_[2-5].pt`              |
|    With template, with ptm |                model_1_ptm<br>model_2_ptm | `params_model_1_ptm.npz`<br>`params_model_2_ptm.npz`                              | `finetuning_ptm_[1-2].pt`          |
|   Without template, no ptm |             model_3<br>model_4<br>model_5 | `parms_model_3.npz`<br>`parms_model_4.npz`<br>`parms_model_5.npz`                 | `finetuning_no_templ_[1-2].pt`     |
| Without template, with ptm | model_3_ptm<br>model_4_ptm<br>model_5_ptm | `parms_model_3_ptm.npz`<br>`parms_model_4_ptm.npz`<br>`parms_model_5_ptm.npz`<br> | `finetuning_no_templ_ptm_1.pt` |

If you use AlphaFold parameters, and the AlphaFold parameters are located in the default parameter directory (e.g. `openfold/resources`) the parameters that match the `--config_preset` will be selected.

The full set of configurations available for all 5 AlphaFold model presets can be viewed in [`config.py`](https://github.com/aqlaboratory/openfold/blob/main/openfold/config.py#L105). The [OpenFold Parameters](OpenFold_Parameters.md) page contains more information about the individual OpenFold parameter files.


#### Model outputs 

The expected output contents are as follows: 
- `alignments`: Directory of alignments. One directory is made per query sequence, and each directory contains alignments against each of the databases used.
- `predictions`: PDB files for predicted structures
- `timings.json`: Json with timings for inference and relaxation, if specified 


### Optional Flags 

Some commonly used command line flags are here. A full list of flags can be viewed from the `--help` menu

- `--config_preset`: Specify a different model configuration. There are 5 available model preset settings, some of which support template modeling, others support template-free modeling. The default is `model_1`. More details can be below in the [[Inference#Template-free modeling]] section 
- `--hmmsearch_binary_path`, `--hmmbuild_binary_path`, etc.  Hmmer, HHsuite, kalign are required to run alignments. `run_pretrained_openfold.py` will search for these packages in the `bin/` directory of your conda environment. If needed, you can specify a different binary directory with these arguments.
- `--openfold_checkpoint_path` : Uses an checkpoint or parameter file. Expected types are Deepspeed checkpoint files or `.pt` files. Make sure your selected checkpoint file matches the configuration setting chosen in `--config_preset`.
- `--data_random_seed`: Specifies a random seed to use.
- `--save_outputs`: Saves a copy of all outputs from the model, e.g. the output of the msa track, ptm heads.
- `--experiment_config_json`: Specify configuration settings using a json file. For example, passing a json with `{globals.relax.max_iterations = 10}` specifies 10 as the maximum number of relaxation iterations. See for  [`openfold/config.py`](https://github.com/aqlaboratory/openfold/blob/main/openfold/config.py#L283) the full dictionary of configuration settings. Any parameters that are not manually set in these configuration settings will refer to the defaults specified by your `config_preset`.


### Advanced Options for Increasing Efficiency

#### Speeding up inference 

The **DeepSpeed DS4Sci_EvoformerAttention kernel** is a memory-efficient attention kernel developed as part of a collaboration between OpenFold and the DeepSpeed4Science initiative. 

If your system supports deepseed, using deepspeed generally leads an inference speedup of 2 - 3x without significant additional memory use. You may specify this option by selecting the `--use_deepspeed_inference` argument. 

If DeepSpeed is unavailable for your system, you may also try using [FlashAttention](https://github.com/HazyResearch/flash-attention) by adding `globals.use_flash = True` to the `--experiment_config_json`. Note that FlashAttention appears to work best for sequences with < 1000 residues.

#### Large-scale batch inference 
For large-scale batch inference, we offer an optional tracing mode, which massively improves runtimes at the cost of a lengthy model compilation process. To enable it, add `--trace_model` to the inference command.

#### Configuring the chunk size for sequence alignments
Note that chunking (as defined in section 1.11.8 of the AlphaFold 2 supplement) is enabled by default in inference mode. To disable it, set `globals.chunk_size` to `None` in the config. If a value is specified, OpenFold will attempt to dynamically tune it, considering the chunk size specified in the config as a minimum. This tuning process automatically ensures consistently fast runtimes regardless of input sequence length, but it also introduces some runtime variability, which may be undesirable for certain users. It is also recommended to disable this feature for very long chains (see below). To do so, set the `tune_chunk_size` option in the config to `False`.

#### Long sequence inference 
To minimize memory usage during inference on long sequences, consider the following changes:

- As noted in the AlphaFold-Multimer paper, the AlphaFold/OpenFold template stack is a major memory bottleneck for inference on long sequences. OpenFold supports two mutually exclusive inference modes to address this issue. One, `average_templates` in the `template` section of the config, is similar to the solution offered by AlphaFold-Multimer, which is simply to average individual template representations. Our version is modified slightly to accommodate weights trained using the standard template algorithm. Using said weights, we notice no significant difference in performance between our averaged template embeddings and the standard ones. The second, `offload_templates`, temporarily offloads individual template embeddings into CPU memory. The former is an approximation while the latter is slightly slower; both are memory-efficient and allow the model to utilize arbitrarily many templates across sequence lengths. Both are disabled by default, and it is up to the user to determine which best suits their needs, if either.
- Inference-time low-memory attention (LMA) can be enabled in the model config. This setting trades off speed for vastly improved memory usage. By default, LMA is run with query and key chunk sizes of 1024 and 4096, respectively. These represent a favorable tradeoff in most memory-constrained cases. Powerusers can choose to tweak these settings in `openfold/model/primitives.py`. For more information on the LMA algorithm, see the aforementioned Staats & Rabe preprint.
- Disable `tune_chunk_size` for long sequences. Past a certain point, it only wastes time.
- As a last resort, consider enabling `offload_inference`. This enables more extensive CPU offloading at various bottlenecks throughout the model.
- Disable FlashAttention, which seems unstable on long sequences.

Using the most conservative settings, we were able to run inference on a 4600-residue complex with a single A100. Compared to AlphaFold's own memory offloading mode, ours is considerably faster; the same complex takes the more efficent AlphaFold-Multimer more than double the time. Use the `long_sequence_inference` config option to enable all of these interventions at once. The `run_pretrained_openfold.py` script can enable this config option with the `--long_sequence_inference` command line option

Input FASTA files containing multiple sequences are treated as complexes. In this case, the inference script runs AlphaFold-Gap, a hack proposed [here](https://twitter.com/minkbaek/status/1417538291709071362?lang=en), using the specified stock AlphaFold/OpenFold parameters (NOT AlphaFold-Multimer).