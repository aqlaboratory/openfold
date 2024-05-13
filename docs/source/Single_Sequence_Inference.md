### Soloseq inference

MSA-free sequence to structure prediction using the [ESM-1b model](https://github.com/facebookresearch/esm) embeddings.

To run inference for a sequence using the SoloSeq single-sequence model, you can either precompute ESM-1b embeddings in bulk, or you can generate them during inference.

For generating ESM-1b embeddings in bulk, use the provided script: [`scripts/precompute_embeddings.py`](https://github.com/aqlaboratory/openfold/blob/main/scripts/precompute_embeddings.py). The script takes a directory of FASTA files (one sequence per file) and generates ESM-1b embeddings in the same format and directory structure as required by SoloSeq. Following is an example command to use the script:

```shell
python scripts/precompute_embeddings.py fasta_dir/ embeddings_output_dir/
```

In the same per-label subdirectories inside `embeddings_output_dir`, you can also place `*.hhr` files (outputs from HHSearch), which can contain the details about the structures that you want to use as templates. If you do not place any such file, templates will not be used and only the ESM-1b embeddings will be used to predict the structure. If you want to use templates, you need to pass the PDB MMCIF dataset to the command.

Then download the SoloSeq model weights, e.g.:

```shell
bash scripts/download_openfold_soloseq_params.sh openfold/resources
```

Now, you are ready to run inference:

```shell
python run_pretrained_openfold.py \
    fasta_dir \
    data/pdb_mmcif/mmcif_files/ \
    --use_precomputed_alignments embeddings_output_dir \
    --output_dir ./ \
    --model_device "cuda:0" \
    --config_preset "seq_model_esm1b_ptm" \
    --openfold_checkpoint_path openfold/resources/openfold_soloseq_params/seq_model_esm1b_ptm.pt
```

For generating the embeddings during inference, skip the `--use_precomputed_alignments` argument. The `*.hhr` files will be generated as well if you pass the paths to the relevant databases and tools, as specified in the command below. If you skip the database and tool arguments, HHSearch will not be used to find templates and only generated ESM-1b embeddings will be used to predict the structure.

```shell
python3 run_pretrained_openfold.py \
    fasta_dir \
    data/pdb_mmcif/mmcif_files/ \
    --output_dir ./ \
    --model_device "cuda:0" \
    --config_preset "seq_model_esm1b_ptm" \
    --openfold_checkpoint_path openfold/resources/openfold_soloseq_params/seq_model_esm1b_ptm.pt \
    --uniref90_database_path data/uniref90/uniref90.fasta \
    --pdb70_database_path data/pdb70/pdb70 \
    --jackhmmer_binary_path lib/conda/envs/openfold_venv/bin/jackhmmer \
    --hhsearch_binary_path lib/conda/envs/openfold_venv/bin/hhsearch \
    --kalign_binary_path lib/conda/envs/openfold_venv/bin/kalign \
```

For generating template information, you will need the UniRef90 and PDB70 databases and the JackHmmer and HHSearch binaries.

SoloSeq allows you to use the same flags and optimizations as the MSA-based OpenFold. For example, you can skip relaxation using `--skip_relaxation`, save all model outputs using `--save_outputs`, and generate output files in MMCIF format using `--cif_output`.

```{note}
Due to the nature of the ESM-1b embeddings, the sequence length for inference using the SoloSeq model is limited to 1022 residues. Sequences longer than that will be truncated.
```