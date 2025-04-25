# Multimer Inference

To run inference on a complex or multiple complexes using a set of DeepMind's pretrained parameters, you will need:

- AlphaFold Multimer v2.3 parameters
- Updated sequence databases, with UniRef and PDB Seqres databases. 


## Upgrade from a previous OpenFold Installation 

If you had previously downloaded OpenFold parameters and or AlphaFold databases, you will need to download updated versions. Here are some instructions for upgrading from an existing openfold installations. 

### Download AlphaFold-Multimer v2.3 Model Parameters 
1. Re-download the alphafold parameters to get the latest
AlphaFold-Multimer v2.3 weights:
    
   ```bash
    bash scripts/download_alphafold_params.sh openfold/resources
   ```

### Download AlphaFold Databases for Multimer 

1. Download the [UniProt](https://www.uniprot.org/uniprotkb/) 
and [PDB SeqRes](https://www.rcsb.org/) databases: 
    
   ```bash
    bash scripts/download_uniprot.sh data/
   ```
    
    The PDB SeqRes and PDB databases must be from the same date to avoid potential 
    errors during template searching. Remove the existing `data/pdb_mmcif` directory 
    and download both databases:
    
   ```bash
    bash scripts/download_pdb_mmcif.sh data/
    bash scripts/download_pdb_seqres.sh data/
   ```

1. Additionally, AlphaFold-Multimer uses upgraded versions of the [MGnify](https://www.ebi.ac.uk/metagenomics) 
and [UniRef30](https://uniclust.mmseqs.com/) (previously UniClust30) databases. To download the upgraded databases, run:
    
   ```bash
    bash scripts/download_uniref30.sh data/
    bash scripts/download_mgnify.sh data/
   ```

```{note}
Multimer inference can also run with the older database versions if desired. 
```

## Running Multimer Inference 

The [`run_pretrained_openfold.py`](https://github.com/aqlaboratory/openfold/blob/main/run_pretrained_openfold.py) script can be used to run multimer inference with the follwoing command.

```bash
python3 run_pretrained_openfold.py \
    fasta_dir \
    data/pdb_mmcif/mmcif_files/ \
    --uniref90_database_path data/uniref90/uniref90.fasta \
    --mgnify_database_path data/mgnify/mgy_clusters_2022_05.fa \
    --pdb_seqres_database_path data/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path data/uniref30/UniRef30_2021_03 \
    --uniprot_database_path data/uniprot/uniprot.fasta \
    --bfd_database_path data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --jackhmmer_binary_path lib/conda/envs/openfold_venv/bin/jackhmmer \
    --hhblits_binary_path lib/conda/envs/openfold_venv/bin/hhblits \
    --hmmsearch_binary_path lib/conda/envs/openfold_venv/bin/hmmsearch \
    --hmmbuild_binary_path lib/conda/envs/openfold_venv/bin/hmmbuild \
    --kalign_binary_path lib/conda/envs/openfold_venv/bin/kalign \
    --config_preset "model_1_multimer_v3" \
    --model_device "cuda:0" \
    --output_dir ./ 
```

**Notes:**
- Template searching in the multimer pipeline uses HMMSearch with the PDB SeqRes database, replacing HHSearch and PDB70 used in the monomer pipeline.
- As with monomer inference, if you've already computed alignments for the query, you can use the `--use_precomputed_alignments` option.
- At this time, only AlphaFold parameter weights are available for multimer mode. 