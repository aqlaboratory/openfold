# Setting up the OpenFold PDB training set from RODA

The multiple sequence alignments of OpenProteinSet and mmCIF structure files required to train OpenFold are freely available at the [Registry of Open Data on AWS (RODA)](https://registry.opendata.aws/openfold/). Additionally, OpenFold requires some postprocessing and [auxiliary files](Aux_seq_files.md) for training that need to be generated from the AWS data manually. This documentation is intended to give a full overview of those steps starting from the data download.

### Pre-Requisites:
- OpenFold conda environment. See [OpenFold Installation](Installation.md) for instructions on how to build this environment. 
- In particular, the [AWS CLI](https://aws.amazon.com/cli/) is used to download data from RODA.
- For this guide, we assume that the OpenFold codebase is located at `$OF_DIR`.

## 1. Downloading alignments and structure files
To fetch all the alignments corresponding to the original PDB training set of OpenFold alongside their mmCIF 3D structures, you can run the following commands:

```bash
mkdir -p alignment_data/alignment_dir_roda
aws s3 cp s3://openfold/pdb/ alignment_data/alignment_dir_roda/ --recursive --no-sign-request

mkdir pdb_data
aws s3 cp s3://openfold/pdb_mmcif.zip pdb_data/ --no-sign-request
aws s3 cp s3://openfold/duplicate_pdb_chains.txt . --no-sign-request
unzip pdb_mmcif.zip -d pdb_data
```

The nested alignment directory structure is not yet exactly what OpenFold expects, so you can run the `flatten_roda.sh` script to convert them to the correct format:

```bash
bash $OF_DIR/scripts/flatten_roda.sh alignment_data/alignment_dir_roda alignment_data/
```

Afterwards, the old directory can be safely removed:

```bash
rm -r alignment_data/alignment_dir_roda
```

## 2. Creating alignment DBs (optional)
As further explained in [Auxiliary Sequence Files in OpenFold](Aux_seq_files.md), OpenFold supports an alternate format for storing alignments that can increase training performance in I/O bottlenecked systems. These so-called `alignment_db` files can be generated with the following script:

```bash
python $OF_DIR/scripts/alignment_db_scripts/create_alignment_db_sharded.py \
    alignment_data/alignments \
    alignment_data/alignment_dbs \
    alignment_db \
    --n_shards 10 \
    --duplicate_chains_file pdb_data/duplicate_pdb_chains.txt
```

We recommend creating 10 total `alignment_db` files (= "shards") for better
filesystem health and fast preprocessing, but note that this script will only run
optimally if the number of CPUs on your machine is at least as big as the number
of shards you are creating.

As an optional check, you can run the following command which should return $634,434$:

```bash
grep "files" alignment_data/alignment_dbs/alignment_db.index | wc -l
```

## 3. Adding duplicate chains to alignments (skip if step 2 was used)
To save space, the OpenProteinSet alignment database is stored without duplicates, meaning that only one representative alignment is stored for all chains with identical sequences in the PDB and duplicate instances are tracked with a [`duplicate_chains.txt`](Aux_seq_files.md#duplicate-pdb-chain-files) file. As OpenFold will select chains during training based on the chains in the alignment directory (or `alignment_db`), we therefore need to add those duplicate chains back in in order to train on the full conformational diversity of chains in the PDB.

If you've followed the optional Step 2, the `.index` file of your `alignment_db` files will have already been adjusted for duplicates and you can proceed to the next step. Otherwise, the standard alignment directory can be expanded to accommodate duplicates by inserting symlinked directories for the duplicate chains that point to their representative alignments:

```bash
python $OF_DIR/scripts/expand_alignment_duplicates.py \
    alignment_data/alignments \
    pdb_data/duplicate_pdb_chains.txt
```

As an optional check, the following command should return $634,434$:

```bash
ls alignment_data/alignments/ | wc -l
```

## 4. Generating cluster-files
The AlphaFold dataloader adjusts the sampling probability of chains by their inverse cluster size, so we need to generate these sequence clusters for our training set.

As a first step, we'll need a `.fasta` file of all sequences in the training set. This can be generated with the following scripts, depending on how you set up your alignment data in the previous steps:

**Use this if you set up the duplicate-expanded alignment directory (faster):**
```bash
python $OF_DIR/scripts/alignment_data_to_fasta.py \
    alignment_data/all-seqs.fasta \
    --alignment_dir alignment_data/alignments
```

**Use this if you set up the `alignment_db` files:**
```bash
python $OF_DIR/scripts/alignment_data_to_fasta.py \
    alignment_data/all-seqs.fasta \
    --alignment_db_index alignment_data/alignment_dbs/alignment_db.index
```

Next, we need to generate a cluster file at 40% sequence identity, which will contain all chains in a particular cluster on the same line. You'll need [MMSeqs2](https://github.com/soedinglab/MMseqs2?tab=readme-ov-file#installation) for this as well, which can be set up either in a conda environment or as a binary.

```bash
python $OF_DIR/scripts/fasta_to_clusterfile.py \
    alignment_data/all-seqs.fasta \
    alignment_data/all-seqs_clusters-40.txt \
    /path/to/mmseqs \
    --seq-id 0.4
```

## 5. Generating cluster-files
As a last step, OpenFold requires ["cache" files](Aux_seq_files.md#chain-cache-files-and-mmcif-cache-files) with metadata information for each chain that are used for choosing templates and samples during training.

The data caches for OpenProteinSet can be downloaded from RODA with the following:

```bash
aws s3 cp s3://openfold/data_caches/ pdb_data/ --recursive --no-sign-request
```
If you wish to create data caches for your own datasets, the steps to generate the cache are as follows:

```bash
mkdir pdb_data/data_caches

python $OF_DIR/scripts/generate_mmcif_cache.py \
    pdb_data/mmcif_files \
    pdb_data/data_caches/mmcif_cache.json \
    --no_workers 16
```

The chain-data-cache is used for filtering training samples and adjusting per-chain sampling probabilities and can be generated with the following script:

```bash
python $OF_DIR/scripts/generate_chain_data_cache.py \
    pdb_data/mmcif_files \
    pdb_data/data_caches/chain_data_cache.json \
    --cluster_file alignment_data/all-seqs_clusters-40.txt \
    --no_workers 16
```
