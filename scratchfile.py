# Determining bond lengths in the predicted PDB
# import openfold.np.protein
# import Bio.PDB
#
# parser = Bio.PDB.PDBParser(QUIET=True)
# structures = parser.get_structure('1IQQ', 'output_dir/epoch_msa_1IQQ_unrelaxed.pdb')
# structure = structures[0]
# residue_list = [_ for _ in structure.get_residues()]
# bb_NCA_bond_lens = []
# bb_CAC_bond_lens = []
# bb_CN_bond_lens = [] # backbone inter-residue CN bond lengths
# bb_CACA_bond_lens = [] # inter-residue CA-CA bond length
#
# for residue1, residue2 in zip(residue_list[:-1], residue_list[1:]):
#     residue1_dict = {a.get_name(): a for a in residue1.get_atoms()}
#     residue2_dict = {a.get_name(): a for a in residue2.get_atoms()}
#     bb_NCA_bond_lens.append(residue1_dict['N']-residue1_dict['CA'])
#     bb_CAC_bond_lens.append(residue1_dict['CA']-residue1_dict['C'])
#     bb_CN_bond_lens.append(residue1_dict['C']-residue2_dict['N'])
#     bb_CACA_bond_lens.append(residue1_dict['CA']-residue2_dict['CA'])
# # *Have not* added the NCA and CAC bond lens from last residue in chain
# print('\n======== NCA =========')
# print(bb_NCA_bond_lens)
# print('\n======== CAC =========')
# print(bb_CAC_bond_lens)
# print('\n======== CN =========')
# print(bb_CN_bond_lens)
# print('\n======== CACA ========')
# print(bb_CACA_bond_lens)

# =================================
# ESM 1-b
#import os

#import torch

# Load ESM-1b model
# import numpy as np
#
# from openfold.data import parsers
#
# model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
# batch_converter = alphabet.get_batch_converter()
# model.eval()  # disables dropout for deterministic results
#
# # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
# with open('../esm_embeddings_pdb/for_esm') as f:
#     seqs, descs = parsers.parse_fasta(f.read())
#
# data = list((desc, seq) for desc, seq in zip(descs, seqs))
# batch_labels, batch_strs, batch_tokens = batch_converter(data)
#
# # Extract per-residue representations (on CPU)
# with torch.no_grad():
#     results = model(batch_tokens, repr_layers=[33], return_contacts=True)
# token_representations = results["representations"][33]
#
# # Generate per-sequence representations via averaging
# # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
# sequence_representations = []
# for i, (id, seq) in enumerate(data):
#     with open(os.path.join('../esm_embeddings_pdb/for_esm_np', id+'.npy'), 'wb') as f:
#         np.save(f, token_representations[i, 1 : len(seq)+1])

#################
# Find the unique PDB sequences in the ESM embeddings sub dirs
import os
fasta_labels = []
fasta_seqs = []
with open("../../pdb_seqres_unique.txt") as infile:
    while True:
            label = infile.readline()
            seq = infile.readline()
            if not seq: break
            fasta_labels.append(label.strip())
            fasta_seqs.append(seq.strip())
esm_base_path = "../../esm_embeddings_pdb/chunk"
esm_file_chunks = [-1]*7
for index in range(7):
    path = esm_base_path + str(index)
    esm_file_chunks[index] = os.listdir(path)
for index in range(len(fasta_labels)):
    fasta_labels[index] = fasta_labels[index][1:].split(' ')[0]

fasta_paths = [-1]*len(fasta_labels)
for fasta_index in range(len(fasta_labels)):
    label = fasta_labels[fasta_index]
    for chunk_index in range(7):
        if label+'.pt' in esm_file_chunks[chunk_index]:
            fasta_paths[fasta_index] = esm_base_path+str(chunk_index)+'/'+label+'.pt'
            break

# ===============================
# To create the dummy alignments for the full Distillation set sequences
# import os

def parse_fasta(fasta_string: str):
    """Parses FASTA string and returns list of strings with amino-acid sequences.

    Arguments:
        fasta_string: The string contents of a FASTA file.

    Returns:
        A tuple of two lists:
        * A list of sequences.
        * A list of sequence descriptions taken from the comment lines. In the
            same order as the sequences.
    """
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append("")
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions


# with open('../esm_embeddings_pdb/fastas/prot_data_cache.fasta') as fa_file:
#     seqs, labels = parse_fasta(fa_file.read())
#
# for seq, label in zip(seqs, labels):
#     seq_folder = os.path.join('../train_data_openfold/alignments', label)
#     if not os.path.exists(seq_folder):
#         os.makedirs(seq_folder)
#     with open(os.path.join(seq_folder,  'bfd_uniclust_hits.a3m'), 'w') as bfd_file:
#         lines = []
#         lines += '>'+label+'\n'
#         lines += str(seq+'\n')
#         bfd_file.writelines(lines)
#
#     with open(os.path.join(seq_folder,  'mgnify_hits.a3m'), 'w') as mgnify_file:
#         lines = []
#         lines += '>'+label+'\n'
#         lines += str(seq+'\n')
#         mgnify_file.writelines(lines)
#
#     with open(os.path.join(seq_folder,  'uniref90_hits.a3m'), 'w') as uniref_file:
#         lines = []
#         lines += '>'+label+'\n'
#         lines += str(seq+'\n')
#         uniref_file.writelines(lines)

### Extract templates out of the databases
# import json
# # Open the database in binary mode
# _alignment_index_path = None
# if(_alignment_index_path is not None):
#     with open(_alignment_index_path, "r") as fp:
#         _alignment_index = json.load(fp)
# all_hits = {}
# fp = open(os.path.join(dir, _alignment_index["db"]))
#
# def read_template(start, size):
#     fp.seek(start)
#     return fp.read(size).decode("utf-8")
#
# ##==== Don't need this, since we are not reading the templates
#     # from the DB directly
# # for (name, start, size) in _alignment_index["files"]:
# #     ext = os.path.splitext(name)[-1]
# #
# #     if(ext == ".hhr"):
# #         hits = parsers.parse_hhr(read_template(start, size))
# #         all_hits[name] = hits
#
# fp.close()

#################
# Find the unique PDB sequences in the ESM embeddings sub dirs
# import os
# fasta_labels = []
# fasta_seqs = []
# with open("esm_embeddings_pdb/chunk0.fasta") as infile:
#     while True:
#             label = infile.readline()
#             seq = infile.readline()
#             if not seq: break
#             fasta_labels.append(label.strip())
#             fasta_seqs.append(seq.strip())
# esm_base_path = "esm_embeddings_pdb/chunk"
# esm_file_chunks = [-1]*7
# for index in range(7):
#     path = esm_base_path + str(index)
#     esm_file_chunks[index] = os.listdir(path)
# for index in range(len(fasta_labels)):
#     fasta_labels[index] = fasta_labels[index][1:].split(' ')[0]
#
# with open('esm_embeddings_pdb/fastas/chunk0_rest.fasta', 'w') as chunk0_rest_fasta:
#     fasta_paths = [-1]*len(fasta_labels)
#     for fasta_index in range(len(fasta_labels)):
#         label = fasta_labels[fasta_index]
#         if label+'.pt' not in esm_file_chunks[0]:
#             chunk0_rest_fasta.write('>'+label+'\n')
#             chunk0_rest_fasta.write(str(fasta_seqs[fasta_index])+'\n')


# Delete proteins from the pdb_esm_embeddings folder which have >1022 residues
# Read the seqs from pdb_seqres.txt
import os
fasta_labels = []
fasta_seqs = []
with open("../../pdb_seqres_unique.txt") as infile:
    while True:
            label = infile.readline()
            seq = infile.readline()
            if not seq: break
            fasta_labels.append(label.strip())
            fasta_seqs.append(seq.strip())
esm_base_path = "../../esm_embeddings_pdb/chunk"
esm_file_chunks = [-1]*7
for index in range(7):
    path = esm_base_path + str(index)
    esm_file_chunks[index] = os.listdir(path)
for index in range(len(fasta_labels)):
    fasta_labels[index] = fasta_labels[index][1:].split(' ')[0]

fasta_paths = [-1]*len(fasta_labels)
for fasta_index in range(len(fasta_labels)):
    label = fasta_labels[fasta_index]
    for chunk_index in range(7):
        if label+'.pt' in esm_file_chunks[chunk_index]:
            fasta_paths[fasta_index] = esm_base_path+str(chunk_index)+'/'+label+'.pt'
            break

#### EXTRACT PROTEIN SEQS FROM CAMEO CIFs
import os
from openfold.data.mmcif_parsing import parse

cameo_chain_dict = {}
cameo_mmcif_dir = '../validation_set_cameo/cameo/mmcif_files'
for filename in os.listdir(cameo_mmcif_dir):
    try:
        mmcif_file = os.path.join(cameo_mmcif_dir, filename)
        prot = os.path.splitext(filename)[0]

        with open(mmcif_file, 'r') as mmcif:
            mmcif_parse = parse(file_id=filename, mmcif_string=mmcif.read())
            mmcif_obj = mmcif_parse.mmcif_object
            if(mmcif_obj is None):
                print(f'Failed to parse {filename}')
                continue

            for chain_id, seq in mmcif_obj.chain_to_seqres.items():
                item = {prot+'_'+chain_id: seq}
                cameo_chain_dict.update(item)
    except Exception as e:
        print(f'Failed to parse {filename}')
        raise e

print(len(cameo_chain_dict))
# Move every file into its own directory
import shutil
cameo_esm_dir = '../validation_set_camep/cameo/cameo_esm_embeddings'
for filename in os.listdir(cameo_esm_dir):
    chain, ext = os.path.splitext(filename)
    os.makedirs(os.path.join(cameo_esm_dir, chain))
    source_file = os.path.join(cameo_esm_dir, filename)
    destination_file = os.path.join(cameo_esm_dir, chain, filename)
    shutil.move(source_file, destination_file)
# Move the ESM embeddings to their respective folders with the templates
cameo_alignment_dir = '../validation_set_cameo/cameo_esm_emeddings'
for filename in os.listdir(cameo_alignment_dir):
    try:
        source_file = os.path.join(cameo_esm_dir, filename, filename+'.pt')
        destination_file = os.path.join(cameo_alignment_dir, filename, filename+'.pt')
        shutil.copy(source_file, destination_file)
    except Exception as e:
        print(f"Exception occured for {filename}")
# Move the ESM embeddings and the template hhr files to new dirs
import os
import shutil
cameo_alignment_dir = '../validation_set_cameo/cameo_esm_emeddings'
cameo_sseq_alignment_dir = '../validation_set_cameo/cameo_sseq_alignments'
for filename in os.listdir(cameo_alignment_dir):
    try:
        source_file = os.path.join(cameo_alignment_dir, filename, filename+'.pt')
        destination_file = os.path.join(cameo_alignment_dir, filename, filename+'.pt')
        shutil.copy(source_file, destination_file)
        source_file = os.path.join(cameo_alignment_dir, filename, 'pdb70_hits.hhr' )
        destination_file = os.path.join(cameo_sseq_alignment_dir, filename, 'pdb70_hits.hhr')
    except Exception as e:
        print(f"Exception occured for {filename}")
# Take out pdb70_hits.hhr files from Flatiron format DBs
# alignment folder
import os
import json

os.chdir('/fsx/distillation_alignments/template_dbs/')
with open('super.index','r') as file:
    super = json.loads(file.read())

for key, value in super.items():
    db = value['db']
    files_list = value['files']
    for file_block in files_list:
        if file_block[0] == 'pdb70_hits.hhr':
            template_file_start = file_block[1]
            template_file_length = file_block[2]
        os.mkdir(os.path.join('alignment_files',key))
        file_name_new = os.path.join('alignment_files', key, 'pdb70_hits.hhr')
        with open(db, 'r') as db_file:
            db_file.seek(template_file_start)
            file_data = db_file.read(template_file_length)
        with open(file_name_new,'w') as template_file:
            template_file.write(file_data)

# Move distillation ESM embeddings to appropriate alignment folders
import os
import shutil
distillation_esm_dir = '../distillation_chunk'
distillation_sseq_alignment_dir = '../distillation_alignments/template_dbs/alignment_files'
for filename in os.listdir(distillation_sseq_alignment_dir):
    try:
        source_file = os.path.join(distillation_esm_dir, filename, filename+'.pt')
        destination_file = os.path.join(cameo_alignment_dir, filename, filename+'.pt')
        shutil.copy(source_file, destination_file)
        source_file = os.path.join(cameo_alignment_dir, filename, 'pdb70_hits.hhr' )
        destination_file = os.path.join(cameo_sseq_alignment_dir, filename, 'pdb70_hits.hhr')
    except Exception as e:
        print(f"Exception occured for {filename}")

# Move PDB ESM embeddings to appropriate alignment folders
import os
import shutil
import progressbar
pdb_esm_dir = '/fsx/pdb_esm_embeddings/'
pdb_sseq_alignment_dir = '/fsx/pdb_alignment_dbs/alignment_dbs/alignment_files/'
widgets = [' [',
         progressbar.Timer(format= 'elapsed time: %(elapsed)s'),
         '] ',
           progressbar.Bar('*'),' (',
           progressbar.ETA(), ') ',
          ]
files_in_pdb_sseq_alignment_dir = os.listdir(pdb_sseq_alignment_dir)
bar = progressbar.ProgressBar(max_value=len(files_in_pdb_sseq_alignment_dir),
                              widgets=widgets).start()
i=0
for filename in os.listdir(pdb_sseq_alignment_dir):
    try:
        for chunk_index in range(7):
            if os.path.exists(os.path.join(pdb_esm_dir, 'chunk' + str(chunk_index), filename + '.pt')):
                source_file = os.path.join(pdb_esm_dir, filename, filename+'.pt')
                break
        else:
            raise FileNotFoundError()
        destination_file = os.path.join(pdb_sseq_alignment_dir, filename, filename+'.pt')
        shutil.copy(source_file, destination_file)
        bar.update(i)
        i += 1
    except Exception as e:
        print(f"Exception occured for {filename}")