# Adapted from https://www.github.com/jproney/AF2Rank/blob/master/test_templates.py

# Copyright 2024 AlQuraishi Laboratory
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import os
import argparse
import traceback
import jax
import jax.numpy as jnp
import numpy as np
import re
import subprocess
import torch
from collections import namedtuple
from copy import deepcopy

parser = argparse.ArgumentParser()
parser.add_argument("name", help="name to save everything under")
parser.add_argument("--target_list", nargs='*', help="List of target names to run")
parser.add_argument("--targets_file", default="", help="File with list of target names to run")
parser.add_argument("--recycles", type=int, default=1, help="Number of recycles when predicting")
parser.add_argument("--model_name", type=str, default="model_1_ptm", help="Which OF model to use")
parser.add_argument("--seed", type=int, default=0, help="RNG Seed")
parser.add_argument("--verbose", action='store_true', help="print extra")
parser.add_argument("--deterministic", action='store_true', help="make all data processing deterministic (no masking, etc.)")
parser.add_argument("--use_native", action='store_true', help="add the native structure as a decoy, and compare outputs against it")
parser.add_argument("--mask_sidechains", action='store_true', help="mask out sidechain atoms except for C-Beta")
parser.add_argument("--mask_sidechains_add_cb", action='store_true', help="mask out sidechain atoms except for C-Beta, and add C-Beta to glycines")
parser.add_argument("--seq_replacement", default='', help="Amino acid residue to fill the decoy sequence with. Default keeps target sequence")
parser.add_argument("--of_dir", default="/home/user/openfold/", help="OpenFold code and weights directory")
parser.add_argument("--esm_dir", help="ESM1b embeddings directory, containing embeddings as *.pt")
parser.add_argument("--decoy_dir", default="/home/user/openfold/rosetta_decoy_set/", help="Rosetta decoy directory")
parser.add_argument("--output_dir", default="/home/user/ofss_ranking_experiment/outputs/", help="Rosetta decoy directory")
parser.add_argument("--openfold_checkpoint_path", help="Path to the OpenFold model checkpoint")
parser.add_argument("--jax_param_path", help="Path to the JAX parameters checkpoint")
parser.add_argument("--model_device", default="cpu", help="Device to run the model on")
parser.add_argument("--tm_exec", default="/home/user/tmscore/TMscore", help="TMScore executable")

args = parser.parse_args()

sys.path.insert(0, args.of_dir)


# openfold imports
from openfold import config

from openfold.data import data_pipeline
from openfold.data import feature_pipeline

from openfold.np import protein
from openfold.np import residue_constants

from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.script_utils import load_models_from_command_line, run_model


# helper functions

"""
Read in a PDB file from a path
"""
def pdb_to_string(pdb_file):
  lines = []
  for line in open(pdb_file,"r"):
    if line[:6] == "HETATM" and line[17:20] == "MSE":
      line = "ATOM  "+line[6:17]+"MET"+line[20:]
    if line[:4] == "ATOM":
      lines.append(line)
  return "".join(lines)

"""
Compute aligned RMSD between two corresponding sets of points
true -- set of reference points. Numpy array of dimension N x 3
pred -- set of predicted points, Numpy array of dimension N x 3
"""
def jnp_rmsd(true, pred):
  def kabsch(P, Q):
    V, S, W = jnp.linalg.svd(P.T @ Q, full_matrices=False)
    flip = jax.nn.sigmoid(-10 * jnp.linalg.det(V) * jnp.linalg.det(W))
    S = flip * S.at[-1].set(-S[-1]) + (1-flip) * S
    V = flip * V.at[:,-1].set(-V[:,-1]) + (1-flip) * V
    return V@W
  p = true - true.mean(0,keepdims=True)
  q = pred - pred.mean(0,keepdims=True)
  p = p @ kabsch(p,q)
  loss = jnp.sqrt(jnp.square(p-q).sum(-1).mean() + 1e-8)
  return float(loss)

"""
Create an OpenFold model runner
name -- The name of the model to get the parameters from. Options: model_[1-5]
"""
def make_model_runner(name, recycles, args):
  cfg = config.model_config(name)      

  cfg.data.common.max_recycling_iters = recycles

  if args.deterministic:
    cfg.data.eval.masked_msa_replace_fraction = 0.0
    cfg.data.predict.masked_msa_replace_fraction = 0.0
  
  model_generator = load_models_from_command_line(cfg, args.model_device, args.openfold_checkpoint_path, args.jax_param_path, args.output_dir)
  model, _ = model_generator.__next__()

  return model, cfg

"""
Make a set of empty features for no-template evaluations
"""
def empty_placeholder_template_features(num_templates, num_res):
  return {
      'template_aatype': np.zeros((num_templates, num_res), dtype=np.int64),
      'template_all_atom_mask': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num),
          dtype=np.float32),
      'template_all_atom_positions': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num, 3),
          dtype=np.float32),
      'template_domain_names': np.zeros([num_templates], dtype=object),
      'template_sequence': np.zeros([num_templates], dtype=object),
      'template_sum_probs': np.zeros([num_templates, 1], dtype=np.float32),
  }

def make_embedding_features(args, label):
  seqemb_features = {}

  path = os.path.join(args.esm_dir, label+'.pt')

  # Load embedding file
  seqemb_data = torch.load(path)
  seqemb_features["seq_embedding"] = seqemb_data["representations"][33]

  return seqemb_features

"""
Create a feature dictionary for input to OpenFold
runner - The model runner being invoked. Returned from `make_model_runner`
sequence - The target sequence being predicted
templates - The template features being added to the inputs
seed - The random seed being used for data processing
"""
def make_processed_feature_dict(cfg, sequence, name="test", templates=None, seed=0):
  feature_dict = {}
  feature_dict.update(data_pipeline.make_sequence_features(sequence, name, len(sequence)))

  msa = [[sequence]]
  deletion_matrix = [[[0 for _ in sequence]]]

  feature_dict.update(data_pipeline.make_msa_features(msa, deletion_matrix))

  if templates is not None:
    feature_dict.update(templates)
  else:
    feature_dict.update(empty_placeholder_template_features(num_templates=0, num_res=len(sequence)))

  feature_dict.update(make_embedding_features(args, name.split('_')[0]))
  
  feature_processor = feature_pipeline.FeaturePipeline(cfg.data)
  processed_feature_dict = feature_processor.process_features(feature_dict, mode='predict')
  processed_feature_dict = {
    k: torch.as_tensor(v, device=args.model_device)
    for k, v in processed_feature_dict.items()
  }

  return processed_feature_dict

"""
Package OpenFold's output into an easy-to-use dictionary
prediction_result - output from running OpenFold on an input dictionary
processed_feature_dict -- The dictionary passed to OpenFold as input. Returned by `make_processed_feature_dict`.
"""
def parse_results(prediction_result, processed_feature_dict):
  b_factors = prediction_result['plddt'][:,None] * prediction_result['final_atom_mask']

  out = {"unrelaxed_protein": protein.from_prediction(processed_feature_dict, prediction_result, b_factors=b_factors),
        "plddt": prediction_result['plddt'],
        "pLDDT": prediction_result['plddt'].mean(),}
        
  out.update({"pTMscore": prediction_result['predicted_tm_score']})

  return out


'''
Function used to add C-Beta to glycine resides
input:  3 coords (a,b,c), (L)ength, (A)ngle, and (D)ihedral
output: 4th coord
'''
def extend(a,b,c, L,A,D):
  N = lambda x: x/np.sqrt(np.square(x).sum(-1,keepdims=True) + 1e-8)
  bc = N(b-c)
  n = N(np.cross(b-a, bc))
  m = [bc,np.cross(n,bc),n]
  d = [L*np.cos(A), L*np.sin(A)*np.cos(D), -L*np.sin(A)*np.sin(D)]
  return c + sum([m*d for m,d in zip(m,d)])

"""
Ingest a decoy protein, pass it to OpenFold as a template, and extract the parsed output
target_seq -- the sequence to be predicted
decoy_prot -- the decoy structure to be injected as a template
model_runner -- the model runner to execute
name -- the name associated with this prediction
"""
def score_decoy(target_seq, decoy_prot, model_runner, name):
  decoy_seq_in = "".join([residue_constants.restypes[x] for x in decoy_prot.aatype]) # the sequence in the decoy PDB file

  mismatch = False
  if decoy_seq_in == target_seq:
    assert jnp.all(decoy_prot.residue_index - 1 == np.arange(len(target_seq)))
  else: # case when template is missing some residues
    if args.verbose:
      print("Sequence mismatch: {}".format(name))
    mismatch=True

    assert "".join(target_seq[i-1] for i in decoy_prot.residue_index) == decoy_seq_in 
  
  # use this to index into the template features
  template_idxs = decoy_prot.residue_index-1
  template_idx_set = set(template_idxs)

  # The sequence associated with the decoy. Always has same length as target sequence.
  decoy_seq = args.seq_replacement*len(target_seq) if len(args.seq_replacement) == 1 else target_seq

  # create empty template features
  pos = np.zeros([1,len(decoy_seq), 37, 3])
  atom_mask = np.zeros([1, len(decoy_seq), 37])

  if args.mask_sidechains_add_cb:
    pos[0, template_idxs, :5] = decoy_prot.atom_positions[:,:5]

    # residues where we have all of the key backbone atoms (N CA C)
    backbone_modelled = np.asarray(jnp.all(decoy_prot.atom_mask[:,[0,1,2]] == 1, axis=1))
    backbone_idx_set = set(decoy_prot.residue_index[backbone_modelled] - 1)

    projected_cb = [i-1 for i,b,m in zip(decoy_prot.residue_index, backbone_modelled, decoy_prot.atom_mask) if m[3] == 0 and b]
    projected_cb_set = set(projected_cb)
    gly_idx = [i for i,a in enumerate(target_seq) if a == "G"]
    assert all([k in projected_cb_set for k in gly_idx if k in template_idx_set and k in backbone_idx_set]) # make sure we are adding CBs to all of the glycines

    cbs = np.array([extend(c,n,ca, 1.522, 1.927, -2.143) for c, n ,ca in zip(pos[0,:,2], pos[0,:,0], pos[0,:,1])])

    pos[0, projected_cb, 3] = cbs[projected_cb]
    atom_mask[0, template_idxs, :5] = decoy_prot.atom_mask[:, :5]
    atom_mask[0, projected_cb, 3] = 1

    template = {"template_aatype":residue_constants.sequence_to_onehot(decoy_seq, residue_constants.HHBLITS_AA_TO_ID)[None],
                "template_all_atom_mask": atom_mask.astype(np.float32),
                "template_all_atom_positions":pos.astype(np.float32),
                "template_domain_names":np.asarray(["None"])}
  elif args.mask_sidechains:
    pos[0, template_idxs, :5] = decoy_prot.atom_positions[:,:5]
    atom_mask[0, template_idxs, :5] = decoy_prot.atom_mask[:,:5]

    template = {"template_aatype":residue_constants.sequence_to_onehot(decoy_seq, residue_constants.HHBLITS_AA_TO_ID)[None],
                "template_all_atom_mask": atom_mask.astype(np.float32),
                "template_all_atom_positions": pos.astype(np.float32),
                "template_domain_names":np.asarray(["None"])}
  else:
    pos[0, template_idxs] = decoy_prot.atom_positions
    atom_mask[0, template_idxs] = decoy_prot.atom_mask

    template = {"template_aatype":residue_constants.sequence_to_onehot(decoy_seq, residue_constants.HHBLITS_AA_TO_ID)[None],
                "template_all_atom_mask":decoy_prot.atom_mask[None].astype(np.float32),
                "template_all_atom_positions":decoy_prot.atom_positions[None].astype(np.float32),
                "template_domain_names":np.asarray(["None"])}
                                    
  features = make_processed_feature_dict(cfg, target_seq, name=name, templates=template, seed=args.seed)
  #with open(os.path.join(args.output_dir, name + '_features.pt'), 'wb') as outfile:
  #  torch.save(features, outfile)
  working_batch = deepcopy(features)
  out, inference_time = run_model(model_runner, working_batch, name, args.output_dir)
  print(f"{name} done. Inference time: ", inference_time)
  working_batch = tensor_tree_map(lambda x: np.array(x[..., -1].cpu()), working_batch)
  out = tensor_tree_map(lambda x: np.array(x.cpu()), out)
  result = parse_results(out, working_batch)
  return result, mismatch


tm_re = re.compile(r'TM-score[\s]*=[\s]*(\d.\d+)')
ref_len_re = re.compile(r'Length=[\s]*(\d+)[\s]*\(by which all scores are normalized\)')
common_re = re.compile(r'Number of residues in common=[\s]*(\d+)')
super_re = re.compile(r'\(":" denotes the residue pairs of distance < 5\.0 Angstrom\)\\n([A-Z\-]+)\\n[" ", :]+\\n([A-Z\-]+)\\n')

"""
Compute TM Scores between two PDBs and parse outputs
pdb_pred -- The path to the predicted PDB
pdb_native -- The path to the native PDB
test_len -- run asserts that the input and output should have the same length
"""
def compute_tmscore(pdb_pred, pdb_native, test_len=True):
  cmd = ([args.tm_exec, pdb_pred, pdb_native])
  tmscore_output = str(subprocess.check_output(cmd))
  try:
    tm_out = float(tm_re.search(tmscore_output).group(1))
    reflen = int(ref_len_re.search(tmscore_output).group(1))
    common = int(common_re.search(tmscore_output).group(1))
    
    seq1 = super_re.search(tmscore_output).group(1)
    seq2 = super_re.search(tmscore_output).group(1)
  except Exception as e:
    print("Failed on: " + " ".join(cmd))
    raise e

  if test_len:
    assert reflen == common, cmd
    assert seq1 == seq2, cmd
    assert len(seq1) == reflen, cmd

  return tm_out

# Simple wrapper for keeping track of the information associated with each decoy. 
decoy_fields_list = ['target', 'decoy_id', 'decoy_path', 'rmsd', 'rosettascore', 'gdt_ts', 'tmscore', 'danscore']
Decoy = namedtuple("Decoy", decoy_fields_list)


# headers for csv outputs
csv_headers = decoy_fields_list + ['output_path', 'rmsd_out', 'tm_diff', 'tm_out', 'plddt', 'ptm']

def write_results(decoy, af_result, prot_native=None, pdb_native=None, mismatch=False):
  plddt = float(af_result['pLDDT'])
  if "pTMscore" not in af_result:
    ptm = -1
  else:
    ptm = float(af_result["pTMscore"])
  if prot_native is None:
    rms_out = -1
  else:
    rms_out = jnp_rmsd(jnp.asarray(prot_native.atom_positions[:,1,:]), jnp.asarray(af_result['unrelaxed_protein'].atom_positions[:,1,:]))

  pdb_lines = protein.to_pdb(af_result["unrelaxed_protein"])
  pdb_out_path = args.output_dir + args.name + "/pdbs/" + decoy.target + "_" + decoy.decoy_id
  with open(pdb_out_path, 'w') as f:
    f.write(pdb_lines)

  if decoy.decoy_id != "none.pdb":
    tm_diff = compute_tmscore(decoy.decoy_path, pdb_out_path, test_len = not mismatch)
  else:
    tm_diff = -1

  if pdb_native is None:
    tm_out = -1
  else:
    tm_out = compute_tmscore(pdb_out_path, pdb_native)

  if not os.path.exists(args.output_dir + args.name + "/results/results_{}.csv".format(decoy.target)):
    with open(args.output_dir + args.name + "/results/results_{}.csv".format(decoy.target), "w") as f:
      f.write(",".join(csv_headers) + "\n")


  with open(args.output_dir + args.name + "/results/results_{}.csv".format(decoy.target), "a") as f:
    result_fields = [str(x) for x in list(decoy) + [pdb_out_path, rms_out, tm_diff, tm_out, plddt, ptm]]
    f.write(",".join(result_fields) + "\n")

    if args.verbose:
      print(",".join([x + "=" + y for x,y in zip(csv_headers, result_fields)]))


# create all of the output directories
os.makedirs(args.output_dir + args.name, exist_ok=True)
os.makedirs(args.output_dir + args.name + "/pdbs", exist_ok=True)
os.makedirs(args.output_dir + args.name + "/results", exist_ok=True)

if len(args.targets_file) > 0:
  natives_list = open(args.targets_file, 'r').read().split("\n")[:-1]
else:
  natives_list = args.target_list


finished_decoys = []
for n in natives_list:
  if os.path.exists(args.output_dir  + args.name + "/results/results_{}.csv".format(n)):
    finished_decoys += [x.split(",")[0] + "_" + x.split(",")[1] for x in open(args.output_dir  + args.name + "/results/results_{}.csv".format(n), "r").readlines()]
finished_decoys = set(finished_decoys)


if os.path.exists(args.output_dir  + args.name + "/finished_targets.txt"):
  finished_targets = set(open(args.output_dir + args.name + "/finished_targets.txt", 'r').read().split("\n")[:-1])
else:
  finished_targets = []


# info of the form "target decoy_id"
decoy_list = [x.split() for x in open(args.decoy_dir + "decoy_list.txt", 'r').read().split("\n")[:-1]] 

# parse all of the information about the decoys
decoy_data = {}
for field in decoy_fields_list[2:]:
  if os.path.exists(args.decoy_dir + field + ".txt"):
    lines = [x.split() for x in open(args.decoy_dir + field + ".txt", 'r').read().split("\n")[:-1]] # form "target decoy_id metric value"

    # make sure everything is in the same order
    for i,l in enumerate(lines):
      assert l[0] == decoy_list[i][0]
      assert l[1] == decoy_list[i][1]

    decoy_data[field] = [l[-1] for l in lines]
  else:
    decoy_data[field] = [-1]*len(decoy_list) # -1 as a placeholder

decoy_dict = {n : [] for n in natives_list if n not in finished_targets} # key = target name, value = list of Decoy objects

for i, d in enumerate(decoy_list):

  decoy = Decoy(target=d[0], decoy_id=d[1], decoy_path=args.decoy_dir + "decoys/" + d[0] + "/" + d[1], 
                  rmsd = decoy_data["rmsd"][i], rosettascore = decoy_data["rosettascore"][i], gdt_ts = decoy_data["gdt_ts"][i], 
                    tmscore=decoy_data["tmscore"][i], danscore = decoy_data["danscore"][i])

  if decoy.target in decoy_dict and decoy.target + "_" + decoy.decoy_id not in finished_decoys:
    decoy_dict[decoy.target].append(decoy)

# add another decoy entry for the native structure
if args.use_native:
  for n in decoy_dict.keys():
    if n + "_native" not in finished_decoys:
      decoy_dict[n].insert(0, Decoy(target=n, decoy_id="native.pdb", decoy_path=args.decoy_dir + "natives/" + n + ".pdb", 
                            rmsd = 0, rosettascore = -1, gdt_ts = 1, tmscore = 1, danscore = -1))

if args.verbose:
  print(finished_decoys)

model_name = args.model_name
results_key = model_name + "_seed_{}".format(args.seed)
for n in natives_list:
  try:
    pdb_native = args.decoy_dir + "natives/" + n + ".pdb"
    prot_native = protein.from_pdb_string(pdb_to_string(pdb_native))
    seq_native = "".join([residue_constants.restypes[x] for x in prot_native.aatype])
    runner, cfg = make_model_runner(model_name, args.recycles, args)

    if n + "_none.pdb" not in finished_decoys:

      # run the model with no templates
      features = make_processed_feature_dict(cfg, seq_native, name=n + "_none", seed=args.seed)
      working_batch = deepcopy(features)
      out, inference_time = run_model(runner, working_batch, n + "_none", args.output_dir)
      print(f"{n}_none done. Inference time: ", inference_time)
      working_batch = tensor_tree_map(lambda x: np.array(x[..., -1].cpu()), working_batch)
      out = tensor_tree_map(lambda x: np.array(x.cpu()), out)
      result = parse_results(out, working_batch)

      dummy_decoy = Decoy(target=n, decoy_id="none.pdb", decoy_path="_", rmsd=-1, rosettascore=-1, gdt_ts=-1, tmscore=-1,danscore=-1)
      write_results(dummy_decoy, result, prot_native=prot_native if args.use_native else None, pdb_native=pdb_native if args.use_native else None)


    # run the model with all of the decoys passed as templates
    for d in decoy_dict[n]:
      prot = protein.from_pdb_string(pdb_to_string(d.decoy_path))
      result, mismatch = score_decoy(seq_native, prot, runner, d.target + "_" + d.decoy_id)
      write_results(d, result, prot_native=prot_native if args.use_native else None, pdb_native=pdb_native if args.use_native else None, mismatch=mismatch)


    with open(args.output_dir + args.name + "/finished_targets.txt", 'a') as f:
      f.write(n + "\n")
  except AssertionError as ae:
    print(f"AssertionError encountered while processing a decoy of native {n}")
    traceback.print_exc()
  except Exception as e:
    print(f"Exception encountered while processing a decoy of native {n}")
    traceback.print_exc()
