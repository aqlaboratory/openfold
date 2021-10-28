# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
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

import os
import datetime
from typing import Mapping, Optional, Sequence, Any

import numpy as np

from openfold.data import templates, parsers, mmcif_parsing
from openfold.data.tools import jackhmmer, hhblits, hhsearch
from openfold.data.tools.utils import to_date 
from openfold.np import residue_constants, protein


FeatureDict = Mapping[str, np.ndarray]

def empty_template_feats(n_res) -> FeatureDict:
    return {
        "template_aatype": np.zeros((0, n_res)).astype(np.int64),
        "template_all_atom_positions": 
            np.zeros((0, n_res, 37, 3)).astype(np.float32),
        "template_sum_probs": np.zeros((0, 1)).astype(np.float32),
        "template_all_atom_mask": np.zeros((0, n_res, 37)).astype(np.float32),
    }


def make_sequence_features(
    sequence: str, description: str, num_res: int
) -> FeatureDict:
    """Construct a feature dict of sequence features."""
    features = {}
    features["aatype"] = residue_constants.sequence_to_onehot(
        sequence=sequence,
        mapping=residue_constants.restype_order_with_x,
        map_unknown_to_x=True,
    )
    features["between_segment_residues"] = np.zeros((num_res,), dtype=np.int32)
    features["domain_name"] = np.array(
        [description.encode("utf-8")], dtype=np.object_
    )
    features["residue_index"] = np.array(range(num_res), dtype=np.int32)
    features["seq_length"] = np.array([num_res] * num_res, dtype=np.int32)
    features["sequence"] = np.array(
        [sequence.encode("utf-8")], dtype=np.object_
    )
    return features


def make_mmcif_features(
    mmcif_object: mmcif_parsing.MmcifObject, chain_id: str
) -> FeatureDict:
    input_sequence = mmcif_object.chain_to_seqres[chain_id]
    description = "_".join([mmcif_object.file_id, chain_id])
    num_res = len(input_sequence)

    mmcif_feats = {}

    mmcif_feats.update(
        make_sequence_features(
            sequence=input_sequence,
            description=description,
            num_res=num_res,
        )
    )

    all_atom_positions, all_atom_mask = mmcif_parsing.get_atom_coords(
        mmcif_object=mmcif_object, chain_id=chain_id
    )
    mmcif_feats["all_atom_positions"] = all_atom_positions
    mmcif_feats["all_atom_mask"] = all_atom_mask

    mmcif_feats["resolution"] = np.array(
        [mmcif_object.header["resolution"]], dtype=np.float32
    )

    mmcif_feats["release_date"] = np.array(
        [mmcif_object.header["release_date"].encode("utf-8")], dtype=np.object_
    )

    mmcif_feats["is_distillation"] = np.array(0., dtype=np.float32)

    return mmcif_feats


def make_pdb_features(
        protein_object: protein.Protein, 
        description: str, 
        confidence_threshold: float = 0.5,
) -> FeatureDict:
    pdb_feats = {}

    pdb_feats.update(
        make_sequence_features(
            sequence=protein_object.aatype,
            description=description,
            num_res=len(protein_object.aatype),
        )
    )

    all_atom_positions = protein_object.atom_positions
    all_atom_mask = protein_object.atom_mask

    high_confidence = protein.b_factors > confidence_threshold
    high_confidence = np.any(high_confidence, axis=-1)
    for i, confident in enumerate(high_confidence):
        if(not confident):
            all_atom_mask[i] = 0

    pdb_feats["all_atom_positions"] = all_atom_positions
    pdb_feats["all_atom_mask"] = all_atom_mask

    pdb_feats["resolution"] = np.array([0.]).astype(np.float32)
    pdb_feats["is_distillation"] = np.array(1.).astype(np.float32)

    return pdb_feats


def make_msa_features(
    msas: Sequence[Sequence[str]],
    deletion_matrices: Sequence[parsers.DeletionMatrix],
) -> FeatureDict:
    """Constructs a feature dict of MSA features."""
    if not msas:
        raise ValueError("At least one MSA must be provided.")

    int_msa = []
    deletion_matrix = []
    seen_sequences = set()
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(
                f"MSA {msa_index} must contain at least one sequence."
            )
        for sequence_index, sequence in enumerate(msa):
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            int_msa.append(
                [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence]
            )
            deletion_matrix.append(deletion_matrices[msa_index][sequence_index])

    num_res = len(msas[0][0])
    num_alignments = len(int_msa)
    features = {}
    features["deletion_matrix_int"] = np.array(deletion_matrix, dtype=np.int32)
    features["msa"] = np.array(int_msa, dtype=np.int32)
    features["num_alignments"] = np.array(
        [num_alignments] * num_res, dtype=np.int32
    )
    return features


class AlignmentRunner:
    """Runs alignment tools and saves the results"""

    def __init__(
        self,
        jackhmmer_binary_path: str,
        hhblits_binary_path: str,
        hhsearch_binary_path: str,
        uniref90_database_path: str,
        mgnify_database_path: str,
        bfd_database_path: Optional[str],
        uniclust30_database_path: Optional[str],
        small_bfd_database_path: Optional[str],
        pdb70_database_path: str,
        use_small_bfd: bool,
        no_cpus: int,
        uniref_max_hits: int = 10000,
        mgnify_max_hits: int = 5000,
    ):
        self._use_small_bfd = use_small_bfd
        self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=uniref90_database_path,
            n_cpu=no_cpus,
        )

        if use_small_bfd:
            self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=small_bfd_database_path,
                n_cpu=no_cpus,
            )
        else:
            self.hhblits_bfd_uniclust_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path, uniclust30_database_path],
                n_cpu=no_cpus,
            )

        self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=mgnify_database_path,
            n_cpu=no_cpus,
        )

        self.hhsearch_pdb70_runner = hhsearch.HHSearch(
            binary_path=hhsearch_binary_path,
            databases=[pdb70_database_path],
            n_cpu=no_cpus,
        )
        self.uniref_max_hits = uniref_max_hits
        self.mgnify_max_hits = mgnify_max_hits

    def run(
        self,
        fasta_path: str,
        output_dir: str,
    ):
        """Runs alignment tools on a sequence"""
        jackhmmer_uniref90_result = self.jackhmmer_uniref90_runner.query(
            fasta_path
        )[0]
        uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(
            jackhmmer_uniref90_result["sto"], max_sequences=self.uniref_max_hits
        )
        uniref90_out_path = os.path.join(output_dir, "uniref90_hits.a3m")
        with open(uniref90_out_path, "w") as f:
            f.write(uniref90_msa_as_a3m)

        jackhmmer_mgnify_result = self.jackhmmer_mgnify_runner.query(
            fasta_path
        )[0]
        mgnify_msa_as_a3m = parsers.convert_stockholm_to_a3m(
            jackhmmer_mgnify_result["sto"], max_sequences=self.mgnify_max_hits
        )
        mgnify_out_path = os.path.join(output_dir, "mgnify_hits.a3m")
        with open(mgnify_out_path, "w") as f:
            f.write(mgnify_msa_as_a3m)

        hhsearch_result = self.hhsearch_pdb70_runner.query(uniref90_msa_as_a3m)
        pdb70_out_path = os.path.join(output_dir, "pdb70_hits.hhr")
        with open(pdb70_out_path, "w") as f:
            f.write(hhsearch_result)

        if self._use_small_bfd:
            jackhmmer_small_bfd_result = self.jackhmmer_small_bfd_runner.query(
                fasta_path
            )[0]
            bfd_out_path = os.path.join(output_dir, "small_bfd_hits.sto")
            with open(bfd_out_path, "w") as f:
                f.write(jackhmmer_small_bfd_result["sto"])
        else:
            hhblits_bfd_uniclust_result = (
                self.hhblits_bfd_uniclust_runner.query(fasta_path)
            )
            if output_dir is not None:
                bfd_out_path = os.path.join(output_dir, "bfd_uniclust_hits.a3m")
                with open(bfd_out_path, "w") as f:
                    f.write(hhblits_bfd_uniclust_result["a3m"])


class DataPipeline:
    """Assembles input features."""

    def __init__(
        self,
        template_featurizer: templates.TemplateHitFeaturizer,
    ):
        self.template_featurizer = template_featurizer

    def _parse_msa_data(
        self,
        alignment_dir: str,
    ) -> Mapping[str, Any]:
        msa_data = {}
        for f in os.listdir(alignment_dir):
            path = os.path.join(alignment_dir, f)
            ext = os.path.splitext(f)[-1]

            if(ext == ".a3m"):
                with open(path, "r") as fp:
                    msa, deletion_matrix = parsers.parse_a3m(fp.read())
                data = {"msa": msa, "deletion_matrix": deletion_matrix}
            elif(ext == ".sto"):
                with open(path, "r") as fp:
                    msa, deletion_matrix, _ = parsers.parse_stockholm(
                        fp.read()
                    )
                data = {"msa": msa, "deletion_matrix": deletion_matrix}
            else:
                continue
            
            msa_data[f] = data

        return msa_data

    def _parse_template_hits(
        self,
        alignment_dir: str,
    ) -> Mapping[str, Any]:
        all_hits = {}
        for f in os.listdir(alignment_dir):
            path = os.path.join(alignment_dir, f)
            ext = os.path.splitext(f)[-1]

            if(ext == ".hhr"):
                with open(path, "r") as fp:
                    hits = parsers.parse_hhr(fp.read())
                all_hits[f] = hits

        return all_hits

    def _process_msa_feats(
        self,
        alignment_dir: str,
    ) -> Mapping[str, Any]:
        msa_data = self._parse_msa_data(alignment_dir)
        msas, deletion_matrices = zip(*[
            (v["msa"], v["deletion_matrix"]) for v in msa_data.values()
        ])
        msa_features = make_msa_features(
            msas=msas,
            deletion_matrices=deletion_matrices,
        )

        return msa_features

    def process_fasta(
        self,
        fasta_path: str,
        alignment_dir: str,
    ) -> FeatureDict:
        """Assembles features for a single sequence in a FASTA file"""
        with open(fasta_path) as f:
            fasta_str = f.read()
        input_seqs, input_descs = parsers.parse_fasta(fasta_str)
        if len(input_seqs) != 1:
            raise ValueError(
                f"More than one input sequence found in {fasta_path}."
            )
        input_sequence = input_seqs[0]
        input_description = input_descs[0]
        num_res = len(input_sequence)

        hits = self._parse_template_hits(alignment_dir)
        hits_cat = sum(hits.values(), [])
        if(len(hits_cat) == 0):
            template_features = empty_template_feats(len(input_sequence))
        else:
            templates_result = self.template_featurizer.get_templates(
                query_sequence=input_sequence,
                query_pdb_code=None,
                query_release_date=None,
                hits=hits_cat,
            )
            template_features = templates_result.features

            # The template featurizer doesn't format empty template features
            # properly. This is a quick fix.
            if(template_features["template_aatype"].shape[0] == 0):
                template_features = empty_template_feats(len(input_sequence))

        sequence_features = make_sequence_features(
            sequence=input_sequence,
            description=input_description,
            num_res=num_res,
        )

        msa_features = self._process_msa_feats(alignment_dir)
        
        return {
            **sequence_features,
            **msa_features, 
            **template_features
        }

    def process_mmcif(
        self,
        mmcif: mmcif_parsing.MmcifObject,  # parsing is expensive, so no path
        alignment_dir: str,
        chain_id: Optional[str] = None,
    ) -> FeatureDict:
        """
            Assembles features for a specific chain in an mmCIF object.

            If chain_id is None, it is assumed that there is only one chain
            in the object. Otherwise, a ValueError is thrown.
        """
        if chain_id is None:
            chains = mmcif.structure.get_chains()
            chain = next(chains, None)
            if chain is None:
                raise ValueError("No chains in mmCIF file")
            chain_id = chain.id

        mmcif_feats = make_mmcif_features(mmcif, chain_id)

        input_sequence = mmcif.chain_to_seqres[chain_id]
        hits = self._parse_template_hits(alignment_dir)
        hits_cat = sum(hits.values(), [])
        print(len(hits_cat))
        if(len(hits_cat) == 0):
            template_features = empty_template_feats(len(input_sequence))
        else:
            templates_result = self.template_featurizer.get_templates(
                query_sequence=input_sequence,
                query_pdb_code=None,
                query_release_date=to_date(mmcif.header["release_date"]),
                hits=hits_cat,
            )
            template_features = templates_result.features

            # The template featurizer doesn't format empty template features
            # properly. This is a quick fix.
            if(template_features["template_aatype"].shape[0] == 0):
                template_features = empty_template_feats(len(input_sequence))


        msa_features = self._process_msa_feats(alignment_dir)

        return {**mmcif_feats, **template_features, **msa_features}

    def process_pdb(
        self,
        pdb_path: str,
        alignment_dir: str,
    ) -> FeatureDict:
        """
            Assembles features for a protein in a PDB file.
        """
        with open(pdb_path, 'r') as f:
            pdb_str = pdb_path

        protein_object = protein.from_pdb_string(pdb_str)
        input_sequence = protein_object.aatype 

        pdb_feats = make_pdb_features(protein_object)

        hits = self._parse_template_hits(alignment_dir)
        hits_cat = sum(hits.values(), [])
        if(len(hits_cat) == 0):
            template_features = empty_template_feats(len(input_sequence))
        else:
            templates_result = self.template_featurizer.get_templates(
                query_sequence=input_sequence,
                query_pdb_code=None,
                query_release_date=None,
                hits=hits_cat,
            )
            template_features = templates_result.features

            # The template featurizer doesn't format empty template features
            # properly. This is a quick fix.
            if(template_features["template_aatype"].shape[0] == 0):
                template_features = empty_template_feats(len(input_sequence))

        msa_features = self._process_msa_feats(alignment_dir)

        return {**pdb_feats, **template_features, **msa_features}
