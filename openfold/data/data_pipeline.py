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
import copy
import collections
import contextlib
import dataclasses
from multiprocessing import cpu_count
import tempfile
from typing import Mapping, Optional, Sequence, Any, MutableMapping, Union
import numpy as np
import torch
from openfold.data import (
    templates,
    parsers,
    mmcif_parsing,
    msa_identifiers,
    msa_pairing,
    feature_processing_multimer,
)
from openfold.data.templates import (
    get_custom_template_features,
    empty_template_feats,
    CustomHitFeaturizer,
)
from openfold.data.tools import jackhmmer, hhblits, hhsearch, hmmsearch
from openfold.np import residue_constants, protein

FeatureDict = MutableMapping[str, np.ndarray]
TemplateSearcher = Union[hhsearch.HHSearch, hmmsearch.Hmmsearch]


def make_template_features(
    input_sequence: str,
    hits: Sequence[Any],
    template_featurizer: Any,
) -> FeatureDict:
    hits_cat = sum(hits.values(), [])
    if template_featurizer is None or (
        len(hits_cat) == 0 and not isinstance(template_featurizer, CustomHitFeaturizer)
    ):
        template_features = empty_template_feats(len(input_sequence))
    else:
        templates_result = template_featurizer.get_templates(
            query_sequence=input_sequence,
            hits=hits_cat,
        )
        template_features = templates_result.features

    return template_features


def unify_template_features(
    template_feature_list: Sequence[FeatureDict]
) -> FeatureDict:
    out_dicts = []
    seq_lens = [fd["template_aatype"].shape[1] for fd in template_feature_list]
    for i, fd in enumerate(template_feature_list):
        out_dict = {}
        n_templates, n_res = fd["template_aatype"].shape[:2]
        for k,v in fd.items():
            seq_keys = [
                "template_aatype",
                "template_all_atom_positions",
                "template_all_atom_mask",
            ]
            if(k in seq_keys):
                new_shape = list(v.shape)
                assert(new_shape[1] == n_res)
                new_shape[1] = sum(seq_lens)
                new_array = np.zeros(new_shape, dtype=v.dtype)

                if(k == "template_aatype"):
                    new_array[..., residue_constants.HHBLITS_AA_TO_ID['-']] = 1

                offset = sum(seq_lens[:i])
                new_array[:, offset:offset + seq_lens[i]] = v
                out_dict[k] = new_array
            else:
                out_dict[k] = v

        chain_indices = np.array(n_templates * [i])
        out_dict["template_chain_index"] = chain_indices

        if(n_templates != 0):
            out_dicts.append(out_dict)

    if(len(out_dicts) > 0):
        out_dict = {
            k: np.concatenate([od[k] for od in out_dicts]) for k in out_dicts[0]
        }
    else:
        out_dict = empty_template_feats(sum(seq_lens))

    return out_dict


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
        [description.encode("utf-8")], dtype=object
    )
    features["residue_index"] = np.array(range(num_res), dtype=np.int32)
    features["seq_length"] = np.array([num_res] * num_res, dtype=np.int32)
    features["sequence"] = np.array(
        [sequence.encode("utf-8")], dtype=object
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
        [mmcif_object.header["release_date"].encode("utf-8")], dtype=object
    )

    mmcif_feats["is_distillation"] = np.array(0., dtype=np.float32)

    return mmcif_feats


def _aatype_to_str_sequence(aatype):
    return ''.join([
        residue_constants.restypes_with_x[aatype[i]]
        for i in range(len(aatype))
    ])


def make_protein_features(
    protein_object: protein.Protein,
    description: str,
    _is_distillation: bool = False,
) -> FeatureDict:
    pdb_feats = {}
    aatype = protein_object.aatype
    sequence = _aatype_to_str_sequence(aatype)
    pdb_feats.update(
        make_sequence_features(
            sequence=sequence,
            description=description,
            num_res=len(protein_object.aatype),
        )
    )

    all_atom_positions = protein_object.atom_positions
    all_atom_mask = protein_object.atom_mask

    pdb_feats["all_atom_positions"] = all_atom_positions.astype(np.float32)
    pdb_feats["all_atom_mask"] = all_atom_mask.astype(np.float32)

    pdb_feats["resolution"] = np.array([0.]).astype(np.float32)
    pdb_feats["is_distillation"] = np.array(
        1. if _is_distillation else 0.
    ).astype(np.float32)

    return pdb_feats


def make_pdb_features(
    protein_object: protein.Protein,
    description: str,
    is_distillation: bool = True,
    confidence_threshold: float = 50.,
) -> FeatureDict:
    pdb_feats = make_protein_features(
        protein_object, description, _is_distillation=True
    )

    if(is_distillation):
        high_confidence = protein_object.b_factors > confidence_threshold
        high_confidence = np.any(high_confidence, axis=-1)
        pdb_feats["all_atom_mask"] *= high_confidence[..., None]

    return pdb_feats


def make_msa_features(msas: Sequence[parsers.Msa]) -> FeatureDict:
    """Constructs a feature dict of MSA features."""
    if not msas:
        raise ValueError("At least one MSA must be provided.")

    int_msa = []
    deletion_matrix = []
    species_ids = []
    seen_sequences = set()
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(
                f"MSA {msa_index} must contain at least one sequence."
            )
        for sequence_index, sequence in enumerate(msa.sequences):
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            int_msa.append(
                [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence]
            )

            deletion_matrix.append(msa.deletion_matrix[sequence_index])
            identifiers = msa_identifiers.get_identifiers(
                msa.descriptions[sequence_index]
            )
            species_ids.append(identifiers.species_id.encode('utf-8'))

    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {}
    features["deletion_matrix_int"] = np.array(deletion_matrix, dtype=np.int32)
    features["msa"] = np.array(int_msa, dtype=np.int32)
    features["num_alignments"] = np.array(
        [num_alignments] * num_res, dtype=np.int32
    )
    features["msa_species_identifiers"] = np.array(species_ids, dtype=object)
    return features


def run_msa_tool(
    msa_runner,
    fasta_path: str,
    msa_out_path: str,
    msa_format: str,
    max_sto_sequences: Optional[int] = None,
) -> Mapping[str, Any]:
    """Runs an MSA tool, checking if output already exists first."""
    if(msa_format == "sto" and max_sto_sequences is not None):
        result = msa_runner.query(fasta_path, max_sto_sequences)[0]
    else:
        result = msa_runner.query(fasta_path)[0]

    assert msa_out_path.split('.')[-1] == msa_format
    with open(msa_out_path, "w") as f:
        f.write(result[msa_format])

    return result


def make_dummy_msa_obj(input_sequence) -> parsers.Msa:
    deletion_matrix = [[0 for _ in input_sequence]]
    return parsers.Msa(sequences=[input_sequence],
                       deletion_matrix=deletion_matrix,
                       descriptions=['dummy'])


# Generate 1-sequence MSA features having only the input sequence
def make_dummy_msa_feats(input_sequence) -> FeatureDict:
    msa_data_obj = make_dummy_msa_obj(input_sequence)
    return make_msa_features([msa_data_obj])


def make_sequence_features_with_custom_template(
        sequence: str,
        mmcif_path: str,
        pdb_id: str,
        chain_id: str,
        kalign_binary_path: str) -> FeatureDict:
    """
    process a single fasta file using features derived from a single template rather than an alignment
    """
    num_res = len(sequence)

    sequence_features = make_sequence_features(
        sequence=sequence,
        description=pdb_id,
        num_res=num_res,
    )

    msa_data = [sequence]
    deletion_matrix = [[0 for _ in sequence]]
    msa_data_obj = parsers.Msa(sequences=msa_data, deletion_matrix=deletion_matrix, descriptions=None)

    msa_features = make_msa_features([msa_data_obj])
    template_features = get_custom_template_features(
        mmcif_path=mmcif_path,
        query_sequence=sequence,
        pdb_id=pdb_id,
        chain_id=chain_id,
        kalign_binary_path=kalign_binary_path
    )

    return {
        **sequence_features,
        **msa_features,
        **template_features.features
    }


class AlignmentRunner:
    """Runs alignment tools and saves the results"""
    def __init__(
        self,
        jackhmmer_binary_path: Optional[str] = None,
        hhblits_binary_path: Optional[str] = None,
        uniref90_database_path: Optional[str] = None,
        mgnify_database_path: Optional[str] = None,
        bfd_database_path: Optional[str] = None,
        uniref30_database_path: Optional[str] = None,
        uniclust30_database_path: Optional[str] = None,
        uniprot_database_path: Optional[str] = None,
        template_searcher: Optional[TemplateSearcher] = None,
        use_small_bfd: Optional[bool] = None,
        no_cpus: Optional[int] = None,
        uniref_max_hits: int = 10000,
        mgnify_max_hits: int = 5000,
        uniprot_max_hits: int = 50000,
    ):
        """
        Args:
            jackhmmer_binary_path:
                Path to jackhmmer binary
            hhblits_binary_path:
                Path to hhblits binary
            uniref90_database_path:
                Path to uniref90 database. If provided, jackhmmer_binary_path
                must also be provided
            mgnify_database_path:
                Path to mgnify database. If provided, jackhmmer_binary_path
                must also be provided
            bfd_database_path:
                Path to BFD database. Depending on the value of use_small_bfd,
                one of hhblits_binary_path or jackhmmer_binary_path must be
                provided.
            uniref30_database_path:
                Path to uniref30. Searched alongside BFD if use_small_bfd is
                false.
            uniclust30_database_path:
                Path to uniclust30. Searched alongside BFD if use_small_bfd is
                false.
            use_small_bfd:
                Whether to search the BFD database alone with jackhmmer or
                in conjunction with uniref30/uniclust30 with hhblits.
            no_cpus:
                The number of CPUs available for alignment. By default, all
                CPUs are used.
            uniref_max_hits:
                Max number of uniref hits
            mgnify_max_hits:
                Max number of mgnify hits
            uniprot_max_hits:
                Max number of uniprot hits
        """
        db_map = {
            "jackhmmer": {
                "binary": jackhmmer_binary_path,
                "dbs": [
                    uniref90_database_path,
                    mgnify_database_path,
                    bfd_database_path if use_small_bfd else None,
                    uniprot_database_path,
                ],
            },
            "hhblits": {
                "binary": hhblits_binary_path,
                "dbs": [
                    bfd_database_path if not use_small_bfd else None,
                ],
            },
        }

        for name, dic in db_map.items():
            binary, dbs = dic["binary"], dic["dbs"]
            if(binary is None and not all([x is None for x in dbs])):
                raise ValueError(
                    f"{name} DBs provided but {name} binary is None"
                )

        self.uniref_max_hits = uniref_max_hits
        self.mgnify_max_hits = mgnify_max_hits
        self.uniprot_max_hits = uniprot_max_hits
        self.use_small_bfd = use_small_bfd

        if(no_cpus is None):
            no_cpus = cpu_count()

        self.jackhmmer_uniref90_runner = None
        if(jackhmmer_binary_path is not None and
            uniref90_database_path is not None
        ):
            self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniref90_database_path,
                n_cpu=no_cpus,
            )

        self.jackhmmer_small_bfd_runner = None
        self.hhblits_bfd_unirefclust_runner = None
        if(bfd_database_path is not None):
            if use_small_bfd:
                self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                    binary_path=jackhmmer_binary_path,
                    database_path=bfd_database_path,
                    n_cpu=no_cpus,
                )
            else:
                dbs = [bfd_database_path]
                if(uniref30_database_path is not None):
                    dbs.append(uniref30_database_path)
                if (uniclust30_database_path is not None):
                    dbs.append(uniclust30_database_path)
                self.hhblits_bfd_unirefclust_runner = hhblits.HHBlits(
                    binary_path=hhblits_binary_path,
                    databases=dbs,
                    n_cpu=no_cpus,
                )

        self.jackhmmer_mgnify_runner = None
        if(mgnify_database_path is not None):
            self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=mgnify_database_path,
                n_cpu=no_cpus,
            )

        self.jackhmmer_uniprot_runner = None
        if(uniprot_database_path is not None):
            self.jackhmmer_uniprot_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniprot_database_path,
                n_cpu=no_cpus
            )

        if(template_searcher is not None and
           self.jackhmmer_uniref90_runner is None
        ):
            raise ValueError(
                "Uniref90 runner must be specified to run template search"
            )

        self.template_searcher = template_searcher

    def run(
        self,
        fasta_path: str,
        output_dir: str,
    ):
        """Runs alignment tools on a sequence"""
        if(self.jackhmmer_uniref90_runner is not None):
            uniref90_out_path = os.path.join(output_dir, "uniref90_hits.sto")

            jackhmmer_uniref90_result = run_msa_tool(
                msa_runner=self.jackhmmer_uniref90_runner,
                fasta_path=fasta_path,
                msa_out_path=uniref90_out_path,
                msa_format='sto',
                max_sto_sequences=self.uniref_max_hits,
            )

            template_msa = jackhmmer_uniref90_result["sto"]
            template_msa = parsers.deduplicate_stockholm_msa(template_msa)
            template_msa = parsers.remove_empty_columns_from_stockholm_msa(
                template_msa
            )

            if(self.template_searcher is not None):
                if(self.template_searcher.input_format == "sto"):
                    pdb_templates_result = self.template_searcher.query(
                        template_msa,
                        output_dir=output_dir
                    )
                elif(self.template_searcher.input_format == "a3m"):
                    uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(
                        template_msa
                    )
                    pdb_templates_result = self.template_searcher.query(
                        uniref90_msa_as_a3m,
                        output_dir=output_dir
                    )
                else:
                    fmt = self.template_searcher.input_format
                    raise ValueError(
                        f"Unrecognized template input format: {fmt}"
                    )

        if(self.jackhmmer_mgnify_runner is not None):
            mgnify_out_path = os.path.join(output_dir, "mgnify_hits.sto")
            jackhmmer_mgnify_result = run_msa_tool(
                msa_runner=self.jackhmmer_mgnify_runner,
                fasta_path=fasta_path,
                msa_out_path=mgnify_out_path,
                msa_format='sto',
                max_sto_sequences=self.mgnify_max_hits
            )

        if(self.use_small_bfd and self.jackhmmer_small_bfd_runner is not None):
            bfd_out_path = os.path.join(output_dir, "small_bfd_hits.sto")
            jackhmmer_small_bfd_result = run_msa_tool(
                msa_runner=self.jackhmmer_small_bfd_runner,
                fasta_path=fasta_path,
                msa_out_path=bfd_out_path,
                msa_format="sto",
            )
        elif(self.hhblits_bfd_unirefclust_runner is not None):
            uni_name = "uni"
            for db_name in self.hhblits_bfd_unirefclust_runner.databases:
                if "uniref" in db_name.lower():
                    uni_name = f"{uni_name}ref"
                elif "uniclust" in db_name.lower():
                    uni_name = f"{uni_name}clust"

            bfd_out_path = os.path.join(output_dir, f"bfd_{uni_name}_hits.a3m")
            hhblits_bfd_unirefclust_result = run_msa_tool(
                msa_runner=self.hhblits_bfd_unirefclust_runner,
                fasta_path=fasta_path,
                msa_out_path=bfd_out_path,
                msa_format="a3m",
            )

        if(self.jackhmmer_uniprot_runner is not None):
            uniprot_out_path = os.path.join(output_dir, 'uniprot_hits.sto')
            result = run_msa_tool(
                self.jackhmmer_uniprot_runner,
                fasta_path=fasta_path,
                msa_out_path=uniprot_out_path,
                msa_format='sto',
                max_sto_sequences=self.uniprot_max_hits,
            )


@dataclasses.dataclass(frozen=True)
class _FastaChain:
    sequence: str
    description: str


def _make_chain_id_map(
    sequences: Sequence[str],
    descriptions: Sequence[str],
) -> Mapping[str, _FastaChain]:
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
      raise ValueError('sequences and descriptions must have equal length. '
                       f'Got {len(sequences)} != {len(descriptions)}.')
    if len(sequences) > protein.PDB_MAX_CHAINS:
      raise ValueError('Cannot process more chains than the PDB format supports. '
                       f'Got {len(sequences)} chains.')
    chain_id_map = {}
    for chain_id, sequence, description in zip(
        protein.PDB_CHAIN_IDS, sequences, descriptions
    ):
        chain_id_map[chain_id] = _FastaChain(
            sequence=sequence, description=description
        )
    return chain_id_map


@contextlib.contextmanager
def temp_fasta_file(fasta_str: str):
    with tempfile.NamedTemporaryFile('w', suffix='.fasta') as fasta_file:
      fasta_file.write(fasta_str)
      fasta_file.seek(0)
      yield fasta_file.name


def convert_monomer_features(
    monomer_features: FeatureDict,
    chain_id: str
) -> FeatureDict:
    """Reshapes and modifies monomer features for multimer models."""
    converted = {}
    converted['auth_chain_id'] = np.asarray(chain_id, dtype=object)
    unnecessary_leading_dim_feats = {
        'sequence', 'domain_name', 'num_alignments', 'seq_length'
    }
    for feature_name, feature in monomer_features.items():
      if feature_name in unnecessary_leading_dim_feats:
        # asarray ensures it's a np.ndarray.
        feature = np.asarray(feature[0], dtype=feature.dtype)
      elif feature_name == 'aatype':
        # The multimer model performs the one-hot operation itself.
        feature = np.argmax(feature, axis=-1).astype(np.int32)
      elif feature_name == 'template_aatype':
        feature = np.argmax(feature, axis=-1).astype(np.int32)
        new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
        feature = np.take(new_order_list, feature.astype(np.int32), axis=0)
      elif feature_name == 'template_all_atom_masks':
        feature_name = 'template_all_atom_mask'
      converted[feature_name] = feature
    return converted


def int_id_to_str_id(num: int) -> str:
    """Encodes a number as a string, using reverse spreadsheet style naming.

    Args:
      num: A positive integer.

    Returns:
      A string that encodes the positive integer using reverse spreadsheet style,
      naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
      usual way to encode chain IDs in mmCIF files.
    """
    if num <= 0:
      raise ValueError(f'Only positive integers allowed, got {num}.')

    num = num - 1  # 1-based indexing.
    output = []
    while num >= 0:
      output.append(chr(num % 26 + ord('A')))
      num = num // 26 - 1
    return ''.join(output)


def add_assembly_features(
    all_chain_features: MutableMapping[str, FeatureDict],
) -> MutableMapping[str, FeatureDict]:
    """Add features to distinguish between chains.

    Args:
      all_chain_features: A dictionary which maps chain_id to a dictionary of
        features for each chain.

    Returns:
      all_chain_features: A dictionary which maps strings of the form
        `<seq_id>_<sym_id>` to the corresponding chain features. E.g. two
        chains from a homodimer would have keys A_1 and A_2. Two chains from a
        heterodimer would have keys A_1 and B_1.
    """
    # Group the chains by sequence
    seq_to_entity_id = {}
    grouped_chains = collections.defaultdict(list)
    for chain_id, chain_features in all_chain_features.items():
      seq = str(chain_features['sequence'])
      if seq not in seq_to_entity_id:
        seq_to_entity_id[seq] = len(seq_to_entity_id) + 1
      grouped_chains[seq_to_entity_id[seq]].append(chain_features)

    new_all_chain_features = {}
    chain_id = 1
    for entity_id, group_chain_features in grouped_chains.items():
      for sym_id, chain_features in enumerate(group_chain_features, start=1):
        new_all_chain_features[
            f'{int_id_to_str_id(entity_id)}_{sym_id}'] = chain_features
        seq_length = chain_features['seq_length']
        chain_features['asym_id'] = (
            chain_id * np.ones(seq_length)
        ).astype(np.int64)
        chain_features['sym_id'] = (
            sym_id * np.ones(seq_length)
        ).astype(np.int64)
        chain_features['entity_id'] = (
            entity_id * np.ones(seq_length)
        ).astype(np.int64)
        chain_id += 1

    return new_all_chain_features


def pad_msa(np_example, min_num_seq):
    np_example = dict(np_example)
    num_seq = np_example['msa'].shape[0]
    if num_seq < min_num_seq:
      for feat in ('msa', 'deletion_matrix', 'bert_mask', 'msa_mask'):
        np_example[feat] = np.pad(
            np_example[feat], ((0, min_num_seq - num_seq), (0, 0)))
      np_example['cluster_bias_mask'] = np.pad(
          np_example['cluster_bias_mask'], ((0, min_num_seq - num_seq),))
    return np_example


class DataPipeline:
    """Assembles input features."""
    def __init__(
        self,
        template_featurizer: Optional[templates.TemplateHitFeaturizer],
    ):
        self.template_featurizer = template_featurizer

    def _parse_msa_data(
        self,
        alignment_dir: str,
        alignment_index: Optional[Any] = None
    ) -> Mapping[str, Any]:
        msa_data = {}
        if alignment_index is not None:
            fp = open(os.path.join(alignment_dir, alignment_index["db"]), "rb")

            def read_msa(start, size):
                fp.seek(start)
                msa = fp.read(size).decode("utf-8")
                return msa

            for (name, start, size) in alignment_index["files"]:
                filename, ext = os.path.splitext(name)

                if ext == ".a3m":
                    msa = parsers.parse_a3m(
                        read_msa(start, size)
                    )
                # The "hmm_output" exception is a crude way to exclude
                # multimer template hits.
                # Multimer "uniprot_hits" processed separately.
                elif ext == ".sto" and filename not in ["uniprot_hits", "hmm_output"]:
                    msa = parsers.parse_stockholm(read_msa(start, size))
                else:
                    continue

                msa_data[name] = msa

            fp.close()
        else:
            for f in os.listdir(alignment_dir):
                path = os.path.join(alignment_dir, f)
                filename, ext = os.path.splitext(f)

                if ext == ".a3m":
                    with open(path, "r") as fp:
                        msa = parsers.parse_a3m(fp.read())
                elif ext == ".sto" and filename not in ["uniprot_hits", "hmm_output"]:
                    with open(path, "r") as fp:
                        msa = parsers.parse_stockholm(
                            fp.read()
                        )
                else:
                    continue

                msa_data[f] = msa

        return msa_data

    def _parse_template_hit_files(
        self,
        alignment_dir: str,
        input_sequence: str,
        alignment_index: Optional[Any] = None
    ) -> Mapping[str, Any]:
        all_hits = {}
        if(alignment_index is not None):
            fp = open(os.path.join(alignment_dir, alignment_index["db"]), 'rb')

            def read_template(start, size):
                fp.seek(start)
                return fp.read(size).decode("utf-8")

            for (name, start, size) in alignment_index["files"]:
                ext = os.path.splitext(name)[-1]

                if(ext == ".hhr"):
                    hits = parsers.parse_hhr(read_template(start, size))
                    all_hits[name] = hits
                elif(name == "hmmsearch_output.sto"):
                    hits = parsers.parse_hmmsearch_sto(
                        read_template(start, size),
                        input_sequence,
                    )
                    all_hits[name] = hits

            fp.close()
        else:
            for f in os.listdir(alignment_dir):
                path = os.path.join(alignment_dir, f)
                ext = os.path.splitext(f)[-1]

                if(ext == ".hhr"):
                    with open(path, "r") as fp:
                        hits = parsers.parse_hhr(fp.read())
                    all_hits[f] = hits
                elif(f == "hmm_output.sto"):
                    with open(path, "r") as fp:
                        hits = parsers.parse_hmmsearch_sto(
                            fp.read(),
                            input_sequence,
                        )
                    all_hits[f] = hits

        return all_hits

    def _get_msas(self,
        alignment_dir: str,
        input_sequence: Optional[str] = None,
        alignment_index: Optional[Any] = None,
    ):
        msa_data = self._parse_msa_data(alignment_dir, alignment_index)
        if(len(msa_data) == 0):
            if(input_sequence is None):
                raise ValueError(
                    """
                    If the alignment dir contains no MSAs, an input sequence 
                    must be provided.
                    """
                )

            msa_data["dummy"] = make_dummy_msa_obj(input_sequence)

        return list(msa_data.values())

    def _process_msa_feats(
        self,
        alignment_dir: str,
        input_sequence: Optional[str] = None,
        alignment_index: Optional[Any] = None
    ) -> Mapping[str, Any]:

        msas = self._get_msas(
            alignment_dir, input_sequence, alignment_index
        )
        msa_features = make_msa_features(
            msas=msas
        )

        return msa_features

    # Load and process sequence embedding features
    def _process_seqemb_features(self,
        alignment_dir: str,
    ) -> Mapping[str, Any]:
        seqemb_features = {}
        for f in os.listdir(alignment_dir):
            path = os.path.join(alignment_dir, f)
            ext = os.path.splitext(f)[-1]

            if (ext == ".pt"):
                # Load embedding file
                seqemb_data = torch.load(path)
                seqemb_features["seq_embedding"] = seqemb_data["representations"][33]

        return seqemb_features

    def process_fasta(
        self,
        fasta_path: str,
        alignment_dir: str,
        alignment_index: Optional[Any] = None,
        seqemb_mode: bool = False,
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

        hits = self._parse_template_hit_files(
            alignment_dir=alignment_dir,
            input_sequence=input_sequence,
            alignment_index=alignment_index,
        )

        template_features = make_template_features(
            input_sequence,
            hits,
            self.template_featurizer,
        )

        sequence_features = make_sequence_features(
            sequence=input_sequence,
            description=input_description,
            num_res=num_res,
        )

        sequence_embedding_features = {}
        # If using seqemb mode, generate a dummy MSA features using just the sequence
        if seqemb_mode:
            msa_features = make_dummy_msa_feats(input_sequence)
            sequence_embedding_features = self._process_seqemb_features(alignment_dir)
        else:
            msa_features = self._process_msa_feats(alignment_dir, input_sequence, alignment_index)
        
        return {
            **sequence_features,
            **msa_features, 
            **template_features,
            **sequence_embedding_features
        }

    def process_mmcif(
        self,
        mmcif: mmcif_parsing.MmcifObject,  # parsing is expensive, so no path
        alignment_dir: str,
        chain_id: Optional[str] = None,
        alignment_index: Optional[Any] = None,
        seqemb_mode: bool = False,
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
        hits = self._parse_template_hit_files(
            alignment_dir=alignment_dir,
            input_sequence=input_sequence,
            alignment_index=alignment_index)

        template_features = make_template_features(
            input_sequence,
            hits,
            self.template_featurizer
        )

        sequence_embedding_features = {}
        # If using seqemb mode, generate a dummy MSA features using just the sequence
        if seqemb_mode:
            msa_features = make_dummy_msa_feats(input_sequence)
            sequence_embedding_features = self._process_seqemb_features(alignment_dir)
        else:
            msa_features = self._process_msa_feats(alignment_dir, input_sequence, alignment_index)

        return {**mmcif_feats, **template_features, **msa_features, **sequence_embedding_features}

    def process_pdb(
        self,
        pdb_path: str,
        alignment_dir: str,
        is_distillation: bool = True,
        chain_id: Optional[str] = None,
        _structure_index: Optional[str] = None,
        alignment_index: Optional[Any] = None,
        seqemb_mode: bool = False,
    ) -> FeatureDict:
        """
            Assembles features for a protein in a PDB file.
        """
        if(_structure_index is not None):
            db_dir = os.path.dirname(pdb_path)
            db = _structure_index["db"]
            db_path = os.path.join(db_dir, db)
            fp = open(db_path, "rb")
            _, offset, length = _structure_index["files"][0]
            fp.seek(offset)
            pdb_str = fp.read(length).decode("utf-8")
            fp.close()
        else:
            with open(pdb_path, 'r') as f:
                pdb_str = f.read()

        protein_object = protein.from_pdb_string(pdb_str, chain_id)
        input_sequence = _aatype_to_str_sequence(protein_object.aatype)
        description = os.path.splitext(os.path.basename(pdb_path))[0].upper()
        pdb_feats = make_pdb_features(
            protein_object,
            description,
            is_distillation=is_distillation
        )

        hits = self._parse_template_hit_files(
            alignment_dir=alignment_dir,
            input_sequence=input_sequence,
            alignment_index=alignment_index,
        )

        template_features = make_template_features(
            input_sequence,
            hits,
            self.template_featurizer,
        )

        sequence_embedding_features = {}
        # If in sequence embedding mode, generate dummy MSA features using just the input sequence
        if seqemb_mode:
            msa_features = make_dummy_msa_feats(input_sequence)
            sequence_embedding_features = self._process_seqemb_features(alignment_dir)
        else:
            msa_features = self._process_msa_feats(alignment_dir, input_sequence, alignment_index)

        return {**pdb_feats, **template_features, **msa_features, **sequence_embedding_features}

    def process_core(
        self,
        core_path: str,
        alignment_dir: str,
        alignment_index: Optional[Any] = None,
        seqemb_mode: bool = False,
    ) -> FeatureDict:
        """
            Assembles features for a protein in a ProteinNet .core file.
        """
        with open(core_path, 'r') as f:
            core_str = f.read()

        protein_object = protein.from_proteinnet_string(core_str)
        input_sequence = _aatype_to_str_sequence(protein_object.aatype)
        description = os.path.splitext(os.path.basename(core_path))[0].upper()
        core_feats = make_protein_features(protein_object, description)

        hits = self._parse_template_hit_files(
            alignment_dir=alignment_dir,
            input_sequence=input_sequence,
            alignment_index=alignment_index,
        )

        template_features = make_template_features(
            input_sequence,
            hits,
            self.template_featurizer,
        )

        sequence_embedding_features = {}
        # If in sequence embedding mode, generate dummy MSA features using just the input sequence
        if seqemb_mode:
            msa_features = make_dummy_msa_feats(input_sequence)
            sequence_embedding_features = self._process_seqemb_features(alignment_dir)
        else:
            msa_features = self._process_msa_feats(alignment_dir, input_sequence)

        return {**core_feats, **template_features, **msa_features, **sequence_embedding_features}

    def process_multiseq_fasta(self,
                               fasta_path: str,
                               super_alignment_dir: str,
                               ri_gap: int = 200,
                               ) -> FeatureDict:
        """
            Assembles features for a multi-sequence FASTA. Uses Minkyung Baek's
            hack from Twitter (a.k.a. AlphaFold-Gap).
        """
        with open(fasta_path, 'r') as f:
            fasta_str = f.read()

        input_seqs, input_descs = parsers.parse_fasta(fasta_str)

        # No whitespace allowed
        input_descs = [i.split()[0] for i in input_descs]

        # Stitch all of the sequences together
        input_sequence = ''.join(input_seqs)
        input_description = '-'.join(input_descs)
        num_res = len(input_sequence)

        sequence_features = make_sequence_features(
            sequence=input_sequence,
            description=input_description,
            num_res=num_res,
        )

        seq_lens = [len(s) for s in input_seqs]
        total_offset = 0
        for sl in seq_lens:
            total_offset += sl
            sequence_features["residue_index"][total_offset:] += ri_gap

        msa_list = []
        deletion_mat_list = []
        for seq, desc in zip(input_seqs, input_descs):
            alignment_dir = os.path.join(
                super_alignment_dir, desc
            )
            msas = self._get_msas(
                alignment_dir, seq, None
            )
            msa_list.append([m.sequences for m in msas])
            deletion_mat_list.append([m.deletion_matrix for m in msas])

        final_msa = []
        final_deletion_mat = []
        final_msa_obj = []
        msa_it = enumerate(zip(msa_list, deletion_mat_list))
        for i, (msas, deletion_mats) in msa_it:
            prec, post = sum(seq_lens[:i]), sum(seq_lens[i + 1:])
            msas = [
                [prec * '-' + seq + post * '-' for seq in msa] for msa in msas
            ]
            deletion_mats = [
                [prec * [0] + dml + post * [0] for dml in deletion_mat]
                for deletion_mat in deletion_mats
            ]

            assert (len(msas[0][-1]) == len(input_sequence))

            final_msa.extend(msas)
            final_deletion_mat.extend(deletion_mats)
            final_msa_obj.extend([parsers.Msa(sequences=msas[k], deletion_matrix=deletion_mats[k], descriptions=None)
                                  for k in range(len(msas))])

        msa_features = make_msa_features(
            msas=final_msa_obj
        )

        template_feature_list = []
        for seq, desc in zip(input_seqs, input_descs):
            alignment_dir = os.path.join(
                super_alignment_dir, desc
            )
            hits = self._parse_template_hit_files(alignment_dir=alignment_dir,
                                                  input_sequence=seq,
                                                  alignment_index=None)

            template_features = make_template_features(
                seq,
                hits,
                self.template_featurizer,
            )
            template_feature_list.append(template_features)

        template_features = unify_template_features(template_feature_list)

        return {
            **sequence_features,
            **msa_features,
            **template_features,
        }


class DataPipelineMultimer:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 monomer_data_pipeline: DataPipeline,
                 ):
        """Initializes the data pipeline.
        Args:
          monomer_data_pipeline: An instance of pipeline.DataPipeline - that runs
            the data pipeline for the monomer AlphaFold system.
          jackhmmer_binary_path: Location of the jackhmmer binary.
          uniprot_database_path: Location of the unclustered uniprot sequences, that
            will be searched with jackhmmer and used for MSA pairing.
          max_uniprot_hits: The maximum number of hits to return from uniprot.
          use_precomputed_msas: Whether to use pre-existing MSAs; see run_alphafold.
        """
        self._monomer_data_pipeline = monomer_data_pipeline

    def _process_single_chain(
            self,
            chain_id: str,
            sequence: str,
            description: str,
            chain_alignment_dir: str,
            chain_alignment_index: Optional[Any],
            is_homomer_or_monomer: bool
    ) -> FeatureDict:
        """Runs the monomer pipeline on a single chain."""
        chain_fasta_str = f'>{chain_id}\n{sequence}\n'

        if chain_alignment_index is None and not os.path.exists(chain_alignment_dir):
            raise ValueError(f"Alignments for {chain_id} not found...")

        with temp_fasta_file(chain_fasta_str) as chain_fasta_path:
            chain_features = self._monomer_data_pipeline.process_fasta(
                fasta_path=chain_fasta_path,
                alignment_dir=chain_alignment_dir,
                alignment_index=chain_alignment_index
            )

            # We only construct the pairing features if there are 2 or more unique
            # sequences.
            if not is_homomer_or_monomer:
                all_seq_msa_features = self._all_seq_msa_features(
                    chain_alignment_dir,
                    chain_alignment_index
                )
                chain_features.update(all_seq_msa_features)
        return chain_features

    @staticmethod
    def _all_seq_msa_features(alignment_dir, alignment_index):
        """Get MSA features for unclustered uniprot, for pairing."""
        if alignment_index is not None:
            fp = open(os.path.join(alignment_dir, alignment_index["db"]), "rb")

            def read_msa(start, size):
                fp.seek(start)
                msa = fp.read(size).decode("utf-8")
                return msa

            start, size = next(iter((start, size) for name, start, size in alignment_index["files"]
                                    if name == 'uniprot_hits.sto'))

            msa = parsers.parse_stockholm(read_msa(start, size))
            fp.close()
        else:
            uniprot_msa_path = os.path.join(alignment_dir, "uniprot_hits.sto")
            if not os.path.exists(uniprot_msa_path):
                chain_id = os.path.basename(os.path.normpath(alignment_dir))
                raise ValueError(f"Missing 'uniprot_hits.sto' for {chain_id}. "
                                 f"This is required for Multimer MSA pairing.")

            with open(uniprot_msa_path, "r") as fp:
                uniprot_msa_string = fp.read()
            msa = parsers.parse_stockholm(uniprot_msa_string)

        all_seq_features = make_msa_features([msa])
        valid_feats = msa_pairing.MSA_FEATURES + (
            'msa_species_identifiers',
        )
        feats = {
            f'{k}_all_seq': v for k, v in all_seq_features.items()
            if k in valid_feats
        }
        return feats

    def process_fasta(self,
                      fasta_path: str,
                      alignment_dir: str,
                      alignment_index: Optional[Any] = None
                      ) -> FeatureDict:
        """Creates features."""
        with open(fasta_path) as f:
            input_fasta_str = f.read()

        input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)

        all_chain_features = {}
        sequence_features = {}
        is_homomer_or_monomer = len(set(input_seqs)) == 1
        for desc, seq in zip(input_descs, input_seqs):
            if seq in sequence_features:
                all_chain_features[desc] = copy.deepcopy(
                    sequence_features[seq]
                )
                continue

            if alignment_index is not None:
                chain_alignment_index = alignment_index.get(desc)
                chain_alignment_dir = alignment_dir
            else:
                chain_alignment_index = None
                chain_alignment_dir = os.path.join(alignment_dir, desc)

            chain_features = self._process_single_chain(
                chain_id=desc,
                sequence=seq,
                description=desc,
                chain_alignment_dir=chain_alignment_dir,
                chain_alignment_index=chain_alignment_index,
                is_homomer_or_monomer=is_homomer_or_monomer
            )

            chain_features = convert_monomer_features(
                chain_features,
                chain_id=desc
            )
            all_chain_features[desc] = chain_features
            sequence_features[seq] = chain_features

        all_chain_features = add_assembly_features(all_chain_features)

        np_example = feature_processing_multimer.pair_and_merge(
            all_chain_features=all_chain_features,
        )

        # Pad MSA to avoid zero-sized extra_msa.
        np_example = pad_msa(np_example, 512)

        return np_example
    
    def get_mmcif_features(
            self, mmcif_object: mmcif_parsing.MmcifObject, chain_id: str
    ) -> FeatureDict:
        mmcif_feats = {}

        all_atom_positions, all_atom_mask = mmcif_parsing.get_atom_coords(
            mmcif_object=mmcif_object, chain_id=chain_id
        )
        mmcif_feats["all_atom_positions"] = all_atom_positions
        mmcif_feats["all_atom_mask"] = all_atom_mask

        mmcif_feats["resolution"] = np.array(
            mmcif_object.header["resolution"], dtype=np.float32
        )

        mmcif_feats["release_date"] = np.array(
            [mmcif_object.header["release_date"].encode("utf-8")], dtype=object
        )

        mmcif_feats["is_distillation"] = np.array(0., dtype=np.float32)

        return mmcif_feats

    def process_mmcif(
            self,
            mmcif: mmcif_parsing.MmcifObject,  # parsing is expensive, so no path
            alignment_dir: str,
            alignment_index: Optional[Any] = None,
    ) -> FeatureDict:

        all_chain_features = {}
        sequence_features = {}
        is_homomer_or_monomer = len(set(list(mmcif.chain_to_seqres.values()))) == 1
        for chain_id, seq in mmcif.chain_to_seqres.items():
            desc= "_".join([mmcif.file_id, chain_id])

            if seq in sequence_features:
                all_chain_features[desc] = copy.deepcopy(
                    sequence_features[seq]
                )
                continue

            if alignment_index is not None:
                chain_alignment_index = alignment_index.get(desc)
                chain_alignment_dir = alignment_dir
            else:
                chain_alignment_index = None
                chain_alignment_dir = os.path.join(alignment_dir, desc)

            chain_features = self._process_single_chain(
                chain_id=desc,
                sequence=seq,
                description=desc,
                chain_alignment_dir=chain_alignment_dir,
                chain_alignment_index=chain_alignment_index,
                is_homomer_or_monomer=is_homomer_or_monomer
            )

            chain_features = convert_monomer_features(
                chain_features,
                chain_id=desc
            )

            mmcif_feats = self.get_mmcif_features(mmcif, chain_id)
            chain_features.update(mmcif_feats)
            all_chain_features[desc] = chain_features
            sequence_features[seq] = chain_features

        all_chain_features = add_assembly_features(all_chain_features)

        np_example = feature_processing_multimer.pair_and_merge(
            all_chain_features=all_chain_features,
        )

        # Pad MSA to avoid zero-sized extra_msa.
        np_example = pad_msa(np_example, 512)

        return np_example
