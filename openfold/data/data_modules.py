import copy
from functools import partial
import json
import logging
import os
import pickle
from typing import Optional, Sequence, List, Any

import ml_collections as mlc
import numpy as np
import pytorch_lightning as pl
import torch
from torch.utils.data import RandomSampler

from openfold.data import (
    data_pipeline,
    feature_pipeline,
    mmcif_parsing,
    templates,
)
from openfold.utils.tensor_utils import tensor_tree_map, dict_multimap


class OpenFoldSingleDataset(torch.utils.data.Dataset):
    def __init__(self,
        data_dir: str,
        alignment_dir: str, 
        template_mmcif_dir: str,
        max_template_date: str,
        config: mlc.ConfigDict,
        kalign_binary_path: str = '/usr/bin/kalign',
        max_template_hits: int = 4,
        obsolete_pdbs_file_path: Optional[str] = None,
        template_release_dates_cache_path: Optional[str] = None,
        shuffle_top_k_prefiltered: Optional[int] = None,
        treat_pdb_as_distillation: bool = True,
        filter_path: Optional[str] = None,
        mode: str = "train", 
        alignment_index: Optional[Any] = None,
        _output_raw: bool = False,
        _structure_index: Optional[Any] = None,
    ):
        """
            Args:
                data_dir:
                    A path to a directory containing mmCIF files (in train
                    mode) or FASTA files (in inference mode).
                alignment_dir:
                    A path to a directory containing only data in the format 
                    output by an AlignmentRunner 
                    (defined in openfold.features.alignment_runner).
                    I.e. a directory of directories named {PDB_ID}_{CHAIN_ID}
                    or simply {PDB_ID}, each containing .a3m, .sto, and .hhr
                    files.
                template_mmcif_dir:
                    Path to a directory containing template mmCIF files.
                config:
                    A dataset config object. See openfold.config
                kalign_binary_path:
                    Path to kalign binary.
                max_template_hits:
                    An upper bound on how many templates are considered. During
                    training, the templates ultimately used are subsampled
                    from this total quantity.
                template_release_dates_cache_path:
                    Path to the output of scripts/generate_mmcif_cache.
                obsolete_pdbs_file_path:
                    Path to the file containing replacements for obsolete PDBs.
                shuffle_top_k_prefiltered:
                    Whether to uniformly shuffle the top k template hits before
                    parsing max_template_hits of them. Can be used to
                    approximate DeepMind's training-time template subsampling
                    scheme much more performantly.
                treat_pdb_as_distillation:
                    Whether to assume that .pdb files in the data_dir are from
                    the self-distillation set (and should be subjected to
                    special distillation set preprocessing steps).
                mode:
                    "train", "val", or "predict"
        """
        super(OpenFoldSingleDataset, self).__init__()
        self.data_dir = data_dir
        self.alignment_dir = alignment_dir
        self.config = config
        self.treat_pdb_as_distillation = treat_pdb_as_distillation
        self.mode = mode
        self.alignment_index = alignment_index
        self._output_raw = _output_raw
        self._structure_index = _structure_index

        self.supported_exts = [".cif", ".core", ".pdb"]

        valid_modes = ["train", "eval", "predict"]
        if(mode not in valid_modes):
            raise ValueError(f'mode must be one of {valid_modes}')

        if(template_release_dates_cache_path is None):
            logging.warning(
                "Template release dates cache does not exist. Remember to run "
                "scripts/generate_mmcif_cache.py before running OpenFold"
            )

        if(alignment_index is not None):
            self._chain_ids = list(alignment_index.keys())
        else:
            self._chain_ids = list(os.listdir(alignment_dir))
        
        if(filter_path is not None):
            with open(filter_path, "r") as f:
                chains_to_include = set([l.strip() for l in f.readlines()])

            self._chain_ids = [c for c in self._chain_ids if c in chains_to_include]
       
        self._chain_id_to_idx_dict = {
            chain: i for i, chain in enumerate(self._chain_ids)
        }

        template_featurizer = templates.TemplateHitFeaturizer(
            mmcif_dir=template_mmcif_dir,
            max_template_date=max_template_date,
            max_hits=max_template_hits,
            kalign_binary_path=kalign_binary_path,
            release_dates_path=template_release_dates_cache_path,
            obsolete_pdbs_path=obsolete_pdbs_file_path,
            _shuffle_top_k_prefiltered=shuffle_top_k_prefiltered,
        )

        self.data_pipeline = data_pipeline.DataPipeline(
            template_featurizer=template_featurizer,
        )

        if(not self._output_raw):
            self.feature_pipeline = feature_pipeline.FeaturePipeline(config) 

    def _parse_mmcif(self, path, file_id, chain_id, alignment_dir, alignment_index):
        with open(path, 'r') as f:
            mmcif_string = f.read()

        mmcif_object = mmcif_parsing.parse(
            file_id=file_id, mmcif_string=mmcif_string
        )

        # Crash if an error is encountered. Any parsing errors should have
        # been dealt with at the alignment stage.
        if(mmcif_object.mmcif_object is None):
            raise list(mmcif_object.errors.values())[0]

        mmcif_object = mmcif_object.mmcif_object

        data = self.data_pipeline.process_mmcif(
            mmcif=mmcif_object,
            alignment_dir=alignment_dir,
            chain_id=chain_id,
            alignment_index=alignment_index
        )

        return data

    def chain_id_to_idx(self, chain_id):
        return self._chain_id_to_idx_dict[chain_id]

    def idx_to_chain_id(self, idx):
        return self._chain_ids[idx]

    def __getitem__(self, idx):
        name = self.idx_to_chain_id(idx)
        alignment_dir = os.path.join(self.alignment_dir, name)

        alignment_index = None
        if(self.alignment_index is not None):
            alignment_dir = self.alignment_dir
            alignment_index = self.alignment_index[name]

        if(self.mode == 'train' or self.mode == 'eval'):
            spl = name.rsplit('_', 1)
            if(len(spl) == 2):
                file_id, chain_id = spl
            else:
                file_id, = spl
                chain_id = None

            path = os.path.join(self.data_dir, file_id)
            structure_index_entry = None
            if(self._structure_index is not None):
                structure_index_entry = self._structure_index[name]
                assert(len(structure_index_entry["files"]) == 1)
                filename, _, _ = structure_index_entry["files"][0]
                ext = os.path.splitext(filename)[1]
            else:
                ext = None
                for e in self.supported_exts:
                    if(os.path.exists(path + e)):
                        ext = e
                        break

                if(ext is None):
                    raise ValueError("Invalid file type")

            path += ext
            if(ext == ".cif"):
                data = self._parse_mmcif(
                    path, file_id, chain_id, alignment_dir, alignment_index,
                )
            elif(ext == ".core"):
                data = self.data_pipeline.process_core(
                    path, alignment_dir, alignment_index,
                )
            elif(ext == ".pdb"):
                structure_index = None
                if(self._structure_index is not None):
                    structure_index = self._structure_index[name]
                data = self.data_pipeline.process_pdb(
                    pdb_path=path,
                    alignment_dir=alignment_dir,
                    is_distillation=self.treat_pdb_as_distillation,
                    chain_id=chain_id,
                    alignment_index=alignment_index,
                    _structure_index=structure_index,
                )
            else:
               raise ValueError("Extension branch missing") 
        else:
            path = os.path.join(name, name + ".fasta")
            data = self.data_pipeline.process_fasta(
                fasta_path=path,
                alignment_dir=alignment_dir,
                alignment_index=alignment_index,
            )

        if(self._output_raw):
            return data

        feats = self.feature_pipeline.process_features(
            data, self.mode 
        )

        feats["batch_idx"] = torch.tensor([idx for _ in range(feats["aatype"].shape[-1])], dtype=torch.int64, device=feats["aatype"].device)

        return feats

    def __len__(self):
        return len(self._chain_ids) 


def deterministic_train_filter(
    chain_data_cache_entry: Any,
    max_resolution: float = 9.,
    max_single_aa_prop: float = 0.8,
) -> bool:
    # Hard filters
    resolution = chain_data_cache_entry.get("resolution", None)
    if(resolution is not None and resolution > max_resolution):
        return False

    seq = chain_data_cache_entry["seq"]
    counts = {}
    for aa in seq:
        counts.setdefault(aa, 0)
        counts[aa] += 1
    largest_aa_count = max(counts.values())
    largest_single_aa_prop = largest_aa_count / len(seq)
    if(largest_single_aa_prop > max_single_aa_prop):
        return False

    return True


def get_stochastic_train_filter_prob(
    chain_data_cache_entry: Any,
) -> List[float]:
    # Stochastic filters
    probabilities = []
    
    cluster_size = chain_data_cache_entry.get("cluster_size", None)
    if(cluster_size is not None and cluster_size > 0):
        probabilities.append(1 / cluster_size)
    
    chain_length = len(chain_data_cache_entry["seq"])
    probabilities.append((1 / 512) * (max(min(chain_length, 512), 256)))

    # Risk of underflow here?
    out = 1
    for p in probabilities:
        out *= p

    return out


class OpenFoldDataset(torch.utils.data.Dataset):
    """
        Implements the stochastic filters applied during AlphaFold's training.
        Because samples are selected from constituent datasets randomly, the
        length of an OpenFoldFilteredDataset is arbitrary. Samples are selected
        and filtered once at initialization.
    """
    def __init__(self,
        datasets: Sequence[OpenFoldSingleDataset],
        probabilities: Sequence[int],
        epoch_len: int,
        chain_data_cache_paths: List[str],
        generator: torch.Generator = None,
        _roll_at_init: bool = True,
    ):
        self.datasets = datasets
        self.probabilities = probabilities
        self.epoch_len = epoch_len
        self.generator = generator
        
        self.chain_data_caches = []
        for path in chain_data_cache_paths:
            with open(path, "r") as fp:
                self.chain_data_caches.append(json.load(fp))

        def looped_shuffled_dataset_idx(dataset_len):
            while True:
                # Uniformly shuffle each dataset's indices
                weights = [1. for _ in range(dataset_len)]
                shuf = torch.multinomial(
                    torch.tensor(weights),
                    num_samples=dataset_len,
                    replacement=False,
                    generator=self.generator,
                )
                for idx in shuf:
                    yield idx

        def looped_samples(dataset_idx):
            max_cache_len = int(epoch_len * probabilities[dataset_idx])
            dataset = self.datasets[dataset_idx]
            idx_iter = looped_shuffled_dataset_idx(len(dataset))
            chain_data_cache = self.chain_data_caches[dataset_idx]
            while True:
                weights = []
                idx = []
                for _ in range(max_cache_len):
                    candidate_idx = next(idx_iter)
                    chain_id = dataset.idx_to_chain_id(candidate_idx)
                    chain_data_cache_entry = chain_data_cache[chain_id]
                    if(not deterministic_train_filter(chain_data_cache_entry)):
                        continue

                    p = get_stochastic_train_filter_prob(
                        chain_data_cache_entry,
                    )
                    weights.append([1. - p, p])
                    idx.append(candidate_idx)

                samples = torch.multinomial(
                    torch.tensor(weights),
                    num_samples=1,
                    generator=self.generator,
                )
                samples = samples.squeeze()

                cache = [i for i, s in zip(idx, samples) if s]

                for datapoint_idx in cache:
                    yield datapoint_idx

        self._samples = [looped_samples(i) for i in range(len(self.datasets))]

        if(_roll_at_init):
            self.reroll()

    def __getitem__(self, idx):
        dataset_idx, datapoint_idx = self.datapoints[idx]
        return self.datasets[dataset_idx][datapoint_idx]

    def __len__(self):
        return self.epoch_len

    def reroll(self):
        dataset_choices = torch.multinomial(
            torch.tensor(self.probabilities),
            num_samples=self.epoch_len,
            replacement=True,
            generator=self.generator,
        )

        self.datapoints = []
        for dataset_idx in dataset_choices:
            samples = self._samples[dataset_idx]
            datapoint_idx = next(samples)
            self.datapoints.append((dataset_idx, datapoint_idx))


class OpenFoldBatchCollator:
    def __call__(self, prots):
        stack_fn = partial(torch.stack, dim=0)
        return dict_multimap(stack_fn, prots) 


class OpenFoldDataLoader(torch.utils.data.DataLoader):
    def __init__(self, *args, config, stage="train", generator=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.config = config
        self.stage = stage    

        if(generator is None):
            generator = torch.Generator()
        
        self.generator = generator
        self._prep_batch_properties_probs()

    def _prep_batch_properties_probs(self):
        keyed_probs = []
        stage_cfg = self.config[self.stage]

        max_iters = self.config.common.max_recycling_iters
        if(stage_cfg.supervised):
            clamp_prob = self.config.supervised.clamp_prob
            keyed_probs.append(
                ("use_clamped_fape", [1 - clamp_prob, clamp_prob])
            )
        
        if(stage_cfg.uniform_recycling):
            recycling_probs = [
                1. / (max_iters + 1) for _ in range(max_iters + 1)
            ]
        else:
            recycling_probs = [
                0. for _ in range(max_iters + 1)
            ]
            recycling_probs[-1] = 1.
        
        keyed_probs.append(
            ("no_recycling_iters", recycling_probs)
        )

        keys, probs = zip(*keyed_probs)
        max_len = max([len(p) for p in probs])
        padding = [[0.] * (max_len - len(p)) for p in probs] 
        
        self.prop_keys = keys
        self.prop_probs_tensor = torch.tensor(
            [p + pad for p, pad in zip(probs, padding)],
            dtype=torch.float32,
        )

    def _add_batch_properties(self, batch):
        samples = torch.multinomial(
            self.prop_probs_tensor,
            num_samples=1, # 1 per row
            replacement=True,
            generator=self.generator
        )

        aatype = batch["aatype"]
        batch_dims = aatype.shape[:-2]
        recycling_dim = aatype.shape[-1]
        no_recycling = recycling_dim
        for i, key in enumerate(self.prop_keys):
            sample = int(samples[i][0])
            sample_tensor = torch.tensor(
                sample, 
                device=aatype.device, 
                requires_grad=False
            )
            orig_shape = sample_tensor.shape
            sample_tensor = sample_tensor.view(
                (1,) * len(batch_dims) + sample_tensor.shape + (1,)
            )
            sample_tensor = sample_tensor.expand(
                batch_dims + orig_shape + (recycling_dim,)
            )
            batch[key] = sample_tensor

            if(key == "no_recycling_iters"):
                no_recycling = sample 
        
        resample_recycling = lambda t: t[..., :no_recycling + 1]
        batch = tensor_tree_map(resample_recycling, batch)

        return batch

    def __iter__(self):
        it = super().__iter__()

        def _batch_prop_gen(iterator):
            for batch in iterator:
                yield self._add_batch_properties(batch)

        return _batch_prop_gen(it)


class OpenFoldDataModule(pl.LightningDataModule):
    def __init__(self,
        config: mlc.ConfigDict,
        template_mmcif_dir: str,
        max_template_date: str,
        train_data_dir: Optional[str] = None,
        train_alignment_dir: Optional[str] = None,
        train_chain_data_cache_path: Optional[str] = None,
        distillation_data_dir: Optional[str] = None,
        distillation_alignment_dir: Optional[str] = None,
        distillation_chain_data_cache_path: Optional[str] = None,
        val_data_dir: Optional[str] = None,
        val_alignment_dir: Optional[str] = None,
        predict_data_dir: Optional[str] = None,
        predict_alignment_dir: Optional[str] = None,
        kalign_binary_path: str = '/usr/bin/kalign',
        train_filter_path: Optional[str] = None,
        distillation_filter_path: Optional[str] = None,
        obsolete_pdbs_file_path: Optional[str] = None,
        template_release_dates_cache_path: Optional[str] = None,
        batch_seed: Optional[int] = None,
        train_epoch_len: int = 50000, 
        _distillation_structure_index_path: Optional[str] = None,
        alignment_index_path: Optional[str] = None,
        distillation_alignment_index_path: Optional[str] = None,
        **kwargs
    ):
        super(OpenFoldDataModule, self).__init__()

        self.config = config
        self.template_mmcif_dir = template_mmcif_dir
        self.max_template_date = max_template_date
        self.train_data_dir = train_data_dir
        self.train_alignment_dir = train_alignment_dir
        self.train_chain_data_cache_path = train_chain_data_cache_path
        self.distillation_data_dir = distillation_data_dir
        self.distillation_alignment_dir = distillation_alignment_dir
        self.distillation_chain_data_cache_path = (
            distillation_chain_data_cache_path
        )
        self.val_data_dir = val_data_dir
        self.val_alignment_dir = val_alignment_dir
        self.predict_data_dir = predict_data_dir
        self.predict_alignment_dir = predict_alignment_dir
        self.kalign_binary_path = kalign_binary_path
        self.train_filter_path = train_filter_path
        self.distillation_filter_path = distillation_filter_path
        self.template_release_dates_cache_path = (
            template_release_dates_cache_path
        )
        self.obsolete_pdbs_file_path = obsolete_pdbs_file_path
        self.batch_seed = batch_seed
        self.train_epoch_len = train_epoch_len

        if(self.train_data_dir is None and self.predict_data_dir is None):
            raise ValueError(
                'At least one of train_data_dir or predict_data_dir must be '
                'specified'
            )

        self.training_mode = self.train_data_dir is not None

        if(self.training_mode and train_alignment_dir is None):
            raise ValueError(
                'In training mode, train_alignment_dir must be specified'
            )
        elif(not self.training_mode and predict_alignment_dir is None):
            raise ValueError(
                'In inference mode, predict_alignment_dir must be specified'
            )      
        elif(val_data_dir is not None and val_alignment_dir is None):
            raise ValueError(
                'If val_data_dir is specified, val_alignment_dir must '
                'be specified as well'
        )

        # An ad-hoc measure for our particular filesystem restrictions
        self._distillation_structure_index = None
        if(_distillation_structure_index_path is not None):
            with open(_distillation_structure_index_path, "r") as fp:
                self._distillation_structure_index = json.load(fp)
        
        self.alignment_index = None
        if(alignment_index_path is not None):
            with open(alignment_index_path, "r") as fp:
                self.alignment_index = json.load(fp)

        self.distillation_alignment_index = None
        if(distillation_alignment_index_path is not None):
            with open(distillation_alignment_index_path, "r") as fp:
                self.distillation_alignment_index = json.load(fp)

    def setup(self):
        # Most of the arguments are the same for the three datasets 
        dataset_gen = partial(OpenFoldSingleDataset,
            template_mmcif_dir=self.template_mmcif_dir,
            max_template_date=self.max_template_date,
            config=self.config,
            kalign_binary_path=self.kalign_binary_path,
            template_release_dates_cache_path=
                self.template_release_dates_cache_path,
            obsolete_pdbs_file_path=
                self.obsolete_pdbs_file_path,
        )

        if(self.training_mode):
            train_dataset = dataset_gen(
                data_dir=self.train_data_dir,
                alignment_dir=self.train_alignment_dir,
                filter_path=self.train_filter_path,
                max_template_hits=self.config.train.max_template_hits,
                shuffle_top_k_prefiltered=
                    self.config.train.shuffle_top_k_prefiltered,
                treat_pdb_as_distillation=False,
                mode="train",
                alignment_index=self.alignment_index,
            )

            distillation_dataset = None
            if(self.distillation_data_dir is not None):
                distillation_dataset = dataset_gen(
                    data_dir=self.distillation_data_dir,
                    alignment_dir=self.distillation_alignment_dir,
                    filter_path=self.distillation_filter_path,
                    max_template_hits=self.config.train.max_template_hits,
                    treat_pdb_as_distillation=True,
                    mode="train",
                    alignment_index=self.distillation_alignment_index,
                    _structure_index=self._distillation_structure_index,
                )

                d_prob = self.config.train.distillation_prob
           
            if(distillation_dataset is not None):
                datasets = [train_dataset, distillation_dataset]
                d_prob = self.config.train.distillation_prob
                probabilities = [1. - d_prob, d_prob]
                chain_data_cache_paths = [
                    self.train_chain_data_cache_path,
                    self.distillation_chain_data_cache_path,
                ]
            else:
                datasets = [train_dataset]
                probabilities = [1.]   
                chain_data_cache_paths = [
                    self.train_chain_data_cache_path,
                ]

            generator = None
            if(self.batch_seed is not None):
                generator = torch.Generator()
                generator = generator.manual_seed(self.batch_seed + 1)
            
            self.train_dataset = OpenFoldDataset(
                datasets=datasets,
                probabilities=probabilities,
                epoch_len=self.train_epoch_len,
                chain_data_cache_paths=chain_data_cache_paths,
                generator=generator,
                _roll_at_init=False,
            )
    
            if(self.val_data_dir is not None):
                self.eval_dataset = dataset_gen(
                    data_dir=self.val_data_dir,
                    alignment_dir=self.val_alignment_dir,
                    filter_path=None,
                    max_template_hits=self.config.eval.max_template_hits,
                    mode="eval",
                )
            else:
                self.eval_dataset = None
        else:           
            self.predict_dataset = dataset_gen(
                data_dir=self.predict_data_dir,
                alignment_dir=self.predict_alignment_dir,
                filter_path=None,
                max_template_hits=self.config.predict.max_template_hits,
                mode="predict",
            )

    def _gen_dataloader(self, stage):
        generator = torch.Generator()
        if(self.batch_seed is not None):
            generator = generator.manual_seed(self.batch_seed)

        dataset = None
        if(stage == "train"):
            dataset = self.train_dataset
            # Filter the dataset, if necessary
            dataset.reroll()
        elif(stage == "eval"):
            dataset = self.eval_dataset
        elif(stage == "predict"):
            dataset = self.predict_dataset
        else:
            raise ValueError("Invalid stage")

        batch_collator = OpenFoldBatchCollator()

        dl = OpenFoldDataLoader(
            dataset,
            config=self.config,
            stage=stage,
            generator=generator,
            batch_size=self.config.data_module.data_loaders.batch_size,
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=batch_collator,
        )

        return dl

    def train_dataloader(self):
        return self._gen_dataloader("train") 

    def val_dataloader(self):
        if(self.eval_dataset is not None):
            return self._gen_dataloader("eval")
        return None

    def predict_dataloader(self):
        return self._gen_dataloader("predict") 


class DummyDataset(torch.utils.data.Dataset):
    def __init__(self, batch_path):
        with open(batch_path, "rb") as f:
            self.batch = pickle.load(f)

    def __getitem__(self, idx):
        return copy.deepcopy(self.batch)

    def __len__(self):
        return 1000


class DummyDataLoader(pl.LightningDataModule):
    def __init__(self, batch_path):
        super().__init__()
        self.dataset = DummyDataset(batch_path)

    def train_dataloader(self):
        return torch.utils.data.DataLoader(self.dataset)
