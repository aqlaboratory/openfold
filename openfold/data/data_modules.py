import copy
from functools import partial
import json
import logging
import os
import pickle
from typing import Optional, Sequence

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
        mapping_path: Optional[str] = None,
        max_template_hits: int = 4,
        obsolete_pdbs_file_path: Optional[str] = None,
        template_release_dates_cache_path: Optional[str] = None,
        shuffle_top_k_prefiltered: Optional[int] = None,
        treat_pdb_as_distillation: bool = True,
        mode: str = "train", 
        _output_raw: bool = False,
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
                mapping_path:
                    A json file containing a mapping from consecutive numerical
                    ids to sample names (matching the directories in data_dir).
                    Samples not in this mapping are ignored. Can be used to 
                    implement the various training-time filters described in
                    the AlphaFold supplement.
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
        self._output_raw = _output_raw

        valid_modes = ["train", "eval", "predict"]
        if(mode not in valid_modes):
            raise ValueError(f'mode must be one of {valid_modes}')

        if(mapping_path is None):
            self.mapping = {
                str(i):os.path.splitext(name)[0] 
                for i, name in enumerate(os.listdir(alignment_dir))
            }
        else:
            with open(mapping_path, 'r') as fp:
                self.mapping = json.load(fp)

        if(template_release_dates_cache_path is None):
            logging.warning(
                "Template release dates cache does not exist. Remember to run "
                "scripts/generate_mmcif_cache.py before running OpenFold"
            )

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

    def _parse_mmcif(self, path, file_id, chain_id, alignment_dir):
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
        )

        return data
    
    def __getitem__(self, idx):
        name = self.mapping[str(idx)]
        alignment_dir = os.path.join(self.alignment_dir, name)

        if(self.mode == 'train' or self.mode == 'eval'):
            spl = name.rsplit('_', 1)
            if(len(spl) == 2):
                file_id, chain_id = spl
            else:
                file_id, = spl
                chain_id = None

            path = os.path.join(self.data_dir, file_id)
            if(os.path.exists(path + ".cif")):
                data = self._parse_mmcif(
                    path + ".cif", file_id, chain_id, alignment_dir
                )
            elif(os.path.exists(path + ".core")):
                data = self.data_pipeline.process_core(
                    path + ".core", alignment_dir
                )
            elif(os.path.exists(path + ".pdb")):
                data = self.data_pipeline.process_pdb(
                    pdb_path=path + ".pdb",
                    alignment_dir=alignment_dir,
                    is_distillation=self.treat_pdb_as_distillation,
                    chain_id=chain_id,
                )
            else:
                raise ValueError("Invalid file type")
        else:
            path = os.path.join(name, name + ".fasta")
            data = self.data_pipeline.process_fasta(
                fasta_path=path,
                alignment_dir=alignment_dir,
            )

        if(self._output_raw):
            return data

        feats = self.feature_pipeline.process_features(
            data, self.mode 
        )

        return feats

    def __len__(self):
        return len(self.mapping.keys()) 


def looped_sequence(sequence):
    while True:
        for x in sequence:
            yield x


class OpenFoldDataset(torch.utils.data.IterableDataset):
    """
        The Dataset is written to accommodate the requirement that proteins are
        sampled from the distillation set with some probability p
        and from the PDB set with probability (1 - p). Proteins are sampled
        from both sets without replacement, and as soon as either set is
        emptied, it is refilled. The Dataset therefore has an arbitrary length.
        Nevertheless, for compatibility with various PyTorch Lightning
        functionalities, it is possible to specify an epoch length. This length
        has no effect on the output of the Dataset.
    """
    def __init__(self,
        datasets: Sequence[OpenFoldSingleDataset],
        probabilities: Sequence[int],
        epoch_len: int,
    ):
        self.datasets = datasets
        self.samplers = [
            looped_sequence(RandomSampler(d)) for d in datasets
        ]
        self.epoch_len = epoch_len

        self.distr = torch.distributions.categorical.Categorical(
            probs=torch.tensor(probabilities),
        )

    def __iter__(self):
        return self

    def __next__(self):
        dataset_idx = self.distr.sample()
        sampler = self.samplers[dataset_idx]
        element_idx = next(sampler)
        return self.datasets[dataset_idx][element_idx] 

    def __len__(self):
        return self.epoch_len


class OpenFoldBatchCollator:
    def __init__(self, config, stage="train"):
        self.stage = stage
        self.feature_pipeline = feature_pipeline.FeaturePipeline(config)

    def __call__(self, raw_prots):
        processed_prots = []
        for prot in raw_prots:
            features = self.feature_pipeline.process_features(
                prot, self.stage
            )
            processed_prots.append(features)

        stack_fn = partial(torch.stack, dim=0)
        return dict_multimap(stack_fn, processed_prots) 


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
        distillation_data_dir: Optional[str] = None,
        distillation_alignment_dir: Optional[str] = None,
        val_data_dir: Optional[str] = None,
        val_alignment_dir: Optional[str] = None,
        predict_data_dir: Optional[str] = None,
        predict_alignment_dir: Optional[str] = None,
        kalign_binary_path: str = '/usr/bin/kalign',
        train_mapping_path: Optional[str] = None,
        distillation_mapping_path: Optional[str] = None,
        obsolete_pdbs_file_path: Optional[str] = None,
        template_release_dates_cache_path: Optional[str] = None,
        batch_seed: Optional[int] = None,
        **kwargs
    ):
        super(OpenFoldDataModule, self).__init__()

        self.config = config
        self.template_mmcif_dir = template_mmcif_dir
        self.max_template_date = max_template_date
        self.train_data_dir = train_data_dir
        self.train_alignment_dir = train_alignment_dir
        self.distillation_data_dir = distillation_data_dir
        self.distillation_alignment_dir = distillation_alignment_dir
        self.val_data_dir = val_data_dir
        self.val_alignment_dir = val_alignment_dir
        self.predict_data_dir = predict_data_dir
        self.predict_alignment_dir = predict_alignment_dir
        self.kalign_binary_path = kalign_binary_path
        self.train_mapping_path = train_mapping_path
        self.distillation_mapping_path = distillation_mapping_path
        self.template_release_dates_cache_path = (
            template_release_dates_cache_path
        )
        self.obsolete_pdbs_file_path = obsolete_pdbs_file_path
        self.batch_seed = batch_seed

        if(self.train_data_dir is None and self.predict_data_dir is None):
            raise ValueError(
                'At least one of train_data_dir or predict_data_dir must be '
                'specified'
            )

        self.training_mode = self.train_data_dir is not None

        if(self.training_mode and self.train_alignment_dir is None):
            raise ValueError(
                'In training mode, train_alignment_dir must be specified'
            )
        elif(not self.training_mode and self.predict_alingment_dir is None):
            raise ValueError(
                'In inference mode, predict_alignment_dir must be specified'
            )      
        elif(val_data_dir is not None and val_alignment_dir is None):
            raise ValueError(
                'If val_data_dir is specified, val_alignment_dir must '
                'be specified as well'
        )

    def setup(self, stage: Optional[str] = None):
        if(stage is None):
            stage = "train"

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
            self.train_dataset = dataset_gen(
                data_dir=self.train_data_dir,
                alignment_dir=self.train_alignment_dir,
                mapping_path=self.train_mapping_path,
                max_template_hits=self.config.train.max_template_hits,
                shuffle_top_k_prefiltered=
                    self.config.train.shuffle_top_k_prefiltered,
                treat_pdb_as_distillation=False,
                mode="train",
                _output_raw=True,
            )

            if(self.distillation_data_dir is not None):
                distillation_dataset = dataset_gen(
                    data_dir=self.distillation_data_dir,
                    alignment_dir=self.distillation_alignment_dir,
                    mapping_path=self.distillation_mapping_path,
                    max_template_hits=self.train.max_template_hits,
                    treat_pdb_as_distillation=True,
                    mode="train",
                    _output_raw=True,
                )

                d_prob = self.config.train.distillation_prob
                self.train_dataset = OpenFoldDataset(
                    datasets=[self.train_dataset, distillation_dataset],
                    probabilities=[1 - d_prob, d_prob],
                    epoch_len=(
                        self.train_dataset.len() + distillation_dataset.len()
                    ),
                )
    
            if(self.val_data_dir is not None):
                self.eval_dataset = dataset_gen(
                    data_dir=self.val_data_dir,
                    alignment_dir=self.val_alignment_dir,
                    mapping_path=None,
                    max_template_hits=self.config.eval.max_template_hits,
                    mode="eval",
                    _output_raw=True,
                )
            else:
                self.eval_dataset = None
        else:           
            self.predict_dataset = dataset_gen(
                data_dir=self.predict_data_dir,
                alignment_dir=self.predict_alignment_dir,
                mapping_path=None,
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
        elif(stage == "eval"):
            dataset = self.eval_dataset
        elif(stage == "predict"):
            dataset = self.predict_dataset
        else:
            raise ValueError("Invalid stage")

        batch_collator = OpenFoldBatchCollator(self.config, stage)

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
