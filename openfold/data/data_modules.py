from functools import partial
import json
import logging
import os
from typing import Optional, Sequence

import ml_collections as mlc
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
        template_release_dates_cache_path: Optional[str] = None,
        use_small_bfd: bool = True,
        output_raw: bool = False,
        mode: str = "train", 
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
                    or simply {PDB_ID}, each containing:
                        * bfd_uniclust_hits.a3m/small_bfd_hits.sto
                        * mgnify_hits.a3m
                        * pdb70_hits.hhr
                        * uniref90_hits.a3m
                config:
                    A dataset config object. See openfold.config
                mapping_path:
                    A json file containing a mapping from consecutive numerical
                    ids to sample names (matching the directories in data_dir).
                    Samples not in this mapping are ignored. Can be used to 
                    implement the various training-time filters described in
                    the AlphaFold supplement
        """
        super(OpenFoldSingleDataset, self).__init__()
        self.data_dir = data_dir
        self.alignment_dir = alignment_dir
        self.config = config
        self.output_raw = output_raw
        self.mode = mode

        valid_modes = ["train", "val", "predict"]
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
                "scripts/generate_mmcif_caches.py before running OpenFold"
            )

        template_featurizer = templates.TemplateHitFeaturizer(
            mmcif_dir=template_mmcif_dir,
            max_template_date=max_template_date,
            max_hits=max_template_hits,
            kalign_binary_path=kalign_binary_path,
            release_dates_path=template_release_dates_cache_path,
            obsolete_pdbs_path=None,
        )

        self.data_pipeline = data_pipeline.DataPipeline(
            template_featurizer=template_featurizer,
            use_small_bfd=use_small_bfd,
        )

        if(not self.output_raw):
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

        if(self.mode == 'train' or self.mode == 'val'):
            spl = name.rsplit('_', 1)
            if(len(spl) == 2):
                file_id, chain_id = spl
            else:
                file_id, = spl
                chain_id = None

            path = os.path.join(self.data_dir, file_id + '.cif')
            if(os.path.exists(path)):
                data = self._parse_mmcif(
                    path, file_id, chain_id, alignment_dir
                )
            else:
                # Try to search for a distillation PDB file instead
                path = os.path.join(self.data_dir, file_id + '.pdb')
                data = self.data_pipeline.process_pdb(
                    pdb_path=path,
                    alignment_dir=alignment_dir
                )
        else:
            path = os.path.join(name, name + ".fasta")
            data = self.data_pipeline.process_fasta(
                fasta_path=feats,
                alignment_dir=alignment_dir,
            )

        if(self.output_raw):
            return data

        feats = self.feature_pipeline.process_features(
            data, self.mode, "unclamped" 
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
        self.batch_size = batch_size
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
    def __init__(self, config, generator, stage="train"):
        self.config = config
        batch_modes = config.common.batch_modes
        batch_mode_names, batch_mode_probs = list(zip(*batch_modes))
        self.batch_mode_names = batch_mode_names
        self.batch_mode_probs = batch_mode_probs
        self.generator = generator
        self.stage = stage
        
        self.batch_mode_probs_tensor = torch.tensor(self.batch_mode_probs)

        self.feature_pipeline = feature_pipeline.FeaturePipeline(self.config)  

    def __call__(self, raw_prots):
        # We use torch.multinomial here rather than Categorical because the
        # latter doesn't accept a generator for some reason
        batch_mode_idx = torch.multinomial(
            self.batch_mode_probs_tensor, 
            1, 
            generator=self.generator
        ).item() 
        batch_mode_name = self.batch_mode_names[batch_mode_idx]

        processed_prots = []
        for prot in raw_prots:
            features = self.feature_pipeline.process_features(
                prot, self.stage, batch_mode_name
            )
            processed_prots.append(features)

        stack_fn = partial(torch.stack, dim=0)
        return dict_multimap(stack_fn, processed_prots) 


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
        template_release_dates_cache_path: Optional[str] = None, 
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

    def setup(self, stage):
        # Most of the arguments are the same for the three datasets 
        dataset_gen = partial(OpenFoldSingleDataset,
            template_mmcif_dir=self.template_mmcif_dir,
            max_template_date=self.max_template_date,
            config=self.config,
            kalign_binary_path=self.kalign_binary_path,
            template_release_dates_cache_path=
                self.template_release_dates_cache_path,
            use_small_bfd=self.config.data_module.use_small_bfd,
        )

        if(self.training_mode):        
            self.train_dataset = dataset_gen(
                data_dir=self.train_data_dir,
                alignment_dir=self.train_alignment_dir,
                mapping_path=self.train_mapping_path,
                max_template_hits=self.config.train.max_template_hits,
                output_raw=True,
                mode="train",
            )

            if(self.distillation_data_dir is not None):
                distillation_dataset = dataset_gen(
                    data_dir=self.distillation_data_dir,
                    alignment_dir=self.distillation_alignment_dir,
                    mapping_path=self.distillation_mapping_path,
                    max_template_hits=self.train.max_template_hits,
                    output_raw=True,
                    mode="train",
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
                self.val_dataset = dataset_gen(
                    data_dir=self.val_data_dir,
                    alignment_dir=self.val_alignment_dir,
                    mapping_path=None,
                    max_template_hits=self.config.eval.max_template_hits,
                    mode="eval",
                )
        else:           
            self.predict_dataset = dataset_gen(
                data_dir=self.predict_data_dir,
                alignment_dir=self.predict_alignment_dir,
                mapping_path=None,
                max_template_hits=self.config.predict.max_template_hits,
                mode="predict",
            )

        self.batch_collation_seed = torch.Generator().seed()

    def _gen_batch_collator(self, stage):
        """ We want each process to use the same batch collation seed """
        generator = torch.Generator()
        generator = generator.manual_seed(self.batch_collation_seed)
        collate_fn = OpenFoldBatchCollator(
            self.config, generator, stage
        )
        return collate_fn

    def train_dataloader(self):
        return torch.utils.data.DataLoader(
            self.train_dataset,
            batch_size=self.config.data_module.data_loaders.batch_size,
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=self._gen_batch_collator("train"),
        )

    def val_dataloader(self):
        return torch.utils.data.DataLoader(
            self.val_dataset,
            batch_size=self.config.data_module.data_loaders.batch_size,
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=self._gen_batch_collator("eval")
        )

    def predict_dataloader(self):
        return torch.utils.data.DataLoader(
            self.predict_dataset,
            batch_size=self.config.data_module.data_loaders.batch_size,
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=self._gen_batch_collator("eval")
        )
