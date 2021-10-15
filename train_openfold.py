import argparse
from functools import partial
import json
import logging
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "4"

import time
from typing import Optional

import ml_collections as mlc
import pytorch_lightning as pl
from pytorch_lightning.plugins.training_type import DeepSpeedPlugin
import torch
from torch.utils.data import RandomSampler

torch.manual_seed(42)

from openfold.config import model_config
from openfold.model.model import AlphaFold
from openfold.features import (
    data_pipeline,
    feature_pipeline,
    mmcif_parsing,
)
from openfold.features import templates
from openfold.features.np.utils import to_date
from openfold.utils.exponential_moving_average import ExponentialMovingAverage
from openfold.utils.loss import AlphaFoldLoss
from openfold.utils.tensor_utils import tensor_tree_map, dict_multimap


class OpenFoldDataset(torch.utils.data.Dataset):
    def __init__(self,
        data_dir: str,
        alignment_dir: str, 
        template_mmcif_dir: str,
        max_template_date: str,
        config: mlc.ConfigDict,
        kalign_binary_path: str = '/usr/bin/kalign',
        mapping_path: Optional[str] = None,
        mmcif_cache_dir: str = 'tmp/', 
        use_small_bfd: bool = True,
        seed: int = 42,
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
        super(OpenFoldDataset, self).__init__()
        self.data_dir = data_dir
        self.alignment_dir = alignment_dir
        self.config = config
        self.seed = seed
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

        template_release_dates_path = os.path.join(
            mmcif_cache_dir, "template_release_dates.json"
        )
        if(not os.path.exists(template_release_dates_path)):
            logging.warning(
                "Template release dates cache does not exist. Remember to run "
                "scripts/generate_mmcif_caches.py before running OpenFold"
            )
            template_release_dates_path = None

        template_featurizer = templates.TemplateHitFeaturizer(
            mmcif_dir=template_mmcif_dir,
            max_template_date=max_template_date,
            max_hits=(20 if (mode == 'train') else 4),
            kalign_binary_path=kalign_binary_path,
            release_dates_path=template_release_dates_path,
            obsolete_pdbs_path=None,
        )

        self.data_pipeline = data_pipeline.DataPipeline(
            template_featurizer=template_featurizer,
            use_small_bfd=use_small_bfd,
        )

        self.feature_pipeline = feature_pipeline.FeaturePipeline(config) 

    def __getitem__(self, idx):
        no_batch_modes = len(self.config.common.batch_modes)
        batch_mode_idx = idx % no_batch_modes
        batch_mode_str = self.config.common.batch_modes[batch_mode_idx][0]
        idx = int(idx / no_batch_modes)

        name = self.mapping[str(idx)]

        if(self.mode == 'train' or self.mode == 'val'):
            spl = name.rsplit('_', 1)
            if(len(spl) == 2):
                file_id, chain_id = spl
            else:
                file_id, = spl
                chain_id = None

            path = os.path.join(self.data_dir, file_id + '.cif')
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

            alignment_dir = os.path.join(self.alignment_dir, name)

            data = self.data_pipeline.process_mmcif(
                mmcif=mmcif_object,
                alignment_dir=alignment_dir,
                chain_id=chain_id,
            )

        else:
            path = os.path.join(name, name + '.fasta')

            data = self.data_pipeline.process_fasta(
                fasta_path = feats,
                alignment_dir = alignment_dir,
            )

        feats = self.feature_pipeline.process_features(
            data, self.mode, batch_mode_str
        )

        return feats

    def __len__(self):
        return len(self.mapping.keys()) 


class OpenFoldBatchSampler(torch.utils.data.BatchSampler):
    """
        A shameful hack.

        In AlphaFold, certain batches are designated for loss clamping. The
        exact method by residue cropping withing that batch is performed 
        depends on that designation.

        In idiomatic PyTorch, such "batch-wide" properties generally do not 
        exist; samples are supposed to be generated independently and only
        later batched. This class and OpenFoldDataset get around this design
        limitation by encoding batch properties in the indices sent to the
        Dataset.

        While this works (and efficiently), it precludes the future use of an
        IterableDataset (such as WebDataset), which doesn't use indices. In
        that case, the same can be accomplished by delaying the feature
        processing step to the collate_fn, an argument of the DataLoader. That
        solution is avoided here because it requires loading an entire batch's
        worth of uncropped features into memory at a time.

        A third option would be to generate two separate Dataset objects, one
        that generates "clamped" batches and another for "unclamped" ones.
        However, this would require parsing the precomputed caches of most
        proteins twice, once for each loader. Given how lopsided the chances of
        drawing a "clamped" batch are, care would also have to be taken not
        to allocate too many resources to the less used DataLoader.
    """
    def __init__(self, config, **kwargs):
        super(OpenFoldBatchSampler, self).__init__(**kwargs)
        self.config = config
        self.no_batch_modes = len(self.config.common.batch_modes)

    def __iter__(self):
        it = super().__iter__()
        distr = torch.distributions.categorical.Categorical(
            torch.tensor(
                [prob for name, prob in self.config.common.batch_modes]
            )
        )
        for sample in it:
            mode_idx = distr.sample().item()
            sample = [s * self.no_batch_modes + mode_idx for s in sample]
            yield sample


class OpenFoldDataModule(pl.LightningDataModule):
    def __init__(self,
        config: mlc.ConfigDict,
        template_mmcif_dir: str,
        max_template_date: str,
        train_data_dir: Optional[str] = None,
        train_alignment_dir: Optional[str] = None,
        val_data_dir: Optional[str] = None,
        val_alignment_dir: Optional[str] = None,
        predict_data_dir: Optional[str] = None,
        predict_alignment_dir: Optional[str] = None,
        kalign_binary_path: str = '/usr/bin/kalign',
        train_mapping_path: Optional[str] = None,
        mmcif_cache_dir: str = 'tmp/', 
        **kwargs
    ):
        super(OpenFoldDataModule, self).__init__()

        self.config = config
        self.template_mmcif_dir = template_mmcif_dir
        self.max_template_date = max_template_date
        self.train_data_dir = train_data_dir
        self.train_alignment_dir = train_alignment_dir
        self.val_data_dir = val_data_dir
        self.val_alignment_dir = val_alignment_dir
        self.predict_data_dir = predict_data_dir
        self.predict_alignment_dir = predict_alignment_dir
        self.kalign_binary_path = kalign_binary_path
        self.train_mapping_path = train_mapping_path
        self.mmcif_cache_dir = mmcif_cache_dir

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
        dataset_gen = partial(OpenFoldDataset,
            template_mmcif_dir=self.template_mmcif_dir,
            max_template_date=self.max_template_date,
            config=self.config,
            kalign_binary_path=self.kalign_binary_path,
            mmcif_cache_dir=self.mmcif_cache_dir,
            use_small_bfd=self.config.data_module.use_small_bfd,
        )

        if(self.training_mode):        
            self.train_dataset = dataset_gen(
                data_dir=self.train_data_dir,
                alignment_dir=self.train_alignment_dir,
                mapping_path=self.train_mapping_path,
                mode='train',
            )
    
            if(self.val_data_dir is not None):
                self.val_dataset = dataset_gen(
                    data_dir=self.val_data_dir,
                    alignment_dir=self.val_alignment_dir,
                    mapping_path=None,
                    mode='val',
                )
        else:           
            self.predict_dataset = dataset_gen(
                data_dir=self.predict_data_dir,
                alignment_dir=self.predict_alignment_dir,
                mapping_path=None,
                mode='predict',
            )

    def train_dataloader(self):
        stack_fn = partial(torch.stack, dim=0)
        stack = lambda l: dict_multimap(stack_fn, l)
        return torch.utils.data.DataLoader(
            self.train_dataset,
            batch_sampler=OpenFoldBatchSampler(
                config=self.config, 
                sampler=RandomSampler(self.train_dataset), 
                batch_size=self.config.data_module.data_loaders.batch_size,
                drop_last=False,
            ),
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=stack,
        )

    def val_dataloader(self):
        stack_fn = partial(torch.stack, dim=0)
        stack = lambda l: dict_multimap(stack_fn, l)
        return torch.utils.data.DataLoader(
            self.val_dataset,
            batch_size=self.config.data_module.data_loaders.batch_size,
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=stack
        )

    def predict_dataloader(self):
        stack_fn = partial(torch.stack, dim=0)
        stack = lambda l: dict_multimap(stack_fn, l)
        return torch.utils.data.DataLoader(
            self.predict_dataset,
            batch_size=self.config.data_module.data_loaders.batch_size,
            num_workers=self.config.data_module.data_loaders.num_workers,
            collate_fn=stack
        )


class OpenFoldWrapper(pl.LightningModule):
    def __init__(self, config):
        super(OpenFoldWrapper, self).__init__()
        self.config = config
        self.model = AlphaFold(config.model)
        self.loss = AlphaFoldLoss(config.loss)
        self.ema = ExponentialMovingAverage(self.model, decay=config.ema.decay)

    def forward(self, batch):
        return self.model(batch)

    def training_step(self, batch, batch_idx):
        # Run the model
        outputs = self(batch)
        
        # Remove the recycling dimension
        batch = tensor_tree_map(lambda t: t[..., -1], batch)

        # Compute loss
        loss = self.loss(outputs, batch)

        return {"loss": loss, "pred": outputs["sm"]["positions"][-1].detach()}

    def training_epoch_end(self, outs):
        out = outs[-1]["pred"].cpu()
        with open("prediction/preds_" + str(time.strftime("%H:%M:%S")) + ".pickle", "wb") as f:
            pickle.dump(out, f, protocol=pickle.HIGHEST_PROTOCOL)

    def configure_optimizers(self, 
        learning_rate: float = 1e-3,
        eps: float = 1e-8
    ) -> torch.optim.Adam:
        # Ignored as long as a DeepSpeed optimizer is configured
        return torch.optim.Adam(
            self.model.parameters(), 
            lr=learning_rate, 
            eps=eps
        )


def main(args):
    config = model_config(
        "model_1", 
        train=True, 
        low_prec=(args.precision == 16)
    )

    plugins = []
    #plugins.append(DeepSpeedPlugin(config="deepspeed_config.json"))

    trainer = pl.Trainer.from_argparse_args(
        args,
        plugins=plugins,
    )

    model_module = OpenFoldWrapper(config) 
    data_module = OpenFoldDataModule(config=config.data, **vars(args))
    
    trainer.fit(model_module, data_module)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "train_data_dir", type=str,
        help="Directory containing training mmCIF files"
    )
    parser.add_argument(
        "train_alignment_dir", type=str,
        help="Directory containing precomputed training alignments"
    )
    parser.add_argument(
        "template_mmcif_dir", type=str,
        help="Directory containing mmCIF files to search for templates"
    )
    parser.add_argument(
        "max_template_date", type=str,
        help="""Cutoff for all templates. In training mode, templates are also 
                filtered by the release date of the target"""
    )
    parser.add_argument(
        "--val_data_dir", type=str, default=None,
        help="Directory containing validation mmCIF files"
    )
    parser.add_argument(
        "--val_alignment_dir", type=str, default=None,
        help="Directory containing precomputed validation alignments"
    )
    parser.add_argument(
        "--kalign_binary_path", type=str, default='/usr/bin/kalign',
        help="Path to the kalign binary"
    )
    parser.add_argument(
        "--train_mapping_path", type=str, default=None,
        help="""Optional path to a .json file containing a mapping from
                consecutive numerical indices to sample names. Used to filter
                the training set"""
    )
    parser.add_argument(
        "--mmcif_cache_dir", type=str, default="tmp/",
        help="Directory containing precomputed mmCIF metadata"
    )
    parser.add_argument(
        "--use_small_bfd", type=bool, default=False,
        help="Whether to use a reduced version of the BFD database"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for the DataModule"
    )
    parser = pl.Trainer.add_argparse_args(parser)
    
    parser.set_defaults(
        num_sanity_val_steps=0,
    )

    args = parser.parse_args()

    # Seed torch
    torch.manual_seed(args.seed)

    main(args)
