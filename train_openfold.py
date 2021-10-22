import argparse
import logging
import os

#os.environ["CUDA_VISIBLE_DEVICES"] = "6"
#os.environ["MASTER_ADDR"]="10.119.81.14"
#os.environ["MASTER_PORT"]="42069"
#os.environ["NODE_RANK"]="0"


import random
import time

import numpy as np
import pytorch_lightning as pl
from pytorch_lightning.plugins import DDPPlugin
from pytorch_lightning.plugins.training_type import DeepSpeedPlugin
import torch

from openfold.config import model_config
from openfold.data.data_modules import (
    OpenFoldDataModule,
    DummyDataLoader,
)
from openfold.model.model import AlphaFold
from openfold.utils.exponential_moving_average import ExponentialMovingAverage
from openfold.utils.loss import AlphaFoldLoss
from openfold.utils.seed import seed_everything
from openfold.utils.tensor_utils import tensor_tree_map


class OpenFoldWrapper(pl.LightningModule):
    def __init__(self, config):
        super(OpenFoldWrapper, self).__init__()
        self.config = config
        self.model = AlphaFold(config)
        self.loss = AlphaFoldLoss(config.loss)
        self.ema = ExponentialMovingAverage(
            model=self.model, decay=config.ema.decay
        )

    def forward(self, batch):
        return self.model(batch)

    def training_step(self, batch, batch_idx):
        if(self.ema.device != batch["aatype"].device):
            self.ema.to(batch["aatype"].device)

        # Run the model
        outputs = self(batch)
        
        # Remove the recycling dimension
        batch = tensor_tree_map(lambda t: t[..., -1], batch)

        # Compute loss
        loss = self.loss(outputs, batch)

        return {"loss": loss}

    def validation_step(self, batch, batch_idx):
        # At the start of validation, load the EMA weights
        if(self.cached_weights is None):
            self.cached_weights = model.state_dict()
            self.model.load_state_dict(self.ema.state_dict()["params"])
        
        # Calculate validation loss
        outputs = self(batch)
        batch = tensor_tree_map(lambda t: t[..., -1], batch)
        loss = self.loss(outputs, batch)
        return {"val_loss": loss}

    def validation_epoch_end(self, _):
        # Restore the model weights to normal
        self.model.load_state_dict(self.cached_weights)
        self.cached_weights = None

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

    def on_before_zero_grad(self, *args, **kwargs):
        self.ema.update(self.model)


def main(args):
    if(args.seed is not None):
        seed_everything(args.seed) 

    config = model_config(
        "model_1", 
        train=True, 
        low_prec=(args.precision == 16)
    ) 

    model_module = OpenFoldWrapper(config) 
    #data_module = DummyDataLoader("batch.pickle")
    data_module = OpenFoldDataModule(
        config=config.data, 
        batch_seed=args.seed,
        **vars(args)
    )
    data_module.prepare_data()
    data_module.setup()
    
    plugins = []
    if(args.deepspeed_config_path is not None):
        plugins.append(DeepSpeedPlugin(config=args.deepspeed_config_path))
    
    trainer = pl.Trainer.from_argparse_args(
        args,
        plugins=plugins,
    )

    trainer.fit(model_module, datamodule=data_module)


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
        "--distillation_data_dir", type=str, default=None,
        help="Directory containing training PDB files"
    )
    parser.add_argument(
        "--distillation_alignment_dir", type=str, default=None,
        help="Directory containing precomputed distillation alignments"
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
        "--distillation_mapping_path", type=str, default=None,
        help="""See --train_mapping_path"""
    )
    parser.add_argument(
        "--template_release_dates_cache_path", type=str, default=None,
        help="""Output of scripts/generate_mmcif_cache.py run on template mmCIF
                files."""
    )
    parser.add_argument(
        "--use_small_bfd", type=bool, default=False,
        help="Whether to use a reduced version of the BFD database"
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed"
    )
    parser.add_argument(
        "--deepspeed_config_path", type=str, default=None,
        help="Path to DeepSpeed config. If not provided, DeepSpeed is disabled"
    )
    parser = pl.Trainer.add_argparse_args(parser)
    
    parser.set_defaults(
        num_sanity_val_steps=0,
    )

    args = parser.parse_args()

    main(args)
