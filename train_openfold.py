import argparse
import logging
import os
import sys
import json

import pytorch_lightning as pl
from pytorch_lightning.callbacks.lr_monitor import LearningRateMonitor
from pytorch_lightning.callbacks import DeviceStatsMonitor
from pytorch_lightning.callbacks.model_checkpoint import ModelCheckpoint
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.strategies import DDPStrategy, DeepSpeedStrategy
from pytorch_lightning.plugins.environments import MPIEnvironment
from pytorch_lightning import seed_everything
import torch
import wandb
from deepspeed.utils import zero_to_fp32 

from openfold.config import model_config
from openfold.data.data_modules import OpenFoldDataModule, OpenFoldMultimerDataModule
from openfold.model.model import AlphaFold
from openfold.model.torchscript import script_preset_
from openfold.np import residue_constants
from openfold.utils.callbacks import (
    EarlyStoppingVerbose,
)
from openfold.utils.exponential_moving_average import ExponentialMovingAverage
from openfold.utils.loss import AlphaFoldLoss, lddt_ca
from openfold.utils.lr_schedulers import AlphaFoldLRScheduler
from openfold.utils.multi_chain_permutation import multi_chain_permutation_align
from openfold.utils.superimposition import superimpose
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.validation_metrics import (
    drmsd,
    gdt_ts,
    gdt_ha,
)
from openfold.utils.import_weights import (
    import_jax_weights_,
    import_openfold_weights_
)
from openfold.utils.logger import PerformanceLoggingCallback


class OpenFoldWrapper(pl.LightningModule):
    def __init__(self, config):
        super(OpenFoldWrapper, self).__init__()
        self.config = config
        self.model = AlphaFold(config)
        self.is_multimer = self.config.globals.is_multimer

        self.loss = AlphaFoldLoss(config.loss)

        self.ema = ExponentialMovingAverage(
            model=self.model, decay=config.ema.decay
        )

        self.cached_weights = None
        self.last_lr_step = -1
        self.save_hyperparameters()

    def forward(self, batch):
        return self.model(batch)

    def _log(self, loss_breakdown, batch, outputs, train=True):
        phase = "train" if train else "val"
        for loss_name, indiv_loss in loss_breakdown.items():
            self.log(
                f"{phase}/{loss_name}", 
                indiv_loss, 
                prog_bar=(loss_name == 'loss'),
                on_step=train, on_epoch=(not train), logger=True, sync_dist=False,
            )

            if (train):
                self.log(
                    f"{phase}/{loss_name}_epoch",
                    indiv_loss,
                    on_step=False, on_epoch=True, logger=True, sync_dist=False,
                )

        with torch.no_grad():
            other_metrics = self._compute_validation_metrics(
                batch,
                outputs,
                superimposition_metrics=(not train)
            )

        for k, v in other_metrics.items():
            self.log(
                f"{phase}/{k}",
                torch.mean(v),
                prog_bar = (k == 'loss'),
                on_step=False, on_epoch=True, logger=True, sync_dist=False,
            )

    def training_step(self, batch, batch_idx):
        if (self.ema.device != batch["aatype"].device):
            self.ema.to(batch["aatype"].device)

        ground_truth = batch.pop('gt_features', None)

        # Run the model
        outputs = self(batch)

        # Remove the recycling dimension
        batch = tensor_tree_map(lambda t: t[..., -1], batch)

        if self.is_multimer:
            batch = multi_chain_permutation_align(out=outputs,
                                                  features=batch,
                                                  ground_truth=ground_truth)

        # Compute loss
        loss, loss_breakdown = self.loss(
            outputs, batch, _return_breakdown=True
        )

        # Log it
        self._log(loss_breakdown, batch, outputs)

        return loss

    def on_before_zero_grad(self, *args, **kwargs):
        self.ema.update(self.model)

    def validation_step(self, batch, batch_idx):
        # At the start of validation, load the EMA weights
        if (self.cached_weights is None):
            # model.state_dict() contains references to model weights rather
            # than copies. Therefore, we need to clone them before calling
            # load_state_dict().
            def clone_param(t): return t.detach().clone()
            self.cached_weights = tensor_tree_map(
                clone_param, self.model.state_dict())
            self.model.load_state_dict(self.ema.state_dict()["params"])

        ground_truth = batch.pop('gt_features', None)

        # Run the model
        outputs = self(batch)
        batch = tensor_tree_map(lambda t: t[..., -1], batch)

        batch["use_clamped_fape"] = 0.

        if self.is_multimer:
            batch = multi_chain_permutation_align(out=outputs,
                                                  features=batch,
                                                  ground_truth=ground_truth)

        # Compute loss and other metrics
        _, loss_breakdown = self.loss(
            outputs, batch, _return_breakdown=True
        )

        self._log(loss_breakdown, batch, outputs, train=False)
        
    def on_validation_epoch_end(self):
        # Restore the model weights to normal
        self.model.load_state_dict(self.cached_weights)
        self.cached_weights = None

    def _compute_validation_metrics(self,
                                    batch,
                                    outputs,
                                    superimposition_metrics=False
                                    ):
        metrics = {}

        gt_coords = batch["all_atom_positions"]
        pred_coords = outputs["final_atom_positions"]
        all_atom_mask = batch["all_atom_mask"]

        # This is super janky for superimposition. Fix later
        gt_coords_masked = gt_coords * all_atom_mask[..., None]
        pred_coords_masked = pred_coords * all_atom_mask[..., None]
        ca_pos = residue_constants.atom_order["CA"]
        gt_coords_masked_ca = gt_coords_masked[..., ca_pos, :]
        pred_coords_masked_ca = pred_coords_masked[..., ca_pos, :]
        all_atom_mask_ca = all_atom_mask[..., ca_pos]

        lddt_ca_score = lddt_ca(
            pred_coords,
            gt_coords,
            all_atom_mask,
            eps=self.config.globals.eps,
            per_residue=False,
        )

        metrics["lddt_ca"] = lddt_ca_score

        drmsd_ca_score = drmsd(
            pred_coords_masked_ca,
            gt_coords_masked_ca,
            mask=all_atom_mask_ca,  # still required here to compute n
        )

        metrics["drmsd_ca"] = drmsd_ca_score

        if (superimposition_metrics):
            superimposed_pred, alignment_rmsd = superimpose(
                gt_coords_masked_ca, pred_coords_masked_ca, all_atom_mask_ca,
            )
            gdt_ts_score = gdt_ts(
                superimposed_pred, gt_coords_masked_ca, all_atom_mask_ca
            )
            gdt_ha_score = gdt_ha(
                superimposed_pred, gt_coords_masked_ca, all_atom_mask_ca
            )

            metrics["alignment_rmsd"] = alignment_rmsd
            metrics["gdt_ts"] = gdt_ts_score
            metrics["gdt_ha"] = gdt_ha_score

        return metrics

    def configure_optimizers(self, 
        learning_rate: float = 1e-3,
        eps: float = 1e-5,
    ) -> torch.optim.Adam:
        # Ignored as long as a DeepSpeed optimizer is configured
        optimizer = torch.optim.Adam(
            self.model.parameters(),
            lr=learning_rate,
            eps=eps
        )

        if self.last_lr_step != -1:
            for group in optimizer.param_groups:
                if 'initial_lr' not in group:
                    group['initial_lr'] = learning_rate

        lr_scheduler = AlphaFoldLRScheduler(
            optimizer,
            last_epoch=self.last_lr_step
        )

        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": lr_scheduler,
                "interval": "step",
                "name": "AlphaFoldLRScheduler",
            }
        }

    def on_load_checkpoint(self, checkpoint):
        ema = checkpoint["ema"]
        if (not self.model.template_config.enabled):
            ema["params"] = {k: v for k,
                             v in ema["params"].items() if not "template" in k}
        self.ema.load_state_dict(ema)

    def on_save_checkpoint(self, checkpoint):
        checkpoint["ema"] = self.ema.state_dict()

    def resume_last_lr_step(self, lr_step):
        self.last_lr_step = lr_step

    def load_from_jax(self, jax_path):
        model_basename = os.path.splitext(
            os.path.basename(
                os.path.normpath(jax_path)
            )
        )[0]
        model_version = "_".join(model_basename.split("_")[1:])
        import_jax_weights_(
            self.model, jax_path, version=model_version
        )

def get_model_state_dict_from_ds_checkpoint(checkpoint_dir):
    latest_path = os.path.join(checkpoint_dir, 'latest')
    if os.path.isfile(latest_path):
        with open(latest_path, 'r') as fd:
            tag = fd.read().strip()
    else:
        raise ValueError(f"Unable to find 'latest' file at {latest_path}")

    ds_checkpoint_dir = os.path.join(checkpoint_dir, tag)
    _DS_CHECKPOINT_VERSION = 2  # based on manual parsing of checkpoint files
    state_file = zero_to_fp32.get_model_state_file(ds_checkpoint_dir, _DS_CHECKPOINT_VERSION)
    return torch.load(state_file)

def main(args):
    if(args.seed is not None):
        seed_everything(args.seed, workers=True) 

    is_low_precision = args.precision in [
        "bf16-mixed", "16", "bf16", "16-true", "16-mixed", "bf16-mixed"]

    config = model_config(
        args.config_preset, 
        train=True, 
        low_prec=is_low_precision,
    ) 
    if args.experiment_config_json: 
        with open(args.experiment_config_json, 'r') as f:
            custom_config_dict = json.load(f)
        config.update_from_flattened_dict(custom_config_dict)

    model_module = OpenFoldWrapper(config)

    if args.resume_from_ckpt:
        if args.resume_model_weights_only:
            # Load the checkpoint
            if os.path.isdir(args.resume_from_ckpt):
                sd = zero_to_fp32.get_fp32_state_dict_from_zero_checkpoint(
                    args.resume_from_ckpt)
            else:
                sd = torch.load(args.resume_from_ckpt)
            # Process the state dict
            if 'module' in sd:
                sd = {k[len('module.'):]: v for k, v in sd['module'].items()}
                import_openfold_weights_(model=model_module, state_dict=sd)
            elif 'state_dict' in sd:
                import_openfold_weights_(
                    model=model_module, state_dict=sd['state_dict'])
            else:
                # Loading from pre-trained model
                sd = {'model.'+k: v for k, v in sd.items()}
                import_openfold_weights_(model=model_module, state_dict=sd)
            logging.info("Successfully loaded model weights...")

        else:  # Loads a checkpoint to start from a specific time step
            if os.path.isdir(args.resume_from_ckpt):
                sd = get_model_state_dict_from_ds_checkpoint(args.resume_from_ckpt)
            else:
                sd = torch.load(args.resume_from_ckpt)
            last_global_step = int(sd['global_step'])
            model_module.resume_last_lr_step(last_global_step)
            logging.info("Successfully loaded last lr step...")

    if args.resume_from_jax_params:
        model_module.load_from_jax(args.resume_from_jax_params)
        logging.info(
            f"Successfully loaded JAX parameters at {args.resume_from_jax_params}...")

    # TorchScript components of the model
    if (args.script_modules):
        script_preset_(model_module)

    if "multimer" in args.config_preset:
        data_module = OpenFoldMultimerDataModule(
            config=config.data,
            batch_seed=args.seed,
            **vars(args)
        )
    else:
        data_module = OpenFoldDataModule(
            config=config.data,
            batch_seed=args.seed,
            **vars(args)
        )

    data_module.prepare_data()
    data_module.setup()

    callbacks = []
    if (args.checkpoint_every_epoch):
        mc = ModelCheckpoint(
            every_n_epochs=1,
            auto_insert_metric_name=False,
            save_top_k=-1,
        )
        callbacks.append(mc)

    if (args.early_stopping):
        es = EarlyStoppingVerbose(
            monitor="val/lddt_ca",
            min_delta=args.min_delta,
            patience=args.patience,
            verbose=False,
            mode="max",
            check_finite=True,
            strict=True,
        )
        callbacks.append(es)

    if (args.log_performance):
        global_batch_size = args.num_nodes * args.gpus
        perf = PerformanceLoggingCallback(
            log_file=os.path.join(args.output_dir, "performance_log.json"),
            global_batch_size=global_batch_size,
        )
        callbacks.append(perf)

    if (args.log_lr):
        lr_monitor = LearningRateMonitor(logging_interval="step")
        callbacks.append(lr_monitor)

    loggers = []
    is_rank_zero = args.mpi_plugin and (int(os.environ.get("PMI_RANK")) == 0)
    if(args.wandb):
        if args.mpi_plugin and is_rank_zero:
            wandb_init_dict = dict(
                name=args.experiment_name,
                project=args.wandb_project,
                id=args.wandb_id,
                dir=args.output_dir,
                resume="allow",
                anonymous=None,
                entity=args.wandb_entity
            )
            wandb.run = wandb.init(**wandb_init_dict)

        wdb_logger = WandbLogger(
            name=args.experiment_name,
            save_dir=args.output_dir,
            id=args.wandb_id,
            project=args.wandb_project,
            **{"entity": args.wandb_entity}
        )
        loggers.append(wdb_logger)

    cluster_environment = MPIEnvironment() if args.mpi_plugin else None
    if(args.deepspeed_config_path is not None):
        strategy = DeepSpeedStrategy(
            config=args.deepspeed_config_path,
            cluster_environment=cluster_environment,
        )
        if(args.wandb and is_rank_zero):
            wdb_logger.experiment.save(args.deepspeed_config_path)
            wdb_logger.experiment.save("openfold/config.py")
    elif (args.gpus is not None and args.gpus > 1) or args.num_nodes > 1:
        strategy = DDPStrategy(find_unused_parameters=False,
                               cluster_environment=cluster_environment)
    else:
        strategy = None
 
    if(args.wandb and is_rank_zero):
        freeze_path = f"{wdb_logger.experiment.dir}/package_versions.txt"
        os.system(f"{sys.executable} -m pip freeze > {freeze_path}")
        wdb_logger.experiment.save(f"{freeze_path}")

    trainer_kws = ['num_nodes', 'precision', 'max_epochs', 'log_every_n_steps',
                   'flush_logs_ever_n_steps', 'num_sanity_val_steps', 'reload_dataloaders_every_n_epochs']
    trainer_args = {k: v for k, v in vars(args).items() if k in trainer_kws}
    trainer_args.update({
        'default_root_dir': args.output_dir,
        'strategy': strategy,
        'callbacks': callbacks,
        'logger': loggers,
    })
    trainer = pl.Trainer(**trainer_args)


    if (args.resume_model_weights_only):
        ckpt_path = None
    else:
        ckpt_path = args.resume_from_ckpt

    trainer.fit(
        model_module,
        datamodule=data_module,
        ckpt_path=ckpt_path,
    )


def bool_type(bool_str: str):
    bool_str_lower = bool_str.lower()
    if bool_str_lower in ('false', 'f', 'no', 'n', '0'):
        return False
    elif bool_str_lower in ('true', 't', 'yes', 'y', '1'):
        return True
    else:
        raise ValueError(f'Cannot interpret {bool_str} as bool')


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
        "output_dir", type=str,
        help='''Directory in which to output checkpoints, logs, etc. Ignored
                if not on rank 0'''
    )
    parser.add_argument(
        "max_template_date", type=str,
        help='''Cutoff for all templates. In training mode, templates are also 
                filtered by the release date of the target'''
    )
    parser.add_argument(
        "--train_mmcif_data_cache_path", type=str, default=None,
        help="Path to the json file which records all the information of mmcif structures used during training"
    )
    parser.add_argument(
        "--use_single_seq_mode", type=str, default=False,
        help="Use single sequence embeddings instead of MSAs."
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
        "--val_mmcif_data_cache_path", type=str, default=None,
        help="path to the json file which records all the information of mmcif structures used during validation"
    )
    parser.add_argument(
        "--kalign_binary_path", type=str, default='/usr/bin/kalign',
        help="Path to the kalign binary"
    )
    parser.add_argument(
        "--train_filter_path", type=str, default=None,
        help='''Optional path to a text file containing names of training
                examples to include, one per line. Used to filter the training 
                set'''
    )
    parser.add_argument(
        "--distillation_filter_path", type=str, default=None,
        help="""See --train_filter_path"""
    )
    parser.add_argument(
        "--obsolete_pdbs_file_path", type=str, default=None,
        help="""Path to obsolete.dat file containing list of obsolete PDBs and 
             their replacements."""
    )
    parser.add_argument(
        "--template_release_dates_cache_path", type=str, default=None,
        help="""Output of scripts/generate_mmcif_cache.py run on template mmCIF
                files."""
    )
    parser.add_argument(
        "--use_small_bfd", type=bool_type, default=False,
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
    parser.add_argument(
        "--checkpoint_every_epoch", action="store_true", default=False,
        help="""Whether to checkpoint at the end of every training epoch"""
    )
    parser.add_argument(
        "--early_stopping", type=bool_type, default=False,
        help="Whether to stop training when validation loss fails to decrease"
    )
    parser.add_argument(
        "--min_delta", type=float, default=0,
        help="""The smallest decrease in validation loss that counts as an 
                improvement for the purposes of early stopping"""
    )
    parser.add_argument(
        "--patience", type=int, default=3,
        help="Early stopping patience"
    )
    parser.add_argument(
        "--resume_from_ckpt", type=str, default=None,
        help="Path to a model checkpoint from which to restore training state"
    )
    parser.add_argument(
        "--resume_model_weights_only", type=bool_type, default=False,
        help="Whether to load just model weights as opposed to training state"
    )
    parser.add_argument(
        "--resume_from_jax_params", type=str, default=None,
        help="""Path to an .npz JAX parameter file with which to initialize the model"""
    )
    parser.add_argument(
        "--log_performance", type=bool_type, default=False,
        help="Measure performance"
    )
    parser.add_argument(
        "--wandb", action="store_true", default=False,
        help="Whether to log metrics to Weights & Biases"
    )
    parser.add_argument(
        "--experiment_name", type=str, default=None,
        help="Name of the current experiment. Used for wandb logging"
    )
    parser.add_argument(
        "--wandb_id", type=str, default=None,
        help="ID of a previous run to be resumed"
    )
    parser.add_argument(
        "--wandb_project", type=str, default=None,
        help="Name of the wandb project to which this run will belong"
    )
    parser.add_argument(
        "--wandb_entity", type=str, default=None,
        help="wandb username or team name to which runs are attributed"
    )
    parser.add_argument(
        "--script_modules", type=bool_type, default=False,
        help="Whether to TorchScript eligible components of them model"
    )
    parser.add_argument(
        "--train_chain_data_cache_path", type=str, default=None,
    )
    parser.add_argument(
        "--distillation_chain_data_cache_path", type=str, default=None,
    )
    parser.add_argument(
        "--train_epoch_len", type=int, default=10000,
        help=(
            "The virtual length of each training epoch. Stochastic filtering "
            "of training data means that training datasets have no "
            "well-defined length. This virtual length affects frequency of "
            "validation & checkpointing (by default, one of each per epoch)."
        )
    )
    parser.add_argument(
        "--log_lr", action="store_true", default=False,
        help="Whether to log the actual learning rate"
    )
    parser.add_argument(
        "--config_preset", type=str, default="initial_training",
        help=(
            'Config setting. Choose e.g. "initial_training", "finetuning", '
            '"model_1", etc. By default, the actual values in the config are '
            'used.'
        )
    )
    parser.add_argument(
        "--_distillation_structure_index_path", type=str, default=None,
    )
    parser.add_argument(
        "--alignment_index_path", type=str, default=None,
        help="Training alignment index. See the README for instructions."
    )
    parser.add_argument(
        "--distillation_alignment_index_path", type=str, default=None,
        help="Distillation alignment index. See the README for instructions."
    )
    parser.add_argument(
        "--experiment_config_json", default="", help="Path to a json file with custom config values to overwrite config setting",
    )
    parser.add_argument(
        "--gpus", type=int, default=1, help='For determining optimal strategy and effective batch size.'
    )
    parser.add_argument("--mpi_plugin", action="store_true", default=False,
                        help="Whether to use MPI for parallele processing")

    trainer_group = parser.add_argument_group(
        'Arguments to pass to PyTorch Lightning Trainer')
    trainer_group.add_argument(
        "--num_nodes", type=int, default=1,
    )
    trainer_group.add_argument(
        "--precision", type=str, default='bf16',
        help='Sets precision, lower precision improves runtime performance.',
    )
    trainer_group.add_argument(
        "--max_epochs", type=int, default=1,
    )
    trainer_group.add_argument(
        "--log_every_n_steps", type=int, default=25,
    )
    trainer_group.add_argument(
        "--flush_logs_every_n_steps", type=int, default=5,
    )
    trainer_group.add_argument(
        "--num_sanity_val_steps", type=int, default=0,
    )
    trainer_group.add_argument(
        "--reload_dataloaders_every_n_epochs", type=int, default=1,
    )
    trainer_group.add_argument(
        "--accumulate_grad_batches", type=int, default=1,
        help="Accumulate gradients over k batches before next optimizer step.")

    args = parser.parse_args()

    if (args.seed is None and
        ((args.gpus is not None and args.gpus > 1) or
         (args.num_nodes is not None and args.num_nodes > 1))):
        raise ValueError("For distributed training, --seed must be specified")

    if (str(args.precision) == "16" and args.deepspeed_config_path is not None):
        raise ValueError("DeepSpeed and FP16 training are not compatible")

    if (args.resume_from_jax_params is not None and args.resume_from_ckpt is not None):
        raise ValueError(
            "Choose between loading pretrained Jax-weights and a checkpoint-path")

    main(args)
