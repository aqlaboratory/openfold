import argparse
import os
import logging
import random

import numpy
import torch
from openfold.config import model_config
from openfold.data import feature_pipeline
from openfold.data.data_pipeline import make_sequence_features_with_custom_template
from openfold.np import protein
from openfold.utils.script_utils import load_models_from_command_line, parse_fasta, run_model, prep_output, \
    relax_protein
from openfold.utils.tensor_utils import (
    tensor_tree_map,
)
from scripts.utils import add_data_args

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

torch_versions = torch.__version__.split(".")
torch_major_version = int(torch_versions[0])
torch_minor_version = int(torch_versions[1])
if(
    torch_major_version > 1 or
    (torch_major_version == 1 and torch_minor_version >= 12)
):
    # Gives a large speedup on Ampere-class GPUs
    torch.set_float32_matmul_precision("high")

torch.set_grad_enabled(False)


def main(args):
    os.makedirs(args.output_dir, exist_ok=True)

    config = model_config(args.config_preset)

    random_seed = args.data_random_seed
    if random_seed is None:
        random_seed = random.randrange(2**32)

    numpy.random.seed(random_seed)
    torch.manual_seed(random_seed + 1)
    feature_processor = feature_pipeline.FeaturePipeline(config.data)

    with open(args.input_fasta) as fasta_file:
        tags, sequences = parse_fasta(fasta_file.read())

    if len(sequences) != 1:
        raise ValueError("the threading script can only process a single sequence")

    query_sequence = sequences[0]
    query_tag = tags[0]
    feature_dict = make_sequence_features_with_custom_template(
        query_sequence,
        args.input_mmcif,
        args.template_id,
        args.chain_id,
        args.kalign_binary_path)
    processed_feature_dict = feature_processor.process_features(
        feature_dict, mode='predict',
    )
    processed_feature_dict = {
        k: torch.as_tensor(v, device=args.model_device)
        for k, v in processed_feature_dict.items()
    }

    model_generator = load_models_from_command_line(
        config,
        args.model_device,
        args.openfold_checkpoint_path,
        args.jax_param_path,
        args.output_dir)
    output_name = f'{query_tag}_{args.config_preset}'
    for model, output_directory in model_generator:
        out = run_model(model, processed_feature_dict, query_tag, args.output_dir)

        # Toss out the recycling dimensions --- we don't need them anymore
        processed_feature_dict = tensor_tree_map(
            lambda x: numpy.array(x[..., -1].cpu()),
            processed_feature_dict
        )
        out = tensor_tree_map(lambda x: numpy.array(x.cpu()), out)


        unrelaxed_protein = prep_output(
            out,
            processed_feature_dict,
            feature_dict,
            feature_processor,
            args.config_preset,
            200, # this is the ri_multimer_gap. There's no multimer sequences here, so it doesnt matter what its set to
            args.subtract_plddt
        )

        unrelaxed_output_path = os.path.join(
            output_directory, f'{output_name}_unrelaxed.pdb'
        )

        with open(unrelaxed_output_path, 'w') as fp:
            fp.write(protein.to_pdb(unrelaxed_protein))

        logger.info(f"Output written to {unrelaxed_output_path}...")

        logger.info(f"Running relaxation on {unrelaxed_output_path}...")
        relax_protein(config, args.model_device, unrelaxed_protein, output_directory, output_name, False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", type=str, help="the path to a fasta file containing a single sequence to thread")
    parser.add_argument("input_mmcif", type=str, help="the path to an mmcif file to thread the sequence on to")

    parser.add_argument("--template_id", type=str, help="a PDB id or other identifier for the template")

    parser.add_argument(
        "--chain_id", type=str,
        help="""The chain ID of the chain in the template to use"""
    )

    parser.add_argument(
        "--model_device", type=str, default="cpu",
        help="""Name of the device on which to run the model. Any valid torch
             device name is accepted (e.g. "cpu", "cuda:0")"""
    )
    parser.add_argument(
        "--config_preset", type=str, default="model_1",
        help="""Name of a model config preset defined in openfold/config.py"""
    )
    parser.add_argument(
        "--jax_param_path", type=str, default=None,
        help="""Path to JAX model parameters. If None, and openfold_checkpoint_path
             is also None, parameters are selected automatically according to 
             the model name from openfold/resources/params"""
    )
    parser.add_argument(
        "--openfold_checkpoint_path", type=str, default=None,
        help="""Path to OpenFold checkpoint. Can be either a DeepSpeed 
             checkpoint directory or a .pt file"""
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction""",
    )
    parser.add_argument(
        "--subtract_plddt", action="store_true", default=False,
        help=""""Whether to output (100 - pLDDT) in the B-factor column instead
                 of the pLDDT itself"""
    )

    parser.add_argument(
        "--data_random_seed", type=str, default=None
    )

    add_data_args(parser)

    args = parser.parse_args()

    if(args.jax_param_path is None and args.openfold_checkpoint_path is None):
        args.jax_param_path = os.path.join(
            "openfold", "resources", "params",
            "params_" + args.config_preset + ".npz"
        )

    if(args.model_device == "cpu" and torch.cuda.is_available()):
        logging.warning(
            """The model is being run on CPU. Consider specifying 
            --model_device for better performance"""
        )

    main(args)