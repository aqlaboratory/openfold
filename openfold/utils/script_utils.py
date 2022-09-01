import json
import logging
import os
import re
import time

import numpy
import torch

from openfold.model.model import AlphaFold
from openfold.np import residue_constants, protein
from openfold.np.relax import relax
from openfold.utils.import_weights import (
    import_jax_weights_,
)

from pytorch_lightning.utilities.deepspeed import (
    convert_zero_checkpoint_to_fp32_state_dict
)

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)


def count_models_to_evaluate(openfold_checkpoint_path, jax_param_path):
    model_count = 0
    if openfold_checkpoint_path:
        model_count += len(openfold_checkpoint_path.split(","))
    if jax_param_path:
        model_count += len(jax_param_path.split(","))
    return model_count


def get_model_basename(model_path):
    return os.path.splitext(
                os.path.basename(
                    os.path.normpath(model_path)
                )
            )[0]


def make_output_directory(output_dir, model_name, multiple_model_mode):
    if multiple_model_mode:
        prediction_dir = os.path.join(output_dir, "predictions", model_name)
    else:
        prediction_dir = os.path.join(output_dir, "predictions")
    os.makedirs(prediction_dir, exist_ok=True)
    return prediction_dir


def load_models_from_command_line(config, model_device, openfold_checkpoint_path, jax_param_path, output_dir):
    # Create the output directory

    multiple_model_mode = count_models_to_evaluate(openfold_checkpoint_path, jax_param_path) > 1
    if multiple_model_mode:
        logger.info(f"evaluating multiple models")

    if jax_param_path:
        for path in jax_param_path.split(","):
            model_basename = get_model_basename(path)
            model_version = "_".join(model_basename.split("_")[1:])
            model = AlphaFold(config)
            model = model.eval()
            import_jax_weights_(
                model, path, version=model_version
            )
            model = model.to(model_device)
            logger.info(
                f"Successfully loaded JAX parameters at {path}..."
            )
            output_directory = make_output_directory(output_dir, model_basename, multiple_model_mode)
            yield model, output_directory

    if openfold_checkpoint_path:
        for path in openfold_checkpoint_path.split(","):
            model = AlphaFold(config)
            model = model.eval()
            checkpoint_basename = get_model_basename(path)
            if os.path.isdir(path):
                # A DeepSpeed checkpoint
                ckpt_path = os.path.join(
                    output_dir,
                    checkpoint_basename + ".pt",
                )

                if not os.path.isfile(ckpt_path):
                    convert_zero_checkpoint_to_fp32_state_dict(
                        path,
                        ckpt_path,
                    )
                d = torch.load(ckpt_path)
                model.load_state_dict(d["ema"]["params"])
            else:
                ckpt_path = path
                d = torch.load(ckpt_path)

                if "ema" in d:
                    # The public weights have had this done to them already
                    d = d["ema"]["params"]
                model.load_state_dict(d)

            model = model.to(model_device)
            logger.info(
                f"Loaded OpenFold parameters at {path}..."
            )
            output_directory = make_output_directory(output_dir, checkpoint_basename, multiple_model_mode)
            yield model, output_directory

    if not jax_param_path and not openfold_checkpoint_path:
        raise ValueError(
            "At least one of jax_param_path or openfold_checkpoint_path must "
            "be specified."
        )


def parse_fasta(data):
    data = re.sub('>$', '', data, flags=re.M)
    lines = [
        l.replace('\n', '')
        for prot in data.split('>') for l in prot.strip().split('\n', 1)
    ][1:]
    tags, seqs = lines[::2], lines[1::2]

    tags = [t.split()[0] for t in tags]

    return tags, seqs


def update_timings(timing_dict, output_file=os.path.join(os.getcwd(), "timings.json")):
    """
    Write dictionary of one or more run step times to a file
    """
    if os.path.exists(output_file):
        with open(output_file, "r") as f:
            try:
                timings = json.load(f)
            except json.JSONDecodeError:
                logger.info(f"Overwriting non-standard JSON in {output_file}.")
                timings = {}
    else:
        timings = {}
    timings.update(timing_dict)
    with open(output_file, "w") as f:
        json.dump(timings, f)
    return output_file


def run_model(model, batch, tag, output_dir):
    with torch.no_grad():
        # Temporarily disable templates if there aren't any in the batch
        template_enabled = model.config.template.enabled
        model.config.template.enabled = template_enabled and any([
            "template_" in k for k in batch
        ])

        logger.info(f"Running inference for {tag}...")
        t = time.perf_counter()
        out = model(batch)
        inference_time = time.perf_counter() - t
        logger.info(f"Inference time: {inference_time}")
        update_timings({"inference": inference_time}, os.path.join(output_dir, "timings.json"))

        model.config.template.enabled = template_enabled

    return out


def prep_output(out, batch, feature_dict, feature_processor, config_preset, multimer_ri_gap, subtract_plddt):
    plddt = out["plddt"]

    plddt_b_factors = numpy.repeat(
        plddt[..., None], residue_constants.atom_type_num, axis=-1
    )

    if subtract_plddt:
        plddt_b_factors = 100 - plddt_b_factors

    # Prep protein metadata
    template_domain_names = []
    template_chain_index = None
    if feature_processor.config.common.use_templates and "template_domain_names" in feature_dict:
        template_domain_names = [
            t.decode("utf-8") for t in feature_dict["template_domain_names"]
        ]

        # This works because templates are not shuffled during inference
        template_domain_names = template_domain_names[
                                :feature_processor.config.predict.max_templates
                                ]

        if "template_chain_index" in feature_dict:
            template_chain_index = feature_dict["template_chain_index"]
            template_chain_index = template_chain_index[
                                   :feature_processor.config.predict.max_templates
                                   ]

    no_recycling = feature_processor.config.common.max_recycling_iters
    remark = ', '.join([
        f"no_recycling={no_recycling}",
        f"max_templates={feature_processor.config.predict.max_templates}",
        f"config_preset={config_preset}",
    ])

    # For multi-chain FASTAs
    ri = feature_dict["residue_index"]
    chain_index = (ri - numpy.arange(ri.shape[0])) / multimer_ri_gap
    chain_index = chain_index.astype(numpy.int64)
    cur_chain = 0
    prev_chain_max = 0
    for i, c in enumerate(chain_index):
        if c != cur_chain:
            cur_chain = c
            prev_chain_max = i + cur_chain * multimer_ri_gap

        batch["residue_index"][i] -= prev_chain_max

    unrelaxed_protein = protein.from_prediction(
        features=batch,
        result=out,
        b_factors=plddt_b_factors,
        chain_index=chain_index,
        remark=remark,
        parents=template_domain_names,
        parents_chain_index=template_chain_index,
    )

    return unrelaxed_protein


def relax_protein(config, model_device, unrelaxed_protein, output_directory, output_name):
    amber_relaxer = relax.AmberRelaxation(
        use_gpu=(model_device != "cpu"),
        **config.relax,
    )

    t = time.perf_counter()
    visible_devices = os.getenv("CUDA_VISIBLE_DEVICES", default="")
    if "cuda" in model_device:
        device_no = model_device.split(":")[-1]
        os.environ["CUDA_VISIBLE_DEVICES"] = device_no
    relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
    os.environ["CUDA_VISIBLE_DEVICES"] = visible_devices
    relaxation_time = time.perf_counter() - t

    logger.info(f"Relaxation time: {relaxation_time}")
    update_timings({"relaxation": relaxation_time}, os.path.join(output_directory, "timings.json"))

    # Save the relaxed PDB.
    relaxed_output_path = os.path.join(
        output_directory, f'{output_name}_relaxed.pdb'
    )
    with open(relaxed_output_path, 'w') as fp:
        fp.write(relaxed_pdb_str)

    logger.info(f"Relaxed output written to {relaxed_output_path}...")