# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import os

import torch

from .tensorrt_lazy_compiler import trt_compile

logger = logging.getLogger("trt_compile")
logger.setLevel(level=logging.INFO)


def instrument_with_trt_compile(model, config):
    if config.trt.mode is None:
        return

    if (
        config.globals.use_cuequivariance_attention
        or config.globals.use_cuequivariance_multiplicative_update
    ):
        from cuequivariance_ops_torch.onnx import op_table
        from cuequivariance_ops_torch.tensorrt import register_plugins

        register_plugins()
    else:
        op_table = None

    engine_dir = config.trt.engine_dir
    os.makedirs(engine_dir, exist_ok=True)
    # Clean the directory if rebuilding
    if config.trt.mode == "build":
        for filename in os.listdir(engine_dir):
            file_path = os.path.join(engine_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)

    # skip_once_registry = trt_compile_make_registry(model.structure_module, "sample")

    S_MIN = 16
    S_MAX = config.trt.max_sequence_len

    # TODO: use config for those numbers
    seq_dim = torch.export.Dim("seq_len", max=S_MAX)

    evoformer_dynamic_shapes = {
        "m": {1: seq_dim},
        "z": {0: seq_dim, 1: seq_dim},
        "msa_mask": {1: seq_dim},
        "pair_mask": {0: seq_dim, 1: seq_dim},
        "chunk_size": None,
        "use_deepspeed_evo_attention": None,
        "use_cuequivariance_attention": None,
        "use_cuequivariance_multiplicative_update": None,
        "use_lma": None,
        "use_flash": None,
        "inplace_safe": None,
        "_mask_trans": None,
    }

    def evoformer_profile_one(min_len, max_len):
        return {
            "m": [[516, min_len, 256], [516, max_len, 256], [516, max_len, 256]],
            "z": [
                [min_len, min_len, 128],
                [max_len, max_len, 128],
                [max_len, max_len, 128],
            ],
            "msa_mask": [[516, min_len], [516, max_len], [516, max_len]],
            "pair_mask": [[min_len, min_len], [max_len, max_len], [max_len, max_len]],
        }

    def msa_profile_one(min_len, max_len):
        return {
            "m": [[5120, min_len, 64], [5120, max_len, 64], [5120, max_len, 64]],
            "z": [
                [min_len, min_len, 128],
                [max_len, max_len, 128],
                [max_len, max_len, 128],
            ],
            "msa_mask": [[5120, min_len], [5120, max_len], [5120, max_len]],
            "pair_mask": [[min_len, min_len], [max_len, max_len], [max_len, max_len]],
        }

    def input_profiles(input_profile_one, num_profiles=1):
        if num_profiles == 4:
            return [
                input_profile_one(S_MIN, S_MAX // 4),
                input_profile_one(S_MAX // 4 + 1, S_MAX // 2),
                input_profile_one(S_MAX // 2 + 1, (S_MAX // 4) * 3),
                input_profile_one((S_MAX // 4) * 3, S_MAX),
            ]
        elif num_profiles == 2:
            return [
                input_profile_one(S_MIN, S_MAX // 2),
                input_profile_one(S_MAX // 2 + 1, S_MAX),
            ]
        else:  # default: num_profiles = 1
            return [
                input_profile_one(S_MIN, S_MAX),
            ]

    evoformer_compile_args = {
        "precision": config.precision,
        "fallback": True,
        "input_profiles": input_profiles(
            evoformer_profile_one, config.trt.num_profiles
        ),
        "export_args": {
            "opset_version": 20,
            "dynamo": True,
            "report": False,
            "dynamic_shapes": evoformer_dynamic_shapes,
        },
        "build_args": {
            "builder_optimization_level": config.trt.optimization_level,
            "precision_constraints": "prefer",
        },
    }

    if op_table is not None:
        evoformer_compile_args["export_args"]["custom_translation_table"] = op_table

    trt_compile(
        model.evoformer,
        f"{engine_dir}/EvoformerStack",
        args=evoformer_compile_args,
        logger=logger,
    )
    logger.info("model.evoformer instrumented")

    """
    msa_dynamic_shapes = copy.deepcopy(evoformer_dynamic_shapes)
    msa_dynamic_shapes.pop("use_flash")

    msa_compile_args = copy.deepcopy(evoformer_compile_args)
    msa_compile_args["export_args"]["dynamic_shapes"] = msa_dynamic_shapes
    msa_compile_args["input_profiles"] = input_profiles(msa_profile_one, 2)

    # Requires too much memory - also need to share context memory if using more than one engine
    if False: # model.extra_msa_config.enabled:
        trt_compile(
            model.extra_msa_stack,
            f"{engine_dir}/ExtraMSAStack",
            args = msa_compile_args,
            logger = logger
        )
    """
