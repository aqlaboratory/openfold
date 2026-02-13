# Copyright 2022 AlQuraishi Laboratory
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
#
# Converts OpenFold .pt checkpoints into AlphaFold .npz ones, which can then be
# used to run inference using DeepMind's JAX code.

import logging
import argparse
import os
import shutil
import torch

from openfold.utils.import_weights import convert_deprecated_v1_keys
from deepspeed.utils.zero_to_fp32 import (
    get_optim_files, parse_optim_states, get_model_state_file
)


def convert_v1_to_v2_weights(args):
    checkpoint_path = args.input_ckpt_path
    is_dir = os.path.isdir(checkpoint_path)
    if is_dir:
        # A DeepSpeed checkpoint
        logging.info(
            'Converting deepspeed checkpoint found at {args.input_checkpoint_path}')
        state_dict_key = 'module'
        latest_path = os.path.join(checkpoint_path, 'latest')
        if os.path.isfile(latest_path):
            with open(latest_path, 'r') as fd:
                tag = fd.read().strip()
        else:
            raise ValueError(f"Unable to find 'latest' file at {latest_path}")

        ds_checkpoint_dir = os.path.join(checkpoint_path, tag)
        model_output_path = os.path.join(args.output_ckpt_path, tag)
        optim_files = get_optim_files(ds_checkpoint_dir)
        zero_stage, _, _ = parse_optim_states(optim_files, ds_checkpoint_dir)
        model_file = get_model_state_file(ds_checkpoint_dir, zero_stage)
    else:
        # A Pytorch Lightning checkpoint
        logging.info(
            'Converting pytorch lightning checkpoint found at {args.input_checkpoint_path}')
        state_dict_key = 'state_dict'
        model_output_path = args.output_ckpt_path
        model_file = checkpoint_path

    model_dict = torch.load(model_file, map_location=torch.device('cpu'))
    model_dict[state_dict_key] = convert_deprecated_v1_keys(
        model_dict[state_dict_key])

    if 'ema' in model_dict:
        ema_state_dict = model_dict['ema']['params']
        model_dict['ema']['params'] = convert_deprecated_v1_keys(
            ema_state_dict)

    if is_dir:
        param_shapes = convert_deprecated_v1_keys(
            model_dict['param_shapes'][0])
        model_dict['param_shapes'] = [param_shapes]

        shutil.copytree(checkpoint_path, args.output_ckpt_path)
        out_fname = os.path.join(
            model_output_path, os.path.basename(model_file))

        for optim_file in optim_files:
            optim_dict = torch.load(optim_file)
            new_optim_dict = optim_dict.copy()
            new_optim_dict['optimizer_state_dict']['param_slice_mappings'][0] = convert_deprecated_v1_keys(
                optim_dict['optimizer_state_dict']['param_slice_mappings'][0])
            out_optim_fname = os.path.join(
                model_output_path, os.path.basename(optim_file))
            torch.save(new_optim_dict, out_optim_fname)
    else:
        out_fname = model_output_path

    torch.save(model_dict, out_fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_ckpt_path", type=str)
    parser.add_argument("output_ckpt_path", type=str)

    args = parser.parse_args()

    convert_v1_to_v2_weights(args)
