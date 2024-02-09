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

import argparse
import os
import shutil
import torch

from openfold.utils.import_weights import convert_deprecated_v1_keys
from zero_to_fp32 import get_optim_files, parse_optim_states, get_model_state_file

def get_latest_checkpoint_dir(checkpoint_dir):
    # Based on zero_to_fp32.get_fp32_state_dict_from_zero_checkpoint
    latest_path = os.path.join(checkpoint_dir, 'latest')
    if os.path.isfile(latest_path):
        with open(latest_path, 'r') as fd:
            tag = fd.read().strip()
    else:
        raise ValueError(f"Unable to find 'latest' file at {latest_path}")

    return os.path.join(checkpoint_dir, tag)


def convert_v1_to_v2_weights(args):
    # TODO can we use zero_to_fp32.get_fp32_state_dict_from_zero_checkpoint here?

    checkpoint_path = args.input_ckpt_path
    is_dir = os.path.isdir(checkpoint_path)
    if is_dir:
        # A DeepSpeed checkpoint
        ds_checkpoint_path = get_latest_checkpoint_dir(checkpoint_path)
        state_dict_key = 'module'
        optim_files = get_optim_files(ds_checkpoint_path)
        zero_stage, _, _ = parse_optim_states(optim_files, ds_checkpoint_path)
        model_file = get_model_state_file(ds_checkpoint_path, zero_stage)
    else:
        # A Pytorch Lightning checkpoint
        state_dict_key = 'state_dict'
        model_file = checkpoint_path

    model_dict = torch.load(model_file, map_location=torch.device('cpu'))
    model_dict[state_dict_key] = convert_deprecated_v1_keys(model_dict[state_dict_key])

    if 'ema' in model_dict:
        ema_state_dict = model_dict['ema']['params']
        model_dict['ema']['params'] = convert_deprecated_v1_keys(ema_state_dict)

    if is_dir:
        shutil.copytree(checkpoint_path, args.output_ckpt_path)
        out_fname = os.path.join(args.output_ckpt_path, os.path.basename(model_file))
    else:
        out_fname = args.output_ckpt_path

    torch.save(model_dict, out_fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_ckpt_path", type=str)
    parser.add_argument("output_ckpt_path", type=str)

    args = parser.parse_args()

    convert_v1_to_v2_weights(args)