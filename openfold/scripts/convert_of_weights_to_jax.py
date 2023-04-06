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

import numpy as np
import torch

from openfold.config import model_config
from openfold.model.model import AlphaFold
from openfold.utils.import_weights import (
    Param, 
    ParamType, 
    generate_translation_dict, 
    process_translation_dict,
)
from openfold.utils.tensor_utils import tree_map


def reshape_fn(of_param, af_weight):
    transformations = {
        ParamType.LinearWeight: lambda w: w.transpose(-1, -2),
        ParamType.LinearWeightMHA: lambda w: w.transpose(-1, -2).reshape(af_weight.shape),
        ParamType.LinearMHAOutputWeight: lambda w: w.transpose(-1, -2).reshape(af_weight.shape),
        ParamType.LinearBiasMHA: lambda w: w.reshape(af_weight.shape),
        ParamType.LinearWeightOPM: lambda w: w.transpose(-1, -2).reshape(af_weight.shape),
        ParamType.Other: lambda w: w,
    }

    if(of_param.stacked):
        of_weight = torch.stack([torch.Tensor(p) for p in of_param.param])
    else:
        of_weight = torch.Tensor(of_param.param)
    
    return transformations[of_param.param_type](of_weight)


def transfer(of_dict, af_weight_template):
    for k in of_dict:
        if(type(of_dict[k]) == dict):
            transfer(of_dict[k], af_weight_template[k])
        else:
            reshaped = reshape_fn(of_dict[k], af_weight_template[k])
            reshaped = reshaped.detach().numpy()
            np.copyto(af_weight_template[k], reshaped)


def main(args):
    d = torch.load(args.of_pt_path)

    config = model_config(args.config_preset)
    model = AlphaFold(config)
    model.load_state_dict(d)
    
    translation = generate_translation_dict(model, args.config_preset)
    translation = process_translation_dict(translation)
    
    af_weight_template = np.load(args.template_npz_path)
    af_weight_template = {k:v for k,v in af_weight_template.items() if k in translation}
    zero = lambda n: n * 0
    af_weight_template = tree_map(zero, af_weight_template, np.ndarray)
    
    transfer(translation, af_weight_template)
    
    np.savez(args.out_path, **af_weight_template)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "of_pt_path", type=str, help="Path to OpenFold .pt checkpoint file"
    )
    parser.add_argument(
        "config_preset", type=str, help="The corresponding config preset"
    )
    parser.add_argument(
        "out_path", type=str, help="Path for output .npz file"
    )
    parser.add_argument(
        "--template_npz_path", 
        type=str, 
        default="openfold/resources/params/params_model_1_ptm.npz",
        help="""Path to an AlphaFold checkpoint w/ a superset of the OF
                checkpoint's parameters. params_model_1_ptm.npz always works.
             """
    )

    args = parser.parse_args()

    main(args)
