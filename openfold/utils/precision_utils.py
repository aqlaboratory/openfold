# Copyright 2022 AlQuraishi Laboratory
# Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
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
import importlib

import torch

def cast_tensor(x, from_dtype, to_dtype):
    return x.to(dtype=to_dtype) if torch.is_tensor(x) and x.dtype == from_dtype else x


def cast_all(x, from_dtype, to_dtype):
    if isinstance(x, torch.Tensor):
        return cast_tensor(x, from_dtype=from_dtype, to_dtype=to_dtype)
    else:
        if isinstance(x, dict):
            new_dict = {}
            for k in x.keys():
                new_dict[k] = cast_all(x[k], from_dtype=from_dtype, to_dtype=to_dtype)
            return new_dict
        elif isinstance(x, tuple):
            return tuple(cast_all(y, from_dtype=from_dtype, to_dtype=to_dtype) for y in x)
        elif isinstance(x, list):
            return list(cast_all(y, from_dtype=from_dtype, to_dtype=to_dtype) for y in x)
        else:
            return x

class PrecisionWrapper(torch.nn.Module):
    def __init__(self, model, precision):
        super().__init__()
        self.precision = precision
        if self.precision == "bf16":
            print(f"Converting {model.__class__} to BF16 ...")
            model = model.bfloat16()
        elif self.precision == "fp16":
            print(f"Converting {model.__class__} to FP16 ...")
            model = model.half()
        self.model = model

    # TODO: generalize!!
    def forward(self, *args, **kwargs):
        if self.precision == "bf16":
            args = cast_all(args, from_dtype=torch.float32, to_dtype=torch.bfloat16)
            kwargs = cast_all(kwargs, from_dtype=torch.float32, to_dtype=torch.bfloat16)
        elif self.precision == "fp16":
            args = cast_all(args, from_dtype=torch.float32, to_dtype=torch.float16)
            kwargs = cast_all(kwargs,  from_dtype=torch.float32, to_dtype=torch.float16)
        out = self.model(*args, **kwargs)
        if self.precision == "bf16":
            out = cast_all(out, from_dtype=torch.bfloat16, to_dtype=torch.float32)
        elif self.precision == "fp16":
            out = cast_all(out, from_dtype=torch.float16, to_dtype=torch.float32)

        return out

def wrap_for_precision(model, precision):
    return PrecisionWrapper(model, precision)

def is_fp16_enabled():
    # Autocast world
    fp16_enabled = torch.get_autocast_gpu_dtype() == torch.float16
    fp16_enabled = fp16_enabled and torch.is_autocast_enabled()

    return fp16_enabled
