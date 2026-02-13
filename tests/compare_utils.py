import os

import importlib
import pkgutil
import sys
import unittest

import numpy as np
import torch

from openfold.config import model_config
from openfold.model.model import AlphaFold
from openfold.utils.import_weights import import_jax_weights_

from tests.config import consts

# Give JAX some GPU memory discipline
# (by default it hogs 90% of GPU memory. This disables that behavior and also
# forces it to proactively free memory that it allocates)
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
os.environ["JAX_PLATFORM_NAME"] = "gpu"


def skip_unless_ds4s_installed():
    deepspeed_is_installed = importlib.util.find_spec("deepspeed") is not None
    ds4s_is_installed = deepspeed_is_installed and importlib.util.find_spec(
        "deepspeed.ops.deepspeed4science") is not None
    return unittest.skipUnless(ds4s_is_installed, "Requires DeepSpeed with version ≥ 0.10.4")


def skip_unless_flash_attn_installed():
    fa_is_installed = importlib.util.find_spec("flash_attn") is not None
    return unittest.skipUnless(fa_is_installed, "Requires Flash Attention")


def alphafold_is_installed():
    return importlib.util.find_spec("alphafold") is not None


def skip_unless_alphafold_installed():
    return unittest.skipUnless(alphafold_is_installed(), "Requires AlphaFold")


def import_alphafold():
    """
    If AlphaFold is installed using the provided setuptools script, this
    is necessary to expose all of AlphaFold's precious insides
    """
    if "alphafold" in sys.modules:
        return sys.modules["alphafold"]
    module = importlib.import_module("alphafold")
    # Forcefully import alphafold's submodules
    submodules = pkgutil.walk_packages(module.__path__, prefix=("alphafold."))
    for submodule_info in submodules:
        importlib.import_module(submodule_info.name)
    sys.modules["alphafold"] = module
    globals()["alphafold"] = module

    return module


def get_alphafold_config():
    config = alphafold.model.config.model_config(consts.model)  # noqa
    config.model.global_config.deterministic = True
    return config


dir_path = os.path.dirname(os.path.realpath(__file__))
_param_path = os.path.join(dir_path, "..", f"openfold/resources/params/params_{consts.model}.npz")
_model = None


def get_global_pretrained_openfold():
    global _model
    if _model is None:
        _model = AlphaFold(model_config(consts.model))
        _model = _model.eval()
        if not os.path.exists(_param_path):
            raise FileNotFoundError(
                """Cannot load pretrained parameters. Make sure to run the 
                installation script before running tests."""
            )
        import_jax_weights_(_model, _param_path, version=consts.model)
        _model = _model.cuda()

    return _model


_orig_weights = None


def _get_orig_weights():
    global _orig_weights
    if _orig_weights is None:
        _orig_weights = np.load(_param_path)

    return _orig_weights


def _remove_key_prefix(d, prefix):
    for k, v in list(d.items()):
        if k.startswith(prefix):
            d.pop(k)
            d[k[len(prefix):]] = v


def fetch_alphafold_module_weights(weight_path):
    orig_weights = _get_orig_weights()
    params = {k: v for k, v in orig_weights.items() if weight_path in k}
    if "/" in weight_path:
        spl = weight_path.split("/")
        spl = spl if len(spl[-1]) != 0 else spl[:-1]
        prefix = "/".join(spl[:-1]) + "/"
        _remove_key_prefix(params, prefix)

    try:
        params = alphafold.model.utils.flat_params_to_haiku(params)  # noqa
    except:
        raise ImportError(
            "Make sure to call import_alphafold before running this function"
        )
    return params


def _assert_abs_diff_small_base(compare_func, expected, actual, eps):
    # Helper function for comparing absolute differences of two torch tensors.
    abs_diff = torch.abs(expected - actual)
    err = compare_func(abs_diff)
    zero_tensor = torch.tensor(0, dtype=err.dtype)
    rtol = 1.6e-2 if err.dtype == torch.bfloat16 else 1.3e-6  
    torch.testing.assert_close(err, zero_tensor, atol=eps, rtol=rtol)


def assert_max_abs_diff_small(expected, actual, eps):
    _assert_abs_diff_small_base(torch.max, expected, actual, eps)


def assert_mean_abs_diff_small(expected, actual, eps):
    _assert_abs_diff_small_base(torch.mean, expected, actual, eps)
