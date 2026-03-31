# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
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

import re
import logging
from enum import Enum
from dataclasses import dataclass
from functools import partial
import numpy as np
import torch
from typing import Union, List


_NPZ_KEY_PREFIX = "alphafold/alphafold_iteration/"


# With Param, a poor man's enum with attributes (Rust-style)
class ParamType(Enum):
    LinearWeight = partial(  # hack: partial prevents fns from becoming methods
        lambda w: w.unsqueeze(-1) if len(w.shape) == 1 else w.transpose(-1, -2)
    )
    LinearWeightMHA = partial(
        lambda w: w.reshape(*w.shape[:-2], -1).transpose(-1, -2)
    )
    LinearMHAOutputWeight = partial(
        lambda w: w.reshape(*w.shape[:-3], -1, w.shape[-1]).transpose(-1, -2)
    )
    LinearBiasMHA = partial(lambda w: w.reshape(*w.shape[:-2], -1))
    LinearWeightOPM = partial(
        lambda w: w.reshape(*w.shape[:-3], -1, w.shape[-1]).transpose(-1, -2)
    )
    LinearWeightMultimer = partial(
        lambda w: w.unsqueeze(-1) if len(w.shape) == 1 else 
            w.reshape(w.shape[0], -1).transpose(-1, -2)
    )
    LinearBiasMultimer = partial(
        lambda w: w.reshape(-1)
    )
    Other = partial(lambda w: w)

    def __init__(self, fn):
        self.transformation = fn


@dataclass
class Param:
    param: Union[torch.Tensor, List[torch.Tensor]]
    param_type: ParamType = ParamType.Other
    stacked: bool = False
    swap: bool = False


def process_translation_dict(d, top_layer=True):
    flat = {}
    for k, v in d.items():
        if type(v) == dict:
            prefix = _NPZ_KEY_PREFIX if top_layer else ""
            sub_flat = {
                (prefix + "/".join([k, k_prime])): v_prime
                for k_prime, v_prime in process_translation_dict(
                    v, top_layer=False
                ).items()
            }
            flat.update(sub_flat)
        else:
            k = "/" + k if not top_layer else k
            flat[k] = v

    return flat


def stacked(param_dict_list, out=None):
    """
    Args:
        param_dict_list:
            A list of (nested) Param dicts to stack. The structure of
            each dict must be the identical (down to the ParamTypes of
            "parallel" Params). There must be at least one dict
            in the list.
    """
    if out is None:
        out = {}
    template = param_dict_list[0]
    for k, _ in template.items():
        v = [d[k] for d in param_dict_list]
        if type(v[0]) is dict:
            out[k] = {}
            stacked(v, out=out[k])
        elif type(v[0]) is Param:
            stacked_param = Param(
                param=[param.param for param in v],
                param_type=v[0].param_type,
                stacked=True,
                swap=v[0].swap
            )

            out[k] = stacked_param

    return out


def assign(translation_dict, orig_weights):
    for k, param in translation_dict.items():
        with torch.no_grad():
            weights = torch.as_tensor(orig_weights[k])
            ref, param_type = param.param, param.param_type
            if param.stacked:
                weights = torch.unbind(weights, 0)
            else:
                weights = [weights]
                ref = [ref]

            try:
                weights = list(map(param_type.transformation, weights))
                for p, w in zip(ref, weights):
                    if param.swap:
                        index = p.shape[0] // 2
                        p[:index].copy_(w[index:])
                        p[index:].copy_(w[:index])
                    else:
                        p.copy_(w)
            except:
                print(k)
                print(ref[0].shape)
                print(weights[0].shape)
                raise


def generate_translation_dict(model, version, is_multimer=False):
    #######################
    # Some templates
    #######################
    LinearWeight = lambda l: (Param(l, param_type=ParamType.LinearWeight))
    LinearBias = lambda l: (Param(l))
    LinearWeightMHA = lambda l: (Param(l, param_type=ParamType.LinearWeightMHA))
    LinearBiasMHA = lambda b: (Param(b, param_type=ParamType.LinearBiasMHA))
    LinearWeightOPM = lambda l: (Param(l, param_type=ParamType.LinearWeightOPM))
    LinearWeightMultimer = lambda l: (
        Param(l, param_type=ParamType.LinearWeightMultimer)
    )
    LinearBiasMultimer = lambda l: (
        Param(l, param_type=ParamType.LinearBiasMultimer)
    )
    LinearWeightSwap = lambda l: (Param(l, param_type=ParamType.LinearWeight, swap=True))
    LinearBiasSwap = lambda l: (Param(l, swap=True))

    LinearParams = lambda l: {
        "weights": LinearWeight(l.weight),
        "bias": LinearBias(l.bias),
    }

    LinearParamsMHA = lambda l: {
        "weights": LinearWeightMHA(l.weight),
        "bias": LinearBiasMHA(l.bias),
    }

    LinearParamsSwap = lambda l: {
        "weights": LinearWeightSwap(l.weight),
        "bias": LinearBiasSwap(l.bias),
    }

    LinearParamsMultimer = lambda l: {
        "weights": LinearWeightMultimer(l.weight),
        "bias": LinearBiasMultimer(l.bias),
    }

    LayerNormParams = lambda l: {
        "scale": Param(l.weight),
        "offset": Param(l.bias),
    }

    AttentionParams = lambda att: {
        "query_w": LinearWeightMHA(att.linear_q.weight),
        "key_w": LinearWeightMHA(att.linear_k.weight),
        "value_w": LinearWeightMHA(att.linear_v.weight),
        "output_w": Param(
            att.linear_o.weight,
            param_type=ParamType.LinearMHAOutputWeight,
        ),
        "output_b": LinearBias(att.linear_o.bias),
    }

    AttentionGatedParams = lambda att: dict(
        **AttentionParams(att),
        **{
            "gating_w": LinearWeightMHA(att.linear_g.weight),
            "gating_b": LinearBiasMHA(att.linear_g.bias),
        },
    )

    GlobalAttentionParams = lambda att: dict(
        AttentionGatedParams(att),
        key_w=LinearWeight(att.linear_k.weight),
        value_w=LinearWeight(att.linear_v.weight),
    )

    TriAttParams = lambda tri_att: {
        "query_norm": LayerNormParams(tri_att.layer_norm),
        "feat_2d_weights": LinearWeight(tri_att.linear.weight),
        "attention": AttentionGatedParams(tri_att.mha),
    }

    def TriMulOutParams(tri_mul, outgoing=True):
        if re.fullmatch("^model_[1-5]_multimer_v3$", version):
            lin_param_type = LinearParams if outgoing else LinearParamsSwap
            d = {
                "left_norm_input": LayerNormParams(tri_mul.layer_norm_in),
                "projection": lin_param_type(tri_mul.linear_ab_p),
                "gate": lin_param_type(tri_mul.linear_ab_g),
                "center_norm": LayerNormParams(tri_mul.layer_norm_out),
            }
        else:
            # see commit b88f8da on the Alphafold repo
            # Alphafold swaps the pseudocode's a and b between the incoming/outcoming
            # iterations of triangle multiplication, which is confusing and not
            # reproduced in our implementation.
            if outgoing:
                left_projection = LinearParams(tri_mul.linear_a_p)
                right_projection = LinearParams(tri_mul.linear_b_p)
                left_gate = LinearParams(tri_mul.linear_a_g)
                right_gate = LinearParams(tri_mul.linear_b_g)
            else:
                left_projection = LinearParams(tri_mul.linear_b_p)
                right_projection = LinearParams(tri_mul.linear_a_p)
                left_gate = LinearParams(tri_mul.linear_b_g)
                right_gate = LinearParams(tri_mul.linear_a_g)

            d = {
                "layer_norm_input": LayerNormParams(tri_mul.layer_norm_in),
                "left_projection": left_projection,
                "right_projection": right_projection,
                "left_gate": left_gate,
                "right_gate": right_gate,
                "center_layer_norm": LayerNormParams(tri_mul.layer_norm_out),
            }

        d.update({
            "output_projection": LinearParams(tri_mul.linear_z),
            "gating_linear": LinearParams(tri_mul.linear_g),
        })

        return d

    TriMulInParams = partial(TriMulOutParams, outgoing=False)

    PairTransitionParams = lambda pt: {
        "input_layer_norm": LayerNormParams(pt.layer_norm),
        "transition1": LinearParams(pt.linear_1),
        "transition2": LinearParams(pt.linear_2),
    }

    MSAAttParams = lambda matt: {
        "query_norm": LayerNormParams(matt.layer_norm_m),
        "attention": AttentionGatedParams(matt.mha),
    }

    MSAColAttParams = lambda matt: {
        "query_norm": LayerNormParams(matt._msa_att.layer_norm_m),
        "attention": AttentionGatedParams(matt._msa_att.mha),
    }

    MSAGlobalAttParams = lambda matt: {
        "query_norm": LayerNormParams(matt.layer_norm_m),
        "attention": GlobalAttentionParams(matt.global_attention),
    }

    MSAAttPairBiasParams = lambda matt: dict(
        **MSAAttParams(matt),
        **{
            "feat_2d_norm": LayerNormParams(matt.layer_norm_z),
            "feat_2d_weights": LinearWeight(matt.linear_z.weight),
        },
    )

    IPAParams = lambda ipa: {
        "q_scalar": LinearParams(ipa.linear_q),
        "kv_scalar": LinearParams(ipa.linear_kv),
        "q_point_local": LinearParams(ipa.linear_q_points.linear),
        "kv_point_local": LinearParams(ipa.linear_kv_points.linear),
        "trainable_point_weights": Param(
            param=ipa.head_weights, param_type=ParamType.Other
        ),
        "attention_2d": LinearParams(ipa.linear_b),
        "output_projection": LinearParams(ipa.linear_out),
    }

    PointProjectionParams = lambda pp: {
        "point_projection": LinearParamsMHA(
            pp.linear,
        ),
    }

    IPAParamsMultimer = lambda ipa: {
        "q_scalar_projection": {
            "weights": LinearWeightMHA(
                ipa.linear_q.weight,
            ),
        },
        "k_scalar_projection": {
            "weights": LinearWeightMHA(
                ipa.linear_k.weight,
            ),
        },
        "v_scalar_projection": {
            "weights": LinearWeightMHA(
                ipa.linear_v.weight,
            ),
        },
        "q_point_projection": PointProjectionParams(
            ipa.linear_q_points
        ),
        "k_point_projection": PointProjectionParams(
            ipa.linear_k_points
        ),
        "v_point_projection": PointProjectionParams(
            ipa.linear_v_points
        ),
        "trainable_point_weights": Param(
            param=ipa.head_weights, param_type=ParamType.Other
        ),
        "attention_2d": LinearParams(ipa.linear_b),
        "output_projection": LinearParams(ipa.linear_out),
    }

    TemplatePairBlockParams = lambda b: {
        "triangle_attention_starting_node": TriAttParams(b.tri_att_start),
        "triangle_attention_ending_node": TriAttParams(b.tri_att_end),
        "triangle_multiplication_outgoing": TriMulOutParams(b.tri_mul_out),
        "triangle_multiplication_incoming": TriMulInParams(b.tri_mul_in),
        "pair_transition": PairTransitionParams(b.pair_transition),
    }

    MSATransitionParams = lambda m: {
        "input_layer_norm": LayerNormParams(m.layer_norm),
        "transition1": LinearParams(m.linear_1),
        "transition2": LinearParams(m.linear_2),
    }

    OuterProductMeanParams = lambda o: {
        "layer_norm_input": LayerNormParams(o.layer_norm),
        "left_projection": LinearParams(o.linear_1),
        "right_projection": LinearParams(o.linear_2),
        "output_w": LinearWeightOPM(o.linear_out.weight),
        "output_b": LinearBias(o.linear_out.bias),
    }

    def EvoformerBlockParams(b, is_extra_msa=False):
        if is_extra_msa:
            col_att_name = "msa_column_global_attention"
            msa_col_att_params = MSAGlobalAttParams(b.msa_att_col)
        else:
            col_att_name = "msa_column_attention"
            msa_col_att_params = MSAColAttParams(b.msa_att_col)

        d = {
            "msa_row_attention_with_pair_bias": MSAAttPairBiasParams(
                b.msa_att_row
            ),
            col_att_name: msa_col_att_params,
            "msa_transition": MSATransitionParams(b.msa_transition),
            "outer_product_mean": 
                OuterProductMeanParams(b.outer_product_mean),
            "triangle_multiplication_outgoing": 
                TriMulOutParams(b.pair_stack.tri_mul_out),
            "triangle_multiplication_incoming": 
                TriMulInParams(b.pair_stack.tri_mul_in),
            "triangle_attention_starting_node": 
                TriAttParams(b.pair_stack.tri_att_start),
            "triangle_attention_ending_node": 
                TriAttParams(b.pair_stack.tri_att_end),
            "pair_transition": 
                PairTransitionParams(b.pair_stack.pair_transition),
        }

        return d

    ExtraMSABlockParams = partial(EvoformerBlockParams, is_extra_msa=True)

    def FoldIterationParams(sm):
        d = {
            "invariant_point_attention": 
                IPAParamsMultimer(sm.ipa) if is_multimer else IPAParams(sm.ipa),
            "attention_layer_norm": LayerNormParams(sm.layer_norm_ipa),
            "transition": LinearParams(sm.transition.layers[0].linear_1),
            "transition_1": LinearParams(sm.transition.layers[0].linear_2),
            "transition_2": LinearParams(sm.transition.layers[0].linear_3),
            "transition_layer_norm": LayerNormParams(sm.transition.layer_norm),
            "affine_update": LinearParams(sm.bb_update.linear),
            "rigid_sidechain": {
                "input_projection": LinearParams(sm.angle_resnet.linear_in),
                "input_projection_1": 
                    LinearParams(sm.angle_resnet.linear_initial),
                "resblock1": LinearParams(sm.angle_resnet.layers[0].linear_1),
                "resblock2": LinearParams(sm.angle_resnet.layers[0].linear_2),
                "resblock1_1": 
                    LinearParams(sm.angle_resnet.layers[1].linear_1),
                "resblock2_1": 
                    LinearParams(sm.angle_resnet.layers[1].linear_2),
                "unnormalized_angles": 
                    LinearParams(sm.angle_resnet.linear_out),
            },
        }

        if(is_multimer):
            d.pop("affine_update")
            d["quat_rigid"] = {
                "rigid": LinearParams(
                   sm.bb_update.linear
                )
            }

        return d

    ############################
    # translations dict overflow
    ############################
    ems_blocks = model.extra_msa_stack.blocks
    ems_blocks_params = stacked([ExtraMSABlockParams(b) for b in ems_blocks])

    evo_blocks = model.evoformer.blocks
    evo_blocks_params = stacked([EvoformerBlockParams(b) for b in evo_blocks])

    if(not is_multimer):
        translations = {
            "evoformer": {
                "preprocess_1d": LinearParams(model.input_embedder.linear_tf_m),
                "preprocess_msa": LinearParams(model.input_embedder.linear_msa_m),
                "left_single": LinearParams(model.input_embedder.linear_tf_z_i),
                "right_single": LinearParams(model.input_embedder.linear_tf_z_j),
                "prev_pos_linear": LinearParams(model.recycling_embedder.linear),
                "prev_msa_first_row_norm": LayerNormParams(
                    model.recycling_embedder.layer_norm_m
                ),
                "prev_pair_norm": LayerNormParams(
                    model.recycling_embedder.layer_norm_z
                ),
                "pair_activiations": LinearParams(
                    model.input_embedder.linear_relpos
                ),
                "extra_msa_activations": LinearParams(
                    model.extra_msa_embedder.linear
                ),
                "extra_msa_stack": ems_blocks_params,
                "evoformer_iteration": evo_blocks_params,
                "single_activations": LinearParams(model.evoformer.linear),
            },
            "structure_module": {
                "single_layer_norm": LayerNormParams(
                    model.structure_module.layer_norm_s
                ),
                "initial_projection": LinearParams(
                    model.structure_module.linear_in
                ),
                "pair_layer_norm": LayerNormParams(
                    model.structure_module.layer_norm_z
                ),
                "fold_iteration": FoldIterationParams(model.structure_module),
            },
            "predicted_lddt_head": {
                "input_layer_norm": LayerNormParams(
                    model.aux_heads.plddt.layer_norm
                ),
                "act_0": LinearParams(model.aux_heads.plddt.linear_1),
                "act_1": LinearParams(model.aux_heads.plddt.linear_2),
                "logits": LinearParams(model.aux_heads.plddt.linear_3),
            },
            "distogram_head": {
                "half_logits": LinearParams(model.aux_heads.distogram.linear),
            },
            "experimentally_resolved_head": {
                "logits": LinearParams(
                    model.aux_heads.experimentally_resolved.linear
                ),
            },
            "masked_msa_head": {
                "logits": LinearParams(model.aux_heads.masked_msa.linear),
            },
        }
    else:
        translations = {
            "evoformer": {
                "preprocess_1d": LinearParams(model.input_embedder.linear_tf_m),
                "preprocess_msa": LinearParams(model.input_embedder.linear_msa_m),
                "left_single": LinearParams(model.input_embedder.linear_tf_z_i),
                "right_single": LinearParams(model.input_embedder.linear_tf_z_j),
                "prev_pos_linear": LinearParams(model.recycling_embedder.linear),
                "prev_msa_first_row_norm": LayerNormParams(
                    model.recycling_embedder.layer_norm_m
                ),
                "prev_pair_norm": LayerNormParams(
                    model.recycling_embedder.layer_norm_z
                ),
                "~_relative_encoding": {
                    "position_activations": LinearParams(
                        model.input_embedder.linear_relpos
                    ),
                },
                "extra_msa_activations": LinearParams(
                    model.extra_msa_embedder.linear
                ),
                "extra_msa_stack": ems_blocks_params,
                "evoformer_iteration": evo_blocks_params,
                "single_activations": LinearParams(model.evoformer.linear),
            },
            "structure_module": {
                "single_layer_norm": LayerNormParams(
                    model.structure_module.layer_norm_s
                ),
                "initial_projection": LinearParams(
                    model.structure_module.linear_in
                ),
                "pair_layer_norm": LayerNormParams(
                    model.structure_module.layer_norm_z
                ),
                "fold_iteration": FoldIterationParams(model.structure_module),
            },
            "predicted_lddt_head": {
                "input_layer_norm": LayerNormParams(
                    model.aux_heads.plddt.layer_norm
                ),
                "act_0": LinearParams(model.aux_heads.plddt.linear_1),
                "act_1": LinearParams(model.aux_heads.plddt.linear_2),
                "logits": LinearParams(model.aux_heads.plddt.linear_3),
            },
            "distogram_head": {
                "half_logits": LinearParams(model.aux_heads.distogram.linear),
            },
            "experimentally_resolved_head": {
                "logits": LinearParams(
                    model.aux_heads.experimentally_resolved.linear
                ),
            },
            "masked_msa_head": {
                "logits": LinearParams(model.aux_heads.masked_msa.linear),
            },
        }

    no_templ = [
        "model_3",
        "model_4",
        "model_5",
        "model_3_ptm",
        "model_4_ptm",
        "model_5_ptm",
    ]

    if version not in no_templ:
        tps_blocks = model.template_embedder.template_pair_stack.blocks
        tps_blocks_params = stacked(
            [TemplatePairBlockParams(b) for b in tps_blocks]
        )
        if (not is_multimer):
            template_param_dict = {
                "template_embedding": {
                    "single_template_embedding": {
                        "embedding2d": LinearParams(
                            model.template_embedder.template_pair_embedder.linear
                        ),
                        "template_pair_stack": {
                            "__layer_stack_no_state": tps_blocks_params,
                        },
                        "output_layer_norm": LayerNormParams(
                            model.template_embedder.template_pair_stack.layer_norm
                        ),
                    },
                    "attention": AttentionParams(model.template_embedder.template_pointwise_att.mha),
                },
                "template_single_embedding": LinearParams(
                    model.template_embedder.template_single_embedder.linear_1
                ),
                "template_projection": LinearParams(
                    model.template_embedder.template_single_embedder.linear_2
                ),
            }
        else:
            temp_embedder = model.template_embedder
            template_param_dict = {
                "template_embedding": {
                    "single_template_embedding": {
                        "query_embedding_norm": LayerNormParams(
                            temp_embedder.template_pair_embedder.query_embedding_layer_norm
                        ),
                        "template_pair_embedding_0": LinearParams(
                            temp_embedder.template_pair_embedder.dgram_linear
                        ),
                        "template_pair_embedding_1": LinearParams(
                            temp_embedder.template_pair_embedder.pseudo_beta_mask_linear
                        ),
                        "template_pair_embedding_2": LinearParams(
                            temp_embedder.template_pair_embedder.aatype_linear_1
                        ),
                        "template_pair_embedding_3": LinearParams(
                            temp_embedder.template_pair_embedder.aatype_linear_2
                        ),
                        "template_pair_embedding_4": LinearParams(
                            temp_embedder.template_pair_embedder.x_linear
                        ),
                        "template_pair_embedding_5": LinearParams(
                            temp_embedder.template_pair_embedder.y_linear
                        ),
                        "template_pair_embedding_6": LinearParams(
                            temp_embedder.template_pair_embedder.z_linear
                        ),
                        "template_pair_embedding_7": LinearParams(
                            temp_embedder.template_pair_embedder.backbone_mask_linear
                        ),
                        "template_pair_embedding_8": LinearParams(
                            temp_embedder.template_pair_embedder.query_embedding_linear
                        ),
                        "template_embedding_iteration": tps_blocks_params,
                        "output_layer_norm": LayerNormParams(
                            temp_embedder.template_pair_stack.layer_norm
                        ),
                    },
                    "output_linear": LinearParams(
                        temp_embedder.linear_t
                    ),
                },
                "template_projection": LinearParams(
                    temp_embedder.template_single_embedder.template_projector,
                ),
                "template_single_embedding": LinearParams(
                    temp_embedder.template_single_embedder.template_single_embedder,
                ),
            }

        translations["evoformer"].update(template_param_dict)

    if is_multimer or "_ptm" in version:
        translations["predicted_aligned_error_head"] = {
            "logits": LinearParams(model.aux_heads.tm.linear)
        }

    return translations


def import_jax_weights_(model, npz_path, version="model_1"):
    data = np.load(npz_path)
    translations = generate_translation_dict(model, version, is_multimer=("multimer" in version))

    # Flatten keys and insert missing key prefixes
    flat = process_translation_dict(translations)

    # Sanity check
    keys = list(data.keys())
    flat_keys = list(flat.keys())
    incorrect = [k for k in flat_keys if k not in keys]
    missing = [k for k in keys if k not in flat_keys]
    # print(f"Incorrect: {incorrect}")
    # print(f"Missing: {missing}")

    assert len(incorrect) == 0
    # assert(sorted(list(flat.keys())) == sorted(list(data.keys())))

    # Set weights
    assign(flat, data)


def convert_deprecated_v1_keys(state_dict):
    """Update older OpenFold model weight names to match the current model code."""

    replacements = {
        'template_angle_embedder': 'template_single_embedder',
        'core.msa_transition': 'msa_transition',
        'core.outer_product_mean': 'outer_product_mean',
        'core.tri_': 'pair_stack.tri_',
        'core.pair_transition': 'pair_stack.pair_transition',
        'ipa.linear_q_points': 'ipa.linear_q_points.linear',
        'ipa.linear_kv_points': 'ipa.linear_kv_points.linear'
    }

    convert_key_re = re.compile("(%s)" % "|".join(map(re.escape, replacements.keys())))
    template_emb_re = re.compile(r"^((module\.)?(model\.)?)(template(?!_embedder).*)") 

    converted_state_dict = {}
    for key, value in state_dict.items():
        # For each match, look-up replacement value in the dictionary
        new_key = convert_key_re.sub(lambda m: replacements[m.group(1)], key)

        # Add prefix for template layers 
        template_match = re.match(template_emb_re, new_key)
        if template_match:
            prefix = template_match.group(1)
            new_key = f'{prefix if prefix else ""}template_embedder.{template_match.group(4)}'

        converted_state_dict[new_key] = value

    return converted_state_dict


def import_openfold_weights_(model, state_dict):
    """
    Import model weights. Several parts of the model were refactored in the process
    of adding support for Multimer. The state dicts of older models are translated
    to match the refactored model code.
    """
    try:
        model.load_state_dict(state_dict)
    except RuntimeError:
        converted_state_dict = convert_deprecated_v1_keys(state_dict)
        model.load_state_dict(converted_state_dict)
