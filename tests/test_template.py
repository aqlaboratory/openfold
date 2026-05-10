# Copyright 2021 AlQuraishi Laboratory
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
import torch
import numpy as np
import unittest
from openfold.config import model_config
from openfold.model.embedders import TemplateEmbedderMultimer
from openfold.model.template import (
    TemplatePointwiseAttention,
    TemplatePairStack,
)
from openfold.np import residue_constants
import tests.compare_utils as compare_utils
from tests.config import consts
from tests.data_utils import random_template_feats

if compare_utils.alphafold_is_installed():
    alphafold = compare_utils.import_alphafold()
    import jax
    import haiku as hk


class TestTemplatePointwiseAttention(unittest.TestCase):
    def test_shape(self):
        batch_size = consts.batch_size
        n_seq = consts.n_seq
        c_t = consts.c_t
        c_z = consts.c_z
        c = 26
        no_heads = 13
        n_res = consts.n_res
        inf = 1e7

        tpa = TemplatePointwiseAttention(
            c_t, c_z, c, no_heads, inf=inf
        )

        t = torch.rand((batch_size, n_seq, n_res, n_res, c_t))
        z = torch.rand((batch_size, n_res, n_res, c_z))

        z_update = tpa(t, z, chunk_size=None)

        self.assertTrue(z_update.shape == z.shape)


class TestTemplatePairStack(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if compare_utils.alphafold_is_installed():
            if consts.is_multimer:
                cls.am_atom = alphafold.model.all_atom_multimer
                cls.am_fold = alphafold.model.folding_multimer
                cls.am_modules = alphafold.model.modules_multimer
                cls.am_rigid = alphafold.model.geometry
            else:
                cls.am_atom = alphafold.model.all_atom
                cls.am_fold = alphafold.model.folding
                cls.am_modules = alphafold.model.modules
                cls.am_rigid = alphafold.model.r3

    def test_shape(self):
        batch_size = consts.batch_size
        c_t = consts.c_t
        c_hidden_tri_att = 7
        c_hidden_tri_mul = 7
        no_blocks = 2
        no_heads = 4
        pt_inner_dim = 15
        dropout = 0.25
        n_templ = consts.n_templ
        n_res = consts.n_res
        tri_mul_first = consts.is_multimer
        fuse_projection_weights = True if re.fullmatch("^model_[1-5]_multimer_v3$", consts.model) else False
        blocks_per_ckpt = None
        chunk_size = 4
        inf = 1e7
        eps = 1e-7

        tpe = TemplatePairStack(
            c_t,
            c_hidden_tri_att=c_hidden_tri_att,
            c_hidden_tri_mul=c_hidden_tri_mul,
            no_blocks=no_blocks,
            no_heads=no_heads,
            pair_transition_n=pt_inner_dim,
            dropout_rate=dropout,
            tri_mul_first=tri_mul_first,
            fuse_projection_weights=fuse_projection_weights,
            blocks_per_ckpt=None,
            inf=inf,
            eps=eps,
        )

        t = torch.rand((batch_size, n_templ, n_res, n_res, c_t))
        mask = torch.randint(0, 2, (batch_size, n_templ, n_res, n_res))
        shape_before = t.shape
        t = tpe(t, mask, chunk_size=chunk_size)
        shape_after = t.shape

        self.assertTrue(shape_before == shape_after)

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def run_template_pair_stack(pair_act, pair_mask):
            config = compare_utils.get_alphafold_config()
            c_ee = config.model.embeddings_and_evoformer

            if consts.is_multimer:
                safe_key = alphafold.model.prng.SafeKey(hk.next_rng_key())
                template_iteration = self.am_modules.TemplateEmbeddingIteration(
                    c_ee.template.template_pair_stack,
                    config.model.global_config,
                    name='template_embedding_iteration')

                def template_iteration_fn(x):
                    act, safe_key = x

                    safe_key, safe_subkey = safe_key.split()
                    act = template_iteration(
                        act=act,
                        pair_mask=pair_mask,
                        is_training=False,
                        safe_key=safe_subkey)
                    return (act, safe_key)

                if config.model.global_config.use_remat:
                    template_iteration_fn = hk.remat(template_iteration_fn)

                safe_key, safe_subkey = safe_key.split()
                template_stack = alphafold.model.layer_stack.layer_stack(
                    c_ee.template.template_pair_stack.num_block)(
                    template_iteration_fn)
                act, _ = template_stack((pair_act, safe_subkey))
            else:
                tps = self.am_modules.TemplatePairStack(
                    c_ee.template.template_pair_stack,
                    config.model.global_config,
                    name="template_pair_stack",
                )
                act = tps(pair_act, pair_mask, is_training=False)
            ln = hk.LayerNorm([-1], True, True, name="output_layer_norm")
            act = ln(act)
            return act

        f = hk.transform(run_template_pair_stack)

        n_res = consts.n_res

        pair_act = np.random.rand(n_res, n_res, consts.c_t).astype(np.float32)
        pair_mask = np.random.randint(
            low=0, high=2, size=(n_res, n_res)
        ).astype(np.float32)

        if consts.is_multimer:
            params = compare_utils.fetch_alphafold_module_weights(
                "alphafold/alphafold_iteration/evoformer/template_embedding/"
                + "single_template_embedding/template_embedding_iteration"
            )
        else:
            params = compare_utils.fetch_alphafold_module_weights(
                "alphafold/alphafold_iteration/evoformer/template_embedding/"
                + "single_template_embedding/template_pair_stack"
            )
        params.update(
            compare_utils.fetch_alphafold_module_weights(
                "alphafold/alphafold_iteration/evoformer/template_embedding/"
                + "single_template_embedding/output_layer_norm"
            )
        )

        out_gt = f.apply(
            params, jax.random.PRNGKey(42), pair_act, pair_mask
        ).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt))

        model = compare_utils.get_global_pretrained_openfold()
        out_repro = model.template_embedder.template_pair_stack(
            torch.as_tensor(pair_act).unsqueeze(-4).cuda(),
            torch.as_tensor(pair_mask).unsqueeze(-3).cuda(),
            chunk_size=None,
            _mask_trans=False,
        ).cpu()

        compare_utils.assert_max_abs_diff_small(out_gt, out_repro, consts.eps)


class TestTemplateEmbedderMultimer(unittest.TestCase):
    def test_batched_template_embedding_shape(self):
        torch.manual_seed(0)
        batch_size = 2
        n_templ = 2
        n_res = 8

        config = model_config("model_1_multimer_v3")
        config.model.template.template_pair_stack.no_blocks = 1
        config.model.template.template_pair_stack.blocks_per_ckpt = None

        embedder = TemplateEmbedderMultimer(config.model.template)
        embedder.eval()

        aatype = torch.full(
            (batch_size, n_templ, n_res),
            residue_constants.restype_order["A"],
            dtype=torch.long,
        )
        atom_positions = torch.zeros(
            batch_size,
            n_templ,
            n_res,
            residue_constants.atom_type_num,
            3,
        )
        atom_mask = torch.zeros(
            batch_size,
            n_templ,
            n_res,
            residue_constants.atom_type_num,
        )

        n_idx = residue_constants.atom_order["N"]
        ca_idx = residue_constants.atom_order["CA"]
        c_idx = residue_constants.atom_order["C"]
        cb_idx = residue_constants.atom_order["CB"]
        residue_positions = torch.arange(n_res, dtype=torch.float32)

        atom_positions[..., n_idx, 0] = residue_positions
        atom_positions[..., ca_idx, 0] = residue_positions
        atom_positions[..., ca_idx, 1] = 1.0
        atom_positions[..., c_idx, 0] = residue_positions
        atom_positions[..., c_idx, 2] = 1.0
        atom_positions[..., cb_idx, 0] = residue_positions
        atom_positions[..., cb_idx, 1] = 1.0
        atom_positions[..., cb_idx, 2] = 1.0
        atom_mask[..., [n_idx, ca_idx, c_idx, cb_idx]] = 1.0

        batch = {
            "template_aatype": aatype,
            "template_all_atom_positions": atom_positions,
            "template_all_atom_mask": atom_mask,
        }
        z = torch.randn(batch_size, n_res, n_res, config.globals.c_z)
        padding_mask_2d = torch.ones(batch_size, n_res, n_res)
        asym_id = torch.tensor(
            [[1] * (n_res // 2) + [2] * (n_res - n_res // 2)]
            * batch_size
        )
        multichain_mask_2d = asym_id[..., None] == asym_id[..., None, :]

        with torch.no_grad():
            out = embedder(
                batch,
                z,
                padding_mask_2d,
                templ_dim=1,
                chunk_size=None,
                multichain_mask_2d=multichain_mask_2d,
                use_lma=False,
                inplace_safe=False,
            )

        self.assertEqual(
            out["template_pair_embedding"].shape,
            (batch_size, n_res, n_res, config.globals.c_z),
        )
        self.assertEqual(
            out["template_single_embedding"].shape,
            (batch_size, n_templ, n_res, config.globals.c_m),
        )
        self.assertEqual(
            out["template_mask"].shape,
            (batch_size, n_templ, n_res),
        )
        self.assertTrue(
            torch.isfinite(out["template_pair_embedding"]).all().item()
        )
        self.assertTrue(
            torch.isfinite(out["template_single_embedding"]).all().item()
        )


class Template(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if compare_utils.alphafold_is_installed():
            if consts.is_multimer:
                cls.am_atom = alphafold.model.all_atom_multimer
                cls.am_fold = alphafold.model.folding_multimer
                cls.am_modules = alphafold.model.modules_multimer
                cls.am_rigid = alphafold.model.geometry
            else:
                cls.am_atom = alphafold.model.all_atom
                cls.am_fold = alphafold.model.folding
                cls.am_modules = alphafold.model.modules
                cls.am_rigid = alphafold.model.r3

    @compare_utils.skip_unless_alphafold_installed()
    def test_compare(self):
        def test_template_embedding(pair, batch, mask_2d, mc_mask_2d):
            config = compare_utils.get_alphafold_config()
            te = self.am_modules.TemplateEmbedding(
                config.model.embeddings_and_evoformer.template,
                config.model.global_config,
            )

            if consts.is_multimer:
                act = te(pair, batch, mask_2d, multichain_mask_2d=mc_mask_2d, is_training=False)
            else:
                act = te(pair, batch, mask_2d, is_training=False)
            return act

        f = hk.transform(test_template_embedding)

        n_res = consts.n_res
        n_templ = consts.n_templ

        pair_act = np.random.rand(n_res, n_res, consts.c_z).astype(np.float32)
        batch = random_template_feats(n_templ, n_res)
        batch["template_all_atom_masks"] = batch["template_all_atom_mask"]

        multichain_mask_2d = None
        if consts.is_multimer:
            asym_id = batch['asym_id'][0]
            multichain_mask_2d = (
                    asym_id[..., None] == asym_id[..., None, :]
            ).astype(np.float32)

        pair_mask = np.random.randint(0, 2, (n_res, n_res)).astype(np.float32)
        # Fetch pretrained parameters (but only from one block)]
        params = compare_utils.fetch_alphafold_module_weights(
            "alphafold/alphafold_iteration/evoformer/template_embedding"
        )

        out_gt = f.apply(
            params, jax.random.PRNGKey(42), pair_act, batch, pair_mask, multichain_mask_2d
        ).block_until_ready()
        out_gt = torch.as_tensor(np.array(out_gt))

        inds = np.random.randint(0, 21, (n_res,))
        batch["target_feat"] = np.eye(22)[inds]

        model = compare_utils.get_global_pretrained_openfold()

        template_feats = {k: torch.as_tensor(v).cuda() for k, v in batch.items()}
        if consts.is_multimer:
            out_repro_all = model.template_embedder(
                template_feats,
                torch.as_tensor(pair_act).cuda(),
                torch.as_tensor(pair_mask).cuda(),
                templ_dim=0,
                chunk_size=consts.chunk_size,
                multichain_mask_2d=torch.as_tensor(multichain_mask_2d).cuda(),
                _mask_trans=False,
                use_lma=False,
                inplace_safe=False
            )
        else:
            out_repro_all = model.template_embedder(
                template_feats,
                torch.as_tensor(pair_act).cuda(),
                torch.as_tensor(pair_mask).cuda(),
                templ_dim=0,
                chunk_size=consts.chunk_size,
                mask_trans=False,
                use_lma=False,
                inplace_safe=False
            )

        out_repro = out_repro_all["template_pair_embedding"]
        out_repro = out_repro.cpu()

        compare_utils.assert_mean_abs_diff_small(out_gt, out_repro, consts.eps)


if __name__ == "__main__":
    unittest.main()
