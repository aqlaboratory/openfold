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

import math
import torch
import unittest

from alphafold.utils.utils import *


X_90_ROT = torch.tensor([
    [1, 0, 0],
    [0, 0,-1],
    [0, 1, 0],
])

X_NEG_90_ROT = torch.tensor([
    [1, 0, 0],
    [0, 0, 1],
    [0,-1, 0],
])


class TestAffineT(unittest.TestCase):
    def test_T_from_3_points_shape(self):
        batch_size = 2
        n_res = 5

        x1 = torch.rand((batch_size, n_res, 3))
        x2 = torch.rand((batch_size, n_res, 3))
        x3 = torch.rand((batch_size, n_res, 3))

        t = T.from_3_points(x1, x2, x3)

        rot, tra = t.rots, t.trans

        self.assertTrue(rot.shape == (batch_size, n_res, 3, 3))
        self.assertTrue(torch.all(tra == x2))

    def test_T_from_4x4(self):
        batch_size = 2
        transf = [
            [1, 0, 0, 1],
            [0, 0,-1, 2],
            [0, 1, 0, 3],
            [0, 0, 0, 1],
        ]
        transf = torch.tensor(transf)

        true_rot = transf[:3, :3]
        true_trans = transf[:3, 3]

        transf = torch.stack(
            [transf for _ in range(batch_size)],
            dim=0
        )

        t = T.from_4x4(transf)

        rot, tra = t.rots, t.trans

        self.assertTrue(torch.all(rot == true_rot.unsqueeze(0)))
        self.assertTrue(torch.all(tra == true_trans.unsqueeze(0)))

    def test_T_shape(self):
        batch_size = 2
        n = 5
        transf = T(
            torch.rand((batch_size, n, 3, 3)),
            torch.rand((batch_size, n, 3))
        )

        self.assertTrue(transf.shape == (batch_size, n))

    def test_T_concat(self):
        batch_size = 2
        n = 5
        transf = T(
            torch.rand((batch_size, n, 3, 3)),
            torch.rand((batch_size, n, 3))
        )

        transf_concat = T.concat([transf, transf], dim=0)
        
        self.assertTrue(transf_concat.rots.shape == (batch_size * 2, n, 3, 3))

        transf_concat = T.concat([transf, transf], dim=1)

        self.assertTrue(transf_concat.rots.shape == (batch_size, n * 2, 3, 3))

        self.assertTrue(torch.all(transf_concat.rots[:, :n] == transf.rots))
        self.assertTrue(torch.all(transf_concat.trans[:, :n] == transf.trans))

    def test_T_compose(self):
        trans_1 = [0, 1, 0]
        trans_2 = [0, 0, 1]

        t1 = T(X_90_ROT, torch.tensor(trans_1))
        t2 = T(X_NEG_90_ROT, torch.tensor(trans_2))

        t3 = t1.compose(t2)

        self.assertTrue(torch.all(t3.rots == torch.eye(3)))
        self.assertTrue(torch.all(t3.trans == 0))

    def test_T_apply(self):
        rots = torch.stack([X_90_ROT, X_NEG_90_ROT], dim=0)
        trans = torch.tensor([1, 1, 1])
        trans = torch.stack([trans, trans], dim=0)

        t = T(rots, trans)

        x = torch.arange(30)
        x = torch.stack([x, x], dim=0)
        x = x.view(2, -1, 3) # [2, 10, 3]

        pts = t[..., None].apply(x)

        # All simple consequences of the two x-axis rotations
        self.assertTrue(torch.all(pts[..., 0] == x[..., 0] + 1))
        self.assertTrue(torch.all(pts[0, :, 1] == x[0, :, 2] * -1 + 1))
        self.assertTrue(torch.all(pts[1, :, 1] == x[1, :, 2] + 1))
        self.assertTrue(torch.all(pts[0, :, 2] == x[0, :, 1] + 1))
        self.assertTrue(torch.all(pts[1, :, 2] == x[1, :, 1] * -1 + 1))

    def test_quat_to_rot(self):
        forty_five = math.pi / 4
        quat = torch.tensor([math.cos(forty_five), math.sin(forty_five), 0, 0])
        rot = quat_to_rot(quat)
        eps = 1e-07
        self.assertTrue(torch.all(torch.abs(rot - X_90_ROT) < eps))

    def test_chunk_layer_tensor(self):
        x = torch.rand(2, 4, 5, 15)
        l = torch.nn.Linear(15, 30)
        chunked = chunk_layer(l, {"input": x}, chunk_size=4, no_batch_dims=3)
        unchunked = l(x)

        self.assertTrue(torch.all(chunked == unchunked))

    def test_chunk_layer_dict(self):
        class LinearDictLayer(torch.nn.Linear):
            def forward(self, input):
                out = super().forward(input)
                return {"out": out, "inner": {"out": out + 1}}

        x = torch.rand(2, 4, 5, 15)
        l = LinearDictLayer(15, 30)

        chunked = chunk_layer(l, {"input": x}, chunk_size=4, no_batch_dims=3)
        unchunked = l(x)

        self.assertTrue(torch.all(chunked["out"] == unchunked["out"]))
        self.assertTrue(
            torch.all(chunked["inner"]["out"] == unchunked["inner"]["out"])
        ) 
