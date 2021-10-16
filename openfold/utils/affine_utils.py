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

import numpy as np
import torch


def rot_matmul(a, b):
    row_1 = torch.stack(
        [
            a[..., 0, 0] * b[..., 0, 0]
            + a[..., 0, 1] * b[..., 1, 0]
            + a[..., 0, 2] * b[..., 2, 0],
            a[..., 0, 0] * b[..., 0, 1]
            + a[..., 0, 1] * b[..., 1, 1]
            + a[..., 0, 2] * b[..., 2, 1],
            a[..., 0, 0] * b[..., 0, 2]
            + a[..., 0, 1] * b[..., 1, 2]
            + a[..., 0, 2] * b[..., 2, 2],
        ],
        dim=-1,
    )
    row_2 = torch.stack(
        [
            a[..., 1, 0] * b[..., 0, 0]
            + a[..., 1, 1] * b[..., 1, 0]
            + a[..., 1, 2] * b[..., 2, 0],
            a[..., 1, 0] * b[..., 0, 1]
            + a[..., 1, 1] * b[..., 1, 1]
            + a[..., 1, 2] * b[..., 2, 1],
            a[..., 1, 0] * b[..., 0, 2]
            + a[..., 1, 1] * b[..., 1, 2]
            + a[..., 1, 2] * b[..., 2, 2],
        ],
        dim=-1,
    )
    row_3 = torch.stack(
        [
            a[..., 2, 0] * b[..., 0, 0]
            + a[..., 2, 1] * b[..., 1, 0]
            + a[..., 2, 2] * b[..., 2, 0],
            a[..., 2, 0] * b[..., 0, 1]
            + a[..., 2, 1] * b[..., 1, 1]
            + a[..., 2, 2] * b[..., 2, 1],
            a[..., 2, 0] * b[..., 0, 2]
            + a[..., 2, 1] * b[..., 1, 2]
            + a[..., 2, 2] * b[..., 2, 2],
        ],
        dim=-1,
    )

    return torch.stack([row_1, row_2, row_3], dim=-2)


def rot_vec_mul(r, t):
    x = t[..., 0]
    y = t[..., 1]
    z = t[..., 2]
    return torch.stack(
        [
            r[..., 0, 0] * x + r[..., 0, 1] * y + r[..., 0, 2] * z,
            r[..., 1, 0] * x + r[..., 1, 1] * y + r[..., 1, 2] * z,
            r[..., 2, 0] * x + r[..., 2, 1] * y + r[..., 2, 2] * z,
        ],
        dim=-1,
    )


class T:
    def __init__(self, rots, trans):
        self.rots = rots
        self.trans = trans

        if self.rots is None and self.trans is None:
            raise ValueError("Only one of rots and trans can be None")
        elif self.rots is None:
            self.rots = T.identity_rot(
                self.trans.shape[:-1],
                self.trans.dtype,
                self.trans.device,
                self.trans.requires_grad,
            )
        elif self.trans is None:
            self.trans = T.identity_trans(
                self.rots.shape[:-2],
                self.rots.dtype,
                self.rots.device,
                self.rots.requires_grad,
            )

        if (
            self.rots.shape[-2:] != (3, 3)
            or self.trans.shape[-1] != 3
            or self.rots.shape[:-2] != self.trans.shape[:-1]
        ):
            raise ValueError("Incorrectly shaped input")

    def __getitem__(self, index):
        if type(index) != tuple:
            index = (index,)
        return T(
            self.rots[index + (slice(None), slice(None))],
            self.trans[index + (slice(None),)],
        )

    def __eq__(self, obj):
        return torch.all(self.rots == obj.rots) and torch.all(
            self.trans == obj.trans
        )

    def __mul__(self, right):
        rots = self.rots * right[..., None, None]
        trans = self.trans * right[..., None]

        return T(rots, trans)

    def __rmul__(self, left):
        return self.__mul__(left)

    @property
    def shape(self):
        s = self.rots.shape[:-2]
        return s if len(s) > 0 else torch.Size([1])

    def get_trans(self):
        return self.trans

    def get_rots(self):
        return self.rots

    def compose(self, t):
        rot_1, trn_1 = self.rots, self.trans
        rot_2, trn_2 = t.rots, t.trans

        rot = rot_matmul(rot_1, rot_2)
        trn = rot_vec_mul(rot_1, trn_2) + trn_1

        return T(rot, trn)

    def apply(self, pts):
        r, t = self.rots, self.trans
        rotated = rot_vec_mul(r, pts)
        return rotated + t

    def invert_apply(self, pts):
        r, t = self.rots, self.trans
        pts = pts - t
        return rot_vec_mul(r.transpose(-1, -2), pts)

    def invert(self):
        rot_inv = self.rots.transpose(-1, -2)
        trn_inv = rot_vec_mul(rot_inv, self.trans)

        return T(rot_inv, -1 * trn_inv)

    def unsqueeze(self, dim):
        if dim >= len(self.shape):
            raise ValueError("Invalid dimension")
        rots = self.rots.unsqueeze(dim if dim >= 0 else dim - 2)
        trans = self.trans.unsqueeze(dim if dim >= 0 else dim - 1)

        return T(rots, trans)

    @staticmethod
    def identity_rot(shape, dtype, device, requires_grad):
        rots = torch.eye(
            3, dtype=dtype, device=device, requires_grad=requires_grad
        )
        rots = rots.view(*((1,) * len(shape)), 3, 3)
        rots = rots.expand(*shape, -1, -1)

        return rots

    @staticmethod
    def identity_trans(shape, dtype, device, requires_grad):
        trans = torch.zeros(
            (*shape, 3), dtype=dtype, device=device, requires_grad=requires_grad
        )
        return trans

    @staticmethod
    def identity(shape, dtype, device, requires_grad=True):
        return T(
            T.identity_rot(shape, dtype, device, requires_grad),
            T.identity_trans(shape, dtype, device, requires_grad),
        )

    @staticmethod
    def from_4x4(t):
        rots = t[..., :3, :3]
        trans = t[..., :3, 3]
        return T(rots, trans)

    def to_4x4(self):
        tensor = self.rots.new_zeros((*self.shape, 4, 4))
        tensor[..., :3, :3] = self.rots
        tensor[..., :3, 3] = self.trans
        tensor[..., 3, 3] = 1
        return tensor

    @staticmethod
    def from_tensor(t):
        return T.from_4x4(t)

    @staticmethod
    def from_3_points(p_neg_x_axis, origin, p_xy_plane, eps=1e-8):
        p_neg_x_axis = torch.unbind(p_neg_x_axis, dim=-1)
        origin = torch.unbind(origin, dim=-1)
        p_xy_plane = torch.unbind(p_xy_plane, dim=-1)

        e0 = [c1 - c2 for c1, c2 in zip(origin, p_neg_x_axis)]
        e1 = [c1 - c2 for c1, c2 in zip(p_xy_plane, origin)]

        denom = torch.sqrt(sum((c * c for c in e0)) + eps)
        e0 = [c / denom for c in e0]
        dot = sum((c1 * c2 for c1, c2 in zip(e0, e1)))
        e1 = [c2 - c1 * dot for c1, c2 in zip(e0, e1)]
        denom = torch.sqrt(sum((c * c for c in e1)) + eps)
        e1 = [c / denom for c in e1]
        e2 = [
            e0[1] * e1[2] - e0[2] * e1[1],
            e0[2] * e1[0] - e0[0] * e1[2],
            e0[0] * e1[1] - e0[1] * e1[0],
        ]

        rots = torch.stack([c for tup in zip(e0, e1, e2) for c in tup], dim=-1)
        rots = rots.reshape(rots.shape[:-1] + (3, 3))

        return T(rots, torch.stack(origin, dim=-1))

    @staticmethod
    def concat(ts, dim):
        rots = torch.cat([t.rots for t in ts], dim=dim if dim >= 0 else dim - 2)
        trans = torch.cat(
            [t.trans for t in ts], dim=dim if dim >= 0 else dim - 1
        )

        return T(rots, trans)

    def map_tensor_fn(self, fn):
        """
        Apply a function that takes a tensor as its only argument to the
        rotations and translations, treating the final two/one
        dimension(s), respectively, as batch dimensions.

        E.g.: Given t, an instance of T of shape [N, M], this function can
        be used to sum out the second dimension thereof as follows:

            t = t.map_tensor_fn(lambda x: torch.sum(x, dim=-1))

        The resulting object has rotations of shape [N, 3, 3] and
        translations of shape [N, 3]
        """
        rots = self.rots.view(*self.rots.shape[:-2], 9)
        rots = torch.stack(list(map(fn, torch.unbind(rots, -1))), dim=-1)
        rots = rots.view(*rots.shape[:-1], 3, 3)

        trans = torch.stack(list(map(fn, torch.unbind(self.trans, -1))), dim=-1)

        return T(rots, trans)

    def stop_rot_gradient(self):
        return T(self.rots.detach(), self.trans)

    def scale_translation(self, factor):
        return T(self.rots, self.trans * factor)

    @staticmethod
    def make_transform_from_reference(n_xyz, ca_xyz, c_xyz, eps=1e-20):
        translation = -1 * c_xyz
        n_xyz = n_xyz + translation
        c_xyz = c_xyz + translation

        c_x, c_y, c_z = [c_xyz[..., i] for i in range(3)]
        norm = torch.sqrt(eps + c_x ** 2 + c_y ** 2)
        sin_c1 = -c_y / norm
        cos_c1 = c_x / norm
        zeros = sin_c1.new_zeros(sin_c1.shape)
        ones = sin_c1.new_ones(sin_c1.shape)

        c1_rots = sin_c1.new_zeros((*sin_c1.shape, 3, 3))
        c1_rots[..., 0, 0] = cos_c1
        c1_rots[..., 0, 1] = -1 * sin_c1
        c1_rots[..., 1, 0] = sin_c1
        c1_rots[..., 1, 1] = cos_c1
        c1_rots[..., 2, 2] = 1

        norm = torch.sqrt(eps + c_x ** 2 + c_y ** 2 + c_z ** 2)
        sin_c2 = c_z / norm
        cos_c2 = torch.sqrt(c_x ** 2 + c_y ** 2) / norm

        c2_rots = sin_c2.new_zeros((*sin_c2.shape, 3, 3))
        c2_rots[..., 0, 0] = cos_c2
        c2_rots[..., 0, 2] = sin_c2
        c2_rots[..., 1, 1] = 1
        c1_rots[..., 2, 0] = -1 * sin_c2
        c1_rots[..., 2, 2] = cos_c2

        c_rots = rot_matmul(c2_rots, c1_rots)
        n_xyz = rot_vec_mul(c_rots, n_xyz)

        _, n_y, n_z = [n_xyz[..., i] for i in range(3)]
        norm = torch.sqrt(eps + n_y ** 2 + n_z ** 2)
        sin_n = -n_z / norm
        cos_n = n_y / norm

        n_rots = sin_c2.new_zeros((*sin_c2.shape, 3, 3))
        n_rots[..., 0, 0] = 1
        n_rots[..., 1, 1] = cos_n
        n_rots[..., 1, 2] = -1 * sin_n
        n_rots[..., 2, 1] = sin_n
        n_rots[..., 2, 2] = cos_n

        rots = rot_matmul(n_rots, c_rots)

        rots = rots.transpose(-1, -2)
        translation = -1 * translation

        return T(rots, translation)

    def cuda(self):
        return T(self.rots.cuda(), self.trans.cuda())


_quat_elements = ["a", "b", "c", "d"]
_qtr_keys = [l1 + l2 for l1 in _quat_elements for l2 in _quat_elements]
_qtr_ind_dict = {key: ind for ind, key in enumerate(_qtr_keys)}


def _to_mat(pairs):
    mat = torch.zeros((4, 4))
    for pair in pairs:
        key, value = pair
        ind = _qtr_ind_dict[key]
        mat[ind // 4][ind % 4] = value

    return mat


_qtr_mat = np.zeros((4, 4, 3, 3))
_qtr_mat[..., 0, 0] = _to_mat([("aa", 1), ("bb", 1), ("cc", -1), ("dd", -1)])
_qtr_mat[..., 0, 1] = _to_mat([("bc", 2), ("ad", -2)])
_qtr_mat[..., 0, 2] = _to_mat([("bd", 2), ("ac", 2)])
_qtr_mat[..., 1, 0] = _to_mat([("bc", 2), ("ad", 2)])
_qtr_mat[..., 1, 1] = _to_mat([("aa", 1), ("bb", -1), ("cc", 1), ("dd", -1)])
_qtr_mat[..., 1, 2] = _to_mat([("cd", 2), ("ab", -2)])
_qtr_mat[..., 2, 0] = _to_mat([("bd", 2), ("ac", -2)])
_qtr_mat[..., 2, 1] = _to_mat([("cd", 2), ("ab", 2)])
_qtr_mat[..., 2, 2] = _to_mat([("aa", 1), ("bb", -1), ("cc", -1), ("dd", 1)])


def quat_to_rot(quat):  # [*, 4]
    # [*, 4, 4]
    quat = quat[..., None] * quat[..., None, :]

    # [4, 4, 3, 3]
    mat = quat.new_tensor(_qtr_mat)

    # [*, 4, 4, 3, 3]
    shaped_qtr_mat = mat.view((1,) * len(quat.shape[:-2]) + mat.shape)
    quat = quat[..., None, None] * shaped_qtr_mat

    # [*, 3, 3]
    return torch.sum(quat, dim=(-3, -4))


def affine_vector_to_4x4(vector):
    quats = vector[..., :4]
    trans = vector[..., 4:]

    four_by_four = vector.new_zeros((*vector.shape[:-1], 4, 4))
    four_by_four[..., :3, :3] = quat_to_rot(quats)
    four_by_four[..., :3, 3] = trans
    four_by_four[..., 3, 3] = 1

    return four_by_four
