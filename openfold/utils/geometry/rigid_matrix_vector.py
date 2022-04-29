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
"""Rigid3Array Transformations represented by a Matrix and a Vector."""

from __future__ import annotations
import dataclasses
from typing import Union, List

import torch

from openfold.utils.geometry import rotation_matrix
from openfold.utils.geometry import struct_of_array
from openfold.utils.geometry import vector


Float = Union[float, torch.Tensor]


@dataclasses.dataclass(frozen=True)
class Rigid3Array:
  """Rigid Transformation, i.e. element of special euclidean group."""

  rotation: rotation_matrix.Rot3Array
  translation: vector.Vec3Array

  def __matmul__(self, other: Rigid3Array) -> Rigid3Array:
    new_rotation = self.rotation @ other.rotation # __matmul__
    new_translation = self.apply_to_point(other.translation)
    return Rigid3Array(new_rotation, new_translation)

  def __getitem__(self, index) -> Rigid3Array:
    return Rigid3Array(
        self.rotation[index],
        self.translation[index],
    )

  def __mul__(self, other: torch.Tensor) -> Rigid3Array:
    return Rigid3Array(
        self.rotation * other,
        self.translation * other,
    )

  def map_tensor_fn(self, fn) -> Rigid3Array:
    return Rigid3Array(
        self.rotation.map_tensor_fn(fn),
        self.translation.map_tensor_fn(fn),
    )

  def inverse(self) -> Rigid3Array:
    """Return Rigid3Array corresponding to inverse transform."""
    inv_rotation = self.rotation.inverse()
    inv_translation = inv_rotation.apply_to_point(-self.translation)
    return Rigid3Array(inv_rotation, inv_translation)

  def apply_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Apply Rigid3Array transform to point."""
    return self.rotation.apply_to_point(point) + self.translation

  def apply_inverse_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Apply inverse Rigid3Array transform to point."""
    new_point = point - self.translation
    return self.rotation.apply_inverse_to_point(new_point)

  def compose_rotation(self, other_rotation):
    rot = self.rotation @ other_rotation
    return Rigid3Array(rot, trans.clone())

  @classmethod
  def identity(cls, shape, device) -> Rigid3Array:
    """Return identity Rigid3Array of given shape."""
    return cls(
        rotation_matrix.Rot3Array.identity(shape, device),
        vector.Vec3Array.zeros(shape, device)
    )

  @classmethod
  def cat(cls, rigids: List[Rigid3Array], dim: int) -> Rigid3Array:
    return cls(
      Rot3Array.cat([r.rotation for r in rigids], dim=dim),
      Vec3Array.cat([r.translation for r in rigids], dim=dim),
    ) 

  def scale_translation(self, factor: Float) -> Rigid3Array:
    """Scale translation in Rigid3Array by 'factor'."""
    return Rigid3Array(self.rotation, self.translation * factor)

  def to_array(self):
    rot_array = self.rotation.to_array()
    vec_array = self.translation.to_array()
    return torch.cat([rot_array, vec_array[..., None]], dim=-1)

  def reshape(self, new_shape) -> Rigid3Array:
    rots = self.rotation.reshape(new_shape)
    trans = self.translation.reshape(new_shape)
    return Rigid3Aray(rots, trans)

  @classmethod
  def from_array(cls, array):
    rot = rotation_matrix.Rot3Array.from_array(array[..., :3])
    vec = vector.Vec3Array.from_array(array[..., -1])
    return cls(rot, vec)

  @classmethod
  def from_tensor_4x4(cls, array):
      return cls.from_array(array)

  @classmethod
  def from_array4x4(cls, array: torch.tensor) -> Rigid3Array:
    """Construct Rigid3Array from homogeneous 4x4 array."""
    rotation = rotation_matrix.Rot3Array(
        array[..., 0, 0], array[..., 0, 1], array[..., 0, 2],
        array[..., 1, 0], array[..., 1, 1], array[..., 1, 2],
        array[..., 2, 0], array[..., 2, 1], array[..., 2, 2]
    )
    translation = vector.Vec3Array(
        array[..., 0, 3], array[..., 1, 3], array[..., 2, 3])
    return cls(rotation, translation)
