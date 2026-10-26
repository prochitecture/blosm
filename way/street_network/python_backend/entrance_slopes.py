"""Conservative longitudinal-slope bounds along complete quad entrance edges.

For a positive bilinear quad, edge gradients are convex combinations of the
two corner gradients. Source offset velocities vary affinely along a profile.
Their unit directions stay in the cone between the two endpoint directions.
Bounding both corner gradients over that cone bounds the whole entrance edge;
checking only each corner's own direction can miss an interior maximum.
"""
from dataclasses import dataclass
import math
import numpy as np
from .mesh import Mesh, cross
from .paths import right
from .slope_projection import corner_matrices


@dataclass(frozen=True)
class EntranceSlopeGuard:
    limit: float = .03
    tolerance: float = 1e-8

    def __post_init__(self):
        if not all(math.isfinite(v) and v > 0. for v in (self.limit, self.tolerance)):
            raise ValueError('Invalid entrance slope limit or tolerance')

    def measure(self, mesh, boundary_ids, velocities, physical_heights):
        """Return an edge-wide bound, independently of display height scaling."""
        velocities = np.asarray(velocities, dtype=float)
        heights = np.asarray(physical_heights, dtype=float)
        if velocities.shape != (len(mesh.xyz), 2) or heights.shape != (len(mesh.xyz),):
            raise ValueError('Entrance fields must address every mesh vertex')
        if not np.all(np.isfinite(velocities)) or not np.all(np.isfinite(heights)):
            raise ValueError('Entrance fields must be finite')
        boundary = np.zeros(len(mesh.xyz), dtype=bool)
        boundary[np.asarray(boundary_ids, dtype=np.intp)] = True
        hit = boundary[mesh.quads]
        face_ids, corners = np.nonzero(hit & np.roll(hit, -1, axis=1))
        if not len(face_ids):
            raise ValueError('Entrance boundary has no complete quad edges')
        unique_faces, inverse = np.unique(face_ids, return_inverse=True)
        selected = Mesh(mesh.xyz, mesh.quads[unique_faces])
        gradients = np.einsum('fcij,fj->fci', corner_matrices(selected), heights[selected.quads])
        next_corners = (corners + 1) % 4
        g = np.stack((gradients[inverse, corners], gradients[inverse, next_corners]), axis=1)
        first = velocities[mesh.quads[face_ids, corners]]
        last = velocities[mesh.quads[face_ids, next_corners]]
        first_length = np.linalg.norm(first, axis=1)
        last_length = np.linalg.norm(last, axis=1)
        if np.any(np.minimum(first_length, last_length) < 1e-12):
            raise ValueError('Entrance has an undefined source offset tangent')
        first, last = first / first_length[:, None], last / last_length[:, None]
        axis = first + last
        length = np.linalg.norm(axis, axis=1)
        if np.any(length < 1e-10):
            raise ValueError('Entrance source offset tangents oppose each other')
        axis /= length[:, None]
        normal = right(axis)
        cosine = np.clip(np.sum(first * axis, axis=1), 0., 1.)[:, None]
        sine = abs(cross(first, axis))[:, None]
        along = abs(np.sum(g * axis[:, None, :], axis=-1))
        across = abs(np.sum(g * normal[:, None, :], axis=-1))
        # The maximum is at a cone endpoint, unless +/-gradient lies inside
        # the cone, where the complete gradient norm is attained.
        bound = np.where(across * cosine <= along * sine,
                         np.linalg.norm(g, axis=-1), along * cosine + across * sine)
        return float(bound.max())

    def validate(self, mesh, boundary_ids, velocities, physical_heights):
        maximum = self.measure(mesh, boundary_ids, velocities, physical_heights)
        if maximum > self.limit + self.tolerance:
            raise ValueError(f'Entrance longitudinal slope bound {100*maximum:.6g}% exceeds {100*self.limit:g}%')
        return maximum
