"""Geometry records and footprint admission shared by building backends."""
from dataclasses import dataclass, field
import math
from typing import Mapping

from .errors import PadlError

EPSILON = 1e-6


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def signed_area(ring):
    return sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(ring, ring[1:] + ring[:1])) / 2


def contains(point, ring):
    inside = False
    x, y = point[:2]
    for a, b in zip(ring, ring[1:] + ring[:1]):
        if (a[1] > y) != (b[1] > y) and x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0]:
            inside = not inside
    return inside


def on_segment(a, b, p):
    return abs(cross(a, b, p)) <= EPSILON and all(min(a[i], b[i]) - EPSILON <= p[i] <= max(a[i], b[i]) + EPSILON for i in (0, 1))


def intersects(a, b, c, d):
    ca, cb, cc, cd = cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b)
    return (ca * cb < 0 and cc * cd < 0) or any((on_segment(a, b, c), on_segment(a, b, d), on_segment(c, d, a), on_segment(c, d, b)))


@dataclass(frozen=True)
class Footprint:
    identity: str
    footprint_id: int
    outer: tuple
    holes: tuple = ()
    attributes: Mapping = field(default_factory=dict)

    def validate(self):
        if not isinstance(self.identity, str) or not self.identity or type(self.footprint_id) is not int or not 0 < self.footprint_id < 2**31:
            raise PadlError('Footprint requires a stable identity and positive signed-32-bit association ID', code='invalid_identity')
        rings = (self.outer, *self.holes)
        if sum(len(ring) for ring in rings) > 256 or len(self.holes) > 16:
            raise PadlError('Initial roof profile admits at most 256 corners and 16 holes', code='geometry_limit')
        edges = []
        for ri, ring in enumerate(rings):
            if len(ring) < 3 or any(len(v) not in (2, 3) or
                    not all(type(x) in (int, float) and math.isfinite(x) and abs(x) <= 100000 for x in v) for v in ring):
                raise PadlError('Invalid footprint coordinates or ring size', code='invalid_geometry')
            if abs(signed_area(ring)) <= EPSILON:
                raise PadlError('Degenerate footprint ring', code='invalid_geometry')
            for i, (a, b) in enumerate(zip(ring, ring[1:] + ring[:1])):
                if math.dist(a[:2], b[:2]) <= EPSILON:
                    raise PadlError('Coincident adjacent footprint corners', code='invalid_geometry')
                edges.append((ri, i, a, b))
        for j, (ri, i, a, b) in enumerate(edges):
            for rj, k, c, d in edges[j + 1:]:
                if ri == rj and (abs(i - k) == 1 or abs(i - k) == len(rings[ri]) - 1):
                    continue
                if intersects(a, b, c, d):
                    raise PadlError('Footprint rings intersect or touch', code='invalid_geometry')
        for i, hole in enumerate(self.holes):
            if not contains(hole[0], self.outer) or any(contains(hole[0], other) or contains(other[0], hole) for other in self.holes[:i]):
                raise PadlError('Holes must be disjoint and strictly inside the outer ring', code='invalid_geometry')
        return self

    def oriented_rings(self):
        """Private winding copies for polygon kernels; source rings stay untouched."""
        return tuple(tuple(ring if (signed_area(ring) > 0) == (i == 0) else reversed(ring))
                     for i, ring in enumerate((self.outer, *self.holes)))

    @property
    def area(self):
        return abs(signed_area(self.outer)) - sum(abs(signed_area(ring)) for ring in self.holes)


@dataclass
class MeshData:
    vertices: list = field(default_factory=list)
    faces: list = field(default_factory=list)
    component_types: list = field(default_factory=list)
    footprint_ids: list = field(default_factory=list)

    def polygon(self, points, component, footprint_id):
        points = list(points)
        if len(points) < 3:
            return
        a = points[0]
        area = sum(math.sqrt(sum(x * x for x in (
            (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]),
            (b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2]),
            (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))))
                   for b, c in zip(points[1:-1], points[2:])) / 2
        if area <= 1e-10:
            return
        first = len(self.vertices)
        self.vertices.extend(tuple(v) for v in points)
        self.faces.append(tuple(range(first, first + len(points))))
        self.component_types.append(component)
        self.footprint_ids.append(footprint_id)

    def extend(self, other):
        offset = len(self.vertices)
        self.vertices.extend(other.vertices)
        self.faces.extend(tuple(i + offset for i in face) for face in other.faces)
        self.component_types.extend(other.component_types)
        self.footprint_ids.extend(other.footprint_ids)

    def to_blender_mesh(self, name='PADL Buildings'):
        import bpy
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(self.vertices, [], self.faces)
        for name, values in (('footprint_id', self.footprint_ids), ('face_type', [4] * len(self.faces)),
                             ('component_type', self.component_types)):
            attr = mesh.attributes.new(name, 'INT', 'FACE')
            attr.data.foreach_set('value', values)
        mesh.update()
        return mesh
