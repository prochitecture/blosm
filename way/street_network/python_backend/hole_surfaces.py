"""Grass surfaces and physical queries from their actual shared curb heights.

A closed pavement hole may be an ordinary road block, a median plateau, or a
median with lowered caps. Its topology alone supplies no elevation policy.
"""
from dataclasses import dataclass
import math
import numpy as np

from .fixed_grass_heights import FixedGrassCornerSlopeError, boundary_mask, reconcile
from .mesh import Mesh
from .publication_slopes import PublicationSlopeGuard
from .profiles import plateau_extension
from .surface_query import QuadHeightField

SOURCE_HOLE_HEIGHT_POLICY = '2-actual-source-hole-fixed-corners'


@dataclass(frozen=True)
class ConstantHoleHeight:
    value: float
    datum: float
    ceiling: float

    def __post_init__(self):
        if (not all(map(math.isfinite, (self.value, self.datum, self.ceiling)))
                or self.ceiling < self.datum or not self.datum-1e-10 <= self.value <= self.ceiling+1e-10):
            raise ValueError('Invalid constant source-hole height bounds')

    def __call__(self, xy):
        xy = np.asarray(xy, dtype=float)
        if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
            raise ValueError('A constant hole query requires finite XY points')
        return np.full(len(xy), self.value)


@dataclass(frozen=True)
class HoleSurface:
    mesh: object
    field: object
    metadata: dict


def build(template, pavement, loop, *, datum, ceiling, exaggeration=1., max_cells=50000):
    options = dict(datum=datum, ceiling=ceiling, exaggeration=exaggeration, max_cells=max_cells)
    try:
        return _build(template, pavement, loop, **options)
    except FixedGrassCornerSlopeError as failure:
        faces = failure.faces
        if len(template.quads) + 4*len(faces) > max_cells:
            raise
        rim = pavement.xyz[pavement.seams[loop]]
        delta = np.roll(rim, -1, axis=0) - rim
        if np.any(abs(delta[:, 2]) > (.1 + PublicationSlopeGuard.tolerance)
                  * np.linalg.norm(delta[:, :2], axis=1)):
            raise
    # A quad with three fixed rim vertices can impose a steep tangent plane
    # even if each rim edge has a feasible grade. An interior ring gives the
    # two edges independent adjacent faces while preserving their exact seam.
    refined = _inset_faces(template, faces)
    result = _build(refined, pavement, loop, **options)
    return HoleSurface(result.mesh, result.field, dict(result.metadata, exact_xy=False, exact_quads=False,
        refined_fixed_corner_faces=tuple(map(int, faces)), original_quads=len(template.quads)))


def _inset_faces(mesh, faces):
    outer = mesh.quads[faces]
    inner = np.arange(len(mesh.xyz), len(mesh.xyz) + outer.size).reshape(outer.shape)

    def inset(values):
        center = values.mean(axis=1, keepdims=True)
        return center + .35 * (values - center)

    def ring(outer, inner):
        return np.stack((outer, np.roll(outer, -1, axis=1), np.roll(inner, -1, axis=1), inner), axis=2)

    xyz = np.vstack((mesh.xyz, inset(mesh.xyz[outer]).reshape(-1, 3)))
    quads = mesh.quads.copy()
    quads[faces] = inner
    quads = np.vstack((quads, ring(outer, inner).reshape(-1, 4)))
    layers = {}
    for name, values in mesh.uv_layers.items():
        base = values.copy()
        base[faces] = inset(values[faces])
        layers[name] = np.concatenate((base, ring(values[faces], base[faces]).reshape(-1, 4, 2)))
    return Mesh(xyz, quads, mesh.material, dict(mesh.seams), layers)


def _build(template, pavement, loop, *, datum, ceiling, exaggeration, max_cells):
    """Build physical Z without changing a planar chart or a shared rim.

The pavement is still in physical metres. Display exaggeration is used only
for stored-coordinate certification here; the caller applies it once later.
"""
    if (not all(map(math.isfinite, (datum, ceiling, exaggeration))) or ceiling < datum
            or exaggeration < 0. or type(max_cells) is not int or max_cells < 1):
        raise ValueError('Invalid source-hole surface bounds')
    template.validate()
    if len(template.quads) > max_cells:
        raise ValueError('Source-hole surface exceeds its cell budget')
    if set(template.seams) != {0} or loop not in pavement.seams:
        raise ValueError('Source-hole surface requires its actual shared pavement rim')
    rim_ids = np.asarray(template.seams[0])
    boundary = pavement.xyz[pavement.seams[loop]][::-1]
    if (len(rim_ids) != len(boundary)
            or not np.array_equal(template.xyz[rim_ids, :2], boundary[:, :2])):
        raise ValueError('Source-hole and pavement physical boundaries do not match exactly')
    fixed = boundary_mask(template)
    if set(np.flatnonzero(fixed)) != set(map(int, rim_ids)):
        raise ValueError('Source-hole rim omits a physical boundary vertex')
    if (not np.isfinite(boundary).all() or np.any(boundary[:, 2] < datum-1e-10)
            or np.any(boundary[:, 2] > ceiling+1e-10)):
        raise ValueError('Source-hole curb heights exceed the incoming elevation bounds')

    xyz = template.xyz.copy()
    if np.all(boundary[:, 2] == boundary[0, 2]):
        value = float(boundary[0, 2])
        xyz[:, 2] = value
        result = template.with_positions(xyz)
        field = ConstantHoleHeight(value, datum, ceiling)
        details = dict(policy='CONSTANT_ACTUAL_SOURCE_CURB', constant_height=value,
                       maximum_physical_slope=0., maximum_stored_physical_slope=0.)
    else:
        height = max(0., float(boundary[:, 2].max())-datum)
        xyz[:, 2], extension = plateau_extension(xyz[:, :2], boundary, datum, height)
        xyz[rim_ids, 2] = boundary[:, 2]
        physical_seed = template.with_positions(xyz)
        display_xyz = xyz.copy()
        display_xyz[:, 2] = datum + exaggeration*(xyz[:, 2]-datum)
        display = template.with_positions(display_xyz)
        _, _, physical, reconciled = reconcile(display, datum=datum, ceiling=ceiling,
            exaggeration=exaggeration, physical_mesh=physical_seed, max_cells=max_cells,
            policy='ACTUAL_SOURCE_CURB_HOLE')
        result = template.with_positions(physical.xyz)
        field = QuadHeightField(result, datum=datum, ceiling=ceiling, max_cells=max_cells,
            metadata=dict(source_hole_heights=True))
        details = dict(policy='BOUNDARY_DRIVEN_SOURCE_CURB', extension=extension,
                       reconciliation=reconciled)

    if (not np.array_equal(result.xyz[:, :2], template.xyz[:, :2])
            or not np.array_equal(result.quads, template.quads)
            or not np.array_equal(result.xyz[rim_ids], boundary)):
        raise ValueError('Source-hole construction changed its fixed chart or curb')
    physical_quality = result.validate()
    stored = PublicationSlopeGuard.physical_buffer(result, datum, exaggeration if exaggeration else 1.)
    stored_quality = stored.validate()
    if max(physical_quality['maximum_edge_ratio'], stored_quality['maximum_edge_ratio']) > 8.:
        raise ValueError('Source-hole surface fails unchanged physical or stored grass quality')
    return HoleSurface(result, field, dict(details, exact_xy=True, exact_quads=True,
        exact_shared_curb=True, physical_quality=physical_quality, stored_quality=stored_quality))
