"""Bounded grass heights on an explicitly selected fixed physical boundary.

The published curbs, free profile, computational exterior connector, XY and
quad topology stay fixed. Callers retain their established successful policy
first. The rounded grass slope limit remains ten percent.
"""
import numpy as np
from .mesh import Mesh
from .drainage import BoundaryDrainage
from .slope_projection import SlopeProjection
from .publication_slopes import PublicationSlopeGuard
from .surface_query import QuadHeightField


class FixedGrassCornerSlopeError(ValueError):
    def __init__(self, faces):
        super().__init__('Fixed grass fixed boundary cannot meet the ten-percent grass slope limit')
        self.faces = faces


def boundary_mask(mesh):
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)),
                             axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    if not len(counts) or counts.max() > 2:
        raise ValueError('Fixed grass needs a manifold quad boundary')
    fixed = np.zeros(len(mesh.xyz), dtype=bool)
    fixed[np.unique(unique[counts == 1])] = True
    if not fixed.any():
        raise ValueError('Fixed grass needs fixed boundary outlets')
    return fixed


def reconcile(mesh, *, datum, ceiling, exaggeration=1., physical_mesh=None,
              max_cells=50000, policy='fixed_boundary_grass'):
    """Return display mesh, matching query and proof from original physical Z.

    A non-unit display requires an explicit, independently constructed physical
    buffer. This includes E0; erased display heights are never inverted.
    """
    original = mesh
    if exaggeration != 1. and physical_mesh is None:
        raise ValueError('Fixed grass requires its original physical height buffer')
    physical = mesh if physical_mesh is None else physical_mesh
    if (not np.array_equal(physical.xyz[:, :2], mesh.xyz[:, :2])
            or not np.array_equal(physical.quads, mesh.quads)):
        raise ValueError('Fixed grass physical reconstruction changed its chart')
    fixed = boundary_mask(physical)
    expected = datum + exaggeration * (physical.xyz[fixed, 2] - datum)
    if not np.allclose(expected, mesh.xyz[fixed, 2], atol=1e-12, rtol=0.):
        raise ValueError('Fixed grass physical reconstruction changed its fixed profiles')
    initial = physical
    floors_by_corner = PublicationSlopeGuard.fixed_gradient_floors(physical, fixed, datum)
    floor = float(floors_by_corner.max())
    if floor > .1 + PublicationSlopeGuard.tolerance:
        raise FixedGrassCornerSlopeError(np.flatnonzero(
            floors_by_corner.max(axis=1) > .1 + PublicationSlopeGuard.tolerance))
    floors = np.full(len(mesh.xyz), datum)
    rounds = []
    target = floor + .99 * (.1-floor)
    for _ in range(8):
        projected, projection = SlopeProjection(max_cells=max_cells).solve(physical, fixed,
            datum=datum, ceiling=ceiling, maximum_slope=target,
            feasibility_slope=.1, minimum_heights=floors)
        filled, drainage = BoundaryDrainage(max_cells=max_cells).fill(projected, fixed)
        rounds.append(dict(projection=projection, drainage=drainage))
        if drainage['maximum_fill'] <= 1e-6:
            break
        raised = filled > projected.xyz[:, 2] + 1e-10
        floors[raised] = np.maximum(floors[raised], filled[raised])
        physical = projected
    else:
        raise ValueError('Fixed grass exhausted its drainage reconciliation budget')
    field = QuadHeightField(projected, datum=datum, ceiling=ceiling, max_cells=max_cells)
    publication_scale = exaggeration if exaggeration else 1.
    physical, field, publication = PublicationSlopeGuard().correct(projected, field, fixed,
        exaggeration=publication_scale, limit=.1, max_cells=max_cells,
        check_joins=lambda candidate: None)
    if (not np.array_equal(physical.xyz[fixed], initial.xyz[fixed])
            or not np.array_equal(physical.xyz[:, :2], initial.xyz[:, :2])
            or not np.array_equal(physical.quads, initial.quads)):
        raise ValueError('Fixed grass reconciliation changed its fixed geometry')
    xyz = physical.xyz.copy()
    xyz[:, 2] = datum + exaggeration * (xyz[:, 2] - datum)
    xyz[fixed] = original.xyz[fixed]
    display = Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams))
    stored = PublicationSlopeGuard.physical_buffer(physical, datum, publication_scale)
    stored_quality = stored.validate()
    if stored_quality['maximum_edge_ratio'] > 8.:
        raise ValueError('Fixed grass failed its unchanged stored quad-quality limit')
    maximum = float(stored.quad_slopes().max())
    final_drainage = BoundaryDrainage(max_cells=max_cells).fill(physical, fixed)[1]
    if maximum > .1 + PublicationSlopeGuard.tolerance or final_drainage['maximum_fill'] > 1e-6:
        raise ValueError('Fixed grass failed its unchanged slope/drainage limits')
    report = dict(policy=policy, slope_limit=.1,
        fixed_gradient_floor=floor, maximum_physical_slope=float(physical.quad_slopes().max()),
        maximum_stored_physical_slope=maximum,
        maximum_interior_height_change=float(abs(physical.xyz[:, 2]-initial.xyz[:, 2]).max()),
        exact_xy=True, exact_quads=True, exact_boundary=True,
        explicit_physical_reconstruction=physical_mesh is not None,
        physical_publication_scale=publication_scale, stored_quality=stored_quality,
        rounds=rounds, publication=publication, drainage=final_drainage)
    query = QuadHeightField(display, datum=datum,
        ceiling=datum + exaggeration * (ceiling-datum), max_cells=max_cells,
        metadata=dict(fixed_grass_reconciliation=report))
    return display, query, physical, report
