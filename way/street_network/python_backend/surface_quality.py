"""Fixed-boundary planar grading followed by authoritative surface evaluation."""
import math
import numpy as np
from .mesh import Mesh
from .quality_correction import QuadShapeCorrection


class ProfiledQuadCorrection:
    """Preserve physical seams and feature guides while changing internal XY.

    Heights are evaluated from the feature's own field, never inferred from a
    flattened mesh or retained at vertices whose XY has changed. This utility
    adds neither vertices nor faces and does not change quad connectivity.
    """
    def __init__(self, *, max_cells=50000, max_iterations=1000):
        self.planar = QuadShapeCorrection(max_cells=max_cells, max_iterations=max_iterations)

    def improve(self, mesh, height_field, *, fixed_vertices=(), datum, ceiling,
                exaggeration, maximum_slope):
        if (not callable(height_field) or not all(math.isfinite(v) for v in
                (datum, ceiling, exaggeration, maximum_slope)) or ceiling < datum
                or exaggeration < 0. or maximum_slope <= 0.):
            raise ValueError('Profiled quad correction needs a finite height envelope and slope limit')
        xyz = mesh.xyz.copy()
        xyz[:, 2] = 0.
        candidate, details = self.planar.improve(Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams)),
                                                fixed_vertices=fixed_vertices)
        if not details['qualified']:
            raise ValueError('Profiled quad correction cannot meet its unchanged planar quality limits')
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        fixed = np.unique(np.r_[unique[counts == 1].ravel(), np.asarray(fixed_vertices, dtype=np.int32)])
        heights = np.asarray(height_field(candidate.xyz[:, :2]), dtype=float)
        if heights.shape != (len(candidate.xyz),) or not np.all(np.isfinite(heights)):
            raise ValueError('Profiled quad correction received invalid surface heights')
        candidate.xyz[:, 2] = heights
        candidate.xyz[fixed] = mesh.xyz[fixed]
        if candidate.xyz[:, 2].min() < datum - 1e-8 or candidate.xyz[:, 2].max() > ceiling + 1e-8:
            raise ValueError('Corrected pavement exceeds its original elevation envelope')
        quality = candidate.validate()
        slopes = candidate.slope_statistics(exaggeration)
        if slopes['maximum'] > maximum_slope + 1e-8:
            raise ValueError('Corrected pavement exceeds its physical slope limit')
        return candidate, dict(planar=details, quality=quality, physical_slopes=slopes,
            measured_slope_limit=maximum_slope, exact_fixed_xyz=bool(np.array_equal(candidate.xyz[fixed], mesh.xyz[fixed])),
            construction='fixed seam/crown guides with authoritative profile field reevaluation')


class ProfiledSlopeCorrection:
    """Repair interpolation amplification caused by nearly collinear quads.

    Keep every physical/cut boundary and all connectivity exact. Only interior
    XY is adjusted; every moved height is evaluated from the original field.
    Stronger angle targets are bounded alternatives, never relaxed slope limits.
    """
    def __init__(self, *, max_cells=50000, max_iterations=500):
        self.max_cells, self.max_iterations = max_cells, max_iterations

    def improve(self, mesh, height_field, *, datum, ceiling, maximum_slope,
                longitudinal_bound, maximum_longitudinal_slope):
        if (not callable(height_field) or not callable(longitudinal_bound)
                or not all(math.isfinite(v) for v in
                    (datum, ceiling, maximum_slope, maximum_longitudinal_slope))
                or ceiling < datum or min(maximum_slope, maximum_longitudinal_slope) <= 0.):
            raise ValueError('Profile slope correction needs finite fields and positive slope limits')
        mesh.validate()
        if len(mesh.quads) > self.max_cells:
            raise ValueError('Profile slope correction exceeds its cell budget')
        before = float(mesh.quad_slopes().max())
        before_along = float(longitudinal_bound(mesh))
        if before <= maximum_slope+1e-8 and before_along <= maximum_longitudinal_slope+1e-8:
            return mesh, dict(applied=False, physical_slope=before, longitudinal_bound=before_along)
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if counts.max() > 2:
            raise ValueError('Profile slope correction requires manifold edge incidence')
        fixed = np.unique(unique[counts == 1])
        flat = mesh.xyz.copy()
        flat[:, 2] = 0.
        failures = []
        for minimum, target in ((.15, .2), (.3, .4)):
            candidate, shape = QuadShapeCorrection(max_cells=self.max_cells,
                max_iterations=self.max_iterations, minimum_sine=minimum,
                target_sine=target).improve(Mesh(flat.copy(), mesh.quads, mesh.material, dict(mesh.seams)))
            heights = np.asarray(height_field(candidate.xyz[:, :2]), dtype=float)
            if heights.shape != (len(candidate.xyz),) or not np.all(np.isfinite(heights)):
                raise ValueError('Profile slope correction received invalid field heights')
            candidate.xyz[:, 2] = heights
            candidate.xyz[fixed] = mesh.xyz[fixed]
            if candidate.xyz[:, 2].min() < datum-1e-8 or candidate.xyz[:, 2].max() > ceiling+1e-8:
                raise ValueError('Profile slope correction exceeds its original elevation envelope')
            quality = candidate.validate()
            maximum = float(candidate.quad_slopes().max())
            along = float(longitudinal_bound(candidate))
            if (shape['qualified'] and quality['maximum_edge_ratio'] <= 8.
                    and maximum <= maximum_slope+1e-8 and along <= maximum_longitudinal_slope+1e-8):
                return candidate, dict(applied=True, before_physical_slope=before,
                    physical_slope=maximum, before_longitudinal_bound=before_along,
                    longitudinal_bound=along, measured_slope_limit=maximum_slope,
                    longitudinal_slope_limit=maximum_longitudinal_slope, quality=quality,
                    planar=shape, exact_fixed_xyz=bool(np.array_equal(candidate.xyz[fixed], mesh.xyz[fixed])),
                    construction='fixed boundary interior angle repair with original profile field')
            failures.append(dict(minimum_sine=minimum, qualified_shape=shape['qualified'],
                                 physical_slope=maximum, longitudinal_bound=along))
        raise ValueError('Profile slope correction cannot meet unchanged physical limits: '+str(failures))
