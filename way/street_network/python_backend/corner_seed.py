"""Alternative planar corner topology for a rejected single-junction chart.

This module supplies no surface heights. Every physical rim sample survives;
the compound boundary-jet, slope, drainage and publication solvers remain
responsible for reconciling the incoming profiles.
"""
import numpy as np
from .intersections import IntersectionProvider
from .mesh import Mesh
from .partition import BoundaryRef, CompoundChartResult
from .quality_correction import QuadShapeCorrection


class StagedCornerBands(IntersectionProvider):
    @staticmethod
    def pattern(context, layout, side):
        first = layout.approaches[side].profile.intervals(-1, layout.radial)
        last = layout.approaches[(side + 1) % len(layout.approaches)].profile.intervals(1, layout.radial)
        if not first or not last:
            return IntersectionProvider.pattern(context, layout, side)
        end_rows = max(1, min(layout.rows // 2 - 1, layout.rows // 10))
        return context.patterns.staged_strip(layout.rows, first, last,
            middle=3 * max(first, last), end_rows=end_rows)


class CornerBandSeed:
    def __init__(self, patterns, *, max_cells=50000, matching_work=4_000_000):
        if (type(max_cells) is not int or max_cells < 4
                or type(matching_work) is not int or matching_work < 1):
            raise ValueError('Corner seed needs positive cell and rim-matching budgets')
        self.patterns, self.max_cells, self.matching_work = patterns, max_cells, matching_work

    def build(self, group, context, partitioner):
        if len(group.members) != 1 or len(group.footprints) != 1 or len(group.loops) != 1 or group.loops[0].hole:
            raise ValueError('Corner seed requires one ordinary junction exterior')
        junction = next(j for j in context.network.intersections if j.key == group.members[0])
        layout = context[junction.layout_key]
        if layout.mesh_strategy != 'bands' or not layout.corners:
            raise ValueError('Corner seed requires a compatible corner-band footprint')
        maximum = max(approach.profile.intervals(side, layout.radial)
                      for approach in layout.approaches for side in (-1, 1))
        estimate = len(layout.approaches) * (layout.rows ** 2 // 4
            + layout.rows * 3 * maximum + 12 * maximum)
        if estimate > self.max_cells:
            raise ValueError('Corner seed exceeds its cell budget')
        raw = StagedCornerBands.planar_chart(context, junction).mesh
        mesh, correction = QuadShapeCorrection(max_cells=self.max_cells, max_iterations=1000).improve(raw)
        if not correction['qualified']:
            raise ValueError('Corner seed cannot meet the unchanged quad-shape limits')
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if counts.max() != 2 or len(mesh.xyz) - len(unique) + len(mesh.quads) != 1:
            raise ValueError('Corner seed must be a conforming quad disk')
        exterior_edges = unique[counts == 1]
        boundary = np.unique(exterior_edges)
        loop = group.loops[0]
        if len(boundary) * len(loop.xy) > self.matching_work:
            raise ValueError('Corner seed exceeds its rim-matching work budget')
        refs = {}
        for vertex in boundary:
            distances = np.linalg.norm(loop.xy - mesh.xyz[vertex, :2], axis=1)
            at = int(np.argmin(distances))
            if distances[at] > 1e-8:
                raise ValueError('Corner seed does not retain the authoritative rim partition')
            refs[int(vertex)] = BoundaryRef(0, at, 0.)
        rim = []
        for point in loop.xy:
            distances = np.linalg.norm(mesh.xyz[boundary, :2] - point, axis=1)
            at = int(np.argmin(distances))
            if distances[at] > 1e-8:
                raise ValueError('Corner seed lost a physical boundary vertex')
            rim.append(int(boundary[at]))
        if len(rim) != len(boundary) or len(set(rim)) != len(rim):
            raise ValueError('Corner seed boundary mapping is not bijective')
        expected = {tuple(sorted((int(a), int(b)))) for a, b in zip(rim, (*rim[1:], rim[0]))}
        if expected != {tuple(edge) for edge in exterior_edges}:
            raise ValueError('Corner seed lost a complete physical boundary edge')
        mesh.seams = {0: np.asarray(rim, dtype=np.int32)}
        partition = partitioner.build(group, context)
        return CompoundChartResult(mesh, {}, partition, None, refs,
            dict(planar_only=True, corner_band_seed=True, exact_original_boundary=True,
                 quality=mesh.validate(), total_quads=len(mesh.quads), holes=0,
                 correction=correction, partition_cuts_used=False,
                 construction='staged corner bands with fixed physical rim and regraded interior'))
