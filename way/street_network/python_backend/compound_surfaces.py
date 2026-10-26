"""Profiled direct-quad surfaces with bounded optional height projection."""
from dataclasses import dataclass
import math
import numpy as np
from .conditions import BoundaryConditions
from .partition import CompoundQuadChart
from .quad_charts import FixedRimQuality
from .refinement import QuadRefiner
from .surface_fields import BoundaryJetField
from .surface_query import QuadHeightField
from .slope_projection import SlopeProjection
from .mesh import Mesh, cross
from .quality_correction import QuadShapeCorrection
from .drainage import BoundaryDrainage
from .publication_slopes import PublicationSlopeGuard


@dataclass
class CompoundSurface:
    mesh: object
    holes: dict
    chart: object
    field: object
    boundary_refs: dict
    metadata: dict


class CompoundSurfaceBuilder:
    def __init__(self, patterns, *, maximum_slope=.10, slope_margin=1.1, max_refinements=3, max_cells=50000,
                 project_heights=True, maximum_join_slope=.03):
        if (not all(math.isfinite(v) and v > 0. for v in (maximum_slope, slope_margin))
                or slope_margin < 1. or type(max_refinements) is not int or max_refinements < 0
                or type(max_cells) is not int or max_cells < 4 or type(project_heights) is not bool
                or not math.isfinite(maximum_join_slope) or maximum_join_slope < 0.):
            raise ValueError('Invalid compound surface limits')
        self.patterns, self.maximum_slope, self.slope_margin = patterns, maximum_slope, slope_margin
        self.max_refinements, self.max_cells = max_refinements, max_cells
        self.project_heights, self.maximum_join_slope = project_heights, maximum_join_slope

    @staticmethod
    def projection_checks(mesh, fixed, conditions):
        edge = np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2)
        a, b = np.r_[edge[:, 0], edge[:, 1]], np.r_[edge[:, 1], edge[:, 0]]
        neighbor = np.full(len(mesh.xyz), np.inf)
        np.minimum.at(neighbor, a, mesh.xyz[b, 2])
        minima = int(np.sum((~fixed) & (mesh.xyz[:, 2] < neighbor - 1e-6)))
        owners = {}
        for fi, face in enumerate(mesh.quads):
            for side in range(4):
                owners.setdefault(tuple(sorted((int(face[side]), int(face[(side + 1) % 4])))), []).append((fi, side))
        boundary_nodes = np.flatnonzero(fixed)
        lookup = {tuple(mesh.xyz[node, :2]): int(node) for node in boundary_nodes}
        maximum = 0.
        for boundary in conditions:
            for port in boundary.ports:
                nodes = []
                for vertex in port.vertices:
                    xy = boundary.xy[vertex]
                    node = lookup.get(tuple(xy))
                    if node is None:
                        at = int(np.argmin(np.linalg.norm(mesh.xyz[boundary_nodes, :2] - xy, axis=1)))
                        node = int(boundary_nodes[at])
                        if np.linalg.norm(mesh.xyz[node, :2] - xy) > 1e-8:
                            raise ValueError('Height projection lost an incoming profile vertex')
                    nodes.append(node)
                for first, last in zip(nodes, nodes[1:]):
                    hits = owners.get(tuple(sorted((first, last))), ())
                    if len(hits) != 1:
                        raise ValueError('Height projection lost an incoming profile edge')
                    fi, side = hits[0]
                    p = mesh.xyz[mesh.quads[fi]]
                    u, v = ((.5, 0.), (1., .5), (.5, 1.), (0., .5))[side]
                    a = (1. - v) * (p[1] - p[0]) + v * (p[2] - p[3])
                    b = (1. - u) * (p[3] - p[0]) + u * (p[2] - p[1])
                    determinant = cross(a[:2], b[:2])
                    gradient = np.array(((a[2] * b[1] - b[2] * a[1]) / determinant,
                                         (a[0] * b[2] - b[0] * a[2]) / determinant))
                    tangent = mesh.xyz[last, :2] - mesh.xyz[first, :2]
                    tangent /= np.linalg.norm(tangent)
                    maximum = max(maximum, abs(float(gradient @ np.array((tangent[1], -tangent[0])))))
        return dict(isolated_interior_minima=minima, maximum_entrance_longitudinal_slope=maximum)

    def correct_refined_grading(self, mesh, field, refs, conditions, limit, max_cells):
        """Repair a refined chart, then reconcile shape, slopes and drainage.

        Only internal XY coordinates and heights may change. The original
        analytical field is sampled again after planar correction. Optional
        height projection uses monotone, bounded spill floors to avoid replacing
        a steep cell with a small basin. All original boundary XYZ stays fixed.
        """
        xyz = mesh.xyz.copy()
        xyz[:, 2] = 0.
        candidate, planar = QuadShapeCorrection(max_cells=max_cells).improve(
            Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams)))
        if not planar['qualified']:
            raise ValueError('Refined compound surface needs better cell grading')
        candidate.xyz[:, 2] = field(candidate.xyz[:, :2])
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[list(refs)] = True
        candidate.xyz[fixed] = mesh.xyz[fixed]
        projector = SlopeProjection(max_cells=max_cells)
        drainage = BoundaryDrainage(max_cells=max_cells)
        floors = np.full(len(mesh.xyz), field.datum)
        rounds = []
        for _ in range(8):
            if np.any(candidate.quad_slopes() > limit + 1e-8):
                if not self.project_heights:
                    raise ValueError('Regraded compound surface requires height projection')
                candidate, projected = projector.solve(candidate, fixed,
                    datum=field.datum, ceiling=field.ceiling, maximum_slope=limit,
                    minimum_heights=floors)
            else:
                projected = None
            filled, spill = drainage.fill(candidate, fixed)
            rounds.append(dict(projection=projected, drainage=spill))
            if spill['maximum_fill'] <= 1e-6:
                break
            if not self.project_heights:
                raise ValueError('Regraded compound surface requires depression correction')
            raised = filled > candidate.xyz[:, 2] + 1e-10
            floors[raised] = np.maximum(floors[raised], filled[raised])
            candidate, projected = projector.solve(candidate, fixed,
                datum=field.datum, ceiling=field.ceiling, maximum_slope=limit,
                minimum_heights=floors)
            rounds[-1]['drainage_projection'] = projected
        else:
            raise ValueError('Regraded compound drainage exceeded its iteration budget')
        checks = self.projection_checks(candidate, fixed, conditions)
        if checks['isolated_interior_minima']:
            raise ValueError('Regraded compound surface introduces an interior depression')
        if checks['maximum_entrance_longitudinal_slope'] > self.maximum_join_slope:
            raise ValueError('Regraded compound surface exceeds the entrance longitudinal slope limit')
        if not np.array_equal(candidate.xyz[fixed], mesh.xyz[fixed]):
            raise ValueError('Regraded compound surface changed its boundary')
        details = dict(planar=planar, rounds=rounds, joins=checks, exact_boundary=True)
        query = QuadHeightField(candidate, datum=field.datum, ceiling=field.ceiling,
            metadata=dict(field.metadata, refined_grading=details), max_cells=max_cells)
        return candidate, query, details

    def build(self, group, context):
        chart = CompoundQuadChart(self.patterns, spacing=context.network.settings.spacing,
                                  max_cells=self.max_cells).build(group, context)
        group = chart.partition.group
        layouts = {p.region.owner: p.layout for p in group.footprints}
        conditions = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
        field = BoundaryJetField(conditions)
        incoming_slope = max(float(np.linalg.norm(boundary.gradients[list(port.vertices)], axis=1).max())
                             for boundary in conditions for port in boundary.ports)
        mesh = Mesh(chart.mesh.xyz.copy(), chart.mesh.quads, chart.mesh.material, dict(chart.mesh.seams))
        refs = dict(chart.boundary_refs)
        mesh.xyz[:, 2] = field(mesh.xyz[:, :2])
        for node, ref in refs.items():
            mesh.xyz[node, 2] = conditions[ref.loop].edge_samples(ref.edge, ref.fraction)[1]
        # Existing profiles remain authoritative, including a narrow street
        # whose crossfall already exceeds the ordinary surface slope target.
        limit = max(self.maximum_slope, self.slope_margin * incoming_slope)
        grass_cells = sum(len(hole.mesh.quads) for hole in chart.holes.values())
        projection, projection_failure = None, None
        if self.project_heights and np.any(mesh.quad_slopes() > limit + 1e-8):
            fixed = np.zeros(len(mesh.xyz), bool)
            fixed[list(refs)] = True
            try:
                adjusted, details = SlopeProjection(max_cells=self.max_cells - grass_cells).solve(mesh, fixed,
                    datum=field.datum, ceiling=field.ceiling, maximum_slope=limit)
                checks = self.projection_checks(adjusted, fixed, conditions)
                if checks['isolated_interior_minima']:
                    raise ValueError('Projected heights introduce an interior depression')
                if checks['maximum_entrance_longitudinal_slope'] > self.maximum_join_slope:
                    raise ValueError('Projected heights exceed the entrance longitudinal slope limit')
                projection = dict(details, **checks)
                query = QuadHeightField(adjusted, datum=field.datum, ceiling=field.ceiling,
                    metadata=dict(field.metadata, height_projection=projection), max_cells=self.max_cells)
                mesh, field = adjusted, query
            except ValueError as exc:
                projection, projection_failure = None, str(exc)
        refinements = []
        for step in range(self.max_refinements + 1):
            slopes = mesh.quad_slopes()
            selected = np.flatnonzero(slopes > limit + 1e-8)
            if not len(selected):
                break
            if step == self.max_refinements:
                detail = f'; height projection: {projection_failure}' if projection_failure else ''
                raise ValueError('Compound surface needs a different chart to meet its physical slope limit' + detail)
            mesh, remap = QuadRefiner(self.patterns, max_cells=self.max_cells - grass_cells).refine(mesh, selected, field)
            refs = {int(remap[node]): ref for node, ref in refs.items()}
            mesh, improvement = FixedRimQuality().improve(mesh)
            interior = np.ones(len(mesh.xyz), bool)
            interior[list(refs)] = False
            mesh.xyz[interior, 2] = field(mesh.xyz[interior, :2])
            refinements.append(dict(selected=len(selected), quads=len(mesh.quads), improvement=improvement))
        quality = mesh.validate()
        refined_grading = None
        if quality['maximum_edge_ratio'] > 8.:
            mesh, field, refined_grading = self.correct_refined_grading(mesh, field, refs,
                conditions, limit, self.max_cells - grass_cells)
            quality = mesh.validate()
            if quality['maximum_edge_ratio'] > 8.:
                raise ValueError('Refined compound surface needs better cell grading')
        if mesh.xyz[:, 2].min() < field.datum - 1e-8 or mesh.xyz[:, 2].max() > field.ceiling + 1e-8:
            raise ValueError('Compound surface exceeded its incoming elevation envelope')
        corner_seed_reconciliation = None
        if chart.metadata.get('corner_band_seed'):
            # A new planar seed must pass drainage as well as shape/slope
            # checks, even when its first sampled field is already below the
            # slope limit. Do this in physical metres before exaggeration.
            mesh, field, corner_seed_reconciliation = self.correct_refined_grading(
                mesh, field, refs, conditions, limit, self.max_cells - grass_cells)
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[list(refs)] = True
        def check_publication_joins(candidate):
            checks = self.projection_checks(candidate, fixed, conditions)
            if checks['isolated_interior_minima']:
                raise ValueError('Publication correction introduces an interior depression')
            if checks['maximum_entrance_longitudinal_slope'] > self.maximum_join_slope:
                raise ValueError('Publication correction exceeds the entrance longitudinal slope limit')
        mesh, field, publication = PublicationSlopeGuard().correct(mesh, field, fixed,
            exaggeration=context.network.settings.exaggeration, limit=limit,
            max_cells=self.max_cells - grass_cells, check_joins=check_publication_joins,
            allow_projection=self.project_heights)
        holes = {}
        for li, hole in chart.holes.items():
            output = Mesh(hole.mesh.xyz.copy(), hole.mesh.quads, hole.mesh.material, dict(hole.mesh.seams))
            output.xyz[:, 2] = conditions[li].ceiling
            holes[li] = output
        # Field construction and all checks use physical metres. Exaggeration
        # affects only the resulting inspection mesh, around the common datum.
        slopes = mesh.slope_statistics()
        exaggeration = context.network.settings.exaggeration
        for output in (mesh, *holes.values()):
            output.xyz[:, 2] = field.datum + exaggeration * (output.xyz[:, 2] - field.datum)
        return CompoundSurface(mesh, holes, chart, field, refs,
            dict(quality=quality, planar_chart=chart.metadata,
                 physical_slopes=slopes, measured_slope_limit=limit,
                 datum=field.datum, ceiling=field.ceiling, field=field.metadata,
                 maximum_incoming_slope=incoming_slope,
                 height_projection=projection, projection_failure=projection_failure,
                 refined_grading=refined_grading,
                 corner_seed_reconciliation=corner_seed_reconciliation,
                 publication_roundoff=publication,
                 refinements=refinements, total_quads=len(mesh.quads) + grass_cells,
                 construction=('bounded boundary jets with constrained quad heights' if projection or refined_grading else
                               'bounded boundary jets on conforming direct quad charts'),
                 default_provider_integration=False))
