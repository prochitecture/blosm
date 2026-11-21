"""Profiled direct-quad surfaces with bounded optional height projection."""
from dataclasses import dataclass,replace
import math
import numpy as np
from .conditions import BoundaryConditions
from .model import Key
from .partition import CompoundQuadChart
from .quad_charts import FixedRimQuality
from .refinement import QuadRefiner
from .surface_fields import BoundaryJetField
from .updates import drain_steps
from .surface_query import QuadHeightField
from .slope_projection import SlopeProjection, ProjectionIterationBudgetError
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
    physical_mesh: object = None
    physical_holes: object = None
    hole_fields: object = None
    conditions: object = None


class PublicationEntranceSlopeError(ValueError):
    """A stored-slope correction disturbed an actual incoming tangent plane."""
    def __init__(self,candidate,checks):
        self.candidate,self.checks=candidate,dict(checks)
        super().__init__('Publication correction exceeds the entrance longitudinal slope limit')


class CompoundSurfaceBuilder:
    def __init__(self, patterns, *, maximum_slope=.10, slope_margin=1.1, max_refinements=3, max_cells=50000,
                 project_heights=True, maximum_join_slope=.03, source_heights=False):
        if (not all(math.isfinite(v) and v > 0. for v in (maximum_slope, slope_margin))
                or slope_margin < 1. or type(max_refinements) is not int or max_refinements < 0
                or type(max_cells) is not int or max_cells < 4 or type(project_heights) is not bool
                or not math.isfinite(maximum_join_slope) or maximum_join_slope < 0.
                or type(source_heights) is not bool):
            raise ValueError('Invalid compound surface limits')
        self.patterns, self.maximum_slope, self.slope_margin = patterns, maximum_slope, slope_margin
        self.max_refinements, self.max_cells = max_refinements, max_cells
        self.project_heights, self.maximum_join_slope = project_heights, maximum_join_slope
        self.source_heights = source_heights

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

    def correct_publication(self,mesh,field,fixed,conditions,*,exaggeration,limit,max_cells,
                            additional_slabs=None):
        return drain_steps(self.correct_publication_steps(mesh,field,fixed,conditions,
            exaggeration=exaggeration,limit=limit,max_cells=max_cells,additional_slabs=additional_slabs))

    def correct_publication_steps(self,mesh,field,fixed,conditions,*,exaggeration,limit,max_cells,
                                  additional_slabs=None):
        """Keep a failed publication projection inside real entrance slabs."""
        def check_joins(candidate):
            checks=self.projection_checks(candidate,fixed,conditions)
            if checks['isolated_interior_minima']:
                raise ValueError('Publication correction introduces an interior depression')
            if checks['maximum_entrance_longitudinal_slope']>self.maximum_join_slope:
                raise PublicationEntranceSlopeError(candidate,checks)
            if additional_slabs is not None and additional_slabs.violation(candidate.xyz[:,2])>1e-8:
                raise ValueError('Publication correction lost original body height constraints')
        guard=PublicationSlopeGuard()
        directional={}
        if additional_slabs is not None:
            from .entrance_height_recovery import combine_slabs,entrance_slabs
            slabs=combine_slabs(entrance_slabs(mesh,fixed,conditions,
                .95*self.maximum_join_slope,at_corners=False),additional_slabs)
            directional=dict(height_slabs=slabs)
        try:
            result,query,publication=yield from guard.correct_steps(mesh,field,fixed,exaggeration=exaggeration,limit=limit,
                max_cells=max_cells,check_joins=check_joins,allow_projection=self.project_heights,**directional)
            if additional_slabs is not None:
                result=mesh.with_positions(result.xyz)
                query=QuadHeightField(result,datum=field.datum,ceiling=field.ceiling,
                    metadata=query.metadata,max_cells=max_cells)
            return result,query,publication
        except PublicationEntranceSlopeError as failure:
            if not self.project_heights:raise
            candidate=failure.candidate
            if (not np.array_equal(candidate.xyz[fixed],mesh.xyz[fixed])
                    or not np.array_equal(candidate.xyz[:,:2],mesh.xyz[:,:2])
                    or not np.array_equal(candidate.quads,mesh.quads)
                    or float(candidate.quad_slopes().max())>limit+guard.tolerance
                    or float(guard.physical_buffer(candidate,field.datum,exaggeration).quad_slopes().max())
                       >limit+guard.tolerance):
                raise ValueError('Publication entrance recovery requires unchanged feasible source geometry') from failure
            from .entrance_height_recovery import recover_steps,entrance_slabs,combine_slabs
            extra={} if additional_slabs is None else dict(additional_slabs=additional_slabs)
            candidate,directional=yield from recover_steps(self,candidate,field,fixed,conditions,limit,max_cells,
                at_corners=False,**extra)
            query=QuadHeightField(candidate,datum=field.datum,ceiling=field.ceiling,
                metadata=dict(field.metadata,publication_entrance_recovery=directional),max_cells=max_cells)
            slabs=combine_slabs(entrance_slabs(candidate,fixed,conditions,
                .95*self.maximum_join_slope,at_corners=False),additional_slabs)
            result,query,publication=yield from guard.correct_steps(candidate,query,fixed,exaggeration=exaggeration,
                limit=limit,max_cells=max_cells,check_joins=check_joins,
                allow_projection=self.project_heights,height_slabs=slabs)
            rounded=guard.physical_buffer(result,field.datum,exaggeration)
            stored_conditions=[replace(boundary,xy=boundary.xy.astype(np.float32).astype(float))
                for boundary in conditions]
            stored_joins=self.projection_checks(rounded,fixed,stored_conditions)
            if (stored_joins['isolated_interior_minima']
                    or stored_joins['maximum_entrance_longitudinal_slope']>self.maximum_join_slope
                    or not np.array_equal(result.xyz[fixed],mesh.xyz[fixed])
                    or not np.array_equal(result.xyz[:,:2],mesh.xyz[:,:2])
                    or not np.array_equal(result.quads,mesh.quads)):
                raise ValueError('Publication entrance recovery failed its stored or fixed-source contract')
            publication=dict(publication or {},directional_entrance_recovery=directional,
                failed_publication_joins=failure.checks,stored_joins=stored_joins,
                initial_published_maximum=float(guard.physical_buffer(mesh,field.datum,exaggeration).quad_slopes().max()),
                published_maximum=float(rounded.quad_slopes().max()),
                maximum_height_change=float(np.abs(result.xyz[:,2]-mesh.xyz[:,2]).max()),
                exact_boundary=True,exact_xy=True,exact_quads=True,
                maximum_reconciliation_attempts=1,exact_fixed_xyz=True,exact_xy_and_quads=True)
            if additional_slabs is not None:
                result=mesh.with_positions(result.xyz)
            query=QuadHeightField(result,datum=field.datum,ceiling=field.ceiling,
                metadata=dict(query.metadata,publication_roundoff=publication),max_cells=max_cells)
            return result,query,publication

    @staticmethod
    def _coupled_projection(projector,mesh,fixed,*,datum,ceiling,limit,floors):
        return drain_steps(CompoundSurfaceBuilder._coupled_projection_steps(projector,mesh,fixed,
            datum=datum,ceiling=ceiling,limit=limit,floors=floors))

    @staticmethod
    def _coupled_projection_steps(projector,mesh,fixed,*,datum,ceiling,limit,floors):
        try:
            return (yield from projector.solve_steps(mesh,fixed,datum=datum,ceiling=ceiling,
                maximum_slope=limit,minimum_heights=floors))
        except ProjectionIterationBudgetError as failure:
            candidate=failure.candidate
            if (not np.array_equal(candidate.xyz[fixed],mesh.xyz[fixed])
                    or not np.array_equal(candidate.xyz[:,:2],mesh.xyz[:,:2])
                    or not np.array_equal(candidate.quads,mesh.quads)
                    or np.any(candidate.xyz[:,2]<floors-1e-10)
                    or np.any(candidate.xyz[:,2]>ceiling+1e-10)):
                raise
            candidate.validate()
            # This is an intermediate feasibility result. The enclosing
            # coupled solver still requires spill floors, drainage, actual
            # incoming joins and the final publication checks before return.
            if float(candidate.quad_slopes().max())<=limit+1e-8:
                return candidate,dict(failure.report,
                    convergence='physical_feasibility_after_iteration_budget',
                    nearest_displacement_converged=False)
            # A nearest-height solve can exhaust its budget just outside the
            # physical constraint. Use its existing strict feasibility mode,
            # with the same solver budgets and fresh immutable-gradient floor.
            floor=PublicationSlopeGuard.fixed_gradient_floor(candidate,fixed,datum)
            target=floor+.75*(limit-floor)
            if not np.isfinite(floor) or not floor<target<limit:
                raise
            result,details=yield from projector.solve_steps(candidate,fixed,datum=datum,ceiling=ceiling,
                maximum_slope=target,feasibility_slope=limit,minimum_heights=floors)
            if (not np.array_equal(result.xyz[fixed],mesh.xyz[fixed])
                    or not np.array_equal(result.xyz[:,:2],mesh.xyz[:,:2])
                    or not np.array_equal(result.quads,mesh.quads)
                    or np.any(result.xyz[:,2]<floors-1e-10)
                    or np.any(result.xyz[:,2]>ceiling+1e-10)
                    or float(result.quad_slopes().max())>limit):
                raise ValueError('Coupled feasible projection failed its unchanged source or physical constraints')
            result.validate()
            return result,dict(details,coupled_feasible_retry=dict(
                failed_nearest_projection=failure.report,fixed_gradient_floor=floor,
                inner_target=target,physical_acceptance_limit=limit,
                maximum_iterations=projector.max_iterations,memory_bytes=projector.memory_bytes,
                exact_fixed_xyz=True,exact_xy_and_quads=True,spill_floors_preserved=True))

    def correct_refined_grading(self, mesh, field, refs, conditions, limit, max_cells, *, entrance_recovery=False):
        return drain_steps(self.correct_refined_grading_steps(mesh, field, refs, conditions, limit, max_cells,
            entrance_recovery=entrance_recovery))

    def correct_refined_grading_steps(self, mesh, field, refs, conditions, limit, max_cells, *, entrance_recovery=False):
        """Repair a refined chart, then reconcile shape, slopes and drainage.

        Only internal XY coordinates and heights may change. The original
        analytical field is sampled again after planar correction. Optional
        height projection uses monotone, bounded spill floors to avoid replacing
        a steep cell with a small basin. All original boundary XYZ stays fixed.
        """
        if type(entrance_recovery) is not bool:
            raise ValueError('Entrance height recovery requires an explicit Boolean policy')
        xyz = mesh.xyz.copy()
        xyz[:, 2] = 0.
        candidate, planar = yield from QuadShapeCorrection(max_cells=max_cells).improve_steps(
            Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams)))
        if not planar['qualified']:
            from .refined_grading_recovery import recover
            try:
                candidate, planar = recover(Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams)),
                    refs, max_cells, planar)
            except ValueError as failure:
                raise ValueError('Refined compound surface needs better cell grading') from failure
        if isinstance(field, BoundaryJetField):
            candidate.xyz[:, 2] = yield from field.evaluate_steps(candidate.xyz[:, :2])
        else:
            candidate.xyz[:, 2] = field(candidate.xyz[:, :2])
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[list(refs)] = True
        candidate.xyz[fixed] = mesh.xyz[fixed]
        projector = SlopeProjection(max_cells=max_cells)
        drainage = BoundaryDrainage(max_cells=max_cells)
        floors = np.full(len(mesh.xyz), field.datum)
        rounds = []
        spill_tolerance = 1e-6
        ordinary_projections = 8
        maximum_projections = 2 * ordinary_projections
        for iteration in range(maximum_projections + 1):
            if np.any(candidate.quad_slopes() > limit + 1e-8):
                if not self.project_heights:
                    raise ValueError('Regraded compound surface requires height projection')
                candidate, projected = yield from self._coupled_projection_steps(projector,candidate,fixed,
                    datum=field.datum,ceiling=field.ceiling,limit=limit,floors=floors)
            else:
                projected = None
            filled, spill = drainage.fill(candidate, fixed)
            rounds.append(dict(projection=projected, drainage=spill))
            if spill['maximum_fill'] <= spill_tolerance:
                break
            if not self.project_heights:
                raise ValueError('Regraded compound surface requires depression correction')
            if iteration >= ordinary_projections:
                # An eighth projection can leave a few micrometres to drain.
                # Continue only while the three latest spill reductions would
                # reach the original tolerance within the remaining budget.
                recent = [row['drainage']['maximum_fill'] for row in rounds[-4:]]
                contraction = max(b/a for a, b in zip(recent, recent[1:]))
                remaining = maximum_projections - iteration
                if remaining == 0 or contraction >= 1. or recent[-1]*contraction**remaining > spill_tolerance:
                    raise ValueError('Regraded compound drainage exceeded its iteration budget')
            yield
            raised = filled > candidate.xyz[:, 2] + 1e-10
            floors[raised] = np.maximum(floors[raised], filled[raised])
            candidate, projected = yield from self._coupled_projection_steps(projector,candidate,fixed,
                datum=field.datum,ceiling=field.ceiling,limit=limit,floors=floors)
            rounds[-1]['drainage_projection'] = projected
        checks = self.projection_checks(candidate, fixed, conditions)
        if checks['isolated_interior_minima']:
            raise ValueError('Regraded compound surface introduces an interior depression')
        directional=None
        if checks['maximum_entrance_longitudinal_slope'] > self.maximum_join_slope:
            if not entrance_recovery or not self.project_heights:
                raise ValueError('Regraded compound surface exceeds the entrance longitudinal slope limit')
            from .entrance_height_recovery import recover_steps
            candidate,directional=yield from recover_steps(self,candidate,field,fixed,conditions,limit,max_cells)
            checks=directional['joins']
        if not np.array_equal(candidate.xyz[fixed], mesh.xyz[fixed]):
            raise ValueError('Regraded compound surface changed its boundary')
        details = dict(planar=planar, rounds=rounds, joins=checks, exact_boundary=True)
        if directional is not None:details['directional_entrance_recovery']=directional
        query = QuadHeightField(candidate, datum=field.datum, ceiling=field.ceiling,
            metadata=dict(field.metadata, refined_grading=details), max_cells=max_cells)
        return candidate, query, details

    @staticmethod
    def retained_crossfall(group, context):
        """Largest physical transverse slope of a Section owned by this junction."""
        maximum = 0.
        for claim in group.claims:
            section = context[Key('path', claim.section_id)].section
            profile = context[Key('profile', section.id)]
            # The parabolic full and half crowns have affine derivatives in q.
            # Their extreme crossfalls are at the edges, with dq/dx = 2/width.
            edge_slopes = profile.derivative(np.array((-1., 1.))) * (2. / section.width)
            maximum = max(maximum, float(np.max(np.abs(edge_slopes))))
        return maximum

    def boundary_conditions(self, group, context):
        """The actual source policy and physical limit, shared with admission."""
        layouts = {p.region.owner: p.layout for p in group.footprints}
        conditions = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
        incoming_slope = max(float(np.linalg.norm(boundary.gradients[list(port.vertices)], axis=1).max())
                             for boundary in conditions for port in boundary.ports)
        limit = max(self.maximum_slope, self.slope_margin * incoming_slope)
        if self.source_heights:
            from .source_curb_conditions import apply as source_conditions
            limit = max(limit, self.slope_margin * self.retained_crossfall(group, context))
            conditions = [source_conditions(context, group, loop, data, maximum_slope=limit)
                          for loop, data in zip(group.loops, conditions)]
        return conditions, incoming_slope, limit

    def build(self, group, context):
        return drain_steps(self.build_steps(group, context))

    def build_steps(self, group, context):
        one_sided = getattr(group.physical_layout, 'one_sided_source_preparation', False)
        if one_sided:
            from .one_sided_paired_junctions import certify
            from .one_sided_charts import OneSidedQuadChart
            certify(group)
            builder = OneSidedQuadChart(self.patterns, spacing=context.network.settings.spacing,
                                        max_cells=self.max_cells)
            chart = yield from builder.build_steps(group, context)
        else:
            chart_options = {}
            if self.source_heights:
                from .source_height_chart_admission import SourceHeightChartAdmission
                chart_options['source_height_admission'] = SourceHeightChartAdmission(self, context)
            builder = CompoundQuadChart(self.patterns, spacing=context.network.settings.spacing,
                                        max_cells=self.max_cells, **chart_options)
            chart = yield from builder.build_steps(group, context)
        yield
        group = chart.partition.group
        if getattr(group.physical_layout,'source_cycle_full_width_proof',None) is not None:
            from .cycle_dimensional_retry import certify_layout
            group.physical_layout.source_cycle_full_width_proof=certify_layout(context,group,
                {claim.section_id for claim in group.claims})
        layouts = {p.region.owner: p.layout for p in group.footprints}
        # Real incoming profile planes determine the physical limit before an
        # opt-in junction bank transition can select its finite source domain.
        conditions, incoming_slope, limit = self.boundary_conditions(group, context)
        yield
        field = BoundaryJetField(conditions)
        mesh = Mesh(chart.mesh.xyz.copy(), chart.mesh.quads, chart.mesh.material, dict(chart.mesh.seams))
        refs = dict(chart.boundary_refs)
        mesh.xyz[:, 2] = yield from field.evaluate_steps(mesh.xyz[:, :2])
        for node, ref in refs.items():
            mesh.xyz[node, 2] = conditions[ref.loop].edge_samples(ref.edge, ref.fraction)[1]
        # Existing profiles remain authoritative, including a narrow street
        # whose crossfall already exceeds the ordinary surface slope target.
        grass_cells = sum(len(hole.mesh.quads) for hole in chart.holes.values())
        projection, projection_failure = None, None
        initial_source_reconciliation = None
        retry_coupled_grading, initial_reconciliation_attempted = False, False
        projection_base = None
        if self.project_heights and np.any(mesh.quad_slopes() > limit + 1e-8):
            projection_base = (Mesh(mesh.xyz.copy(), mesh.quads, mesh.material, dict(mesh.seams)), field, dict(refs))
            fixed = np.zeros(len(mesh.xyz), bool)
            fixed[list(refs)] = True
            try:
                adjusted, details = yield from SlopeProjection(max_cells=self.max_cells - grass_cells).solve_steps(mesh, fixed,
                    datum=field.datum, ceiling=field.ceiling, maximum_slope=limit)
                checks = self.projection_checks(adjusted, fixed, conditions)
                if checks['isolated_interior_minima']:
                    retry_coupled_grading = True
                    raise ValueError('Projected heights introduce an interior depression')
                if checks['maximum_entrance_longitudinal_slope'] > self.maximum_join_slope:
                    retry_coupled_grading = True
                    raise ValueError('Projected heights exceed the entrance longitudinal slope limit')
                projection = dict(details, **checks)
                query = QuadHeightField(adjusted, datum=field.datum, ceiling=field.ceiling,
                    metadata=dict(field.metadata, height_projection=projection), max_cells=self.max_cells)
                mesh, field = adjusted, query
            except ValueError as exc:
                projection, projection_failure = None, str(exc)
                if isinstance(exc, ProjectionIterationBudgetError):
                    retry_coupled_grading = True
                if (chart.metadata.get('budget_reflex_fallback') or chart.metadata.get('piecewise_curb_sampling')
                        or chart.metadata.get('radial_core_reflex_fallback') or chart.metadata.get('fixed_rim_compound_recovery')
                        or (group.physical_layout is not None and
                        (hasattr(group.physical_layout, 'source_cycle_proofs') or
                         getattr(group.physical_layout, 'source_region_union', False))) or one_sided):
                    # A single slope projection can create a basin. Source
                    # pools already require coupled spill-floor/drainage
                    # reconciliation; try it before refining the unchanged
                    # analytical field and exhausting that unrelated budget.
                    initial_reconciliation_attempted = True
                    try:
                        adjusted, query, details = yield from self.correct_refined_grading_steps(
                            mesh, field, refs, conditions, limit, self.max_cells - grass_cells)
                        mesh, field = adjusted, query
                        initial_source_reconciliation = details
                    except ValueError as coupled_failure:
                        projection_failure += '; initial source reconciliation: ' + str(coupled_failure)
        refinements = []
        deferred_projection_reconciliation = None
        for step in range(self.max_refinements + 1):
            yield
            slopes = mesh.quad_slopes()
            selected = np.flatnonzero(slopes > limit + 1e-8)
            if not len(selected):
                break
            if step == self.max_refinements:
                if (self.project_heights and retry_coupled_grading
                        and not initial_reconciliation_attempted):
                    # Preserve every previously successful refinement first.
                    # An exhausted nearest-height solve can still admit a
                    # feasible surface. Reconcile from the original chart
                    # and certify the spill floors, slopes and entrance joins.
                    base, base_field, base_refs = projection_base
                    try:
                        mesh, field, details = yield from self.correct_refined_grading_steps(base, base_field, base_refs,
                            conditions, limit, self.max_cells-grass_cells,entrance_recovery=True)
                    except ValueError as failure:
                        projection_failure += '; deferred coupled reconciliation: '+str(failure)
                    else:
                        deferred_projection_reconciliation = dict(details,
                            initial_projection_failure=projection_failure,
                            discarded_refinement_attempts=refinements)
                        refs, refinements = base_refs, []
                        break
                detail = f'; height projection: {projection_failure}' if projection_failure else ''
                raise ValueError('Compound surface needs a different chart to meet its physical slope limit' + detail)
            refiner = QuadRefiner(self.patterns, max_cells=self.max_cells - grass_cells)
            mesh, remap = yield from refiner.refine_steps(mesh, selected, field)
            refs = {int(remap[node]): ref for node, ref in refs.items()}
            mesh, improvement = yield from FixedRimQuality().improve_steps(mesh)
            interior = np.ones(len(mesh.xyz), bool)
            interior[list(refs)] = False
            if isinstance(field, BoundaryJetField):
                mesh.xyz[interior, 2] = yield from field.evaluate_steps(mesh.xyz[interior, :2])
            else:
                mesh.xyz[interior, 2] = field(mesh.xyz[interior, :2])
            refinements.append(dict(selected=len(selected), quads=len(mesh.quads), improvement=improvement))
        yield
        quality = mesh.validate()
        refined_grading = None
        if quality['maximum_edge_ratio'] > 8.:
            mesh, field, refined_grading = yield from self.correct_refined_grading_steps(mesh, field, refs,
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
            mesh, field, corner_seed_reconciliation = yield from self.correct_refined_grading_steps(
                mesh, field, refs, conditions, limit, self.max_cells - grass_cells)
        finite_port_reconciliation = None
        if any(getattr(layout, 'finite_scope_ports', ()) or getattr(layout, 'finite_free_pads', ())
               for layout in layouts.values()):
            # A finite source pool introduces a new physical crown or owns a
            # whole terminal pad. Its drainage and incoming joins are required
            # even when the first boundary field is already below the slope
            # limit. Keep all original rim XYZ fixed during reconciliation.
            mesh, field, finite_port_reconciliation = yield from self.correct_refined_grading_steps(
                mesh, field, refs, conditions, limit, self.max_cells - grass_cells)
        source_cycle_reconciliation = None
        # Advanced ordinary mouths own complete source prefixes too. Require
        # the same height/drainage/join reconciliation even below the slope cap.
        ordinary_clearance = (getattr(group.physical_layout,'ordinary_complete_mouth_clearance',None)
            or any(getattr(source.layout,'ordinary_complete_mouth_clearance',None) for source in group.footprints))
        if ((group.physical_layout is not None and (hasattr(group.physical_layout,'source_cycle_proofs')
                or getattr(group.physical_layout,'source_region_union',False)))
                or ordinary_clearance
                or chart.metadata.get('source_width_visibility')
                or chart.metadata.get('budget_reflex_fallback')
                or chart.metadata.get('radial_core_reflex_fallback')
                or chart.metadata.get('fixed_rim_compound_recovery')
                or chart.metadata.get('piecewise_curb_sampling')
                or chart.metadata.get('positive_assembled_repair')
                or chart.metadata.get('banded_grass_recovery') or one_sided):
            # A closed source pool owns complete original roads and may retain
            # plateau holes. Check drainage and every incoming tangent plane
            # even when its initial field already meets the slope bound.
            mesh, field, source_cycle_reconciliation = yield from self.correct_refined_grading_steps(
                mesh, field, refs, conditions, limit, self.max_cells - grass_cells)
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[list(refs)] = True
        body_guides, body_reconciliation = None, None
        if self.source_heights and getattr(group.physical_layout,'source_cycle_full_width_proof',None) is not None:
            from .source_body_heights import prepare_steps as prepare_source_body_steps
            from .source_body_references import BodyGuideLimits
            from .entrance_height_recovery import recover_steps
            body_guides = yield from prepare_source_body_steps(context,group,mesh,refs,
                limits=BodyGuideLimits(maximum_mesh_cells=min(50000,self.max_cells-grass_cells)))
            from .junction_bank_transitions import certify_body_domains
            certify_body_domains(body_guides,conditions)
            if body_guides.height_slabs is not None:
                if not self.project_heights:
                    raise ValueError('Source-body height constraints require height projection')
                try:
                    mesh, body_reconciliation = yield from recover_steps(self,body_guides.seed,field,fixed,conditions,
                        limit,self.max_cells-grass_cells,at_corners=False,
                        additional_slabs=body_guides.height_slabs, body_guides=body_guides)
                except ProjectionIterationBudgetError as failure:
                    if not len(body_guides.core_nodes):
                        raise
                    from .native_body_stencil_refinement import recover_steps as refine_body_steps
                    mesh, refs, fixed, body_guides, body_reconciliation = yield from refine_body_steps(
                        self, context, group, mesh, refs, conditions, field, body_guides, limit,
                        self.max_cells-grass_cells, previous_failure=failure)
            body_field = body_guides.physical_field(mesh,datum=field.datum,ceiling=field.ceiling)
            body_field.metadata = dict(field.metadata,**body_field.metadata)
            field = body_field
        publication_options = (dict(additional_slabs=body_guides.height_slabs)
            if body_guides is not None and body_guides.height_slabs is not None else {})
        yield
        mesh, field, publication = yield from self.correct_publication_steps(mesh,field,fixed,conditions,
            exaggeration=context.network.settings.exaggeration,limit=limit,
            max_cells=self.max_cells-grass_cells,**publication_options)
        if body_guides is not None:
            body_guides.assert_unchanged(mesh)
            body_field = body_guides.physical_field(mesh,datum=field.datum,ceiling=field.ceiling)
            body_field.metadata = dict(field.metadata,**body_field.metadata)
            field = body_field
        one_sided_quality = None
        if one_sided:
            from .four_side_charts import quality as pavement_quality
            from .mesh import publication_coordinates
            physical_quality = pavement_quality(mesh)
            stored_quality = pavement_quality(mesh.with_positions(
                publication_coordinates(mesh.xyz, 1., physical=True)))
            if any(q['maximum_edge_ratio'] > 8. or q['minimum_corner_sine'] < .03
                   for q in (physical_quality, stored_quality)):
                raise ValueError('One-sided physical surface fails unchanged pavement quality')
            one_sided_quality = dict(physical=physical_quality, stored=stored_quality)
        from .hole_surfaces import build as build_hole_surface
        holes, hole_fields, hole_height_reports = {}, {}, {}
        for li, hole in chart.holes.items():
            available = self.max_cells-len(mesh.quads)-(grass_cells-len(hole.mesh.quads))
            output = build_hole_surface(hole.mesh, mesh, li, datum=field.datum,
                ceiling=field.ceiling, exaggeration=context.network.settings.exaggeration,
                max_cells=available)
            holes[li], hole_fields[li], hole_height_reports[li] = output.mesh, output.field, output.metadata
            grass_cells += len(output.mesh.quads) - len(hole.mesh.quads)
        # Field construction and all checks use physical metres. Exaggeration
        # affects only the resulting inspection mesh, around the common datum.
        slopes = mesh.slope_statistics()
        exaggeration = context.network.settings.exaggeration
        # The exact physical construction remains available to Bundle/query
        # consumers even when an inspection exaggeration of zero erases Z.
        # At E1 these views share the normal output buffers.
        physical_mesh = mesh if exaggeration == 1. else mesh.with_positions(mesh.xyz.copy())
        physical_holes = {li:(output if exaggeration == 1. else output.with_positions(output.xyz.copy()))
                          for li,output in holes.items()}
        for output in (mesh, *holes.values()):
            output.xyz[:, 2] = field.datum + exaggeration * (output.xyz[:, 2] - field.datum)
        return CompoundSurface(mesh, holes, chart, field, refs,
            dict(quality=quality, planar_chart=chart.metadata,
                 physical_slopes=slopes, measured_slope_limit=limit,
                 datum=field.datum, ceiling=field.ceiling, field=field.metadata,
                 maximum_incoming_slope=incoming_slope,
                 **(dict(maximum_retained_crossfall=self.retained_crossfall(group, context))
                    if self.source_heights else {}),
                 height_projection=projection, projection_failure=projection_failure,
                 initial_source_reconciliation=initial_source_reconciliation,
                 **(dict(deferred_projection_reconciliation=deferred_projection_reconciliation)
                    if deferred_projection_reconciliation is not None else {}),
                 source_union_fallback=getattr(group.physical_layout, 'source_union_fallback', None),
                 refined_grading=refined_grading,
                 corner_seed_reconciliation=corner_seed_reconciliation,
                 finite_port_reconciliation=finite_port_reconciliation,
                 source_cycle_reconciliation=source_cycle_reconciliation,
                 **(dict(source_body_heights=dict(preparation=body_guides.ledger,
                     reconciliation=body_reconciliation)) if body_guides is not None else {}),
                 **(dict(source_height_conditions=True) if self.source_heights else {}),
                 **(dict(source_bank_transitions=[dict(loop=li,**row)
                      for li,boundary in enumerate(conditions)
                      for row in getattr(boundary,'source_transitions',())])
                    if any(getattr(boundary,'source_transitions',()) for boundary in conditions) else {}),
                 publication_roundoff=publication,
                 **(dict(one_sided_source_quality=one_sided_quality) if one_sided else {}),
                 refinements=refinements, total_quads=len(mesh.quads) + grass_cells,
                 construction=('bounded boundary jets with constrained quad heights' if projection or refined_grading else
                               'bounded boundary jets on conforming direct quad charts'),
                 default_provider_integration=False, source_hole_heights=hole_height_reports),
            physical_mesh, physical_holes, hole_fields,tuple(conditions))
