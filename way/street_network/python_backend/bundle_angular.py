"""Curved common-fibre Bundle union fallback with explicit source ownership.

The angular provider covers two complete members with compatible explicit
Section seams, real free ends and shared crowns. Existing axial and attached
providers stay first. Eligibility does not imply successful construction:
all physical boundary, quad, slope, height and allocation guards still apply.
"""
from dataclasses import replace
import numpy as np
from .bundle_axis import AxialChartUnavailable, _simplify
from .bundle_angular_curves import AngularCoordinates, AngularCrown
from .bundle_angular_islands import AngularIslandBuilder
from .bundle_angular_routes import AngularMemberRoute
from .bundle_angular_surfaces import AngularSurfaceBuilder
from .drainage import BoundaryDrainage
from .mesh import Assembler, Mesh
from .model import Key
from .paths import Path
from .providers import Part
from .quad_subdivision import UniformQuadSubdivision
from .regions import Region, RegionUnion
from .surface_query import QuadHeightField
from .surface_quality import ProfiledQuadCorrection
from .validation import validate_overlaps


class AngularBundlePlan:
    """A source-derived union plan; no generated vertex acts as an input ID."""
    mode = 'union'
    transition_owners = ()
    max_cells = 100000

    def __init__(self, context, bundle):
        if len(bundle.members) != 2:
            raise AxialChartUnavailable('Angular Bundle construction needs two lateral members')
        self.bundle, self.settings = bundle, context.network.settings
        self.source_ids = tuple(t.section_id for m in bundle.members for t in m.path)
        self.paths, self.profiles = [], []
        for role, member in enumerate(bundle.members):
            if len(member.path)>1:
                route = AngularMemberRoute(context, member, role)
                self.paths.append(route)
                self.profiles.append(route.profile)
                continue
            visit = member.path[0]
            path = context[Key('path', visit.section_id)]
            if not visit.forward:
                source = replace(path.section, points=path.section.points[::-1], offset=-path.section.offset)
                path = Path(source, self.settings.curve_tolerance)
            profile = context[Key('profile', visit.section_id)].oriented(0 if visit.forward else 1)
            if profile.half != (1 if role == 0 else -1):
                raise AxialChartUnavailable('Angular members require inward-rising shared half crowns')
            self.paths.append(path)
            self.profiles.append(profile)
        self.datum, self.rise = self.profiles[0].datum, self.profiles[0].rise
        if any(abs(p.datum-self.datum)>1e-7 or abs(p.rise-self.rise)>1e-7 for p in self.profiles):
            raise AxialChartUnavailable('Angular members need a common crown datum and rise')
        self.frame = AngularCoordinates(self.paths)
        self.angular_proofs = tuple(path.certify(self.frame) if isinstance(path, AngularMemberRoute)
                                   else self.frame.certify(path) for path in self.paths)
        self.smooth = all(path.smooth if isinstance(path, AngularMemberRoute) else not path.exact
                          for path in self.paths)
        curves, regions, counts = [], [], []
        owners, certificates, exact_corners = set(), [], []
        for role, path in enumerate(self.paths):
            if isinstance(path, AngularMemberRoute):
                if path.transitions:
                    from .bundle_angular_boundaries import AngularMemberBoundary
                    boundary = AngularMemberBoundary(context, bundle.members[role], self.frame,
                                                     join_tolerance=path.join_tolerance)
                    left, right = (boundary.points(side) for side in (-1.,1.))
                    path.maximum_curb_discrepancy = max(path.maximum_curb_discrepancy,boundary.maximum_join_error)
                    owners.update(boundary.owners)
                    certificates.extend(boundary.transition_certificates)
                    exact_corners.extend(boundary.exact_corners)
                else:
                    left, right = path.curbs(self.settings)
                    if path.exact:
                        for piece in path.pieces:
                            original = piece.path
                            for side in (-1.,1.):
                                q = side if piece.forward else -side
                                exact_corners.extend(original.point(original.knots,
                                    original.section.offset+q*.5*original.section.width)[1:-1])
            else:
                samples = path.samples(0., path.length, self.settings.spacing,
                                       tolerance=.1*self.settings.curve_tolerance)
                left, right = (path.point(samples, path.section.offset+side*.5*path.section.width)
                               for side in (-1., 1.))
                if path.exact:
                    for side in (-1.,1.):
                        exact_corners.extend(path.point(path.knots,
                            path.section.offset+side*.5*path.section.width)[1:-1])
            curves.extend((left, right))
            regions.append(Region(path.section.key, np.vstack((right, left[::-1]))))
            counts.append((len(right),len(left)))
        self.transition_owners = tuple(sorted(owners))
        self.transition_certificates = tuple(certificates)
        self.exact_corners = np.asarray(exact_corners).reshape(-1,2) if owners else None
        self.field = AngularCrown(self.frame, curves, datum=self.datum, rise=self.rise)
        common_first = max(curve.stations[0] for curve in self.field.curves)
        common_last = min(curve.stations[-1] for curve in self.field.curves)
        if common_last <= common_first+self.settings.spacing:
            raise AxialChartUnavailable('Angular members have no positive common body')
        knots = np.unique(np.concatenate((np.array((common_first, common_last)),
            *(c.stations[(c.stations>common_first)&(c.stations<common_last)] for c in self.field.curves))))
        ao, ai, bi, bo = self.field.edges(knots)
        if np.min(ao-ai)<=1e-8 or np.min(bi-bo)<=1e-8 or np.min(ao+ai-bi-bo)<=1e-8:
            raise AxialChartUnavailable('Angular members disagree with declared physical lateral order')
        gaps = np.maximum(ai, bo)-np.minimum(ao, bi)
        self.minimum_gap, self.maximum_gap = float(gaps.min()), float(gaps.max())
        if self.minimum_gap >= -1e-8:
            raise AxialChartUnavailable('No overlapping pavement needs an angular union')
        self.loops = RegionUnion().build(regions)
        outer = [loop for loop in self.loops if not loop.hole]
        holes = [loop for loop in self.loops if loop.hole]
        if len(outer) != 1:
            raise AxialChartUnavailable('Angular union needs another multi-region decomposition')
        self.islands = tuple(holes)
        owners = {path.section.key:(role, counts[role]) for role, path in enumerate(self.paths)}
        labels = []
        for span in outer[0].spans:
            role, (first_count,last_count) = owners[span.owner]
            side = (1 if span.segment<first_count-1 else
                    -1 if first_count<=span.segment<first_count+last_count-1 else 0)
            if side == (-1 if role==0 else 1):
                labels.append(('outer', role))
            else:
                labels.append(('cap' if side==0 else 'inner', role))
        self.outer, self.outer_labels = outer[0].xy, labels
        incoming = max(2.*p.rise/path.section.width for p, path in zip(self.profiles, self.paths))
        if self.transition_owners:
            incoming = max(2.*context[Key('profile',sid)].rise/context.network.section_index[sid].width
                           for sid in self.source_ids)
        self.maximum_slope = max(.10, 1.1*incoming)
        builder_type = AngularSurfaceBuilder
        if len(self.islands) > 1:
            from .bundle_angular_island_chain import MultiIslandSurfaceBuilder
            builder_type = MultiIslandSurfaceBuilder
        self.builder = builder_type(context.patterns, spacing=self.settings.cell_size,
            curve_tolerance=self.settings.curve_tolerance, maximum_slope=self.maximum_slope,
            max_cells=self.max_cells, cap_scale=self.settings.cell_size,
            mandatory_points=self.exact_corners if self.exact_corners is not None and len(self.exact_corners) else None)
        self._pavement = self._separator = self.query = None
        # Build before publishing union ownership or exposing furniture queries.
        self._construct(context)

    def _construct(self, context):
        high, low = self.builder.split_outer(self.outer, self.outer_labels, self.field)
        tolerance = .6*self.settings.curve_tolerance if self.smooth else 1e-9
        high = self.builder.simplify_sample(*high, tolerance, 'outer high')
        low = self.builder.simplify_sample(*low, tolerance, 'outer low')
        grass = None
        grass_metadata = None
        if self.islands:
            upper, lower = self.builder.island_halves(self.islands[0], self.frame)
            if not self.smooth:
                upper, lower = (_simplify(p, 1e-9) for p in (upper, lower))
            exact_sources = all(path.exact for path in self.paths)
            join_error = max((path.maximum_curb_discrepancy for path in self.paths
                              if isinstance(path,AngularMemberRoute)),default=0.)
            # Exact offset segments are linear, hence sampling adds no error.
            # Reserve the actual route-join bound, union weld and simplifier;
            # Smooth curb sampling retains its independent 10% tolerance.
            source_error = 0. if exact_sources else .1*self.settings.curve_tolerance
            if self.transition_owners:
                # Path sampling plus simplification each reserve 10%; quintic
                # transition sampling reserves 10% even for Exact centerlines.
                source_error = max(source_error, (.1 if exact_sources else .2)*self.settings.curve_tolerance)
            boundary_error = source_error+join_error+1e-7+1e-9
            island_builder = AngularIslandBuilder(self.frame, context.patterns, self.bundle.key,
                spacing=self.settings.cell_size, curve_tolerance=self.settings.curve_tolerance,
                height=self.field.ceiling, max_cells=self.max_cells,boundary_error=boundary_error)
            if len(self.islands) > 1:
                from .bundle_angular_island_chain import AngularIslandChain
                grass, seam_high, seam_low, grass_metadata = AngularIslandChain(
                    self.builder, island_builder).build(self.islands, self.field,
                        high[0][0], high[0][-1], smooth=self.smooth)
            else:
                grass, rims, grass_metadata = island_builder.build(upper, lower, smooth=self.smooth)
                head = self.builder.ridge(self.frame, self.field, high[0][0], upper[0], 'ridge head')
                tail = self.builder.ridge(self.frame, self.field, upper[-1], high[0][-1], 'ridge tail')
                head[1][-1], tail[1][0] = rims[0][1][0], rims[0][1][-1]
                seam_high = self.builder.join_chains((head, rims[0], tail))
                seam_low = self.builder.join_chains((head, rims[1], tail))
        else:
            seam_high = seam_low = self.builder.ridge(self.frame, self.field,
                high[0][0], high[0][-1], 'ridge')
        assembly, reports = Assembler(), []
        for name, outer, seam in (('high', high, seam_high), ('low', low, seam_low)):
            mesh, keys, report = self.builder.direct_patch(self.frame, self.field, outer, seam, name)
            assembly.add(mesh.xyz, mesh.quads, keys)
            reports.append(report)
        mesh = assembly.mesh()
        coarse_cells = len(mesh.quads)+(len(grass.quads) if grass is not None else 0)
        if self.settings.resolution > 1:
            subdivision = UniformQuadSubdivision(context.patterns, max_cells=self.max_cells)
            mesh = subdivision.build(mesh, self.settings.resolution, self.field.values,
                datum=self.datum, ceiling=self.field.ceiling)
            if grass is not None:
                grass = subdivision.build(grass, self.settings.resolution,
                    lambda xy:np.full(len(xy), self.field.ceiling),
                    datum=self.datum, ceiling=self.field.ceiling)
        refinement_correction = None
        separator_refinement_correction = None
        if mesh.validate()['maximum_edge_ratio'] > 8.:
            mesh, refinement_correction = ProfiledQuadCorrection(
                max_cells=self.max_cells, max_iterations=500).improve(mesh, self.field.values,
                fixed_vertices=mesh.seams.get('crown_guides', ()), datum=self.datum,
                ceiling=self.field.ceiling, exaggeration=1., maximum_slope=self.maximum_slope)
        if grass is not None and grass.validate()['maximum_edge_ratio'] > 8.:
            grass, separator_refinement_correction = ProfiledQuadCorrection(
                max_cells=self.max_cells, max_iterations=500).improve(grass,
                lambda xy:np.full(len(xy), self.field.ceiling), datum=self.datum,
                ceiling=self.field.ceiling, exaggeration=1., maximum_slope=self.maximum_slope)
            if grass.validate()['maximum_edge_ratio'] > 8.:
                raise ValueError('Refined angular separator exceeds the world quad edge-ratio limit')
        total_cells = len(mesh.quads)+(len(grass.quads) if grass is not None else 0)
        if total_cells > self.max_cells:
            raise ValueError('Angular Bundle exceeds its complete surface cell budget')
        boundary, contract = self.builder.boundary_contract(mesh, high[0], low[0], grass)
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[boundary] = True
        original = mesh.xyz[fixed].copy()
        filled, drainage = BoundaryDrainage(max_cells=self.max_cells).fill(mesh, fixed)
        mesh.xyz[:, 2] = filled
        if not np.array_equal(mesh.xyz[fixed], original):
            raise ValueError('Angular drainage changed a physical boundary profile')
        if mesh.xyz[:, 2].min()<self.datum-1e-10 or mesh.xyz[:, 2].max()>self.field.ceiling+1e-10:
            raise ValueError('Angular pavement exceeded its shared crown envelope')
        physical_slopes = mesh.slope_statistics()
        longitudinal = self.frame.longitudinal_bound(mesh)
        if physical_slopes['maximum']>self.maximum_slope+1e-8 or longitudinal>.03+1e-8:
            raise ValueError('Angular pavement exceeded its physical slope limits')
        validate_overlaps({Key('bundle_pavement', self.bundle.id):mesh,
            **({Key('separator', self.bundle.id):grass} if grass is not None else {})},
            cell_size=max(1., 2.*self.settings.spacing))
        self.query = QuadHeightField(mesh, datum=self.datum, ceiling=self.field.ceiling,
            max_cells=self.max_cells)
        self.grass_query = None if grass is None else QuadHeightField(grass,
            datum=self.datum, ceiling=self.field.ceiling, max_cells=self.max_cells)
        metadata = dict(construction='angular common-fibre pavement union', islands=len(self.islands),
            source_sections=list(self.source_ids), minimum_gap=self.minimum_gap, maximum_gap=self.maximum_gap,
            angular_certificate_lower=list(self.angular_proofs), physical_boundary_contract=contract,
            physical_slopes=physical_slopes, measured_slope_limit=self.maximum_slope,
            angular_longitudinal_bound=longitudinal, crown_ceiling=self.field.ceiling,
            drainage=drainage, patches=reports, all_direct_quads=True, quality=mesh.validate(),
            resolution=self.settings.resolution, coarse_cells=coarse_cells,
            refinement_correction=refinement_correction,
            refinement='uniform cached quad subdivision after qualified cap selection')
        if any(isinstance(path, AngularMemberRoute) for path in self.paths):
            metadata['member_routes'] = [path.source_spans() if isinstance(path, AngularMemberRoute)
                else [dict(section=path.section.id, forward=self.bundle.members[role].path[0].forward,
                           first=0., last=path.length)] for role, path in enumerate(self.paths)]
            metadata['member_route_curb_discrepancy'] = [path.maximum_curb_discrepancy
                if isinstance(path,AngularMemberRoute) else 0. for path in self.paths]
            metadata['member_route_curb_tolerance'] = [path.join_tolerance
                if isinstance(path,AngularMemberRoute) else 0. for path in self.paths]
        if self.transition_owners:
            metadata['dimensional_transition_owners'] = [key.text() for key in self.transition_owners]
            metadata['angular_transition_certificates'] = list(self.transition_certificates)
            metadata['physical_boundary_sampling'] = list(self.builder.boundary_sampling)
        if metadata['quality']['maximum_edge_ratio'] > 8.:
            raise ValueError('Refined angular pavement exceeds the world quad edge-ratio limit')
        mesh.xyz[:, 2] = self.datum+self.settings.exaggeration*(mesh.xyz[:, 2]-self.datum)
        self._pavement = Part(self.bundle.key, mesh, metadata)
        if grass is not None:
            grass.xyz[:, 2] = self.datum+self.settings.exaggeration*(grass.xyz[:, 2]-self.datum)
            separator_metadata = dict(construction='angular finite-tip direct quad separator', islands=len(self.islands),
                                      quality=grass.validate(), chart=grass_metadata)
            if separator_refinement_correction is not None:
                separator_metadata['refinement_correction'] = separator_refinement_correction
            self._separator = Part(self.bundle.key, grass, separator_metadata)

    def pavement(self, context):
        return self._pavement

    def separator(self, context):
        return self._separator

    def height(self, xy):
        points = np.asarray(xy, dtype=float)
        shape = points.shape[:-1]
        points = points.reshape(-1, 2)
        physical = self.query.evaluate(points, default=self.datum)
        if self.grass_query is not None:
            physical = np.maximum(physical, self.grass_query.evaluate(points, default=self.datum))
        return (self.datum+self.settings.exaggeration*(physical-self.datum)).reshape(shape)

    def release_geometry(self):
        # Published Part buffers are released by the backend. Retain only
        # immutable bounded quad queries for source-relative attachments.
        self.outer = None
        self.loops = ()
