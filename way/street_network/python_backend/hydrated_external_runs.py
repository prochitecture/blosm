"""Finite hydrated source leaves from a merged cut to actual separated ports.

This consumer has no runtime discovery hook. Original visits and published
physical profiles are supplied explicitly; a common fibre already transported
to the independently certified run chart is the only merged-side interface.
"""
from types import SimpleNamespace
import numpy as np

from .bundle_angular_curves import AngularCrown
from .bundle_angular_internal import SampledAngularPort
from .bundle_merged_runs import ExternalRunDomain, MergedRunWithExternalCap
from .bundle_visible_ports import ProfileSamplingPatterns, VisibleCaps
from .bundle_visible_runs import BundleRun, MemberRun
from .bundles import BundleProvider, CenterlineRange, certify_disjoint
from .drainage import BoundaryDrainage
from .four_side_charts import quality
from .mesh import Mesh
from .model import End, Key
from .ownership import PLAN, surface_key
from .paired_source_body import PairedSourceBodyPiece
from .source_mesh_containment import certify_mesh
from .validation import validate_overlaps


class HydratedExternalDomain(ExternalRunDomain):
    def __init__(self, context, bundle, view, spans, mouth, merged_cut):
        self.context, self.bundle = context, bundle
        self.frame, self.field, self.rays = view.frame_chart, view.field, view.field.curves
        self.end = 1-mouth.approach.end_index
        boundaries = ('hydrated run head', 'hydrated run tail')
        self.run = BundleRun(tuple(MemberRun(span.visits, (span.head.endpoint, span.tail.endpoint), boundaries)
                                   for span in spans), boundaries)
        self.ports = {}
        self.physical_owners = set()
        for role, span in enumerate(spans):
            endpoint = span.head if self.end == 0 else span.tail
            port = endpoint.port
            visit = span.visits[0 if self.end == 0 else -1]
            original = End(visit.section_id, int(not visit.forward) if self.end == 0 else int(visit.forward))
            if (endpoint.kind != 'PHYSICAL' or endpoint.endpoint != original or endpoint.registry_key is not None
                    or endpoint.prepared_mouth is not None or not isinstance(port, SampledAngularPort)
                    or (port.role, port.end, port.forward) != (role, self.end, visit.forward)
                    or port.path is not context[Key('path', visit.section_id)]
                    or port.mouth.endpoint != original or port.mouth.owner != endpoint.owner
                    or port.station != endpoint.source_station):
                raise ValueError('Hydrated external run needs two actual separated original sampled profiles')
            visible = context[Key('separator_ports', 0)].resolve(original)
            plan = context.products.get(PLAN)
            owner = ((plan.owner(visible.site.key) if plan is not None else None) or visible.site.key) if visible else None
            if owner != endpoint.owner or context.products.get(surface_key(owner)) is None:
                raise ValueError('Hydrated external profile lost its accepted original physical owner')
            profile = context[Key('profile', visit.section_id)].oriented(int(not visit.forward))
            if (profile != port.profile or profile.half != (1 if role == 0 else -1)
                    or abs(profile.datum-self.field.datum)+abs(profile.rise-self.field.rise) > 1e-9):
                raise ValueError('Hydrated external run changed an original half crown')
            row = span.source_intervals[0 if self.end == 0 else -1]
            station = row['first'] if (self.end == 0) == row['forward'] else row['last']
            if row['section'] != original.section_id or abs(station-port.station) > 1e-8:
                raise ValueError('Hydrated external source interval does not reach its complete physical profile')
            self.ports[role, self.end] = port
            self.physical_owners.add(owner)
        values = np.concatenate([self.frame.local(port.mouth.physical_xyz[:, :2])[:, 0]
                                 for port in self.ports.values()])
        local = self.frame.local(merged_cut[:, :2])
        if np.ptp(local[:, 0]) > 1e-8:
            raise ValueError('Hydrated external run needs a complete transported common-fibre cut')
        merged_station = float(local[:, 0].mean())
        if self.end == 0: self.head, self.tail = float(values.max()), merged_station
        else: self.head, self.tail = merged_station, float(values.min())
        if self.tail-self.head <= 4.*context.network.settings.spacing:
            raise ValueError('Hydrated external profiles leave no finite source body interval')

    def patterns_for(self, port):
        return ProfileSamplingPatterns(self.context.patterns, port.mouth.parameters)


class HydratedExternalRunWithCap(MergedRunWithExternalCap):
    """Reuse the accepted external-cap algorithm on explicit hydrated leaves."""
    def make_domain(self, context, bundle, approach):
        return HydratedExternalDomain(context, bundle, approach, self._spans, self._mouth, self._merged_cut)

    @staticmethod
    def cap(context, members, end, sampling):
        actual = VisibleCaps.common_compound(context, *members, end)
        return actual if actual is not None else BundleProvider.cap(context, *members, end, rim_source=sampling)

    def build(self, context, bundle, spans, mouth, merged_cut, *, frame, namespace):
        try: hash(namespace)
        except TypeError: raise ValueError('Hydrated external assembly needs a stable hashable namespace')
        merged_end = mouth.approach.end_index
        root = PairedSourceBodyPiece._mouth(context, mouth, merged_end)
        curves, sections, progress = PairedSourceBodyPiece._sources(context, bundle, spans, frame)
        incident = {end: site for site in (*context.network.intersections, *context.network.connections) for end in site.ends}
        for span in spans:
            for previous, following in zip(span.visits, span.visits[1:]):
                site = incident[End(previous.section_id, int(previous.forward))]
                if len(site.ends) != 2:
                    raise ValueError('Hydrated external source leaves cannot chain a physical multi-arm junction')
        aliases = mouth.approach.source_mouths(mouth.station)
        for role, (span, alias) in enumerate(zip(spans, aliases)):
            end = span.head if merged_end == 0 else span.tail
            if (end.prepared_mouth is None or end.prepared_mouth.approach is not mouth.approach
                    or end.owner != mouth.owner or end.source_station != alias.station):
                raise ValueError('Hydrated external span lost its actual prepared merged mouth')
            if alias.original_incidence:
                if end.endpoint != alias.endpoint or end.registry_key is not None:
                    raise ValueError('Hydrated external span changed its original merged incidence')
            elif (end.endpoint is not None or end.registry_key is None
                    or (end.registry_key.owner, end.registry_key.section_id, end.registry_key.end)
                    != (mouth.owner, alias.endpoint.section_id, alias.endpoint.end)):
                raise ValueError('Hydrated external span lost its owner-scoped derived merged alias')
            row = span.source_intervals[0 if merged_end == 0 else -1]
            station = row['first'] if (merged_end == 0) == row['forward'] else row['last']
            if row['section'] != alias.endpoint.section_id or abs(station-alias.station) > 1e-8:
                raise ValueError('Hydrated external source interval does not reach its actual merged mouth')
        active = tuple(tuple(CenterlineRange.build(context[Key('path', row['section'])], row['first'], row['last'])
                             for row in span.source_intervals) for span in spans)
        separation = certify_disjoint(*active)
        profile = mouth.approach.field
        field = AngularCrown(frame, curves, datum=profile.datum, rise=profile.rise)
        merged_cut = np.asarray(merged_cut, float)
        if (merged_cut.ndim != 2 or merged_cut.shape[1] != 3 or len(merged_cut) < 3
                or not np.isfinite(merged_cut).all()
                or np.max(abs(field.values(merged_cut[:, :2])-merged_cut[:, 2])) > 1e-8):
            raise ValueError('Hydrated external cut changed the complete target-frame merged crown')
        self._spans, self._mouth, self._merged_cut = tuple(spans), mouth, merged_cut
        view = SimpleNamespace(field=field, frame_chart=frame,
            paths=tuple(context[Key('path', sid)] for sid in sections), sign=1. if merged_end == 0 else -1.)
        pavement, grass, legacy = super().build(context, bundle, view, merged_cut)
        limit = self.body.max_cells
        checks = []
        for mesh in (pavement, grass):
            stored = Mesh(mesh.xyz.astype(np.float32).astype(float), mesh.quads, mesh.material)
            double, published = quality(mesh), quality(stored)
            if (max(double['maximum_edge_ratio'], published['maximum_edge_ratio']) > 8.
                    or (mesh.material != 'grass'
                        and min(double['minimum_corner_sine'], published['minimum_corner_sine']) < .03)):
                raise ValueError('Hydrated external run fails complete stored quad quality')
            slopes, along = stored.slope_statistics()['maximum'], frame.longitudinal_bound(stored)
            if slopes > .1+1e-6 or (mesh.material != 'grass' and along > .03+1e-8):
                raise ValueError('Hydrated external run exceeds its physical slope bounds: '+
                                 str(dict(material=mesh.material, maximum=slopes, longitudinal=along)))
            fixed = np.zeros(len(mesh.xyz), bool)
            fixed[np.unique(self.body._boundary_vertices(mesh))] = True
            drainage = BoundaryDrainage(max_cells=limit).fill(stored, fixed)[1]
            if drainage['maximum_fill'] > 1e-6:
                raise ValueError('Hydrated external run contains an interior drainage depression')
            checks.append(dict(material=mesh.material, quality=published, double_quality=double,
                maximum_slope=slopes, longitudinal_bound=along,
                minimum_corner_sine_limit=.03 if mesh.material != 'grass' else None,
                carriageway_longitudinal_limit=.03 if mesh.material != 'grass' else None, drainage=drainage))
        containment = certify_mesh(pavement, curves, source_error=.1*context.network.settings.curve_tolerance,
            mesh_error=.6*context.network.settings.curve_tolerance,
            maximum_error=.7*context.network.settings.curve_tolerance+1e-7, max_cells=limit)
        stored_containment = certify_mesh(Mesh(pavement.xyz.astype(np.float32).astype(float), pavement.quads),
            curves, source_error=.1*context.network.settings.curve_tolerance,
            mesh_error=.6*context.network.settings.curve_tolerance,
            maximum_error=.7*context.network.settings.curve_tolerance+1e-7, max_cells=limit)
        old = self.pavement_keys
        self.pavement_keys = [(namespace, key) for key in old]
        self.grass_keys = [(namespace, 'external grass', i) for i in range(len(grass.xyz))]
        self.merged_cut_keys = [(namespace, ('merged cut', i)) for i in range(len(merged_cut))]
        self.cut_keys = {key: [(namespace, name) for name in names] for key, names in self.body.cut_keys.items()}
        lookup = {key: i for i, key in enumerate(self.pavement_keys)}
        self.merged_cut_ids = np.asarray([lookup[key] for key in self.merged_cut_keys], dtype=np.int32)
        if not np.array_equal(pavement.xyz[self.merged_cut_ids], merged_cut):
            raise ValueError('Hydrated external run changed its transported complete merged cut')
        self.external_profile_ids = {}
        for (role, end), port in self.domain.ports.items():
            nodes = []
            for point in port.mouth.physical_xyz:
                matches = np.flatnonzero(np.all(pavement.xyz == point, axis=1))
                if len(matches) != 1:
                    raise ValueError('Hydrated external run lost an exact complete published source profile')
                nodes.append(int(matches[0]))
            self.external_profile_ids[role, end] = np.asarray(nodes, dtype=np.int32)
        payloads = {Key('hydrated_external_pavement', 0): pavement, Key('hydrated_external_grass', 0): grass}
        for owner in self.domain.physical_owners | {mouth.owner}:
            product = context[surface_key(owner)]
            if product.pavement.mesh is None:
                raise ValueError('Hydrated external run needs accepted hydrated original junction products')
            payloads[Key('physical_junction_pavement', owner.id)] = product.pavement.mesh
            if product.grass is not None: payloads[Key('physical_junction_grass', owner.id)] = product.grass.mesh
        overlap = validate_overlaps(payloads, cell_size=2.)
        self.frame, self.field = frame, field
        return pavement, grass, dict(construction='hydrated original source run with actual external rounded cap',
            source_sections=list(sections), source_intervals=[list(span.source_intervals) for span in spans],
            original_run_progress_lower=progress, centerline_separation=separation,
            stored_checks=checks, source_containment=containment, stored_source_containment=stored_containment,
            overlap=overlap,
            actual_merged_cut=True, actual_separated_source_profiles=True, legacy=legacy,
            total_quads=len(pavement.quads)+len(grass.quads), runtime_discovery_enabled=False)
