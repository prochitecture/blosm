"""Finite merged terminal prefixes with the complete original road remainder.

Evaluation planes cut original Sections before their next physical junction.
They create no End, junction, source alias or connection. Actual side junctions
and their inner curb chains remain in the established separated construction.
"""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np

from .bundle_axis import AxialChartUnavailable, SeparatedPlan
from .bundle_angular_curves import AngularCrown
from .bundle_angular_open_islands import gap_intervals
from .bundle_angular_ports import AngularPort
from .bundle_angular_routes import AngularMemberRoute
from .bundle_cap_sampling import SourceRimSampling
from .bundle_compounds import CompoundBundleBuilder
from .bundle_merged_bridge import MergedProfileBridge
from .bundle_merged_entrances import PairedAngularEntrance
from .bundle_merged_plan import physical_context
from .bundle_merged_runs import ExternalRunDomain, MergedRunWithExternalCap
from .bundles import BundleProvider, Cap, CenterlineRange, certify_disjoint
from .drainage import BoundaryDrainage
from .four_side_charts import quality
from .intersections import half_intervals
from .merged_run_inventory import inventory_frame, sample_sources, _SourceBoundary
from .mesh import Assembler, Mesh, publication_coordinates
from .model import BundleMember, End, Key
from .ownership import PLAN, CompoundProduct, surface_key
from .prepared_merged_profiles import PreparedMergedMouth
from .providers import CarriagewayProvider, Part, RoadSpan
from .routes import RouteTopology
from .separator_ports import PORTS
from .source_mesh_containment import certify_mesh
from .surface_query import QuadHeightField
from .updates import BuildContext
from .validation import validate_overlaps


def terminal_mouths(context, bundle):
    """Resolve paired terminals solely through original incidence registries."""
    if len(bundle.members) != 2 or not all(member.complete for member in bundle.members):
        return ()
    result = []
    registry = context[PORTS]
    plan = context[PLAN]
    for end in (0, 1):
        visits = tuple(member.path[0 if end == 0 else -1] for member in bundle.members)
        ends = tuple(End(v.section_id, int(not v.forward) if end == 0 else int(v.forward)) for v in visits)
        visible = tuple(registry.resolve(endpoint) for endpoint in ends)
        if not all(visible): continue
        from .separator_ports import bind_physical_port
        bindings = tuple(bind_physical_port(context, port, endpoint) for port, endpoint in zip(visible, ends))
        owners = tuple(binding.owner for binding in bindings)
        if owners[0] != owners[1] or bindings[0].approach_index != bindings[1].approach_index:
            continue
        product = context.products.get(surface_key(owners[0]))
        if not isinstance(product, CompoundProduct) or product.group.physical_layout is None: continue
        index = bindings[0].approach_index
        approach = product.group.physical_layout.approaches[index]
        if not isinstance(approach, PairedAngularEntrance): continue
        if (approach.end_index != end or tuple(approach.original_ends) != ends
                or not hasattr(approach, 'run')
                or tuple(tuple(m.visits) for m in approach.run.members)
                    != tuple(tuple(m.path) for m in bundle.members)):
            raise ValueError('Merged terminal lost its complete original ordered Bundle run')
        result.append((owners[0], index, end))
    return tuple(result)


def _terminal_support(context, bundle, role, end):
    """Stop at the first original physical junction, retaining real C1 seams."""
    visits = tuple(bundle.members[role].path)
    by_end = {e: j for j in (*context.network.intersections, *context.network.connections) for e in j.ends}
    ordered = visits if end == 0 else visits[::-1]
    count = 1
    for previous, following in zip(ordered, ordered[1:]):
        exit = End(previous.section_id, int(previous.forward) if end == 0 else int(not previous.forward))
        entry = End(following.section_id, int(not following.forward) if end == 0 else int(following.forward))
        site = by_end.get(exit)
        if site is None or entry not in site.ends:
            raise ValueError('Terminal support lacks its original source incidence')
        if len(site.ends) != 2: break
        a, b = context[Key('path', previous.section_id)], context[Key('path', following.section_id)]
        offsets = (a.section.offset*(1 if previous.forward else -1),
                   b.section.offset*(1 if following.forward else -1))
        if abs(a.section.width-b.section.width) > 1e-8 or abs(offsets[0]-offsets[1]) > 1e-8:
            break
        candidate = ordered[:count+1] if end == 0 else tuple(reversed(ordered[:count+1]))
        try:
            route = AngularMemberRoute(context, BundleMember(tuple(candidate)), role)
        except AxialChartUnavailable:
            break
        if route.transitions: break  # Its real seam surface remains independently owned.
        count += 1
    return tuple(ordered[:count]) if end == 0 else tuple(reversed(ordered[:count]))


def _source_parameter(path, frame, target):
    """Invert certified progress on the original finite centerline."""
    first, last = 0., path.length
    for _ in range(64):
        middle = .5 * (first + last)
        value = float(frame.local(path.frame(middle)[0])[0])
        if value < target: first = middle
        else: last = middle
    return .5 * (first + last)


class _EvaluationDomain(ExternalRunDomain):
    def __init__(self, context, bundle, view, entrance, cuts, merged_cut):
        self.context, self.bundle = context, bundle
        self.frame, self.field, self.rays = view.frame_chart, view.field, view.field.curves
        self.run, self.end = entrance.run, 1-entrance.end_index
        self.ports = {}
        for role, (visit, station) in enumerate(cuts):
            path = context[Key('path', visit.section_id)]
            profile = context[Key('profile', visit.section_id)].oriented(int(not visit.forward))
            if profile.half != (1 if role == 0 else -1):
                raise ValueError('Terminal evaluation plane changed an original half-profile')
            self.ports[role, self.end] = AngularPort(role, self.end, path, visit.forward,
                station, 2*half_intervals(context.network.settings), profile)
        values = np.concatenate([self.frame.local(port.coordinates(np.array((0., 1.)))[0])[:, 0]
                                 for port in self.ports.values()])
        local = self.frame.local(merged_cut[:, :2])
        if np.ptp(local[:, 0]) > 1e-8:
            raise ValueError('Terminal prefix needs its transported complete common-fibre cut')
        merged = float(local[:, 0].mean())
        self.head, self.tail = ((float(values.max()), merged) if self.end == 0 else
                                (merged, float(values.min())))


class _EvaluationPrefix(MergedRunWithExternalCap):
    def __init__(self, entrance, cuts, merged_cut):
        self.entrance, self.evaluation_cuts, self.merged_cut = entrance, cuts, merged_cut

    def make_domain(self, context, bundle, view):
        self.evaluation_domain = _EvaluationDomain(context, bundle, view, self.entrance,
                                                   self.evaluation_cuts, self.merged_cut)
        return self.evaluation_domain

    def entrance_margin(self, context):
        return self.margin_factor*context.network.settings.spacing

    def build(self, context, bundle, view, merged_cut):
        # Keep the established margin first. Only computational inward planes
        # vary; the complete physical source profiles are fixed for all tries.
        from .bundle_angular_attached import AttachedAngularPlan
        failures = []
        for factor in AttachedAngularPlan.margin_factors:
            self.margin_factor = factor
            try:
                pavement, grass, report = super().build(context, bundle, view, merged_cut)
            except ValueError as error:
                failures.append(dict(margin_factor=factor, reason=str(error)))
                continue
            report['evaluation_margin_selection'] = dict(margin_factor=factor, prior_failures=failures)
            return pavement, grass, report
        raise ValueError('No qualified finite source evaluation bridge: '+str(failures))

    def cap(self, context, members, end, sampling):
        domain = self.evaluation_domain
        points = [domain.ports[role, end].coordinates(np.array(1.))[0] for role in (0, 1)]
        count = max(4, 2*int(np.ceil(np.linalg.norm(points[0]-points[1]) /
                                   (2.*context.network.settings.spacing))))
        return Cap.free(np.r_[points[0], domain.field.ceiling],
                        np.r_[points[1], domain.field.ceiling], count)


def _boundary_edges(mesh):
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    if len(counts) and counts.max() > 2: raise ValueError('Terminal union has a nonmanifold edge')
    return unique[counts == 1]


def _edge_coordinates(mesh):
    return {tuple(sorted((tuple(mesh.xyz[a]), tuple(mesh.xyz[b])))) for a, b in _boundary_edges(mesh)}


def _declared_nodes(mesh, points):
    """Restore only an explicitly supplied interface, never infer its owner."""
    nodes = []
    for point in points:
        matches = np.flatnonzero(np.linalg.norm(mesh.xyz-point, axis=1) <= 1e-8)
        if len(matches) != 1:
            raise ValueError('A declared original evaluation profile is absent or ambiguous')
        nodes.append(int(matches[0]))
    if len(set(nodes)) != len(nodes): raise ValueError('An evaluation profile repeats a physical vertex')
    mesh.xyz[nodes] = points
    return np.asarray(nodes, np.int32)


class TerminalMergedBundlePlan:
    mode = 'union'
    transition_owners = ()
    requires_separator_query = True

    @classmethod
    def build(cls, context, bundle):
        terminals = terminal_mouths(context, bundle)
        if not terminals: return None
        if len(terminals) != 1:
            raise AxialChartUnavailable('Two merged Bundle terminals need a multiple-prefix partition')
        return cls(context, bundle, terminals[0])

    def __init__(self, context, bundle, terminal):
        self.bundle, self.settings = bundle, context.network.settings
        self.source_ids = tuple(v.section_id for member in bundle.members for v in member.path)
        attachments = RouteTopology(context.network).attachments(context.network)
        plan = context[PLAN]
        owners = {plan.owner(site.key) or site.key for sid in self.source_ids
                  for site, end in attachments[sid] if len(site.ends) > 2}
        physical = physical_context(context, required_owners=owners)
        owner, index, end = terminal
        product = physical[surface_key(owner)]
        mouth = PreparedMergedMouth.from_product(product, index)
        entrance = mouth.approach
        frame = inventory_frame(physical, bundle, entrance)
        support_visits = tuple(_terminal_support(physical, bundle, role, end) for role in (0, 1))
        support = tuple((SimpleNamespace(visits=visits),) for visits in support_visits)
        samples = sample_sources(physical, support, entrance)
        curves, routes = [], []
        for role, visits in enumerate(support_visits):
            route = AngularMemberRoute(physical, BundleMember(visits), role)
            route.certify(frame)
            routes.append(route)
            boundary = _SourceBoundary(physical, BundleMember(visits), frame, samples,
                                       join_tolerance=route.join_tolerance)
            if boundary.owners:
                raise AxialChartUnavailable('Terminal prefix dimensional seams need explicit transition body ownership')
            curves.extend(boundary.points(side) for side in (-1., 1.))
        field = AngularCrown(frame, curves, datum=entrance.field.datum, rise=entrance.field.rise)
        domain = SimpleNamespace(frame=frame, field=field, rays=field.curves)
        lower = max(ray.stations[0] for ray in field.curves)
        upper = min(ray.stations[-1] for ray in field.curves)
        regions, roots = gap_intervals(domain, lower, upper)
        if not regions or regions[0 if end == 0 else -1][2]:
            raise ValueError('The accepted merged terminal needs a connected original overlap prefix')
        separated = [region for region in regions if region[2]]
        if not separated:
            raise AxialChartUnavailable('Merged terminal has no separated evaluation cut before its next original junction')
        margin = max(2.*self.settings.spacing,
                     2.*max(piece.section.width for route in routes for piece in route.pieces))
        event = (max(b for a, b, positive in regions if not positive) if end == 0 else
                 min(a for a, b, positive in regions if not positive))
        target = float(event + (margin if end == 0 else -margin))
        if not lower < target < upper:
            raise AxialChartUnavailable('Terminal overlap leaves no bounded source evaluation reserve')
        cuts, fully_consumed = [], set()
        for role, route in enumerate(routes):
            at = _source_parameter(route, frame, target)
            position = min(int(np.searchsorted(route.stations, at, side='right')-1), len(route.pieces)-1)
            visit = support_visits[role][position]
            piece = route.pieces[position]
            offset = at-route.stations[position]
            station = float(offset if visit.forward else piece.path.length-offset)
            path = physical[Key('path', visit.section_id)]
            metadata = CarriagewayProvider.span_metadata(physical, path.section, attachments[visit.section_id])
            if not metadata['start'] < station < metadata['end']:
                raise AxialChartUnavailable('Terminal evaluation plane would enter another original physical junction')
            profile = physical[Key('profile', visit.section_id)].oriented(int(not visit.forward))
            port = AngularPort(role, 1-end, path, visit.forward, station, 2*half_intervals(self.settings), profile)
            progress = frame.local(port.coordinates(np.array((0., 1.)))[0])[:, 0]
            if ((end == 0 and progress.min() <= event+2.*self.settings.spacing)
                    or (end == 1 and progress.max() >= event-2.*self.settings.spacing)):
                raise ValueError('Complete evaluation profile needs finite separated source reserve')
            cuts.append((visit, float(station)))
            consumed = support_visits[role][:position] if end == 0 else support_visits[role][position+1:]
            fully_consumed.update(v.section_id for v in consumed)
        visits = tuple(v for v, station in cuts)
        target_inventory = SimpleNamespace(frame=frame, field=field, entrance=entrance,
            source_sections=entrance.run.source_sections,
            source_visits=tuple(member.visits for member in entrance.run.members))
        bridge = MergedProfileBridge()
        connector, connector_report = bridge.build(physical, entrance, mouth.station, mouth.xyz,
                                                    target_inventory=target_inventory)
        view = SimpleNamespace(frame_chart=frame, field=field,
            paths=tuple(piece.path for route in routes for piece in route.pieces), sign=1. if end == 0 else -1.)
        prefix = _EvaluationPrefix(entrance, tuple(cuts), bridge.body_cut)
        # The new prefix's physical width proof uses a two-millimetre budget.
        # Refine only its mesh approximation; original Paths, prepared source
        # chords, actual mouth and field remain authoritative and unchanged.
        precise = BuildContext(replace(physical.network, settings=replace(physical.network.settings,
            curve_tolerance=min(.001, self.settings.curve_tolerance))), physical.products,
            physical.patterns, physical.task_digests, physical.source_digests)
        prefix_pavement, prefix_grass, prefix_report = prefix.build(precise, bundle, view, bridge.body_cut)
        prefix_report['actual_external_rounded_cap'] = False
        prefix_report['original_source_evaluation_planes'] = True

        products = dict(physical.products)
        products[Key('bundle_plan', bundle.id)] = SeparatedPlan()
        ordinary = BuildContext(physical.network, products, physical.patterns,
                                physical.task_digests, physical.source_digests)
        cut_by_sid = {visit.section_id: (visit, station) for visit, station in cuts}
        roads, coverage, accepted_metadata = {}, [], {}
        for sid in self.source_ids:
            section = ordinary.network.section_index[sid]
            metadata = CarriagewayProvider.span_metadata(ordinary, section, attachments[sid])
            accepted_metadata[sid] = metadata
            interval = None
            if sid in cut_by_sid:
                visit, station = cut_by_sid[sid]
                interval = ((station, metadata['end']) if (end == 0) == visit.forward else
                            (metadata['start'], station))
            if sid in fully_consumed:
                road = RoadSpan(dict(metadata, end=metadata['start'],
                    terminal_union_prefix=bundle.key.text(), computational_source_interval=True))
            else:
                road = CarriagewayProvider.build(ordinary, section, attachments[sid], source_interval=interval)
            ordinary.products[Key('road', sid)] = road
            roads[sid] = road
            prefix_interval = None
            if sid in fully_consumed:
                prefix_interval = (metadata['start'], metadata['end'])
            elif sid in cut_by_sid:
                root_end = 0 if (end == 0) == cut_by_sid[sid][0].forward else 1
                station = cut_by_sid[sid][1]
                prefix_interval = ((metadata['start'], station) if root_end == 0 else
                                   (station, metadata['end']))
            pieces = [(road.metadata['start'], road.metadata['end'])]
            if prefix_interval is not None: pieces.append(prefix_interval)
            source_owned = metadata['compound_spans']+metadata['junction_spans']+metadata['transition_spans']
            pieces.extend((span['start'], span['end']) for span in source_owned)
            covered = 0.
            for a, b in sorted(pieces):
                if b <= a+1e-9: continue
                if a > covered+1e-7: raise ValueError('Terminal Bundle drops an original source interval')
                covered = max(covered, b)
            if abs(covered-physical[Key('path', sid)].length) > 1e-7:
                raise ValueError('Terminal Bundle does not cover its complete original source extent')
            coverage.append(dict(section=sid, visible_first=metadata['start'], visible_last=metadata['end'],
                remainder_first=road.metadata['start'], remainder_last=road.metadata['end'],
                terminal_union_interval=prefix_interval, actual_owned_prefixes=source_owned,
                complete_original_extent=True))
        runs = [CompoundBundleBuilder.runs(ordinary, member, role) for role, member in enumerate(bundle.members)]
        regions, continuations = CompoundBundleBuilder.regions(ordinary, bundle, runs)
        if len(regions) != 1 or continuations:
            raise AxialChartUnavailable('Terminal prefix remainder needs one complete actual separated curb pair')
        first, second = regions[0]
        if any(run.terminals[end].endpoint.section_id != visits[role].section_id
               for role, run in enumerate((first, second))):
            raise ValueError('Terminal remainder lost its original evaluation source Section')
        evaluation_cap = prefix.cap(ordinary, (), 1-end, None)
        caps = [None, None]
        caps[end] = evaluation_cap
        other = 1-end
        caps[other] = BundleProvider.cap(ordinary, first.member, second.member, other,
            (first.coordinates[other], second.coordinates[other]),
            rim_source=SourceRimSampling(attachments))
        source_ranges = tuple(tuple(CenterlineRange.build(physical[Key('path', sid)],
            accepted_metadata[sid]['start'], accepted_metadata[sid]['end'])
            for sid in (visit.section_id for visit in member.path)
            if accepted_metadata[sid]['end'] > accepted_metadata[sid]['start']+1e-8)
            for member in bundle.members)
        separation = certify_disjoint(*source_ranges)
        remainder = BundleProvider.region(ordinary, bundle, (first.member, second.member), caps,
                                          separation, correct_positive=True)

        namespace = (bundle.key, 'terminal prefix', end)
        paving = Assembler()
        connector_keys = [(namespace, 'connector', i) for i in range(len(connector.xyz))]
        prefix_keys = [(namespace, key) for key in prefix.pavement_keys]
        for row, node in enumerate(bridge.back_ids): connector_keys[int(node)] = (namespace, ('merged cut', row))
        profile_keys = {}
        for role in (0, 1):
            port = prefix.domain.ports[role, other]
            xy, _, z = port.coordinates(np.linspace(0., 1., port.intervals+1))
            xyz = np.c_[xy, z]
            nodes = _declared_nodes(prefix_pavement, xyz)
            names = [(namespace, 'evaluation profile', role, i) for i in range(len(xyz))]
            for node, name in zip(nodes, names): prefix_keys[int(node)] = name
            profile_keys[role] = (xyz, names)
        paving.add(connector.xyz, connector.quads, connector_keys)
        paving.add(prefix_pavement.xyz, prefix_pavement.quads, prefix_keys)
        for sid, road in roads.items():
            if getattr(road, 'mesh', None) is None: continue
            names = [(namespace, 'original road', sid, i) for i in range(len(road.mesh.xyz))]
            if sid in cut_by_sid:
                role = next(i for i, v in enumerate(visits) if v.section_id == sid)
                visit = visits[role]
                original_end = entrance.original_ends[role].end
                nodes = road.mesh.seams[original_end]
                if (1 if role == 0 else -1)*(1 if visit.forward else -1) < 0.: nodes = nodes[::-1]
                xyz, declared = profile_keys[role]
                if not np.array_equal(road.mesh.xyz[nodes], xyz):
                    raise ValueError('Terminal evaluation profile differs from its original complete road profile')
                for node, name in zip(nodes, declared): names[int(node)] = name
            paving.add(road.mesh.xyz, road.mesh.quads, names)
        pavement = paving.mesh()
        all_edges = np.sort(np.stack((pavement.quads, np.roll(pavement.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(all_edges, axis=0, return_counts=True)
        incidence = {tuple(edge): int(count) for edge, count in zip(unique, counts)}
        for xyz, names in profile_keys.values():
            nodes = [paving.lookup[name] for name in names]
            if any(incidence.get(tuple(sorted((a, b)))) != 2 for a, b in zip(nodes, nodes[1:])):
                raise ValueError('Terminal prefix and source remainder lack whole shared profile edges')
        nodes = [paving.lookup[(namespace, ('merged cut', row))] for row in range(len(bridge.back_ids))]
        if any(incidence.get(tuple(sorted((a, b)))) != 2 for a, b in zip(nodes, nodes[1:])):
            raise ValueError('Merged terminal connector and prefix lack whole shared profile edges')
        grass = Assembler()
        declared_cap = evaluation_cap.rim
        cap_names = [(namespace, 'evaluation grass', i) for i in range(len(declared_cap))]
        for label, mesh in (('prefix grass', prefix_grass), ('remainder grass', remainder.mesh)):
            names = [(namespace, label, i) for i in range(len(mesh.xyz))]
            nodes = _declared_nodes(mesh, declared_cap)
            for node, name in zip(nodes, cap_names): names[int(node)] = name
            grass.add(mesh.xyz, mesh.quads, names)
        separator = grass.mesh('grass')
        edges = np.sort(np.stack((separator.quads, np.roll(separator.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        incidence = {tuple(edge): int(count) for edge, count in zip(unique, counts)}
        nodes = [grass.lookup[name] for name in cap_names]
        if any(incidence.get(tuple(sorted((a, b)))) != 2 for a, b in zip(nodes, nodes[1:])):
            raise ValueError('Terminal prefix and separator remainder lack whole shared computational edges')
        limit = 50000*self.settings.resolution
        if len(pavement.quads)+len(separator.quads) > limit:
            raise ValueError('Complete terminal union and original remainder exceed their finite quad budget')
        checks = []
        for mesh in (pavement, separator):
            stored = mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))
            measured, published = quality(mesh), quality(stored)
            if (max(measured['maximum_edge_ratio'], published['maximum_edge_ratio']) > 8.
                    or mesh.material != 'grass' and min(measured['minimum_corner_sine'], published['minimum_corner_sine']) < .03):
                raise ValueError('Complete terminal union fails physical or stored quad grading')
            slopes = max(mesh.slope_statistics()['maximum'], stored.slope_statistics()['maximum'])
            along = max(frame.longitudinal_bound(mesh), frame.longitudinal_bound(stored))
            if slopes > .1+1e-6 or mesh.material != 'grass' and along > .03+1e-8:
                raise ValueError('Complete terminal union exceeds physical slope limits')
            fixed = np.zeros(len(mesh.xyz), bool); fixed[np.unique(_boundary_edges(mesh))] = True
            drainage = BoundaryDrainage(max_cells=limit).fill(stored, fixed)[1]
            if drainage['maximum_fill'] > 1e-6: raise ValueError('Terminal union contains an interior drainage depression')
            checks.append(dict(material=mesh.material, physical_quality=measured, stored_quality=published,
                               maximum_slope=slopes, longitudinal_bound=along, drainage=drainage))
        source_error = .1*self.settings.curve_tolerance
        maximum_error = min(.002, self.settings.curve_tolerance)
        if source_error >= maximum_error: raise ValueError('Terminal prefix source sampling needs its finer containment budget')
        containment = certify_mesh(prefix_pavement, curves, source_error=source_error,
            mesh_error=maximum_error-source_error, maximum_error=maximum_error+1e-7,
            max_cells=prefix.body.max_cells)
        declared = _edge_coordinates(pavement)
        for owned in owners: declared.update(_edge_coordinates(physical[surface_key(owned)].pavement.mesh))
        for cap in caps:
            if cap is not None:
                declared.update(tuple(sorted((tuple(a), tuple(b)))) for a, b in zip(cap.rim, cap.rim[1:]))
        if not _edge_coordinates(separator) <= declared:
            raise ValueError('Terminal union grass lacks a complete actual pavement or original outer-cap edge')
        payloads = {Key('terminal_union_pavement', bundle.id): pavement,
                    Key('terminal_union_grass', bundle.id): separator}
        query = Assembler(); query.add(pavement.xyz, pavement.quads); query.add(separator.xyz, separator.quads)
        ceilings = [field.ceiling]
        for owned in sorted(owners):
            actual = physical[surface_key(owned)]
            if abs(actual.field.datum-field.datum) > 1e-9:
                raise ValueError('Terminal Bundle and physical source owners need a common minimum datum')
            ceilings.append(actual.field.ceiling)
            for kind, part in (('pavement', actual.pavement), ('grass', actual.grass)):
                if part is None: continue
                payloads[Key('terminal_owned_'+kind, owned.id)] = part.mesh
                query.add(part.mesh.xyz, part.mesh.quads)
        overlap = validate_overlaps(payloads, cell_size=2.)
        self.datum = field.datum
        self.query = QuadHeightField(query.mesh(), datum=self.datum, ceiling=max(ceilings), max_cells=200000)
        self.field, self.frame = field, frame
        self.evaluation_profiles = profile_keys
        self.pavement_keys = [None]*len(pavement.xyz)
        for key, at in paving.lookup.items(): self.pavement_keys[at] = key
        self.grass_keys = [None]*len(separator.xyz)
        for key, at in grass.lookup.items(): self.grass_keys[at] = key
        metadata = dict(construction='finite merged original terminal prefix and complete source road remainder',
            source_sections=list(self.source_ids), terminal_owner=owner.text(), terminal_end=end,
            evaluation_planes=[dict(section=v.section_id, station=s, original_forward=v.forward) for v, s in cuts],
            computational_interfaces_only=True, physical_aliases_created=False,
            complete_selected_source_coverage=True, exact_whole_evaluation_profile_edges=True,
            original_source_coverage=coverage, centerline_proof_checks=separation,
            prefix_complete_sections=sorted(fully_consumed),
            merged_connector=connector_report, prefix=prefix_report, remainder=remainder.metadata,
            source_containment=containment, stored_checks=checks, overlap=overlap,
            query_only_junction_owners=[owned.text() for owned in sorted(owners)],
            junction_geometry_republished=False, complete_quads=len(pavement.quads)+len(separator.quads))
        for mesh in (pavement, separator):
            mesh.xyz[:, 2] = self.datum+self.settings.exaggeration*(mesh.xyz[:, 2]-self.datum)
        self._pavement = Part(bundle.key, pavement, dict(metadata, quality=pavement.validate()))
        self._separator = Part(bundle.key, separator, dict(metadata, quality=separator.validate()))

    def pavement(self, context): return self._pavement

    def separator(self, context): return self._separator

    def height(self, xy):
        xy = np.asarray(xy)
        result = self.query.evaluate(xy.reshape(-1, 2), default=self.datum)
        return (self.datum+self.settings.exaggeration*(result-self.datum)).reshape(xy.shape[:-1])

    def release_geometry(self):
        for part in (self._pavement, self._separator):
            if part.mesh is not None: part.metadata['bounds'] = part.mesh.bounds; part.mesh = None
