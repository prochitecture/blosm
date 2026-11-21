"""Explicit paired-source assembly; discovery remains the caller's policy.

Physical owners and approach indices pair leaves. Original member terminals
identify separated caps. Coordinate proximity never establishes a connection.
"""
from types import SimpleNamespace
from collections.abc import Mapping
import numpy as np

from .bundle_angular_curves import AngularCrown
from .bundle_merged_bridge import MergedProfileBridge
from .bundle_visible_ports import VisibleMouth
from .hydrated_external_runs import HydratedExternalRunWithCap
from .mesh import Assembler, Mesh
from .model import End, Key
from .ownership import CompoundProduct, surface_key
from .separator_ports import PORTS
from .paired_source_body import PairedSourceBodyPiece
from .paired_source_sampling import verify_shared_source_phase
from .prepared_merged_profiles import PreparedMergedMouth
from .source_mouths import source_mouths
from .surface_query import QuadHeightField
from .validation import validate_overlaps


def interval_inventory(spans_by_role):
    """Capture the complete canonical inventory before any selection occurs."""
    return tuple(tuple((row['section'], row['forward'], float(row['first']), float(row['last']))
                       for span in spans for row in span.source_intervals) for spans in spans_by_role)


def _ordered(context, bundle, inventory):
    if isinstance(inventory, Mapping):
        if set(inventory) != {0, 1}:
            raise ValueError('Hydrated source inventory needs both original lateral roles')
        inventory = tuple(inventory[role] for role in (0, 1))
    if len(inventory) != 2:
        raise ValueError('Hydrated run needs two ordered original interval inventories')
    result = []
    for role, rows in enumerate(inventory):
        member = bundle.members[role]
        offsets, distance = {}, 0.
        for visit in member.path:
            if visit.section_id in offsets:
                raise ValueError('Hydrated run repeats an original member Section')
            path = context[Key('path', visit.section_id)]
            offsets[visit.section_id] = (distance, path.length, visit.forward)
            distance += path.length
        previous = -np.inf
        canonical = []
        for sid, forward, first, last in rows:
            if sid not in offsets:
                raise ValueError('Hydrated run inventory invents an original member Section')
            base, length, expected_forward = offsets[sid]
            if forward != expected_forward or not 0. <= first < last <= length:
                raise ValueError('Hydrated run inventory changes original direction or finite source bounds')
            low, high = (first, last) if forward else (length-last, length-first)
            if base+low < previous-1e-9:
                raise ValueError('Hydrated run inventory repeats or reorders an original source interval')
            previous = base+high
            canonical.append((sid, bool(forward), float(first), float(last)))
        result.append(tuple(canonical))
    return tuple(result)


def _token(context, bundle, span, role, end):
    endpoint = span.head if end == 0 else span.tail
    mouth = endpoint.prepared_mouth
    if mouth is not None:
        if (not isinstance(mouth, PreparedMergedMouth) or endpoint.owner != mouth.owner
                or endpoint.approach_index != mouth.approach_index or mouth.approach.end_index != end):
            raise ValueError('Hydrated leaf lost its authoritative merged boundary token')
        PairedSourceBodyPiece._mouth(context, mouth, end)
        return 'MERGED', mouth.owner, mouth.approach_index
    member = bundle.members[role]
    visit = member.path[0 if end == 0 else -1]
    original = End(visit.section_id, int(not visit.forward) if end == 0 else int(visit.forward))
    authored = endpoint.original_terminal or endpoint.endpoint
    if (endpoint.kind != 'PHYSICAL' or authored != original or endpoint.registry_key is not None
            or endpoint.port is None):
        raise ValueError('An unmatched physical leaf boundary is not an actual original Bundle terminal')
    if endpoint.original_terminal is not None:
        visits = span.original_visits
        visit = visits[0 if end == 0 else -1] if visits else None
        if (visit is None or original != End(visit.section_id, int(not visit.forward) if end == 0 else int(visit.forward))
                or not any(claim['section'] == original.section_id and claim['first'] == 0.
                    and claim['last'] == context[Key('path', original.section_id)].length
                    and claim['owner'] == endpoint.owner
                    and claim['source_path'] is context[Key('path', original.section_id)]
                    for claim in span.owned_source_intervals)):
            raise ValueError('Moved outer terminal lost its original complete owner-consumption proof')
    return 'OUTER', end


def _claims(context, bundle, spans_by_role, paired_owners):
    membership = {visit.section_id: (role, visit.forward) for role, member in enumerate(bundle.members)
                  for visit in member.path}
    claims = []
    for owner in sorted(paired_owners):
        product = context[surface_key(owner)]
        layout = product.group.physical_layout
        false = {}
        for index, (approach, station) in enumerate(zip(layout.approaches, layout.outward_stations)):
            if not hasattr(approach, 'source_mouths'): continue
            aliases = source_mouths(approach, station)
            relevant = [alias for alias in aliases if alias.endpoint.section_id in membership]
            if not relevant: continue
            if not hasattr(approach, 'end_index') or len(aliases) != 2:
                if any(not alias.original_incidence for alias in relevant):
                    raise ValueError('An owned interior source cut needs its paired physical mouth policy')
                continue  # True outer profiles are checked explicitly below.
            if len(relevant) != 2:
                raise ValueError('Owned source claim needs both original roles of this paired Bundle')
            if {membership[alias.endpoint.section_id][0] for alias in relevant} != {0, 1}:
                raise ValueError('A paired physical claim changed its original lateral roles')
            PreparedMergedMouth.from_product(product, index)  # Full accepted physical rim is required.
            for alias in aliases:
                sid = alias.endpoint.section_id
                if sid not in membership: continue
                role, forward = membership[sid]
                length = context[Key('path', sid)].length
                if not 0. <= alias.station <= length:
                    raise ValueError('Owned paired source alias leaves its original finite Section')
                if alias.original_incidence:
                    first, last = (0., alias.station) if alias.endpoint.end == 0 else (alias.station, length)
                    claims.append((role, sid, forward, first, last, owner, 'original_profile_prefix'))
                else:
                    false.setdefault((role, sid, forward), []).append(alias.station)
        for (role, sid, forward), stations in false.items():
            if len(stations) != 2 or stations[0] == stations[1]:
                raise ValueError('Owned opposite source interval needs both distinct authoritative paired cuts')
            claims.append((role, sid, forward, min(stations), max(stations), owner, 'derived_paired_interval'))
        for claim in product.group.claims:
            if claim.section_id in membership:
                role, forward = membership[claim.section_id]
                claims.append((role, claim.section_id, forward, 0.,
                               context[Key('path', claim.section_id)].length, owner, 'complete_source_claim'))
    for role, spans in enumerate(spans_by_role):
        for span in spans:
            for claim in span.owned_source_intervals:
                sid = claim['section']
                if (sid not in membership or claim['source_path'] is not context[Key('path', sid)]
                        or claim['proof'] != 'complete_original_approach_route_prefix'):
                    raise ValueError('Owned source prefix lost its original controlled Path proof')
                claims.append((role, sid, membership[sid][1], claim['first'], claim['last'],
                               claim['owner'], claim['proof']))
            for endpoint in (span.head, span.tail):
                if endpoint.prepared_mouth is not None: continue
                sid, end = endpoint.endpoint.section_id, endpoint.endpoint.end
                port = endpoint.port
                if (port is None or port.mouth.endpoint != endpoint.endpoint or port.mouth.owner != endpoint.owner
                        or port.station != endpoint.source_station):
                    raise ValueError('Original terminal ownership needs its actual complete sampled profile')
                product = context.products.get(surface_key(endpoint.owner))
                visible = context[PORTS].resolve(endpoint.endpoint)
                if product is None or visible is None or product.pavement.mesh is None:
                    raise ValueError('Original terminal prefix lost its accepted physical product')
                if isinstance(product, CompoundProduct):
                    actual = VisibleMouth.from_compound(context, endpoint.endpoint, role, membership[sid][1])
                    if (not np.array_equal(actual.parameters, port.mouth.parameters)
                            or not np.array_equal(actual.physical_xyz, port.mouth.physical_xyz)):
                        raise ValueError('Original terminal ownership changed its complete accepted profile')
                else:
                    points = port.mouth.physical_xyz.copy()
                    points[:, 2] = port.profile.datum+product.exaggeration*(points[:, 2]-port.profile.datum)
                    if any(np.count_nonzero(np.all(product.pavement.mesh.xyz == point, axis=1)) != 1 for point in points):
                        raise ValueError('Original terminal profile is not one actual accepted physical boundary')
                length = context[Key('path', sid)].length
                first, last = (0., endpoint.source_station) if end == 0 else (endpoint.source_station, length)
                claims.append((role, sid, membership[sid][1], first, last, endpoint.owner, 'original_terminal_profile'))
    return tuple(dict.fromkeys(row for row in claims if row[4]-row[3] > 1e-9))


def _scope_audit(context, bundle, inventory, claims, scope):
    if scope is None:
        return dict(complete_supplied_inventory=True, complete_selected_scope=False)
    scope = _ordered(context, bundle, scope)
    for role, rows in enumerate(inventory):
        for sid, forward, first, last in rows:
            if not any(sid == s and forward == f and a <= first+1e-9 and last <= b+1e-9
                       for s, f, a, b in scope[role]):
                raise ValueError('A hydrated body interval lies outside the declared selected source scope')
    count = 0
    for role, rows in enumerate(scope):
        for sid, forward, low, high in rows:
            parts = [(max(low, a), min(high, b), ('body', i))
                     for i, (s, f, a, b) in enumerate(inventory[role]) if s == sid and f == forward]
            parts.extend((max(low, a), min(high, b), ('owner', owner))
                         for r, s, f, a, b, owner, _ in claims if (r, s, f) == (role, sid, forward))
            by_owner = {}
            for a, b, owner in parts:
                if b > a+1e-9: by_owner.setdefault(owner, []).append((a, b))
            merged = []
            for owner, intervals in by_owner.items():
                for a, b in sorted(intervals):
                    if merged and merged[-1][2] == owner and a <= merged[-1][1]+1e-9:
                        merged[-1] = (merged[-1][0], max(b, merged[-1][1]), owner)
                    else: merged.append((a, b, owner))
            cursor = low
            for a, b, owner in sorted(merged):
                if abs(a-cursor) > 1e-9:
                    raise ValueError('Body and actual junction ownership do not partition the complete selected source scope')
                cursor = b
            if abs(cursor-high) > 1e-9:
                raise ValueError('Selected original source scope has an unowned interval')
            count += 1
    return dict(complete_supplied_inventory=True, complete_selected_scope=True,
                scope_intervals=count, original_source_station_tolerance=1e-9)


class HydratedRunAssembly:
    max_cells = 50000

    def build(self, context, bundle, spans_by_role, *, frame, namespace, expected_intervals,
              scope_intervals=None):
        try: hash(namespace)
        except TypeError: raise ValueError('Hydrated run assembly needs a stable hashable namespace')
        if context.network.settings.exaggeration != 1.:
            raise ValueError('Hydrated run assembly requires an authoritative physical product context')
        expected = _ordered(context, bundle, expected_intervals)
        observed = _ordered(context, bundle, interval_inventory(spans_by_role))
        if observed != expected or not all(observed):
            raise ValueError('Hydrated run assembly must retain the complete supplied original interval inventory')
        indexed, mouths = [], {}
        incident = {end: site for site in (*context.network.intersections, *context.network.connections) for end in site.ends}
        for role, spans in enumerate(spans_by_role):
            leaves = {}
            for span in spans:
                if span.role != role: raise ValueError('Hydrated leaf changed its original lateral role')
                for first, last in zip(span.visits, span.visits[1:]):
                    site = incident.get(End(first.section_id, int(first.forward)))
                    if (site is None or len(site.ends) != 2
                            or End(last.section_id, int(not last.forward)) not in site.ends):
                        raise ValueError('Hydrated source leaves need original two-ended consecutive seams')
                tokens = tuple(_token(context, bundle, span, role, end) for end in (0, 1))
                if tokens in leaves: raise ValueError('Hydrated run repeats one owned source leaf')
                leaves[tokens] = span
                for end, token in enumerate(tokens):
                    if token[0] == 'MERGED':
                        mouth = (span.head if end == 0 else span.tail).prepared_mouth
                        if token in mouths and (mouth.approach is not mouths[token].approach
                                                or not np.array_equal(mouth.xyz, mouths[token].xyz)):
                            raise ValueError('Paired source roles disagree on an authoritative physical mouth')
                        mouths[token] = mouth
            indexed.append(leaves)
        if set(indexed[0]) != set(indexed[1]):
            raise ValueError('Hydrated source leaves have missing or unmatched physical ownership boundaries')
        paired_owners = {mouth.owner for mouth in mouths.values()}
        physical_owners = paired_owners | {end.owner for own in spans_by_role for span in own
                                           for end in (span.head, span.tail) if end.prepared_mouth is None}
        phase = verify_shared_source_phase(tuple({id(m.approach): m.approach for m in mouths.values()}.values()))
        claims = _claims(context, bundle, spans_by_role, paired_owners)
        scope = _scope_audit(context, bundle, observed, claims, scope_intervals)
        paving, grassy = Assembler(), Assembler()
        reports, pieces, total = [], [], 0
        limit = self.max_cells*context.network.settings.resolution
        for number, tokens in enumerate(indexed[0]):
            spans = tuple(indexed[role][tokens] for role in (0, 1))
            name = (namespace, 'source piece', number)
            merged = [i for i, token in enumerate(tokens) if token[0] == 'MERGED']
            if len(merged) == 2:
                piece = PairedSourceBodyPiece()
                selected = tuple(mouths[token] for token in tokens)
                pavement, grass, report = piece.build(context, bundle, spans, selected, frame=frame, namespace=name)
                keys = list(piece.pavement_keys)
                fronts = {end: piece.front_ids[end] for end in (0, 1)}
            elif len(merged) == 1:
                end = merged[0]; mouth = mouths[tokens[end]]
                curves, ids, _ = PairedSourceBodyPiece._sources(context, bundle, spans, frame)
                field = AngularCrown(frame, curves, datum=mouth.approach.field.datum, rise=mouth.approach.field.rise)
                target = SimpleNamespace(frame=frame, field=field, source_sections=ids,
                                         source_visits=tuple(span.visits for span in spans), entrance=mouth.approach)
                bridge = MergedProfileBridge()
                connector, connector_report = bridge.build(context, mouth.approach, mouth.station, mouth.xyz,
                                                           target_inventory=target)
                piece = HydratedExternalRunWithCap()
                pavement, grass, report = piece.build(context, bundle, spans, mouth, bridge.body_cut,
                                                       frame=frame, namespace=name)
                joined = Assembler()
                connector_keys = [(name, 'merged connector', i) for i in range(len(connector.xyz))]
                for node, key in zip(bridge.back_ids, piece.merged_cut_keys): connector_keys[int(node)] = key
                joined.add(connector.xyz, connector.quads, connector_keys)
                joined.add(pavement.xyz, pavement.quads, piece.pavement_keys)
                pavement = joined.mesh()
                keys = [None]*len(pavement.xyz)
                for key, at in joined.lookup.items(): keys[at] = key
                fronts = {end: np.asarray([joined.lookup[connector_keys[int(i)]] for i in bridge.front_ids])}
                report = dict(report, actual_merged_connector=connector_report)
            else:
                raise ValueError('A separated source run retains its established provider')
            for end, ids in fronts.items():
                mouth = mouths[tokens[end]]
                if not np.array_equal(pavement.xyz[ids], mouth.xyz):
                    raise ValueError('Hydrated assembly moved an authoritative physical mouth')
                for node, reference in zip(ids, mouth.rim_references):
                    keys[int(node)] = (namespace, 'physical mouth', mouth.owner, mouth.approach_index, *reference)
            paving.add(pavement.xyz, pavement.quads, keys)
            if grass is not None: grassy.add(grass.xyz, grass.quads, piece.grass_keys)
            total += len(pavement.quads)+(len(grass.quads) if grass is not None else 0)
            if total > limit: raise ValueError('Hydrated source assembly exceeds its complete run quad budget')
            reports.append(report);pieces.append(piece)
        pavement = paving.mesh()
        grass = grassy.mesh('grass') if grassy.quads else None
        payloads = {Key('hydrated_run_pavement', 0): pavement}
        if grass is not None: payloads[Key('hydrated_run_grass', 0)] = grass
        query = Assembler()
        query.add(pavement.xyz, pavement.quads)
        if grass is not None: query.add(grass.xyz, grass.quads)
        self.body_query_mesh = query.mesh()
        first = next(iter(mouths.values())).approach.field
        ceiling = first.ceiling
        for owner in sorted(physical_owners):
            product = context[surface_key(owner)]
            for material, part in (('pavement', product.pavement), ('grass', product.grass)):
                if part is None: continue
                if part.mesh is None: raise ValueError('Owned paired junction query needs its authoritative physical meshes')
                payloads[Key('owned_junction_'+material, owner.id)] = part.mesh
                if owner in paired_owners:
                    query.add(part.mesh.xyz, part.mesh.quads)
            if owner in paired_owners:
                if abs(product.field.datum-first.datum) > 1e-8:
                    raise ValueError('Owned paired junction queries need the original common source datum')
                ceiling = max(ceiling, product.field.ceiling)
        overlap = validate_overlaps(payloads, cell_size=2.)
        self.pavement_keys = [None]*len(pavement.xyz)
        for key, at in paving.lookup.items(): self.pavement_keys[at] = key
        self.grass_keys = [None]*len(grass.xyz) if grass is not None else []
        for key, at in grassy.lookup.items(): self.grass_keys[at] = key
        self.query_mesh = query.mesh()
        self.query = QuadHeightField(self.query_mesh, datum=first.datum, ceiling=ceiling)
        self.pieces, self.paired_owners = tuple(pieces), tuple(sorted(paired_owners))
        self.source_intervals = observed
        self.transition_owners = tuple(sorted({owner for own in spans_by_role for span in own
                                               for owner in span.transition_owners}))
        return pavement, grass, dict(construction='qualified hydrated original source run assembly',
            original_interval_inventory=observed, coverage=scope, source_phase=phase,
            junction_owned_intervals=[dict(role=r, section=s, forward=f, first=a, last=b,
                                            owner=owner.text(), proof=proof)
                                      for r, s, f, a, b, owner, proof in claims],
            pieces=reports, total_quads=total, quad_budget=limit, overlap=overlap,
            query_only_junction_owners=[owner.text() for owner in self.paired_owners],
            original_terminal_owners=[owner.text() for owner in sorted(physical_owners-paired_owners)],
            query_scope='published body plus paired-owner meshes; original terminal queries retain their JunctionSpan',
            junction_geometry_republished=False, runtime_discovery_enabled=False,
            exact_original_mouth_fronts=True, all_direct_quads=True)

    def height(self, xy):
        xy = np.asarray(xy, dtype=float)
        if xy.shape[-1] != 2: raise ValueError('Hydrated source height queries need XY coordinates')
        return self.query.evaluate(xy.reshape(-1, 2), default=self.query.datum).reshape(xy.shape[:-1])
