"""Explicit prepared dimensional seams for original physical Street surfaces."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class PhysicalDimensionTransition:
    street_id: int
    before_section: int
    after_section: int
    source_kind: str
    source_id: int
    minimum_half_length: float
    maximum_edge_slope: float


class PreparedPhysicalDimensions:
    def read(self, snapshot, street_ids):
        """Retain prepared transition identity; never invent a missing taper.

        The caller assigns disposable native table indices. Stable ownership is
        the original Street, prepared item kind/ID and adjacent Section IDs.
        Geometry, dimensions and station calculations remain native inputs.
        """
        return self._read(snapshot, street_ids)

    def read_spans(self, snapshot, spans, *, source_plan=None):
        """Read tapers inside original spans, excluding their junction joins.

        Width changes across a native junction are owned by that junction.
        A deferred attachment does not cut a span and therefore cannot conceal
        a missing taper. Check the ranges against the original Street items so
        a caller cannot invent shorter spans to bypass a dimensional change.
        """
        requested = tuple(spans)
        if source_plan is not None and any(span not in source_plan.spans for span in requested):
            raise ValueError('Dimensional scope must retain its original port partition')
        if len({s.identity for s in requested}) != len(requested):
            raise ValueError('Physical dimension selection repeats an original span')
        scopes = {}
        used = set()
        for span in requested:
            sid = span.street_id
            source = snapshot.streets.get(sid)
            if source is None or not (source_plan.source_supported(source) if source_plan else source.get('status') == 'SUPPORTED'):
                raise ValueError(f'Physical dimension source Street {sid} is missing or unsupported')
            sections = [s for s in source['items'] if s['kind'] == 'Section']
            ids = tuple(s['id'] for s in sections)
            if not ids or len(set(ids)) != len(ids):
                raise ValueError(f'Street {sid} has missing or multiply owned Sections')
            controls = [0]
            for section in sections:
                if len(section['points']) < 2:
                    raise ValueError('A dimensional source Section needs original controls')
                controls.append(controls[-1] + len(section['points']) - 1)
            cuts = {0: None, len(ids): None}
            category_cuts = {}
            selection_cuts = set()
            if source_plan is not None:
                for scope in source_plan.section_scopes:
                    if scope.street_id != sid:continue
                    for boundary in scope.boundaries:
                        cut = ids.index(boundary.after_section)
                        category_cuts[cut] = boundary
                        cuts[cut] = None
                for boundary in source_plan.selection_boundaries:
                    if boundary.endpoint.street_id != sid:continue
                    peer = source_plan.selection_peer(boundary.endpoint)
                    if peer is None:continue
                    cut = ids.index(boundary.endpoint.section_id)+boundary.endpoint.end
                    if cut in cuts and cut not in selection_cuts:
                        raise ValueError('Interior selection cut conflicts with another span boundary')
                    if not 0 < cut < len(ids):
                        raise ValueError('Paired selection cut must retain an interior original seam')
                    expected = {(ids[cut-1],1),(ids[cut],0)}
                    if {(boundary.endpoint.section_id,boundary.endpoint.end),(peer.section_id,peer.end)} != expected:
                        raise ValueError('Paired selection cut lost its consecutive original Sections')
                    selection_cuts.add(cut);cuts[cut] = None
            for item in source['items']:
                if item.get('surface_treatment') != 'NATIVE_JUNCTION':
                    continue
                ends = tuple(tuple(e) for e in item.get('through', ()))
                matches = [i for i in range(1, len(ids))
                           if len(ends) == 2 and set(ends) == {(ids[i-1], 1), (ids[i], 0)}]
                if len(matches) != 1 or matches[0] in cuts:
                    raise ValueError('Dimensional spans need unique original native-junction cuts')
                cuts[matches[0]] = item['id']
            pairs = list(zip(sorted(cuts), sorted(cuts)[1:]))
            matches = [(a, b) for a, b in pairs if ids[a:b] == span.section_ids]
            if len(matches) != 1:
                raise ValueError('Dimensional span differs from original owned Section range')
            a, b = matches[0]
            if (span.first_control, span.last_control, span.source_control_count) != (
                controls[a], controls[b], controls[-1] + 1
            ):
                raise ValueError('Dimensional span has stale original control bindings')
            for end, section_id, orientation, cut in (
                (span.start, ids[a], 0, a), (span.end, ids[b-1], 1, b)
            ):
                if (end.endpoint.street_id, end.endpoint.section_id, end.endpoint.end) != (sid, section_id, orientation):
                    raise ValueError('Dimensional span has a different original endpoint')
                if cuts[cut] is not None and (end.owner_kind, end.junction_id) != ('INTERIOR', cuts[cut]):
                    raise ValueError('Dimensional span lost its original native-junction owner')
                if cut in category_cuts and (end.owner_kind, end.junction_id) != ('CATEGORY_BOUNDARY', None):
                    raise ValueError('Dimensional span lost its original category boundary owner')
                if cut in selection_cuts and (end.owner_kind, end.junction_id) != ('SELECTION_BOUNDARY', None):
                    raise ValueError('Dimensional span lost its paired selection boundary owner')
            if used.intersection(span.section_ids):
                raise ValueError('Dimensional spans have overlapping original Section owners')
            used.update(span.section_ids)
            pairs = zip(span.section_ids, span.section_ids[1:]) if span.in_scope else ()
            scopes.setdefault(sid, set()).update(pairs)
        return self._read(snapshot, tuple(scopes), scopes=scopes, source_plan=source_plan)

    def _read(self, snapshot, street_ids, *, scopes=None, source_plan=None):
        requested = tuple(street_ids)
        if len(set(requested)) != len(requested):
            raise ValueError('Physical dimension selection contains duplicate Street owners')
        result = []
        for sid in requested:
            street = snapshot.streets.get(sid)
            if street is None or not (source_plan.source_supported(street) if source_plan else street.get('status') == 'SUPPORTED'):
                raise ValueError(f'Physical dimension source Street {sid} is missing or unsupported')
            items = street['items']
            sections = [item for item in items if item['kind'] == 'Section']
            if not sections or len({item['id'] for item in sections}) != len(sections):
                raise ValueError(f'Street {sid} has missing or multiply owned Sections')
            adjacent = {(a['id'], b['id']): (a, b) for a, b in zip(sections, sections[1:])}
            # Actual neighboring Sections determine whether an item belongs to
            # this scope; a corrupt declaration must not select its own scope.
            neighbors = {}
            previous = None
            pending = []
            for index, item in enumerate(items):
                if item['kind'] == 'Section':
                    for at in pending:
                        neighbors[at] = (previous, item['id'])
                    pending.clear()
                    previous = item['id']
                else:
                    pending.append(index)
            declared = set(); source_items = set()
            for index, item in enumerate(items):
                deferred = item['kind'] == 'Intersection' and item.get('surface_treatment') == 'DEFERRED_ATTACHMENT'
                if item['kind'] not in ('SymLane', 'SideLane', 'AuthoredTransition') and not deferred:
                    continue
                if scopes is not None and neighbors.get(index) not in scopes[sid]:
                    continue
                if item.get('surface_treatment') == 'COMPATIBLE_SEAM':
                    continue
                record = item.get('transition')
                if deferred and not record:
                    continue  # A changed seam still rejects in the coverage check below.
                if (not deferred and item.get('surface_treatment') != 'WIDTH_TRANSITION') or not record or item.get('points') or item.get('transition_reason'):
                    raise ValueError(f'Street {sid} has an unsupported prepared transition extent')
                pair = (record['before_section'], record['after_section'])
                if pair not in adjacent or pair in declared or not 0 < index < len(items)-1:
                    raise ValueError(f'Street {sid} transition must have unique consecutive Section owners')
                a, b = adjacent[pair]
                if items[index-1] is not a or items[index+1] is not b:
                    raise ValueError(f'Street {sid} transition does not occupy its declared original seam')
                if record['source_kind'] != item['kind'] or record['source_id'] != item['id']:
                    raise ValueError(f'Street {sid} transition lost its prepared source identity')
                if deferred:
                    ends = tuple(tuple(end) for end in item.get('through', ()))
                    location = item.get('location', ())
                    if (len(ends) != 2 or set(ends) != {(a['id'], 1), (b['id'], 0)}
                            or not item.get('attachments') or len(location) < 2
                            or not all(math.isfinite(v) for v in location)
                            or not a['points'] or math.dist(location[:2], a['points'][-1][:2]) > 1e-6):
                        raise ValueError(f'Street {sid} transition lost its original deferred attachment')
                identity = (item['kind'], item['id'])
                if identity in source_items:
                    raise ValueError(f'Street {sid} has multiply owned prepared transition identities')
                source_items.add(identity)
                half, slope = float(record['minimum_half_length']), float(record['maximum_edge_slope'])
                if not math.isfinite(half) or half <= 0 or not math.isfinite(slope) or not 0 < slope <= 1:
                    raise ValueError(f'Street {sid} transition has invalid extent settings')
                declared.add(pair)
                result.append(PhysicalDimensionTransition(sid, *pair, item['kind'], item['id'], half, slope))
            for pair, (a, b) in adjacent.items():
                if not a['points'] or not b['points'] or math.dist(a['points'][-1], b['points'][0]) > 1e-6:
                    raise ValueError(f'Street {sid} has disconnected original Section endpoints')
                if (scopes is None or pair in scopes[sid]) and (a['width'], a['offset']) != (b['width'], b['offset']) and pair not in declared:
                    raise ValueError(f'Street {sid} has an undeclared dimension change at Sections {pair}')
        return tuple(result)
