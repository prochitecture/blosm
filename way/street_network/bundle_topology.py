"""Explicit prepared Bundle extents, separate from surface-backend capability.

The source Street remains the editing/ownership unit. Spans are derived ranges
of its original Sections; they never allocate substitute Street or Section IDs.
Prepared head/tail locations identify those ranges once, before authoring.
Coordinate edits must subsequently use the retained references, not repeat a
nearest-point search. A managed topology edit must rebuild affected references.

This module has no Blender dependency. It does not claim that a mapped extent
has a valid offset surface, and does not infer connections from proximity.
"""
from dataclasses import asdict, dataclass, replace
import math


@dataclass(frozen=True)
class SectionAnchor:
    section_id: int
    segment_index: int
    factor: float


@dataclass(frozen=True)
class SectionRange:
    start: SectionAnchor
    end: SectionAnchor

    @property
    def forward(self):
        return (self.end.segment_index + self.end.factor
                > self.start.segment_index + self.start.factor)


@dataclass(frozen=True)
class SpanBoundary:
    bundle_end: str | None = None
    junction_id: int | None = None


@dataclass(frozen=True)
class MemberSpan:
    ordinal: int
    start: SpanBoundary
    end: SpanBoundary
    sections: tuple[SectionRange, ...]


@dataclass(frozen=True)
class MemberExtent:
    street_id: int
    forward: bool
    head_order: int
    tail_order: int
    complete: bool
    spans: tuple[MemberSpan, ...]


@dataclass(frozen=True)
class Attachment:
    street_id: int
    section_id: int
    end: int
    connector_id: int
    category: str
    carriageway: bool
    side: str | None = None  # Left/right in the authored through-Street direction.


@dataclass(frozen=True)
class MemberJunction:
    junction_id: int
    street_id: int
    through: tuple[tuple[int, int], ...]
    attachments: tuple[Attachment, ...]


@dataclass(frozen=True)
class ChainJoin:
    junction_id: int
    # Original (Street, Section, end) identities in Bundle head-to-tail order.
    before: tuple[int, int, int]
    after: tuple[int, int, int]


@dataclass(frozen=True)
class BundleTopology:
    bundle_id: int
    members: tuple[MemberExtent, ...]
    junctions: tuple[MemberJunction, ...]
    interior_streets: tuple[int, ...]
    required_features: tuple[str, ...]
    chains: tuple[tuple[int, ...], ...] = ()
    chain_joins: tuple[ChainJoin, ...] = ()
    exterior_streets: tuple[int, ...] = ()

    def report(self):
        result = dict(status='MAPPED', **asdict(self))
        if not self.chains:
            result.pop('chains'); result.pop('chain_joins')
        if not self.exterior_streets:
            result.pop('exterior_streets')
        return result


class BundleTopologyError(ValueError):
    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


class _StreetPath:
    """Control-edge rank is an addressing coordinate, never a metric station."""
    def __init__(self, street, tolerance):
        self.street = street
        self.tolerance = tolerance
        self.sections = []
        self.by_id = {}
        self.length = 0  # Number of authored edges, NOT distance in metres.
        previous = None
        for item in street['items']:
            if item['kind'] != 'Section':
                continue
            sid, points = item['id'], item['points']
            if sid in self.by_id:
                raise BundleTopologyError('DUPLICATE_SECTION', f'Section {sid} occurs twice in Street {street["id"]}')
            if len(points) < 2 or any(len(p) < 2 or not all(math.isfinite(v) for v in p[:2]) for p in points):
                raise BundleTopologyError('INVALID_SECTION', f'Section {sid} has no finite path')
            if previous is not None and math.dist(previous[:2], points[0][:2]) > tolerance:
                raise BundleTopologyError('DISCONNECTED_CHAIN', f'Section {sid} does not meet its predecessor')
            if any(math.dist(a[:2], b[:2]) <= 2*tolerance for a, b in zip(points, points[1:])):
                raise BundleTopologyError('DEGENERATE_EDGE', f'Section {sid} has a collapsed control edge')
            entry = (item, self.length, self.length + len(points)-1)
            self.sections.append(entry)
            self.by_id[sid] = entry
            self.length = entry[2]
            previous = points[-1]
        if not self.sections:
            raise BundleTopologyError('EMPTY_MEMBER', f'Street {street["id"]} contains no Section')

    def locate(self, point):
        if len(point) < 2 or not all(math.isfinite(v) for v in point[:2]):
            raise BundleTopologyError('INVALID_LOCATION', 'Bundle endpoint is not a finite planar location')
        found = set()
        for item, first, _last in self.sections:
            for i, (a, b) in enumerate(zip(item['points'], item['points'][1:])):
                d = (b[0]-a[0], b[1]-a[1])
                t = ((point[0]-a[0])*d[0] + (point[1]-a[1])*d[1]) / (d[0]*d[0]+d[1]*d[1])
                t = min(1., max(0., t))
                if math.dist(point[:2], a[:2]) <= self.tolerance:
                    t = 0.
                elif math.dist(point[:2], b[:2]) <= self.tolerance:
                    t = 1.
                if math.dist(point[:2], (a[0]+t*d[0], a[1]+t*d[1])) <= self.tolerance:
                    found.add(first+i+t)
        if not found:
            raise BundleTopologyError('UNLOCATED_ENDPOINT', f'Bundle endpoint is not on Street {self.street["id"]}')
        if len(found) != 1:
            raise BundleTopologyError('AMBIGUOUS_ENDPOINT', f'Bundle endpoint has several incidences on Street {self.street["id"]}')
        return found.pop()

    def endpoint(self, sid, end):
        if sid not in self.by_id or end not in (0, 1):
            raise BundleTopologyError('INVALID_THROUGH_REFERENCE', 'Minor junction references an absent Section or endpoint')
        _item, start, finish = self.by_id[sid]
        return start if end == 0 else finish

    def ranges(self, start, end):
        lo, hi = sorted((start, end))
        forward = start < end
        result = []
        for item, first, last in self.sections:
            a, b = max(lo, first), min(hi, last)
            if b <= a:
                continue

            def anchor(value):
                q = value-first
                segment = min(len(item['points'])-2, int(math.floor(q)))
                return SectionAnchor(item['id'], segment, q-segment)

            result.append(SectionRange(anchor(a if forward else b), anchor(b if forward else a)))
        return tuple(result if forward else reversed(result))


class _BundleOwnership:
    """Separate broad preparation ownership from physical median interiors.

    Bundle construction also owns attached exterior Streets. Membership in that
    inventory alone cannot mean that a road crosses the separator. Only prove
    an exterior role when the finite authored path stays outside a simple
    two-member extent polygon and every known attachment faces outwards.
    Ambiguous, coincident, crossing and genuinely interior runs retain their
    previous interior role for explicit downstream validation.
    """
    def __init__(self, members, streets, events, tolerance):
        self.tolerance = tolerance
        self.events = events
        self.members = {member.street_id: member for member in members}
        self.rim = ()
        if len(members) != 2:
            return
        paths = [self._points(member, streets[member.street_id]) for member in members]
        rim=[]
        for point in paths[0]+list(reversed(paths[1])):
            if not rim or math.dist(rim[-1],point)>tolerance:rim.append(point)
        if len(rim)>1 and math.dist(rim[0],rim[-1])<=tolerance:rim.pop()
        if len(rim)<3:return
        if any(math.dist(a, b) <= tolerance for a, b in zip(rim, rim[1:] + rim[:1])):
            return
        self.rim = tuple(rim)
        edges = list(zip(rim, rim[1:] + rim[:1]))
        # A self-intersecting extent has no reliable inside/outside. This is a
        # classification guard, not a substitute for native Bundle validation.
        for i, (a, b) in enumerate(edges):
            for j in range(i + 1, len(edges)):
                if j == i + 1 or (i == 0 and j == len(edges) - 1):
                    continue
                if self._cuts(a, b, *edges[j]) is not None:
                    self.rim = ()
                    return
        origin = rim[0]
        area = sum(self._cross((a[0]-origin[0], a[1]-origin[1]),
                               (b[0]-origin[0], b[1]-origin[1])) for a, b in edges)
        if abs(area) <= tolerance * tolerance:
            self.rim = ()
            return
        # The first path follows the rim, the second is reversed in it.
        self.inside = {members[0].street_id: 'LEFT' if area > 0 else 'RIGHT',
                       members[1].street_id: 'RIGHT' if area > 0 else 'LEFT'}

    @staticmethod
    def _cross(a, b):
        return a[0]*b[1] - a[1]*b[0]

    @staticmethod
    def _points(member, street):
        sections = {item['id']: item for item in street['items'] if item['kind'] == 'Section'}
        result = []
        for span in member.spans:
            for extent in span.sections:
                points = sections[extent.start.section_id]['points']
                first = extent.start.segment_index + extent.start.factor
                last = extent.end.segment_index + extent.end.factor
                def sample(parameter):
                    index = min(int(math.floor(parameter)), len(points)-2)
                    t = parameter-index
                    return tuple(points[index][k] + t*(points[index+1][k]-points[index][k]) for k in (0, 1))
                lo, hi = sorted((first, last))
                controls = [float(i) for i in range(int(math.floor(lo))+1, int(math.ceil(hi))) if lo < i < hi]
                if first > last:
                    controls.reverse()
                piece = [sample(value) for value in (first, *controls, last)]
                result.extend(piece[1:] if result else piece)
        return result

    def _cuts(self, a, b, c, d):
        """None means disjoint; () means collinear overlap/touch; else cuts."""
        u, v = (b[0]-a[0], b[1]-a[1]), (d[0]-c[0], d[1]-c[1])
        w = c[0]-a[0], c[1]-a[1]
        lu, lv = math.hypot(*u), math.hypot(*v)
        determinant = self._cross(u, v)
        if abs(determinant) <= self.tolerance * max(lu, lv):
            if abs(self._cross(w, u)) > self.tolerance * lu:
                return None
            parameters = tuple(((p[0]-a[0])*u[0] + (p[1]-a[1])*u[1])/(lu*lu) for p in (c, d))
            return () if max(min(parameters), 0.) <= min(max(parameters), 1.) + self.tolerance/lu else None
        t, s = self._cross(w, v)/determinant, self._cross(w, u)/determinant
        if (-self.tolerance/lu <= t <= 1.+self.tolerance/lu
                and -self.tolerance/lv <= s <= 1.+self.tolerance/lv):
            return (min(1., max(0., t)),)
        return None

    def _location(self, point):
        inside = False
        for a, b in zip(self.rim, self.rim[1:] + self.rim[:1]):
            u = b[0]-a[0], b[1]-a[1]
            length2 = u[0]*u[0] + u[1]*u[1]
            t = min(1., max(0., ((point[0]-a[0])*u[0] + (point[1]-a[1])*u[1])/length2))
            if math.dist(point, (a[0]+t*u[0], a[1]+t*u[1])) <= self.tolerance:
                return 'BOUNDARY'
            if (a[1] > point[1]) != (b[1] > point[1]):
                if point[0] < a[0] + (point[1]-a[1])*u[0]/u[1]:
                    inside = not inside
        return 'INSIDE' if inside else 'OUTSIDE'

    def exterior(self, street):
        if not self.rim:
            return False
        for event in self.events:
            member = self.members[event.street_id]
            inside = self.inside[event.street_id]
            if not member.forward:
                inside = 'RIGHT' if inside == 'LEFT' else 'LEFT'
            for attachment in event.attachments:
                if attachment.street_id == street['id'] and attachment.side in (None, inside):
                    return False
        segments = [pair for item in street['items'] if item['kind'] == 'Section'
                    for pair in zip(item['points'], item['points'][1:])]
        if not segments:
            return False
        for a, b in segments:
            if any(len(point) < 2 or not all(math.isfinite(v) for v in point[:2]) for point in (a, b)):
                return False
            if math.dist(a[:2], b[:2]) <= self.tolerance:
                return False
            cuts = {0., 1.}
            for c, d in zip(self.rim, self.rim[1:] + self.rim[:1]):
                found = self._cuts(a, b, c, d)
                if found == ():
                    return False
                if found is not None:
                    cuts.update(found)
            ordered = sorted(cuts)
            for first, last in zip(ordered, ordered[1:]):
                t = .5*(first+last)
                if self._location(tuple(a[k]+t*(b[k]-a[k]) for k in (0, 1))) != 'OUTSIDE':
                    return False
        return True


class PreparedBundleTopology:
    """Map declared head/tail incidences and internal junction cuts exactly.

    An ordered member can contain consecutive original Streets when original
    physical incidences prove every join. Ambiguous branches are not paired by
    proximity. Interior owned Streets remain separate from lateral members.
    Native surface support is reviewed independently by PreparedBundleReview.
    """
    def __init__(self, categories=(), tolerance=1e-5):
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('Bundle addressing tolerance must be positive and finite')
        self.categories = frozenset(categories)
        self.tolerance = tolerance

    def _junctions(self, path, streets):
        result = []
        for item in path.street['items']:
            if item.get('surface_treatment') != 'NATIVE_JUNCTION':
                continue
            refs = tuple(tuple(ref) for ref in item.get('through', ()))
            if len(refs) != 2 or refs[0] == refs[1]:
                raise BundleTopologyError('INVALID_THROUGH_REFERENCE', 'Internal junction needs two distinct through Section ends')
            stations = [path.endpoint(*ref) for ref in refs]
            if stations[0] != stations[1] or not 0 < stations[0] < path.length:
                raise BundleTopologyError('INVALID_THROUGH_REFERENCE', 'Internal junction must join consecutive Section endpoints')
            if abs(path.locate(item['location'])-stations[0]) > 1e-9:
                raise BundleTopologyError('JUNCTION_LOCATION_MISMATCH', 'Junction location differs from its through Section endpoints')
            attachments = []
            seen = set()
            for a in item.get('attachments', ()):
                if a.get('side') not in (None, 'LEFT', 'RIGHT'):
                    raise BundleTopologyError('INVALID_ATTACHMENT_SIDE', 'Attachment side must be declared LEFT or RIGHT')
                if a['street_id'] not in streets or a['end'] not in (0, 1):
                    raise BundleTopologyError('INVALID_ATTACHMENT',
                        f'Junction {item["id"]}: attachment Street {a["street_id"]}, Section {a["section_id"]}, end {a["end"]} is absent or invalid')
                sections = [s for s in streets[a['street_id']]['items'] if s['kind'] == 'Section' and s['id'] == a['section_id']]
                if len(sections) != 1 or len(sections[0]['points']) < 2:
                    raise BundleTopologyError('INVALID_ATTACHMENT',
                        f'Junction {item["id"]}: Section {a["section_id"]} is absent or ambiguous in attachment Street {a["street_id"]}')
                point = sections[0]['points'][0 if a['end'] == 0 else -1]
                if math.dist(point[:2], item['location'][:2]) > self.tolerance:
                    raise BundleTopologyError('ATTACHMENT_LOCATION_MISMATCH', 'Attachment endpoint does not meet its declared junction')
                key = (a['section_id'], a['end'])
                if key in seen:
                    raise BundleTopologyError('DUPLICATE_ATTACHMENT', 'A Section end occurs twice at an internal junction')
                seen.add(key)
                attachments.append(Attachment(a['street_id'], a['section_id'], a['end'], a['connector_id'],
                    a['category'], a['category'] in self.categories, a.get('side')))
            if not attachments:
                raise BundleTopologyError('MISSING_ATTACHMENT', 'Internal junction has no declared attachment')
            result.append((stations[0], MemberJunction(item['id'], path.street['id'], refs, tuple(attachments))))
        if len({station for station, _event in result}) != len(result):
            raise BundleTopologyError('DUPLICATE_JUNCTION', 'Several junction records occupy one member boundary')
        return result

    def _chained(self, bundle, streets, origins):
        """Trace whole original Streets between the declared lateral terminals.

        No source is concatenated or renumbered. A shared physical site must
        have exactly two incident Streets of this Bundle at a continuation;
        side Streets of other Bundles retain their own identities and roles.
        """
        head, tail = tuple(bundle['head']), tuple(bundle['tail'])
        owned = tuple(bundle.get('owned_streets', ()))
        if (len(head) != len(tail) or not owned or len(set(owned)) != len(owned)
                or not set(head + tail) <= set(owned)):
            raise BundleTopologyError('UNRESOLVED_MEMBERSHIP', 'Chain terminals need complete original Bundle ownership')
        paths = {}
        for sid in owned:
            if sid not in streets or streets[sid].get('bundle_id') != bundle['id']:
                raise BundleTopologyError('MISSING_MEMBER', 'Chain member lacks its original Bundle ownership')
            paths[sid] = _StreetPath(streets[sid], self.tolerance)
        def terminal(sid, location):
            path = paths[sid]; station = path.locate(location)
            if station not in (0, path.length):
                raise BundleTopologyError('PARTIAL_CHAIN_MEMBER', 'A longitudinal chain needs complete original Streets')
            section = path.sections[0 if station == 0 else -1][0]
            return sid, section['id'], int(station != 0)
        starts = tuple(terminal(sid, p) for sid, p in zip(head, bundle['head_locations']))
        stops = tuple(terminal(sid, p) for sid, p in zip(tail, bundle['tail_locations']))
        if set(starts).intersection(stops):
            raise BundleTopologyError('EMPTY_EXTENT', 'A Bundle terminal cannot be both a chain head and tail')
        incidence, sites = {}, {}
        for site in origins.get('physical', ()):
            jid = site['junction_id']
            if jid in sites:
                raise BundleTopologyError('AMBIGUOUS_CHAIN_JOIN', 'Duplicate original physical site')
            entries = []
            for connector in site['ring']:
                for endpoint in connector['endpoints']:
                    sid = endpoint['street_id']
                    if sid not in paths:
                        continue
                    if connector.get('item_kind') != 'Street' or len(connector['endpoints']) != 1:
                        raise BundleTopologyError('INVALID_CHAIN_JOIN', 'Chain joins require individual original Street incidences')
                    end = endpoint['end']; path = paths[sid]
                    if end not in (0, 1) or endpoint['terminal_kind'] != 'Section':
                        raise BundleTopologyError('INVALID_CHAIN_JOIN', 'Chain join needs original terminal Section ends')
                    section = path.sections[-1 if end else 0][0]
                    key = (sid, section['id'], end)
                    if (endpoint['terminal_id'] != section['id'] or key in incidence or
                            math.dist(section['points'][-1 if end else 0][:2], site['location'][:2]) > self.tolerance):
                        raise BundleTopologyError('INVALID_CHAIN_JOIN', 'Physical incidence differs from the original chain endpoint')
                    incidence[key] = jid; entries.append(key)
            sites[jid] = tuple(entries)
        chains, joins, used = [], [], set()
        for order, start in enumerate(starts):
            entry = start; chain = []
            while True:
                sid, _, end = entry
                if sid in used:
                    raise BundleTopologyError('REPEATED_CHAIN_MEMBER', 'Longitudinal chains revisit or share an original Street')
                used.add(sid); chain.append((sid, end == 0))
                path = paths[sid]; section = path.sections[-1 if end == 0 else 0][0]
                exit_port = (sid, section['id'], 1 - end)
                if exit_port in stops:
                    if exit_port != stops[order]:
                        raise BundleTopologyError('MEMBER_ORDER_CHANGE', 'Longitudinal chain order changes between Bundle ends')
                    break
                jid = incidence.get(exit_port)
                ring = sites.get(jid, ())
                if len(ring) != 2 or exit_port not in ring:
                    raise BundleTopologyError('AMBIGUOUS_CHAIN_JOIN',
                        f'Chain continuation needs exactly two original Bundle incidences; '
                        f'Bundle {bundle["id"]}, junction {jid}, exit {exit_port}, incidences {ring}')
                next_port = next(p for p in ring if p != exit_port)
                if next_port in starts or next_port in stops or next_port[2] == exit_port[2]:
                    raise BundleTopologyError('INVALID_CHAIN_JOIN', 'Chain join conflicts with terminal or directed Street incidences')
                joins.append(ChainJoin(jid, exit_port, next_port)); entry = next_port
            chains.append(tuple(chain))
        if used != set(owned):
            raise BundleTopologyError('UNRESOLVED_MEMBERSHIP', 'Owned Streets outside the terminal chains need explicit interior roles')
        members, events, features, sections = [], [], {'LONGITUDINAL_CHAINS'}, set()
        if len(chains) != 2:
            features.add('MULTIPLE_MEMBERS')
        before = {j.before[0]: j.junction_id for j in joins}
        after = {j.after[0]: j.junction_id for j in joins}
        for order, chain in enumerate(chains):
            for i, (sid, forward) in enumerate(chain):
                path = paths[sid]
                if sections.intersection(path.by_id):
                    raise BundleTopologyError('DUPLICATE_SECTION', 'Chain members share a semantic Section ID')
                sections.update(path.by_id)
                junctions = self._junctions(path, streets); events.extend(e for _, e in junctions)
                if junctions: features.add('INTERNAL_JUNCTIONS')
                start, end = (0, path.length) if forward else (path.length, 0)
                cuts = {station: SpanBoundary(junction_id=e.junction_id) for station, e in junctions}
                cuts[start] = SpanBoundary('HEAD' if i == 0 else None, after.get(sid))
                cuts[end] = SpanBoundary('TAIL' if i == len(chain)-1 else None, before.get(sid))
                ordered = sorted(cuts, reverse=not forward)
                spans = tuple(MemberSpan(k, cuts[a], cuts[b], path.ranges(a, b))
                              for k, (a, b) in enumerate(zip(ordered, ordered[1:])))
                members.append(MemberExtent(sid, forward, order, order, True, spans))
        return BundleTopology(bundle['id'], tuple(members), tuple(events), (), tuple(sorted(features)),
                              tuple(tuple(sid for sid, _ in chain) for chain in chains), tuple(joins))

    def _main_chains(self,bundle,streets,origins):
        from .bundle_main_chains import normalize_main_chains
        prepared=normalize_main_chains(bundle,streets,origins,self.tolerance)
        if prepared is None:return None
        normalized,others=prepared
        result=self._chained(normalized,streets,origins)
        by_id={member.street_id:member for member in result.members}
        proxies=[];sources={}
        for role,chain in enumerate(result.chains):
            first=chain[0]
            proxies.append(MemberExtent(first,True,role,role,True,
                tuple(span for sid in chain for span in by_id[sid].spans)))
            sources[first]=dict(id=first,items=[item for sid in chain for item in streets[sid]['items']
                                              if item['kind']=='Section'])
        ownership=_BundleOwnership(tuple(proxies),sources,(),self.tolerance)
        if not ownership.rim:return None
        # The native grouping may tie its end ordering at a common node. A
        # simple complete control rim establishes the physical lateral order.
        if ownership.inside[proxies[0].street_id]=='LEFT':
            normalized=dict(normalized,head=normalized['head'][::-1],tail=normalized['tail'][::-1],
                head_locations=normalized['head_locations'][::-1],tail_locations=normalized['tail_locations'][::-1])
            result=self._chained(normalized,streets,origins)
        exterior=tuple(sid for sid in others if ownership.exterior(streets[sid]))
        interior=tuple(sid for sid in others if sid not in exterior)
        features={*result.required_features,'SOURCE_CARRIAGEWAY_CHAINS'}
        if normalized.get('source_extent_derived'):features.add('FULL_NAMED_SOURCE_EXTENT')
        if interior:features.add('INTERIOR_STREETS')
        return replace(result,interior_streets=interior,exterior_streets=exterior,
                       required_features=tuple(sorted(features)))

    def read(self, bundle, streets, *, origins=None):
        head, tail = tuple(bundle['head']), tuple(bundle['tail'])
        if (not head or len(set(head)) != len(head) or len(set(tail)) != len(tail)):
            raise BundleTopologyError('DUPLICATE_MEMBER', 'Bundle endpoint member lists must be nonempty and unique')
        if len(head) < 2:
            raise BundleTopologyError('INSUFFICIENT_MEMBERS', 'A Bundle needs at least two distinct member Streets')
        if len(bundle['head_locations']) != len(head) or len(bundle['tail_locations']) != len(tail):
            raise BundleTopologyError('ENDPOINT_COUNT', 'Each Bundle endpoint member needs its own prepared location')
        if set(head) != set(tail):
            if origins is not None:
                try:return self._chained(bundle, streets, origins)
                except BundleTopologyError:
                    repaired=self._main_chains(bundle,streets,origins)
                    if repaired is not None:return repaired
                    raise
            raise BundleTopologyError('UNRESOLVED_MEMBERSHIP', 'Members entering or leaving inside a Bundle need explicit incidence records')
        if head != tail:
            if origins is not None:
                from .bundle_source_chains import paired_source_chains,shared_terminal_order
                normalized=shared_terminal_order(bundle,streets,origins,self.tolerance)
                if normalized is not None:
                    result=self.read(normalized,streets,origins=origins)
                    return replace(result,required_features=tuple(sorted(
                        {*result.required_features,'SHARED_TERMINAL_ORDER'})))
                normalized = paired_source_chains(bundle, streets, origins, self.tolerance)
                if normalized is not None:
                    result = self._chained(normalized, streets, origins)
                    return replace(result, required_features=tuple(sorted(
                        {*result.required_features, 'SHARED_TERMINALS'})))
            raise BundleTopologyError('MEMBER_ORDER_CHANGE',
                'Prepared member order changes between ends; review ordering and the no-crossing requirement before constructing a surface')
        members, events, features, used_sections = [], [], set(), set()
        if len(head) != 2:
            features.add('MULTIPLE_MEMBERS')
        for order, sid in enumerate(head):
            if sid not in streets:
                raise BundleTopologyError('MISSING_MEMBER', f'Bundle references absent Street {sid}')
            path = _StreetPath(streets[sid], self.tolerance)
            if used_sections.intersection(path.by_id):
                raise BundleTopologyError('DUPLICATE_SECTION', 'Bundle members share a semantic Section ID')
            used_sections.update(path.by_id)
            start = path.locate(bundle['head_locations'][order])
            end = path.locate(bundle['tail_locations'][order])
            if start == end:
                raise BundleTopologyError('EMPTY_EXTENT', f'Bundle member Street {sid} has coincident extent ends')
            complete = min(start, end) == 0 and max(start, end) == path.length
            if not complete:
                features.add('PARTIAL_EXTENTS')
            junctions = [(station, event) for station, event in self._junctions(path, streets)
                         if min(start, end) <= station <= max(start, end)]
            events.extend(event for _station, event in junctions)
            if junctions:
                features.add('INTERNAL_JUNCTIONS')
            cuts = {station: SpanBoundary(junction_id=event.junction_id) for station, event in junctions}
            cuts[start] = SpanBoundary('HEAD', cuts.get(start, SpanBoundary()).junction_id)
            cuts[end] = SpanBoundary('TAIL', cuts.get(end, SpanBoundary()).junction_id)
            ordered = sorted(cuts, reverse=start > end)
            spans = tuple(MemberSpan(i, cuts[a], cuts[b], path.ranges(a, b))
                          for i, (a, b) in enumerate(zip(ordered, ordered[1:])))
            members.append(MemberExtent(sid, start < end, order, order, complete, spans))
        owned = tuple(sorted(set(bundle.get('owned_streets', ())) - set(head)))
        if any(sid not in streets for sid in owned):
            raise BundleTopologyError('MISSING_INTERIOR_STREET', 'Bundle owns an absent interior Street')
        ownership = _BundleOwnership(members, streets, events, self.tolerance)
        exterior = tuple(sid for sid in owned if ownership.exterior(streets[sid]))
        interior = tuple(sid for sid in owned if sid not in exterior)
        if interior:
            features.add('INTERIOR_STREETS')
        return BundleTopology(bundle['id'], tuple(members), tuple(events), interior, tuple(sorted(features)),
                              exterior_streets=exterior)

    def apply(self, snapshot):
        from .physical_sources import PhysicalStreetInventory
        sources=PhysicalStreetInventory.view(snapshot).streets
        for bundle in snapshot.bundles.values():
            try:
                bundle['topology'] = self.read(bundle, sources,
                    origins=getattr(snapshot, 'junction_origins', None)).report()
            except BundleTopologyError as exc:
                bundle['topology'] = dict(status='UNRESOLVED', code=exc.code, reason=str(exc))
                snapshot.issue('Bundle', bundle['id'], 'BUNDLE_TOPOLOGY_'+exc.code, str(exc))
