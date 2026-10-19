"""Original Section-end ports and once-owned spans for physical surfaces.

A Street remains a complete editable source. Junctions inside it reference
Section ends; spans partition its Sections without allocating replacement
Street/Section IDs. This module prepares topology only, never surface geometry
or sampled curves. Geometry edits retain these references; managed topology
edits must rebuild the affected plan.
"""
from dataclasses import dataclass
import math


class GroundSurfacePolicy:
    """Scope of the current ground carriageway provider, without changing tags.

    A building passage can use the ordinary pavement construction only when
    it has no indicated vertical displacement. Buildings and their openings
    remain the responsibility of their own feature provider.
    """

    @staticmethod
    def reason(tags, *, building_passages=True):
        unsupported = 'Grade/elevation treatment is not yet implemented'
        if tags.get('bridge') not in (None, 'no') or str(tags.get('layer', '0')) != '0':
            return unsupported
        tunnel = tags.get('tunnel')
        if tunnel in (None, 'no'):
            return None
        if (building_passages and tunnel == 'building_passage' and str(tags.get('level', '0')) == '0'
                and tags.get('indoor') in (None, 'no')
                and tags.get('location') in (None, 'surface')):
            return None
        return unsupported


def _identifier(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(label+' must be a nonnegative integer')
    return value


@dataclass(frozen=True, order=True)
class SectionEnd:
    street_id: int
    section_id: int
    end: int

    def __post_init__(self):
        _identifier(self.street_id, 'Street ID')
        _identifier(self.section_id, 'Section ID')
        if type(self.end) is not int or self.end not in (0, 1):
            raise ValueError('A Section end must be the integer zero or one')


@dataclass(frozen=True)
class PhysicalPort:
    junction_id: int
    endpoint: SectionEnd
    source_kind: str
    connector_id: int | None


@dataclass(frozen=True)
class PhysicalPortSite:
    junction_id: int
    location: tuple
    ports: tuple[PhysicalPort, ...]


@dataclass(frozen=True)
class SpanEnd:
    endpoint: SectionEnd
    owner_kind: str
    junction_id: int | None


@dataclass(frozen=True)
class SelectionBoundary:
    """A declared metric trim on an original Section end, not a free-end alias."""
    endpoint: SectionEnd
    trim_m: float

    def __post_init__(self):
        if type(self.trim_m) not in (int,float) or not math.isfinite(self.trim_m) or self.trim_m <= 0:
            raise ValueError('Selection boundary trim must be positive and finite')


@dataclass(frozen=True)
class PhysicalSpan:
    street_id: int
    section_ids: tuple[int, ...]
    start: SpanEnd
    end: SpanEnd
    # Inclusive private knot indices in the complete authored Street. These
    # are bindings, not stable IDs; managed topology edits rebuild them.
    first_control: int
    last_control: int
    source_control_count: int
    in_scope: bool = True

    @property
    def identity(self):
        return self.street_id, self.section_ids[0], self.section_ids[-1]


@dataclass(frozen=True)
class DeferredPort:
    junction_id: int
    endpoint: SectionEnd
    connector_id: int
    category: str


@dataclass(frozen=True)
class CategoryBoundary:
    street_id: int
    before_section: int
    after_section: int
    source_kind: str
    source_id: int | None

    @property
    def endpoints(self):
        return (SectionEnd(self.street_id, self.before_section, 1),
                SectionEnd(self.street_id, self.after_section, 0))


@dataclass(frozen=True)
class CategoryConnection:
    """A feature boundary at an original two-Street node, not a road patch.

    The excluded feature remains in the source inventory. Its identity is a
    handoff declaration, not permission to construct its pavement or to use its
    tangent to reshape the in-scope road.
    """
    source_id: int
    endpoints: tuple[SectionEnd, SectionEnd]
    categories: tuple[str, str]
    source_kind: str = 'OriginalStreetSeam'


@dataclass(frozen=True)
class StreetSectionScope:
    """Feature ownership is independent of the complete editable Street."""
    street_id: int
    categories: tuple[str, ...]
    section_ids: tuple[int, ...]
    included_sections: tuple[int, ...]
    boundaries: tuple[CategoryBoundary, ...]

    @property
    def curve_groups(self):
        """Maximal consecutive feature runs, independent of junction spans.

        A carriageway must not be bent or rejected by an excluded footway,
        while consecutive road Sections still share one original cubic.
        """
        result=[];group=0;previous=None
        for section in self.section_ids:
            included=section in self.included_sections
            if previous is not None and included!=previous:group+=1
            result.append(group);previous=included
        return tuple(result)

    @classmethod
    def prepare(cls, source, categories):
        categories = tuple(sorted(categories))
        items = source['items']
        sections = [s for s in items if s['kind'] == 'Section']
        included = tuple(s['id'] for s in sections if s['category'] in categories)
        if not sections or not included or len(included) == len(sections):
            raise ValueError('Section scope needs both included and excluded original Sections')
        if source.get('status') != 'SUPPORTED' and not (source.get('reason') or '').startswith('Reserved component/category:'):
            raise ValueError('Section scope cannot bypass another unsupported source treatment')
        # A category reason can mask a later invalid road or Street item in the
        # legacy first-failure report. Recheck the entire source deliberately.
        for item in items:
            if item['kind'] != 'Section':
                if (item.get('transition_reason') or item.get('points') or
                        item.get('surface_treatment') not in
                        ('DEFERRED_ATTACHMENT', 'NATIVE_JUNCTION', 'COMPATIBLE_SEAM', 'WIDTH_TRANSITION')):
                    raise ValueError('Section scope has an unsupported original Street item')
                continue
            if not isinstance(item['category'], str) or not item['category']:
                raise ValueError('Section scope needs an explicit category')
            if item['id'] in included:
                reason = GroundSurfacePolicy.reason(item.get('tags', {}))
                if reason:
                    raise ValueError(reason)
            if (len(item['points']) < 2 or item['points'][0] == item['points'][-1]
                    or not math.isfinite(float(item['width'])) or item['width'] <= 0
                    or not math.isfinite(float(item.get('offset') or 0.))):
                raise ValueError('Section scope needs open controls and valid physical dimensions')
        if sections[0]['points'][0] == sections[-1]['points'][-1]:
            raise ValueError('Closed Street requires a dedicated path treatment')
        boundaries = []
        positions = [i for i, item in enumerate(items) if item['kind'] == 'Section']
        for a, b, ia, ib in zip(sections, sections[1:], positions, positions[1:]):
            if (a['id'] in included) == (b['id'] in included):
                continue
            between = items[ia+1:ib]
            if len(between) > 1 or any(s.get('surface_treatment') == 'NATIVE_JUNCTION' for s in between):
                raise ValueError('Category boundary with junction ownership needs a mixed-mode junction contract')
            kind, identity = 'SectionSeam', None
            if between:
                seam = between[0]; kind, identity = seam['kind'], seam['id']
                _identifier(identity, 'Category boundary source ID')
                if kind not in ('SymLane', 'SideLane'):
                    raise ValueError('Category boundary attachment needs an explicit shared feature contract')
                record = seam.get('transition')
                if seam.get('surface_treatment') == 'WIDTH_TRANSITION' and not record:
                    raise ValueError('Category boundary lost its prepared transition declaration')
                if record and (record.get('before_section'), record.get('after_section'),
                               record.get('source_kind'), record.get('source_id')) != (a['id'], b['id'], kind, identity):
                    raise ValueError('Category boundary lost its prepared transition identity')
                if record:
                    half, slope = float(record['minimum_half_length']), float(record['maximum_edge_slope'])
                    if not math.isfinite(half) or half <= 0 or not math.isfinite(slope) or not 0 < slope <= 1:
                        raise ValueError('Category boundary has invalid prepared transition settings')
            boundaries.append(CategoryBoundary(source['id'], a['id'], b['id'], kind, identity))
        return cls(source['id'], categories, tuple(s['id'] for s in sections), included, tuple(boundaries))


@dataclass(frozen=True)
class PhysicalPortPlan:
    sites: tuple[PhysicalPortSite, ...]
    spans: tuple[PhysicalSpan, ...]
    source_streets: tuple[int, ...]
    deferred: tuple[DeferredPort, ...]
    section_scopes: tuple[StreetSectionScope, ...] = ()
    category_connections: tuple[CategoryConnection, ...] = ()
    selection_boundaries: tuple[SelectionBoundary, ...] = ()

    def selection_boundary(self, endpoint):
        rows = [b for b in self.selection_boundaries if b.endpoint == endpoint]
        if len(rows) != 1:
            raise ValueError('Selection boundary needs one exact original endpoint declaration')
        return rows[0]

    def selection_peer(self, endpoint):
        """Return the other side of a declared interior cut, or None at a terminal.

        References use original span/control ownership, never point proximity.
        Both sides remain part of the same complete editable Street.
        """
        self.selection_boundary(endpoint)
        own = [s for s in self.spans if endpoint in (s.start.endpoint, s.end.endpoint)]
        if len(own) != 1:
            raise ValueError('Selection boundary requires exactly one original span')
        span, = own
        knot = span.first_control if endpoint.end == 0 else span.last_control
        if knot == (0 if endpoint.end == 0 else span.source_control_count-1):
            return None
        peers = [s.end if endpoint.end == 0 else s.start for s in self.spans
                 if s.street_id == endpoint.street_id and s != span
                 and (s.last_control if endpoint.end == 0 else s.first_control) == knot]
        if len(peers) != 1 or peers[0].owner_kind != 'SELECTION_BOUNDARY':
            raise ValueError('Interior selection boundary requires its paired span owner')
        peer = peers[0].endpoint
        self.selection_boundary(peer)
        return peer

    def source_supported(self, source):
        scopes = [s for s in self.section_scopes if s.street_id == source['id']]
        if scopes:
            return len(scopes) == 1 and StreetSectionScope.prepare(source, scopes[0].categories) == scopes[0]
        return source.get('status') == 'SUPPORTED'

    def category_boundary(self, endpoint):
        found = [b for s in self.section_scopes for b in s.boundaries if endpoint in b.endpoints]
        found.extend(b for b in self.category_connections if endpoint in b.endpoints)
        if len(found) != 1:
            raise ValueError('Category endpoint needs exactly one original Section seam owner')
        return found[0]

    def span_for(self, port):
        found = [span for span in self.spans
                 if port.endpoint in (span.start.endpoint, span.end.endpoint)]
        if len(found) != 1:
            raise ValueError('A port needs exactly one original source span')
        return found[0]


class PreparedPhysicalPorts:
    def __init__(self, categories, *, tolerance=1e-5, defer_physical_categories=True, section_scopes=False,
                 seam_connections=False, category_connections=False):
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('Port-address tolerance must be positive and finite')
        self.categories = frozenset(categories)
        if not self.categories:
            raise ValueError('Declare the categories owned by this surface provider')
        self.tolerance = tolerance
        self.defer_physical_categories = defer_physical_categories
        self.section_scopes = section_scopes
        self.seam_connections = seam_connections
        self.category_connections = category_connections
        if category_connections and not (section_scopes and seam_connections):
            raise ValueError('Cross-Street category boundaries need Section scopes and original seam support')

    def category_seam_ids(self, snapshot):
        """Classify seeds only; prepare still validates complete provenance.

        A feature boundary has no carriageway junction patch. Malformed or
        missing references are not silently admitted by this classification.
        """
        if not self.category_connections:return frozenset()
        result=set()
        for seam in snapshot.junction_origins.get('seams', ()):
            categories=[]
            for r in seam['endpoints']:
                source=snapshot.streets.get(r['street_id'], {})
                items=[s for s in source.get('items', ()) if s['kind']=='Section'
                       and s['id']==r['terminal_id'] and r['terminal_kind']=='Section']
                if len(items)!=1:break
                category=items[0].get('category')
                if not isinstance(category,str) or not category:break
                categories.append(category in self.categories)
            if len(categories)==2 and set(categories)=={False,True}:
                result.add(seam['junction_id'])
        return frozenset(result)

    @staticmethod
    def _selection(values, label):
        values = tuple(values)
        if any(type(v) is not int or v < 0 for v in values) or len(set(values)) != len(values):
            raise ValueError(label+' IDs must be distinct nonnegative integers')
        return values

    def prepare(self, snapshot, *, physical_ids=(), interior_ids=(), street_ids=(), seam_ids=()):
        physical_ids = self._selection(physical_ids, 'Physical junction')
        interior_ids = self._selection(interior_ids, 'Interior junction')
        street_ids = self._selection(street_ids, 'Original Street')
        seam_ids = self._selection(seam_ids, 'Original seam')
        if seam_ids and not self.seam_connections:
            raise ValueError('Original two-ended Street seam requires a native connection provider')
        if not physical_ids and not interior_ids and not street_ids and not seam_ids:
            raise ValueError('Select at least one junction or original Street')
        if set(physical_ids) & set(interior_ids):
            raise ValueError('A junction cannot have two source representations')
        if set(seam_ids) & (set(physical_ids) | set(interior_ids)):
            raise ValueError('A seam cannot also have an ordinary junction representation')
        physical = {}
        for site in snapshot.junction_origins.get('physical', ()):
            physical.setdefault(site['junction_id'], []).append(site)
        interiors = {}
        for sid, street in snapshot.streets.items():
            for item in street['items']:
                if item.get('surface_treatment') == 'NATIVE_JUNCTION':
                    interiors.setdefault(item['id'], []).append((sid, item))
        cache = {}; scopes = {}

        def street_sections(sid, *, require_supported=True):
            source = snapshot.streets.get(sid)
            if source is None:
                raise ValueError(f'Port source Street {sid} is missing from the original inventory')
            if require_supported and source.get('status') != 'SUPPORTED':
                flags = {s['category'] in self.categories for s in source['items'] if s['kind'] == 'Section'}
                if self.section_scopes and flags == {False, True}:
                    scopes[sid] = StreetSectionScope.prepare(source, self.categories)
                else:
                    reason = source.get('reason') or 'No supported source treatment is declared'
                    raise ValueError(f'Port source Street {sid} is unsupported: {reason}')
            if require_supported and self.section_scopes and source.get('status') == 'SUPPORTED':
                flags = {s['category'] in self.categories for s in source['items'] if s['kind'] == 'Section'}
                if flags == {False, True}:scopes[sid] = StreetSectionScope.prepare(source, self.categories)
            if sid in cache:
                return cache[sid]
            if _identifier(source['id'], 'Street ID') != sid:
                raise ValueError('Source Street ID differs from its inventory key')
            sections = tuple(i for i in source['items'] if i['kind'] == 'Section')
            if not sections or len({s['id'] for s in sections}) != len(sections):
                raise ValueError(f'Street {sid} has missing or duplicate Sections')
            previous = None
            for section in sections:
                _identifier(section['id'], 'Section ID')
                points = section['points']
                if len(points) < 2 or any(len(p) < 2 or not all(math.isfinite(v) for v in p) for p in points):
                    raise ValueError(f'Section {section["id"]} needs finite original controls')
                if previous is not None and math.dist(previous[:2], points[0][:2]) > self.tolerance:
                    raise ValueError(f'Street {sid} has disconnected original Sections')
                previous = points[-1]
            cache[sid] = sections
            return sections

        def section_at(endpoint, *, terminal=False, require_supported=True):
            if endpoint.end not in (0, 1):
                raise ValueError('A Section end must be zero or one')
            sections = street_sections(endpoint.street_id, require_supported=require_supported)
            found = [s for s in sections if s['id'] == endpoint.section_id]
            if len(found) != 1:
                raise ValueError('A port references an absent original Section')
            if terminal and sections[0 if endpoint.end == 0 else -1]['id'] != endpoint.section_id:
                raise ValueError('An attached endpoint is not its original terminal Section')
            return found[0]

        def at_location(endpoint, location, *, terminal=False, allow_deferred=False):
            if len(location) < 2 or not all(math.isfinite(v) for v in location):
                raise ValueError('A junction needs a finite location')
            section = section_at(endpoint, terminal=terminal, require_supported=not allow_deferred)
            if not isinstance(section['category'], str) or not section['category']:
                raise ValueError('A port needs an explicit source category')
            if section['category'] not in self.categories:
                if not allow_deferred:
                    raise ValueError('Selected port category is outside this provider scope')
            elif allow_deferred:
                # Only category exclusions may bypass surface support. A road
                # rejected for grade, dimensions or source topology still fails.
                street_sections(endpoint.street_id)
            point = section['points'][0 if endpoint.end == 0 else -1]
            if math.dist(point[:2], location[:2]) > self.tolerance:
                raise ValueError('Section endpoint differs from its recorded junction location')
            return section

        def through(sid, item):
            sections = street_sections(sid)
            refs = tuple(tuple(r) for r in item.get('through', ()))
            if len(refs) != 2 or len(set(refs)) != 2:
                raise ValueError('An interior junction needs two distinct through Section ends')
            for ref in refs:
                if len(ref) != 2:
                    raise ValueError('A through reference needs a Section ID and end')
                SectionEnd(sid, *ref)
            pairs = [i for i in range(1, len(sections))
                     if set(refs) == {(sections[i-1]['id'], 1), (sections[i]['id'], 0)}]
            if len(pairs) != 1:
                raise ValueError('Through ports must join consecutive original Sections')
            split = pairs[0]
            ends = (SectionEnd(sid, sections[split-1]['id'], 1), SectionEnd(sid, sections[split]['id'], 0))
            for endpoint in ends:
                at_location(endpoint, item['location'])
            # Imported editable inputs may already assign Connector identities
            # to both through arms. Preserve those IDs; ordinary OSM preparation
            # still omits them when it never created through Connectors.
            connectors = {}
            if 'through_connectors' in item:
                for record in item['through_connectors']:
                    endpoint = SectionEnd(sid, record['section_id'], record['end'])
                    if endpoint not in ends or endpoint in connectors:
                        raise ValueError('Through Connectors must address each original through end once')
                    connectors[endpoint] = _identifier(record['connector_id'], 'Through Connector ID')
                if set(connectors) != set(ends) or len(set(connectors.values())) != 2:
                    raise ValueError('Through Connectors need two distinct IDs on the original through ends')
            return split, ends, tuple(connectors.get(end) for end in ends)

        sites, deferred, used = [], [], set()
        def add(site, *, seam=False):
            if (len(site.ports) != 2 if seam else len(site.ports) < 3):
                raise ValueError('A generated junction needs at least three in-scope ports')
            if seam and any(p.source_kind != 'SEAM' or p.connector_id is not None for p in site.ports):
                raise ValueError('A two-ended seam retains original endpoints without invented Connectors')
            connectors = [p.connector_id for p in site.ports if p.connector_id is not None]
            if len(set(connectors)) != len(connectors):
                raise ValueError('A selected junction repeats an original Connector ID')
            for port in site.ports:
                if port.connector_id is not None:
                    _identifier(port.connector_id, 'Connector ID')
                if port.endpoint in used:
                    raise ValueError('An original Section end has multiple surface owners')
                used.add(port.endpoint)
            sites.append(site)

        for jid in physical_ids:
            matches = physical.get(jid, ())
            if len(matches) != 1:
                raise ValueError('Physical junction needs one original ring record')
            site = matches[0]
            ports = []
            ring_connectors, ring_endpoints = set(), set()
            for connector in site['ring']:
                if connector['item_kind'] != 'Street' or len(connector['endpoints']) != 1:
                    raise ValueError('Physical ring must retain individual original incidences')
                record = connector['endpoints'][0]
                if record['terminal_kind'] != 'Section':
                    raise ValueError('A physical port must address an original Section')
                endpoint = SectionEnd(record['street_id'], record['terminal_id'], record['end'])
                cid = _identifier(connector['connector_id'], 'Connector ID')
                if cid in ring_connectors:
                    raise ValueError('A selected junction repeats an original Connector ID')
                if endpoint in ring_endpoints:
                    raise ValueError('A physical ring repeats an original Section end')
                ring_connectors.add(cid)
                ring_endpoints.add(endpoint)
                section = at_location(endpoint, site['location'], terminal=True,
                                      allow_deferred=self.defer_physical_categories)
                if section['category'] not in self.categories:
                    # Preserve excluded incidences separately and retain the
                    # original cyclic order of the carriageway subsequence.
                    # The source itself is neither rewritten nor generated.
                    deferred.append(DeferredPort(jid, endpoint, cid, section['category']))
                    continue
                ports.append(PhysicalPort(jid, endpoint, 'PHYSICAL', cid))
            add(PhysicalPortSite(jid, tuple(site['location']), tuple(ports)))

        for jid in interior_ids:
            matches = interiors.get(jid, ())
            if len(matches) != 1:
                raise ValueError('Interior junction needs one original Street item')
            sid, item = matches[0]
            if item.get('attachment_error'):
                raise ValueError(item['attachment_error'])
            _split, ends, connectors = through(sid, item)
            ports = [PhysicalPort(jid, end, 'THROUGH', cid) for end, cid in zip(ends, connectors)]
            for a in item.get('attachments', ()):
                endpoint = SectionEnd(a['street_id'], a['section_id'], a['end'])
                _identifier(a['connector_id'], 'Connector ID')
                if a['category'] not in self.categories:
                    deferred.append(DeferredPort(jid, endpoint, a['connector_id'], a['category']))
                    continue
                at_location(endpoint, item['location'], terminal=True)
                ports.append(PhysicalPort(jid, endpoint, 'ATTACHMENT', a['connector_id']))
            # Infer cyclic order once, as in the ordinary minor-junction adapter.
            # Endpoint identities remain independent of geometry and row indices.
            def angle(port):
                points = section_at(port.endpoint)['points']
                a, b = (points[0], points[1]) if port.endpoint.end == 0 else (points[-1], points[-2])
                if math.dist(a[:2], b[:2]) <= self.tolerance:
                    raise ValueError('A port has a degenerate original direction')
                return math.atan2(b[1]-a[1], b[0]-a[0]) % math.tau
            arms = sorted((angle(p), p.endpoint, p) for p in ports)
            angles = [a[0] for a in arms]
            if any((b-a) % math.tau < 1e-7 for a, b in zip(angles, angles[1:]+angles[:1])):
                raise ValueError('Interior ports have coincident directions')
            ports = [a[2] for a in arms]
            first = min(range(len(ports)), key=lambda i: ports[i].endpoint)
            ports = ports[first:]+ports[:first]
            add(PhysicalPortSite(jid, tuple(item['location']), tuple(ports)))

        for jid in seam_ids:
            records = [s for s in snapshot.junction_origins.get('seams', ()) if s['junction_id'] == jid]
            if len(records) != 1 or len(records[0]['endpoints']) != 2:
                raise ValueError('Selected seam needs one original two-incidence record')
            seam = records[0];ports = []
            for r in seam['endpoints']:
                if r['terminal_kind'] != 'Section':
                    raise ValueError('Selected seam must address original terminal Sections')
                endpoint = SectionEnd(r['street_id'], r['terminal_id'], r['end'])
                at_location(endpoint, seam['location'], terminal=True)
                ports.append(PhysicalPort(jid, endpoint, 'SEAM', None))
            add(PhysicalPortSite(jid, tuple(seam['location']), tuple(ports)), seam=True)

        source_ids = tuple(sorted({p.endpoint.street_id for site in sites for p in site.ports} | set(street_ids)))
        # An explicitly selected whole Street can begin with an excluded
        # Section. Establish its validated category scope before examining
        # original endpoint records; evaluating allow_deferred before this
        # lazy lookup made admission depend on the authored Street direction.
        for sid in source_ids:
            street_sections(sid)
        # Every original Section is represented once, including spans outside the
        # selection. Their missing junction owners cannot become fake free ends.
        owners = {}
        def owner(endpoint, kind, jid):
            _identifier(jid, 'Junction ID')
            value = kind, jid
            if endpoint in owners:
                raise ValueError('An original Section end has duplicate or conflicting junction records')
            owners[endpoint] = value
        for records in physical.values():
            for site in records:
                for connector in site['ring']:
                    for r in connector['endpoints']:
                        if r['street_id'] in source_ids and r['terminal_kind'] == 'Section':
                            endpoint = SectionEnd(r['street_id'], r['terminal_id'], r['end'])
                            at_location(endpoint, site['location'], terminal=True,
                                        allow_deferred=r['street_id'] in scopes)
                            owner(endpoint, 'PHYSICAL', site['junction_id'])
        recorded_seam_ids = set(); category_connections = []
        for seam in snapshot.junction_origins.get('seams', ()):
            jid = _identifier(seam['junction_id'], 'Original seam junction ID')
            if jid in recorded_seam_ids:
                raise ValueError('Original seam junction IDs must be unique')
            recorded_seam_ids.add(jid)
            records = seam['endpoints']
            if len(records) != 2:
                raise ValueError('Original seam needs exactly two endpoint incidences')
            if not any(r['street_id'] in source_ids for r in records):
                continue
            original = physical.get(jid, ())
            if len(original) != 1 or original[0]['ring'] or tuple(original[0]['location']) != tuple(seam['location']):
                raise ValueError('Original seam needs its retained connector-free physical node')
            references = []; categories = []
            for r in records:
                if r['terminal_kind'] != 'Section':
                    raise ValueError('Original seam must retain its terminal Sections')
                endpoint = SectionEnd(r['street_id'], r['terminal_id'], r['end'])
                # Validate the partner's identity even if its category or grade
                # is deferred. Such a partner must not make this end FREE.
                section = section_at(endpoint, terminal=True, require_supported=False)
                if not isinstance(section['category'], str) or not section['category']:
                    raise ValueError('Original seam needs an explicit source category')
                categories.append(section['category'])
                point = section['points'][0 if endpoint.end == 0 else -1]
                location = seam['location']
                if (len(location) < 2 or not all(math.isfinite(v) for v in location)
                        or math.dist(point[:2], location[:2]) > self.tolerance):
                    raise ValueError('Original seam terminal differs from its recorded node location')
                if endpoint in references:
                    raise ValueError('Original seam repeats an endpoint incidence')
                references.append(endpoint)
                if endpoint.street_id in source_ids:
                    owner(endpoint, 'SEAM', jid)
            if self.category_connections and {c in self.categories for c in categories} == {False, True}:
                if references[0].street_id == references[1].street_id:
                    raise ValueError('Cross-Street category boundary needs two distinct original Streets')
                category_connections.append(CategoryConnection(jid, tuple(references), tuple(categories)))
                for endpoint in references:
                    if endpoint.street_id in source_ids:
                        # The original node remains in the declaration. There
                        # is no road junction surface to traverse or generate.
                        owners[endpoint] = ('CATEGORY_BOUNDARY', None)
        # Logical Bundle aggregation can attach a member's original free end
        # to the Bundle connector at another member's junction. Keep the source
        # boundary from the immutable pre-aggregation record, never from a gap
        # in the physical-ring inventory or proximity to a junction.
        for record in snapshot.junction_origins.get('free_endpoints', ()):
            r = record['endpoint']
            if r['street_id'] not in source_ids:
                continue
            if r['terminal_kind'] != 'Section':
                raise ValueError('An original free endpoint must retain its terminal Section')
            endpoint = SectionEnd(r['street_id'], r['terminal_id'], r['end'])
            at_location(endpoint, record['location'], terminal=True,
                        allow_deferred=r['street_id'] in scopes)
            if endpoint in owners:
                raise ValueError('An original free endpoint has duplicate or conflicting ownership')
            owners[endpoint] = ('FREE', None)
        selection_boundaries = []
        for record in snapshot.junction_origins.get('selection_boundaries', ()):
            endpoint = SectionEnd(**record['endpoint'])
            if endpoint.street_id not in source_ids:
                continue
            section_at(endpoint)
            boundary = SelectionBoundary(endpoint, record['trim_m'])
            source = snapshot.streets[endpoint.street_id]
            sections = street_sections(endpoint.street_id)
            terminal = sections[0 if endpoint.end == 0 else -1]['id'] == endpoint.section_id
            if endpoint in owners or (terminal and source['endpoints'][endpoint.end] is not None):
                raise ValueError('Selection boundary conflicts with an original endpoint owner')
            owners[endpoint] = ('SELECTION_BOUNDARY', None)
            selection_boundaries.append(boundary)
        selection_splits = {}
        declared = {b.endpoint for b in selection_boundaries}
        for endpoint in declared:
            sections = street_sections(endpoint.street_id)
            index = next(i for i, s in enumerate(sections) if s['id'] == endpoint.section_id)
            split = index if endpoint.end == 0 else index+1
            if split in (0, len(sections)):
                continue
            left, right = sections[split-1:split+1]
            pair = {SectionEnd(endpoint.street_id, left['id'], 1),
                    SectionEnd(endpoint.street_id, right['id'], 0)}
            if not pair <= declared:
                raise ValueError('Interior selection cut requires paired declarations on consecutive Sections')
            items = snapshot.streets[endpoint.street_id]['items']
            if items.index(right) != items.index(left)+1:
                raise ValueError('Interior selection cut conflicts with a Street seam item')
            selection_splits.setdefault(endpoint.street_id, set()).add(split)
        cuts = {}
        for sid in source_ids:
            sections = street_sections(sid)
            cuts[sid] = {0, len(sections)} | selection_splits.get(sid, set())
            boundaries = scopes[sid].boundaries if sid in scopes else ()
            for boundary in boundaries:
                split = next(i for i, s in enumerate(sections) if s['id'] == boundary.after_section)
                if split in cuts[sid]:
                    raise ValueError('Category seam conflicts with an interior selection cut')
                cuts[sid].add(split)
                for endpoint in boundary.endpoints:
                    if endpoint in owners:
                        raise ValueError('Category seam has conflicting junction ownership')
                    owners[endpoint] = ('CATEGORY_BOUNDARY', None)
            for item in snapshot.streets[sid]['items']:
                if item.get('surface_treatment') != 'NATIVE_JUNCTION':
                    continue
                split, ends, _connectors = through(sid, item)
                if split in cuts[sid]:
                    raise ValueError('Several junctions occupy one original Section seam')
                cuts[sid].add(split)
                for end in ends:
                    owner(end, 'INTERIOR', item['id'])
        for records in interiors.values():
            for _sid, item in records:
                for a in item.get('attachments', ()):
                    if a['street_id'] in source_ids and a['category'] in self.categories:
                        endpoint = SectionEnd(a['street_id'], a['section_id'], a['end'])
                        at_location(endpoint, item['location'], terminal=True)
                        owner(endpoint, 'INTERIOR', item['id'])
        spans = []
        for sid in source_ids:
            sections = street_sections(sid)
            source = snapshot.streets[sid]
            controls = [0]
            for section in sections:
                controls.append(controls[-1]+len(section['points'])-1)
            for a, b in zip(sorted(cuts[sid]), sorted(cuts[sid])[1:]):
                first, last = SectionEnd(sid, sections[a]['id'], 0), SectionEnd(sid, sections[b-1]['id'], 1)
                ends = []
                for endpoint in (first, last):
                    value = owners.get(endpoint)
                    if value is None:
                        outer = endpoint == SectionEnd(sid, sections[0]['id'], 0) or endpoint == SectionEnd(sid, sections[-1]['id'], 1)
                        if not outer or len(source.get('endpoints', ())) != 2:
                            raise ValueError('An unresolved span boundary needs an explicit owner')
                        jid = source['endpoints'][endpoint.end]
                        if jid is not None:
                            _identifier(jid, 'Remote junction ID')
                        value = ('FREE', None) if jid is None else ('REMOTE_LOGICAL', jid)
                    ends.append(SpanEnd(endpoint, *value))
                included = all(s['category'] in self.categories for s in sections[a:b])
                spans.append(PhysicalSpan(sid, tuple(s['id'] for s in sections[a:b]), *ends,
                                          controls[a], controls[b], controls[-1]+1, included))
        plan = PhysicalPortPlan(tuple(sites), tuple(spans), source_ids, tuple(deferred),
                                tuple(scopes[sid] for sid in source_ids if sid in scopes), tuple(category_connections),
                                tuple(selection_boundaries))
        for site in plan.sites:
            for port in site.ports:
                span = plan.span_for(port)
                end = span.start if span.start.endpoint == port.endpoint else span.end
                if end.junction_id != site.junction_id or end.owner_kind not in ('PHYSICAL', 'INTERIOR', 'SEAM'):
                    raise ValueError('Selected port disagrees with its original span owner')
                if (end.owner_kind == 'SEAM') != (port.source_kind == 'SEAM'):
                    raise ValueError('Selected seam incidence differs from its original owner kind')
        all_sections = [sid for span in spans for sid in span.section_ids]
        if len(set(all_sections)) != len(all_sections):
            raise ValueError('An original Section has multiple source-span owners')
        return plan
