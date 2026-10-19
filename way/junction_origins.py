"""Retain physical connector incidences before Bundle junction aggregation.

Logical Bundle junctions can share a physical site. These immutable records
retain its original ring and authored endpoints without changing source IDs,
coordinates, dimensions, or the legacy linked lists. They are provenance and
ownership inputs, not a certificate that a surface can be generated.
"""
from collections import defaultdict
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class StreetEndpoint:
    street_id: int
    end: int
    terminal_id: int | None
    terminal_kind: str

    @property
    def key(self):
        return self.street_id, self.end


@dataclass(frozen=True)
class ConnectorOrigin:
    connector_id: int
    item_id: int
    item_kind: str
    leaving: bool
    endpoints: tuple[StreetEndpoint, ...]


@dataclass(frozen=True)
class JunctionOrigin:
    junction_id: int
    location: tuple
    ring: tuple[ConnectorOrigin, ...]


@dataclass(frozen=True)
class EndpointOrigin:
    logical_connector_id: int
    endpoint: StreetEndpoint
    physical_junction_id: int | None
    physical_connector_id: int | None
    status: str


@dataclass(frozen=True)
class FreeEndpointOrigin:
    endpoint: StreetEndpoint
    location: tuple


@dataclass(frozen=True)
class StreetSeamOrigin:
    """Two retained incidences of an original node, without invented connectors.

    Directed chain assembly cannot concatenate two starts or two tails without
    reversing a source. Their shared node still owns both ends. A seam is not
    a three-way junction ring and is not two independent free ends.
    """
    junction_id: int
    location: tuple
    endpoints: tuple[StreetEndpoint, StreetEndpoint]


@dataclass(frozen=True)
class LogicalJunctionOrigins:
    junction: JunctionOrigin
    endpoints: tuple[EndpointOrigin, ...]

    @property
    def physical_ids(self):
        return tuple(sorted({e.physical_junction_id for e in self.endpoints
                             if e.status == 'EXACT'}))


@dataclass(frozen=True)
class JunctionOrigins:
    physical: tuple[JunctionOrigin, ...]
    logical: tuple[LogicalJunctionOrigins, ...]
    free_endpoints: tuple[FreeEndpointOrigin, ...] = ()
    seams: tuple[StreetSeamOrigin, ...] = ()

    def components(self):
        """Connected logical records requiring joint physical-site ownership.

        Only recorded endpoint identity establishes a shared site. Equal or
        nearby coordinates never connect two otherwise unrelated junctions.
        Unresolved incidences remain in the component report and prevent a
        planner from claiming complete endpoint ownership.
        """
        by_id = {row.junction.junction_id: row for row in self.logical}
        by_site = defaultdict(set)
        for jid, row in by_id.items():
            for pid in row.physical_ids:
                by_site[pid].add(jid)
        remaining = set(by_id); result = []
        while remaining:
            pending = [min(remaining)]; members = set(); sites = set()
            while pending:
                jid = pending.pop()
                if jid in members:
                    continue
                members.add(jid); remaining.discard(jid)
                for pid in by_id[jid].physical_ids:
                    sites.add(pid); pending.extend(by_site[pid]-members)
            unresolved = tuple((jid, e) for jid in sorted(members)
                               for e in by_id[jid].endpoints if e.status != 'EXACT')
            result.append(dict(logical_junction_ids=tuple(sorted(members)),
                physical_junction_ids=tuple(sorted(sites)),
                shared_physical_junction_ids=tuple(sorted(p for p in sites if len(by_site[p]) > 1)),
                status='MAPPED' if not unresolved else 'UNRESOLVED',
                unresolved=[dict(logical_junction_id=jid, **asdict(e)) for jid, e in unresolved]))
        return tuple(result)

    def report(self):
        return dict(stage='createBundleIntersections', physical=[asdict(j) for j in self.physical],
                    logical=[asdict(j) for j in self.logical], components=self.components(),
                    free_endpoints=[asdict(e) for e in self.free_endpoints],
                    seams=[asdict(s) for s in self.seams])


class JunctionOriginRecorder:
    @staticmethod
    def retain_physical_sources(streets, physical, candidates=()):
        """Keep removed Bundle splitters still referenced by physical rings.

        Bundle grouping may remove an interior splitter from the logical Street
        dictionary before intersection aggregation. Its original connector ring
        still needs that exact source. Preserve the object in the physical view,
        without restoring it to the logical dictionary or inventing a free end.
        Later inventory review checks terminal and Section ownership as usual.
        """
        result = {}
        required = {e.street_id for j in physical for c in j.ring for e in c.endpoints}
        for sources, referenced_only in ((streets, False), (candidates, True)):
            for street in sources:
                if referenced_only and street.id not in required:
                    continue
                previous = result.get(street.id)
                if previous is not None and previous is not street:
                    raise ValueError('Physical Street ID refers to different source objects')
                result[street.id] = street
        return result

    @staticmethod
    def _endpoint(street, end):
        terminal = street.head if end == 0 else street.tail
        return StreetEndpoint(street.id, end, getattr(terminal, 'id', None), type(terminal).__name__)

    def _endpoints(self, connector):
        item = connector.item
        if not hasattr(item, 'streetsHead'):
            return (self._endpoint(item, int(not connector.leaving)),)
        streets = item.streetsHead if connector.leaving else item.streetsTail
        points = item.headLocs if connector.leaving else item.tailLocs
        if len(streets) != len(points):
            raise ValueError('Bundle endpoint locations must match its member inventory')
        result = []
        for street, point in zip(streets, points):
            ends = [end for end, value in enumerate((street.src, street.dst)) if value == point]
            if len(ends) != 1:
                # Partial or ambiguous membership needs its own retained range;
                # it must not be guessed from whichever endpoint is nearest.
                result.append(StreetEndpoint(street.id, -1, None, 'UNRESOLVED_BUNDLE_ENDPOINT'))
            else:
                result.append(self._endpoint(street, ends[0]))
        return tuple(result)

    def capture(self, intersections):
        result = []
        for junction in intersections:
            ring = []; seen = set()
            for c in junction:
                if id(c) in seen:
                    raise ValueError('Physical connector ring repeats before closing')
                seen.add(id(c))
                ring.append(ConnectorOrigin(c.id, c.item.id, type(c.item).__name__,
                    bool(c.leaving), self._endpoints(c)))
            result.append(JunctionOrigin(junction.id, tuple(junction.location), tuple(ring)))
        if len({j.junction_id for j in result}) != len(result):
            raise ValueError('Physical junction IDs must be unique')
        return tuple(result)

    def capture_seams(self, intersections, streets):
        """Use retained Section incidences, never proximity or stale backrefs.

        A LeavingWay may still refer to a Street replaced by directed assembly.
        Its Section and direction identify the actual terminal in the current
        source inventory. No source object, list, direction or ID is changed.
        """
        terminals = {}
        for street in streets:
            for end, terminal in enumerate((street.head, street.tail)):
                key = id(terminal), end
                if key in terminals:
                    raise ValueError('Original terminal Section incidence must have one Street owner')
                terminals[key] = street
        result = []; used = set(); junction_ids = set()
        for junction in intersections:
            ways = getattr(junction, 'leaveWays', ())
            if len(ways) != 2 or getattr(junction, 'isMinor', False) or any(True for _ in junction):
                continue
            if junction.id in junction_ids:
                raise ValueError('Original seam junction IDs must be unique')
            junction_ids.add(junction.id)
            endpoints = []
            for way in ways:
                if type(way.leaving) is not bool:
                    raise ValueError('Original seam incidence needs an authored direction')
                end = int(not way.leaving)
                street = terminals.get((id(way.section), end))
                if street is None:
                    raise ValueError('Original seam incidence has no current terminal Section owner')
                endpoint = self._endpoint(street, end)
                if endpoint.terminal_kind != 'Section':
                    raise ValueError('Original seam incidence must address a Section')
                if tuple(street.src if end == 0 else street.dst) != tuple(junction.location):
                    raise ValueError('Original seam terminal differs from its recorded node location')
                if (street.pred if end == 0 else street.succ) is not None:
                    raise ValueError('Original seam endpoint already has a connector owner')
                if endpoint.key in used:
                    raise ValueError('Original seam endpoint has duplicate ownership')
                used.add(endpoint.key); endpoints.append(endpoint)
            result.append(StreetSeamOrigin(junction.id, tuple(junction.location), tuple(endpoints)))
        return tuple(result)

    def capture_free_endpoints(self, streets, *, seams=()):
        """Record actual unattached ends before logical aggregation mutates links.

        A missing physical-ring incidence is not evidence of a free end. Only
        an explicit None connector on the original Street, with no retained
        degree-two seam incidence, establishes it.
        """
        attached = {e.key for seam in seams for e in seam.endpoints}
        result = []; seen = set()
        for street in streets:
            for end, connector, point in ((0, street.pred, street.src), (1, street.succ, street.dst)):
                endpoint = self._endpoint(street, end)
                if endpoint.key in seen:
                    raise ValueError('Original Street endpoint identities must be unique')
                seen.add(endpoint.key)
                if connector is None and endpoint.key not in attached:
                    result.append(FreeEndpointOrigin(endpoint, tuple(point)))
        if not attached.issubset(seen):
            raise ValueError('Original seam references an absent Street endpoint')
        return tuple(result)

    def resolve(self, physical, intersections, *, free_endpoints=(), seams=()):
        free = {r.endpoint.key: r for r in free_endpoints}
        if len(free) != len(free_endpoints):
            raise ValueError('Original free endpoint identities must be unique')
        owners = defaultdict(list)
        for junction in physical:
            for connector in junction.ring:
                for endpoint in connector.endpoints:
                    if endpoint.end in (0, 1):
                        owners[endpoint.key].append((junction, connector, endpoint))
        seam_ends = {}
        if len({s.junction_id for s in seams}) != len(seams):
            raise ValueError('Original seam junction IDs must be unique')
        for seam in seams:
            if len(seam.endpoints) != 2:
                raise ValueError('Original seam needs exactly two endpoint incidences')
            for endpoint in seam.endpoints:
                if endpoint.key in seam_ends or endpoint.key in owners or endpoint.key in free:
                    raise ValueError('Original seam endpoint has duplicate or conflicting ownership')
                seam_ends[endpoint.key] = (seam, endpoint)
        logical = []
        for junction in self.capture(intersections):
            if not any(c.item_kind == 'Bundle' for c in junction.ring):
                continue
            endpoints = []
            for connector in junction.ring:
                for endpoint in connector.endpoints:
                    candidates = owners.get(endpoint.key, ())
                    status = ('EXACT' if len(candidates) == 1 else
                              'MISSING_ORIGIN' if not candidates else 'AMBIGUOUS_ORIGIN')
                    owner, original, old_endpoint = candidates[0] if status == 'EXACT' else (None, None, None)
                    if old_endpoint is not None and (old_endpoint.terminal_id, old_endpoint.terminal_kind) != (
                            endpoint.terminal_id, endpoint.terminal_kind):
                        status = 'CHANGED_TERMINAL'
                    original_free = free.get(endpoint.key)
                    if original_free is not None:
                        if original_free.endpoint != endpoint:
                            status = 'CHANGED_TERMINAL'
                        elif candidates:
                            status = 'CONFLICTING_FREE_ORIGIN'
                        else:
                            # Known source boundary, not a physical junction.
                            # A logical patch must support this explicitly; it
                            # must not become MAPPED merely because it is known.
                            status = 'ORIGINAL_FREE'
                    original_seam = seam_ends.get(endpoint.key)
                    if original_seam is not None:
                        status = 'ORIGINAL_SEAM' if original_seam[1] == endpoint else 'CHANGED_TERMINAL'
                    endpoints.append(EndpointOrigin(connector.connector_id, endpoint,
                        owner.junction_id if owner else original_seam[0].junction_id if original_seam else None,
                        original.connector_id if original else None, status))
            logical.append(LogicalJunctionOrigins(junction, tuple(endpoints)))
        return JunctionOrigins(tuple(physical), tuple(logical), tuple(free_endpoints), tuple(seams))
