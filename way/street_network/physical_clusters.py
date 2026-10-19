"""Join recorded physical sites into one explicit compound boundary contract.

This opt-in construction contracts a tree of original internal Streets in the
recorded ribbon graph. It never sorts endpoints by angle or substitutes a line
for a curved internal Street. Complete internal geometry is returned separately
and must pass native footprint coverage before any surface can be published.
"""
from dataclasses import dataclass, replace


@dataclass(frozen=True)
class PhysicalClusterInputs:
    owner_id: int
    logical_junction_ids: tuple
    physical_junction_ids: tuple
    junction: object
    internal_streets: tuple
    # Street ID and two (physical site, connector) incidences, authored order.
    internal_incidents: tuple
    # Every original external incidence, in contracted circular order.
    external_ring: tuple
    # Logical connector ID -> its original physical endpoint incidences.
    connector_origins: tuple


class PreparedPhysicalClusters:
    def __init__(self, physical_reader):
        self.reader = physical_reader

    def prepare(self, snapshot, logical_junction_id):
        original = self.reader.prepare(snapshot, logical_junction_id)
        k = self.reader.contract
        sites = {j.source_id: j for j in original.junctions}
        if len(sites) < 2:
            raise ValueError('A compound physical cluster requires multiple sites')

        following, records, by_street, sources = {}, {}, {}, {}
        for site in original.junctions:
            ring = [(site.source_id, c.source_id) for c in site.connectors]
            if len(set(ring)) != len(ring):
                raise ValueError('Physical connector identity is not unique')
            for at, connector in enumerate(site.connectors):
                key = ring[at]
                member, = connector.members
                street = site.streets[member.street_index]
                old = sources.setdefault(street.source_id, street)
                if old != street:
                    raise ValueError('Physical sites disagree on a complete original Street')
                records[key] = (street, member, connector)
                by_street.setdefault(street.source_id, []).append(key)
                following[key] = ring[(at+1) % len(ring)]

        crossing, internal, adjacency = {}, [], {pid: set() for pid in sites}
        bundle_by_incidence = {(pid, cid): bid for pid, cid, bid in original.endpoint_bundles}
        for sid, incidences in sorted(by_street.items()):
            if len(incidences) == 1:
                continue
            if len(incidences) != 2 or {records[i][1].leaving for i in incidences} != {False, True}:
                raise ValueError('Internal Streets require two distinct authored endpoint incidences')
            start = next(i for i in incidences if records[i][1].leaving)
            end = next(i for i in incidences if not records[i][1].leaving)
            if start[0] == end[0]:
                raise ValueError('A returning loop requires a separate surface with an island')
            if any(bundle_by_incidence[i] >= 0 for i in incidences):
                raise ValueError('An internal Bundle member needs explicit separator ownership')
            crossing[start], crossing[end] = end, start
            internal.append((sid, (start, end)))
            adjacency[start[0]].add(end[0]); adjacency[end[0]].add(start[0])

        visited, pending = set(), [min(sites)]
        while pending:
            pid = pending.pop()
            if pid not in visited:
                visited.add(pid); pending.extend(adjacency[pid]-visited)
        if visited != set(sites) or len(internal) != len(sites)-1:
            raise ValueError('The single-boundary compound branch requires a connected tree; cycles retain islands')

        external = set(records)-set(crossing)
        if not external:
            raise ValueError('A compound cluster requires external approaches')
        ring, seen, current = [], set(), min(external)
        while current not in seen:
            if current not in external:
                raise ValueError('Contracted traversal did not reach an external endpoint')
            ring.append(current); seen.add(current)
            current = following[current]
            skipped = set()
            while current in crossing:
                if current in skipped:
                    raise ValueError('Internal traversal closes without reaching the exterior')
                skipped.add(current)
                current = following[crossing[current]]
        if current != ring[0] or set(ring) != external:
            raise ValueError('Physical contraction has more than one external boundary')

        # Existing logical connector identity defines pairing. Positions and
        # angles never establish that two members belong to the same Bundle.
        logical = {row['junction']['junction_id']: row for row in snapshot.junction_origins['logical']}
        groups, group_of, observed = {}, {}, set()
        for jid in original.logical_junction_ids:
            row = logical[jid]
            connectors = {c['connector_id']: c for c in row['junction']['ring']}
            for endpoint in row['endpoints']:
                if endpoint['status'] != 'EXACT':
                    raise ValueError('Compound construction requires exact physical provenance')
                key = (endpoint['physical_junction_id'], endpoint['physical_connector_id'])
                if key not in records:
                    raise ValueError('Logical provenance references a site outside its component')
                if key in observed:
                    raise ValueError('A physical endpoint has more than one logical owner')
                observed.add(key)
                street, member, _ = records[key]
                source = endpoint['endpoint']
                terminal = street.sections[0 if member.leaving else -1]
                if (source['street_id'], source['end'], source['terminal_kind'], source['terminal_id']) != (
                        street.source_id, int(not member.leaving), 'Section', terminal.source_id):
                    raise ValueError('Logical provenance disagrees with its original terminal')
                if key in crossing:
                    continue
                cid = endpoint['logical_connector_id']
                if cid not in connectors or key in group_of:
                    raise ValueError('An external physical endpoint has ambiguous logical ownership')
                group_of[key] = cid
                entry = groups.setdefault(cid, (connectors[cid], []))
                if entry[0] != connectors[cid]:
                    raise ValueError('Logical connector identities collide')
                entry[1].append(key)
        if set(group_of) != external:
            raise ValueError('Every external endpoint requires its original logical owner')
        if observed != set(records):
            raise ValueError('Every internal endpoint also requires its exact logical provenance')

        # Rotate only at an existing group boundary, preserving cyclic order.
        cut = next((i for i, key in enumerate(ring) if group_of[key] != group_of[ring[i-1]]), None)
        if cut is None:
            raise ValueError('A compound patch requires multiple external connectors')
        ring = ring[cut:]+ring[:cut]
        ordered, emitted = [], set()
        for key in ring:
            cid = group_of[key]
            if not ordered or ordered[-1][0] != cid:
                if cid in emitted:
                    raise ValueError('Members of a logical connector are interleaved in the physical boundary')
                ordered.append((cid, [])); emitted.add(cid)
            ordered[-1][1].append(key)

        streets, connectors, profiles, origins = [], [], [], []
        bundle_indices = {}
        for cid, keys in ordered:
            source, endpoints = groups[cid]
            if set(keys) != set(endpoints):
                raise ValueError('A logical connector is missing a physical member')
            is_bundle = source['item_kind'] == 'Bundle'
            if source['item_kind'] not in ('Street', 'Bundle') or len(keys) != (2 if is_bundle else 1):
                raise ValueError('Compound mouths require ordinary Streets or complete two-member Bundles')
            source_ends = {(e['street_id'], e['end']) for e in source['endpoints']}
            if source_ends != {(records[key][0].source_id, int(not records[key][1].leaving)) for key in keys}:
                raise ValueError('A logical connector does not retain all original endpoints')
            # Keep Blosm's member order, separately from the exterior ring.
            by_endpoint = {(records[key][0].source_id, int(not records[key][1].leaving)): key for key in keys}
            members = []
            bid = source['item_id'] if is_bundle else -1
            bindex = bundle_indices.setdefault(bid, len(bundle_indices)) if is_bundle else -1
            for order, endpoint in enumerate(source['endpoints']):
                key = by_endpoint[endpoint['street_id'], endpoint['end']]
                street, member, _ = records[key]
                if bundle_by_incidence[key] != bid:
                    raise ValueError('Physical and logical Bundle membership disagree')
                if is_bundle:
                    role = next(m for m in snapshot.bundles[bid]['topology']['members'] if m['street_id'] == street.source_id)
                    head, tail = role['head_order'], role['tail_order']
                else:
                    head = tail = -1
                members.append(k.MemberInput(len(streets), member.leaving, order, bindex, head, tail))
                streets.append(street)
            profiles.append(k.ProfileInput(self.reader.rise, separator=is_bundle))
            # The final original incidence precedes the next boundary corner.
            # Reuse its declared radius, without averaging or shrinking it.
            radius = records[keys[-1]][2].radius
            connectors.append(k.ConnectorInput(cid, tuple(members), len(profiles)-1,
                radius, k.BUNDLE_ITEM if is_bundle else k.STREET_ITEM, source['leaving'], bid))
            origins.append((cid, tuple(keys)))

        owner = min(original.logical_junction_ids)
        # Reference-frame origin only; every original endpoint stays unchanged.
        location = logical[owner]['junction']['location']
        junction = k.PreparedJunction(owner, tuple(location), tuple(streets), tuple(connectors),
                                     tuple(profiles), cell_size=self.reader.cell_size)
        k.validate(junction)
        return PhysicalClusterInputs(owner, original.logical_junction_ids, tuple(sorted(sites)),
            junction, tuple(sources[sid] for sid, _ in internal), tuple(internal),
            tuple(ring), tuple(origins))


@dataclass(frozen=True)
class PhysicalCorridorInputs:
    cluster: PhysicalClusterInputs
    junctions: tuple
    # Original internal Street ID -> crown side in its authored frame.
    internal_profiles: tuple


class PreparedPhysicalCorridors:
    """Reconcile median boundaries inside an explicitly owned compound region.

    External incoming profiles remain fixed. Only Streets whose two ends are
    internal to the region receive a derived transverse profile. They are not
    reclassified as Bundle members. Coordinate, width, offset, radius and source
    identity stay unchanged. The resulting native patches still require normal
    geometric, curvature, seam and coverage qualification.
    """
    def __init__(self, cluster_reader):
        self.reader = cluster_reader

    @staticmethod
    def _profiles(junctions, internal_ids):
        domains = {sid: {-1, 0, 1} for sid in internal_ids}
        constraints = []
        for junction in junctions:
            endpoints = []
            for connector in junction.connectors:
                member, = connector.members
                sid = junction.streets[member.street_index].source_id
                half = junction.profiles[connector.profile_index].half_profile
                endpoints.append((sid if sid in domains else None,
                                  1 if member.leaving else -1, half))
            constraints.extend(zip(endpoints, endpoints[1:]+endpoints[:1]))

        def values(endpoint):
            sid, sign, fixed = endpoint
            return {(value, value*sign) for value in domains[sid]} if sid is not None else {(fixed, fixed)}

        def allowed(a, b):
            return (a == -1) == (b == 1)

        # Arc consistency over original incidences. No nearest-neighbor or
        # angular inference is used to continue a separator boundary.
        changed = True
        while changed:
            changed = False
            for first, second in constraints:
                left, right = values(first), values(second)
                for endpoint, possible, other, invert in ((first, left, right, False), (second, right, left, True)):
                    sid = endpoint[0]
                    supported = {value for value, half in possible if any(
                        allowed(other_half, half) if invert else allowed(half, other_half) for _, other_half in other)}
                    if not supported:
                        raise ValueError('External profiles have incompatible median continuations through the compound region')
                    if sid is not None and supported != domains[sid]:
                        domains[sid] = supported; changed = True
        chosen = {}
        for sid, possible in domains.items():
            if 0 in possible: chosen[sid] = 0
            elif len(possible) == 1: chosen[sid] = next(iter(possible))
            else: raise ValueError('Internal profile continuation is ambiguous')
        for a, b in constraints:
            value = lambda endpoint: chosen[endpoint[0]]*endpoint[1] if endpoint[0] is not None else endpoint[2]
            if not allowed(value(a), value(b)):
                raise ValueError('Keeping unconstrained full crowns needs a more general internal profile design')
        return chosen

    def prepare(self, snapshot, logical_junction_id):
        cluster = self.reader.prepare(snapshot, logical_junction_id)
        physical = self.reader.reader.prepare(snapshot, logical_junction_id)
        chosen = self._profiles(physical.junctions, {s.source_id for s in cluster.internal_streets})
        junctions = []
        for junction in physical.junctions:
            profiles = list(junction.profiles)
            for connector in junction.connectors:
                member, = connector.members
                sid = junction.streets[member.street_index].source_id
                if sid in chosen:
                    index = connector.profile_index
                    profiles[index] = replace(profiles[index], half_profile=chosen[sid]*(1 if member.leaving else -1))
            junctions.append(replace(junction, profiles=tuple(profiles)))
        return PhysicalCorridorInputs(cluster, tuple(junctions), tuple(sorted(chosen.items())))
