"""Original logical ownership for a connected cluster of physical port rings.

This prepares topology and provenance, not a surface or a hole count. Internal
Bundle spans remain explicit sources. A physical cycle inside one logical
intersection must not be reclassified as a traffic island by this adapter.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class LogicalPortGroup:
    connector_id: int
    item_kind: str
    item_id: int
    leaving: bool
    # Original Blosm member order, independently of the physical boundary walk.
    ports: tuple
    bundle_member_orders: tuple


@dataclass(frozen=True)
class SharedLogicalPortGroup(LogicalPortGroup):
    # A jointly owned physical patch retains each original logical ring ID.
    logical_junction_id: int


@dataclass(frozen=True)
class LogicalPortBoundary:
    start_port: object
    end_port: object
    # A corner is addressed by its first original port in the physical ring.
    original_corners: tuple
    # Original span identity and the endpoint from which it is traversed.
    traversed_spans: tuple


@dataclass(frozen=True)
class LogicalPortCluster:
    logical_junction_id: int
    physical_junction_ids: tuple
    internal_spans: tuple
    groups: tuple[LogicalPortGroup, ...]
    external_ring: tuple
    boundaries: tuple[LogicalPortBoundary, ...]
    source_cycle_rank: int
    # Logical aggregation can retain both terminal connectors of a Street
    # wholly inside this physical component. Keep those incidences explicitly;
    # they are not exterior mouths and never disappear from source ownership.
    internal_groups: tuple[LogicalPortGroup, ...] = ()
    footprint_spans: tuple = ()
    @property
    def owner_ids(self):
        return (self.logical_junction_id,)


@dataclass(frozen=True)
class SharedLogicalPortCluster(LogicalPortCluster):
    # The inherited scalar is the deterministic declaration key, not a
    # replacement for these original junction IDs.
    logical_junction_ids: tuple = ()

    @property
    def owner_ids(self):
        return self.logical_junction_ids


def _end(endpoint):
    return endpoint.street_id, endpoint.section_id, endpoint.end


def _record_end(endpoint):
    if endpoint['terminal_kind'] != 'Section':
        raise ValueError('Logical cluster endpoints must name original terminal Sections')
    return endpoint['street_id'], endpoint['terminal_id'], endpoint['end']


def _same_cycle(actual, expected):
    return len(actual) == len(expected) and any(
        actual == expected[i:]+expected[:i] for i in range(len(expected)))


class PreparedLogicalPortClusters:
    """Keep one fully mapped logical intersection and its original physical rings.

    Shared ownership is opt-in and requires the complete connected component.
    The input may include a larger network, but must contain this whole cluster.
    Geometry, profiles, separator shape and coverage are native responsibilities.
    """
    def prepare(self, snapshot, plan, logical_junction_id, *, allow_internal_connectors=False,
                allow_shared_owners=False):
        if type(logical_junction_id) is not int or logical_junction_id < 0:
            raise ValueError('Logical junction ID must be a nonnegative integer')
        origins = snapshot.junction_origins
        found = [c for c in origins.get('components', ()) if logical_junction_id in c['logical_junction_ids']]
        if len(found) != 1:
            raise ValueError('A logical junction needs one original physical component')
        component, = found
        owners = tuple(sorted(component['logical_junction_ids']))
        shared = len(owners) > 1
        if (component.get('status') != 'MAPPED' or component.get('unresolved') or
                len(set(owners)) != len(owners) or
                any(type(jid) is not int or jid < 0 for jid in owners) or
                (not allow_shared_owners and (component.get('shared_physical_junction_ids') or
                                             owners != (logical_junction_id,)))):
            raise ValueError('Cluster requires one fully mapped logical owner without shared physical sites')
        ids = tuple(component['physical_junction_ids'])
        if len(ids) < 2 or len(set(ids)) != len(ids):
            raise ValueError('Cluster needs multiple distinct original physical sites')
        sites = {s.junction_id: s for s in plan.sites}
        if len(sites) != len(plan.sites) or not set(ids).issubset(sites):
            raise ValueError('The complete original physical component must be selected')
        rows = [r for r in origins.get('logical', ()) if r['junction']['junction_id'] in owners]
        if len(rows) != len(owners) or {r['junction']['junction_id'] for r in rows} != set(owners):
            raise ValueError('Logical junction provenance must be unique')
        rows.sort(key=lambda r:r['junction']['junction_id'])
        if any(r['junction']['junction_id'] not in owners and
               any(e['status']=='EXACT' and e['physical_junction_id'] in ids for e in r['endpoints'])
               for r in origins.get('logical', ())):
            raise ValueError('Component omits an original logical owner of a selected physical site')
        # Reconstruct component membership from original incidences, including
        # its shared sites. Do not trust a convenient but incomplete declaration.
        memberships = {jid: set() for jid in ids}
        for logical in rows:
            jid = logical['junction']['junction_id']
            for record in logical['endpoints']:
                pid = record['physical_junction_id']
                if record['status'] != 'EXACT' or pid not in memberships:
                    raise ValueError('Logical cluster requires exact complete physical provenance')
                memberships[pid].add(jid)
        if any(not members for members in memberships.values()):
            raise ValueError('Declared physical component contains a site with no logical provenance')
        if set(component.get('shared_physical_junction_ids', ())) != {p for p,m in memberships.items() if len(m)>1}:
            raise ValueError('Shared physical sites disagree with original endpoint provenance')
        reached_owners = {owners[0]}
        while True:
            expanded = reached_owners.union(*(m for m in memberships.values() if m & reached_owners))
            if expanded == reached_owners:break
            reached_owners = expanded
        if reached_owners != set(owners):
            raise ValueError('Logical owners are not connected by shared original physical sites')
        ports, following, by_connector = {}, {}, {}
        for sid in ids:
            ring = sites[sid].ports
            if len(ring) < 3:
                raise ValueError('Each physical site requires its complete original ring')
            records = [r for r in origins.get('physical', ()) if r['junction_id'] == sid]
            if len(records) != 1:
                raise ValueError('Each physical site needs unique original ring provenance')
            expected = []
            for connector in records[0]['ring']:
                if connector['item_kind'] != 'Street' or len(connector['endpoints']) != 1:
                    raise ValueError('Physical ring must retain individual Street incidences')
                expected.append((connector['connector_id'], _record_end(connector['endpoints'][0])))
            actual = [(p.connector_id, _end(p.endpoint)) for p in ring]
            if not _same_cycle(actual, expected):
                raise ValueError('Selected physical ring differs from its original cyclic order')
            for at, port in enumerate(ring):
                key = _end(port.endpoint); connector = (sid, port.connector_id)
                if port.junction_id != sid or port.source_kind != 'PHYSICAL' or key in ports or connector in by_connector:
                    raise ValueError('Physical port incidence has conflicting original ownership')
                ports[key] = port; by_connector[connector] = port
                following[key] = _end(ring[(at+1) % len(ring)].endpoint)
        internal = []; twins = {}; internal_at = {}; reached = {sid: set() for sid in ids}
        span_ids = set(); seen = set()
        for span in plan.spans:
            if span.identity in span_ids:
                raise ValueError('Original source spans must have unique identities')
            span_ids.add(span.identity)
            ends = (span.start, span.end)
            for end in ends:
                if end.junction_id in reached:
                    key = _end(end.endpoint)
                    if end.owner_kind != 'PHYSICAL' or key not in ports or ports[key].junction_id != end.junction_id or key in seen:
                        raise ValueError('Every selected endpoint needs one original span incidence')
                    seen.add(key)
            if all(end.junction_id in reached for end in ends):
                a, b = (_end(end.endpoint) for end in ends)
                if a == b:
                    raise ValueError('An internal span cannot repeat one endpoint')
                twins[a], twins[b] = b, a; internal.append(span)
                internal_at[a] = internal_at[b] = span.identity
                reached[span.start.junction_id].add(span.end.junction_id)
                reached[span.end.junction_id].add(span.start.junction_id)
        if seen != set(ports):
            raise ValueError('Cluster contains an unowned original port')
        connected = set(); pending = [ids[0]]
        while pending:
            at = pending.pop()
            if at not in connected:
                connected.add(at); pending.extend(reached[at]-connected)
        if connected != set(ids):
            raise ValueError('Original internal spans must connect the physical cluster')
        external = set(ports)-set(twins)
        if not external:
            raise ValueError('A logical intersection needs external approaches')
        provenance = {}
        for jid, row in ((r['junction']['junction_id'], e) for r in rows for e in r['endpoints']):
            if row['status'] != 'EXACT':
                raise ValueError('Logical cluster requires exact endpoint provenance')
            key = _record_end(row['endpoint'])
            port = by_connector.get((row['physical_junction_id'], row['physical_connector_id']))
            if (port is None or key != _end(port.endpoint) or key in provenance or
                    (key not in external and not allow_internal_connectors)):
                raise ValueError('Logical endpoint disagrees with its original external port')
            provenance[key] = (jid, row['logical_connector_id'])
        if not external.issubset(provenance):
            raise ValueError('Every external port needs exactly one logical connector')
        declared_internal = set(provenance)-external
        footprint_spans = []
        for span in internal:
            ends = {_end(span.start.endpoint), _end(span.end.endpoint)}
            if not ends.intersection(declared_internal):
                continue
            if not ends.issubset(declared_internal):
                raise ValueError('An internal logical Street needs both original terminal incidences')
            footprint_spans.append(span.identity)
        groups = []; internal_groups = []; assigned = {}; connector_ids = set()
        for jid, connector in ((r['junction']['junction_id'], con) for r in rows for con in r['junction']['ring']):
            cid = connector['connector_id']; kind = connector['item_kind']; bid = connector['item_id']
            token = (jid, cid)
            if token in connector_ids or kind not in ('Street', 'Bundle'):
                raise ValueError('Logical connectors must have distinct IDs and supported kinds')
            connector_ids.add(token)
            keys = tuple(_record_end(e) for e in connector['endpoints'])
            if len(keys) != (2 if kind == 'Bundle' else 1) or len(set(keys)) != len(keys):
                raise ValueError('A logical approach needs its complete one- or two-member group')
            if any(key not in provenance or provenance[key] != token or key in assigned for key in keys):
                raise ValueError('Logical connector membership disagrees with exact endpoint provenance')
            inside = any(key in declared_internal for key in keys)
            if inside and (kind != 'Street' or any(key not in declared_internal for key in keys)):
                raise ValueError('Internal logical connectors require a complete individual Street')
            if inside and connector['leaving'] != (keys[0][2] == 0):
                raise ValueError('Internal logical connector direction disagrees with its original endpoint')
            orders = ()
            if kind == 'Bundle':
                bundle = snapshot.bundles.get(bid)
                # Overall Bundle status includes surface-provider support. An
                # unsupported remote patch must not erase known topology here;
                # this preparation does not claim that Bundle has been rendered.
                topology = bundle.get('topology', {}) if bundle is not None else {}
                if topology.get('status') != 'MAPPED':
                    raise ValueError('A logical Bundle needs its retained mapped source topology')
                members = topology.get('members', ())
                chains = topology.get('chains', ())
                if chains:
                    terminal_ids = [chain[0 if connector['leaving'] else -1] for chain in chains if chain]
                    members = [m for m in members if m['street_id'] in terminal_ids]
                    if len(chains) != 2 or len(terminal_ids) != 2:
                        raise ValueError('Logical Bundle approach needs two declared chain terminals')
                if len(members) != 2 or {m['street_id'] for m in members} != {key[0] for key in keys}:
                    raise ValueError('Logical Bundle membership does not match its original member Streets')
                roles = {m['street_id']: m for m in members}
                if chains and any((key[2] == 0) != (roles[key[0]]['forward'] == connector['leaving']) for key in keys):
                    raise ValueError('Logical chain terminal uses the wrong original Street end')
                orders = tuple((roles[key[0]]['head_order'], roles[key[0]]['tail_order']) for key in keys)
                if any(not roles[key[0]].get('complete') or snapshot.streets.get(key[0], {}).get('bundle_id') != bid for key in keys):
                    raise ValueError('Logical Bundle members need complete original ownership')
                if any(sorted(order[i] for order in orders) != [0, 1] for i in (0, 1)):
                    raise ValueError('Bundle endpoint orders must be permutations of its two members')
            elif keys[0][0] != bid:
                raise ValueError('A Street connector must retain its original Street identity')
            assigned.update((key, token) for key in keys)
            target = internal_groups if inside else groups
            args=(cid, kind, bid, connector['leaving'], tuple(ports[k] for k in keys), orders)
            target.append(SharedLogicalPortGroup(*args,jid) if shared else LogicalPortGroup(*args))
        if set(assigned) != set(provenance):
            raise ValueError('The logical ring omits an original port')
        # Contract only original twin incidences, never nearest points. Keep
        # the exterior cyclic order; internal cycles remain counted as source
        # topology, and are not turned into inferred grass/island boundaries.
        current = min(external); ring = []; boundaries = []; visited = set()
        while current not in visited:
            if current not in external:
                raise ValueError('Physical contraction did not reach an external port')
            start = current
            ring.append(current); visited.add(current)
            corners = [ports[current]]; crossed = []; current = following[current]
            skipped = set()
            while current in twins:
                if current in skipped:
                    raise ValueError('Internal traversal does not reach the exterior')
                skipped.add(current)
                crossed.append((internal_at[current], ports[current].endpoint))
                corners.append(ports[twins[current]])
                current = following[twins[current]]
            boundaries.append(LogicalPortBoundary(ports[start], ports[current], tuple(corners), tuple(crossed)))
        if current != ring[0] or visited != external:
            raise ValueError('Logical cluster has multiple external boundary walks')
        physical_groups = [assigned[k] for i,k in enumerate(ring) if assigned[k] != assigned[ring[i-1]]]
        expected = [(g.logical_junction_id if shared else owners[0], g.connector_id) for g in groups]
        if len(physical_groups) != len(expected) or set(physical_groups) != set(expected):
            raise ValueError('Logical approaches must be complete contiguous exterior groups')
        for jid in owners:
            actual = [token for token in physical_groups if token[0] == jid]
            original = [token for token in expected if token[0] == jid]
            if not original or not _same_cycle(actual, original):
                raise ValueError('Logical group order disagrees with the original physical boundary')
        if shared:
            by_token = dict(zip(expected, groups))
            groups = [by_token[token] for token in physical_groups]
        args=(min(owners), tuple(sorted(ids)),
            tuple(sorted(internal, key=lambda s:s.identity)), tuple(groups),
            tuple(ports[k] for k in ring), tuple(boundaries), len(internal)-len(ids)+1,
            tuple(internal_groups), tuple(sorted(footprint_spans)))
        return SharedLogicalPortCluster(*args,owners) if shared else LogicalPortCluster(*args)
