"""Separator identities across several original Streets in one lateral chain.

Coordinates never determine connectivity here. The prepared Bundle chain and
original generated endpoint incidences must agree at every join. This opt-in
preparation is separate from native asset acceptance and surface qualification.
"""
from .original_span_separator_chains import ChainedBundleSeparator, FreeChainedBundleSeparator


class LongitudinalSeparatorPreparation:
    def __init__(self, planner):
        self.planner = planner

    def candidate(self, bundle, snapshot, component, owners=None, *, free_ends=False):
        p = self.planner
        owners = p.owners(component) if owners is None else owners
        topology = bundle.get('topology', {})
        chains = tuple(tuple(chain) for chain in topology.get('chains', ()))
        members = topology.get('members', ())
        ids = tuple(sid for chain in chains for sid in chain)
        if (topology.get('status') != 'MAPPED' or len(chains) != 2 or any(not chain for chain in chains)
                or len(ids) != len(set(ids)) or len(members) != len(ids)
                or {m['street_id'] for m in members} != set(ids)):
            raise ValueError('Longitudinal separator needs two distinct complete mapped chains')
        allowed = {'LONGITUDINAL_CHAINS', 'INTERNAL_JUNCTIONS', 'INTERIOR_STREETS'}
        if set(topology.get('required_features', ())) - allowed:
            raise ValueError('Longitudinal separator has unsupported additional topology requirements')
        p.deferred_interiors(bundle, snapshot, component)
        if not set(ids) <= set(component.ports.source_streets):
            raise ValueError('Every original Street in both chains must belong to the component')
        by_id = {m['street_id']: m for m in members}
        declarations = {}
        for junction in topology.get('junctions', ()):
            key = ('SECTION', junction['street_id'], junction['junction_id'])
            if key in declarations:
                raise ValueError('Duplicate intermediate Bundle junction declaration')
            declarations[key] = junction
        for junction in topology.get('chain_joins', ()):
            before, after = tuple(junction['before']), tuple(junction['after'])
            key = ('STREET', before, after, junction['junction_id'])
            if len(before) != 3 or len(after) != 3 or key in declarations:
                raise ValueError('Duplicate or malformed original Street-chain join')
            declarations[key] = junction
        def endpoint_key(end):
            point = end.endpoint
            return point.street_id, point.section_id, point.end
        all_routes, all_middle, all_ends, all_terminals, used = [], [], [], [], set()
        for side, chain in enumerate(chains):
            ordered = []
            for piece, sid in enumerate(chain):
                member = by_id[sid]
                if (member['head_order'] != side or member['tail_order'] != side or
                        not member['complete'] or type(member['forward']) is not bool):
                    raise ValueError('Original chain member has incomplete or inconsistent lateral ownership')
                source = snapshot.streets.get(sid)
                if source is None or source.get('bundle_id') != bundle['id']:
                    raise ValueError('Original chain Street has another Bundle owner')
                routes = sorted((r for r in component.network.routes if r.span.street_id == sid),
                                key=lambda r: r.span.first_control)
                kinds = {'CORRIDOR', 'APPROACH'} if free_ends else {'CORRIDOR'}
                if not routes or any(r.kind not in kinds or not r.span.in_scope for r in routes):
                    raise ValueError('Every chain member needs accepted in-scope original routes')
                sections = tuple(s['id'] for s in source['items'] if s['kind'] == 'Section')
                actual = tuple(s for r in routes for s in r.span.section_ids)
                spans = [r.span for r in routes]; count = spans[0].source_control_count
                if (actual != sections or len(set(actual)) != len(actual) or spans[0].first_control != 0
                        or spans[-1].last_control != count-1 or
                        any(s.source_control_count != count or not 0 <= s.first_control < s.last_control < count for s in spans)
                        or any(a.last_control != b.first_control for a, b in zip(spans, spans[1:]))):
                    raise ValueError('Longitudinal routes must cover each complete original Street without gaps or repeats')
                directed = []
                for route in routes:
                    ends = tuple(p.route_end(route, at, allow_free=free_ends) for at in (0, 1))
                    directed.append((route, ends if member['forward'] else ends[::-1]))
                if not member['forward']:
                    directed.reverse()
                extents = member['spans']
                if len(extents) != len(directed):
                    raise ValueError('Original member extents and generated routes disagree')
                for ordinal, (extent, (route, ends)) in enumerate(zip(extents, directed)):
                    expected = route.span.section_ids if member['forward'] else route.span.section_ids[::-1]
                    rows = extent.get('sections', ())
                    if (extent.get('ordinal') != ordinal or tuple(s['start']['section_id'] for s in rows) != expected
                            or any(s['start']['section_id'] != s['end']['section_id'] for s in rows)):
                        raise ValueError('Mapped extent differs from original directed Sections')
                    for at, label in enumerate(('start', 'end')):
                        edge = (piece == 0 and ordinal == 0 and at == 0 or
                                piece == len(chain)-1 and ordinal == len(directed)-1 and at == 1)
                        expected_label = ('HEAD' if at == 0 else 'TAIL') if edge else None
                        expected_site = None if edge else ends[at].junction_id
                        if (extent[label].get('bundle_end') != expected_label or
                                extent[label].get('junction_id') != expected_site):
                            raise ValueError('Bundle extent does not retain its exact chain boundary')
                ordered.extend(directed)
            free_positions = [(i, at) for i, (_, ends) in enumerate(ordered)
                              for at, end in enumerate(ends) if end.owner_kind == 'FREE']
            if free_ends:
                if len(free_positions) != 1 or free_positions[0] not in ((0, 0), (len(ordered)-1, 1)):
                    raise ValueError('A complete free chain needs exactly one free outer endpoint')
                free_route = free_positions[0][0]
                if any(r.kind != ('APPROACH' if i == free_route else 'CORRIDOR') for i, (r, _) in enumerate(ordered)):
                    raise ValueError('Only the complete chain terminal may be an approach')
            elif free_positions:
                raise ValueError('A closed chain cannot contain a free endpoint')
            middle = []
            for (_, prior), (_, following) in zip(ordered, ordered[1:]):
                before, after = prior[1], following[0]
                a, b = endpoint_key(before), endpoint_key(after)
                jid = before.junction_id
                if jid is None or jid != after.junction_id or jid not in owners:
                    raise ValueError('Adjacent routes need their exact original generated junction')
                if a[0] == b[0]:
                    key = ('SECTION', a[0], jid)
                    declaration = declarations.get(key, {})
                    through = tuple(tuple(ref) for ref in declaration.get('through', ()))
                    if len(through) != 2 or set(through) != {a[1:], b[1:]}:
                        raise ValueError('Original through Sections do not match the longitudinal chain')
                else:
                    key = ('STREET', a, b, jid)
                    if key not in declarations:
                        raise ValueError('Consecutive Streets lack their explicit original chain join')
                if key in used:
                    raise ValueError('An original chain junction cannot be consumed twice')
                used.add(key); middle.append(owners[jid])
            ends = (ordered[0][1][0], ordered[-1][1][1]); terminal_owners = []
            for end in ends:
                if free_ends and end.owner_kind == 'FREE':
                    terminal_owners.append(())
                elif end.junction_id in owners:
                    terminal_owners.append(owners[end.junction_id])
                else:
                    raise ValueError('A chained separator end has no accepted surface owner')
            all_routes.append(tuple((r.span.street_id, r.span.section_ids[0], r.span.section_ids[-1]) for r, _ in ordered))
            all_middle.append(tuple(middle)); all_ends.append(tuple(e.endpoint for e in ends))
            all_terminals.append(tuple(terminal_owners))
        if not used or used != set(declarations):
            raise ValueError('Longitudinal separator must consume every declared intermediate join')
        if all_terminals[0] != all_terminals[1] or all_terminals[0][0] == all_terminals[0][1]:
            raise ValueError('Lateral chains must agree on attached owners and free ends')
        middle = [owner for rows in all_middle for owner in rows]
        if len(middle) != len(set(middle)) or set(middle).intersection(all_terminals[0]):
            raise ValueError('Shared or repeated intermediate owners need multiple separator regions')
        args = (bundle['id'], tuple(chain[0] for chain in chains), *all_terminals[0], tuple(all_routes),
                tuple(ends[0] for ends in all_ends), tuple(ends[1] for ends in all_ends), tuple(all_middle))
        return (FreeChainedBundleSeparator(*args, bool(all_terminals[0][0]), member_chains=chains) if free_ends else
                ChainedBundleSeparator(*args, member_chains=chains))
