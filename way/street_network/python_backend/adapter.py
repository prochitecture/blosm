"""Transport prepared Blosm records into the numerical backend without bpy.

Original physical incidence, retained Sections and explicit Street item order
remain authoritative. Coordinate proximity never invents a connection.
"""
from dataclasses import dataclass, replace
import math
from .model import Section, End, Intersection, Connection, Bundle, BundleMember, Traversal, Network, ground_extent_id


@dataclass(frozen=True)
class Adaptation:
    network: Network
    issues: tuple
    blocked_sections: frozenset


class PreparedPythonAdapter:
    def __init__(self, categories, ground_policy, category_properties):
        self.categories = frozenset(categories)
        self.ground_policy = ground_policy
        self.category_properties = category_properties

    def read(self, snapshot, network_uuid, settings, *, rise=.08, interpolation=3):
        sections, issues, blocked = {}, [], set()
        originals = {i['id']: i for s in snapshot.streets.values()
                     for i in s['items'] if i['kind'] == 'Section'}
        excluded_grade = {}
        categories = {i['id']: i.get('category') for s in snapshot.streets.values()
                      for i in s['items'] if i['kind'] == 'Section'}
        for sid, street in sorted(snapshot.streets.items()):
            for item in street['items']:
                if item['kind'] != 'Section' or item.get('category') not in self.categories:
                    continue
                reason = self.ground_policy(item.get('tags', {}))
                if reason:
                    issues.append(dict(kind='Section', id=item['id'], reason=reason,
                        code='OUTSIDE_GROUND_SCOPE', source_way_ids=item.get('source_way_ids', ())))
                    excluded_grade[item['id']] = reason
                    continue
                identity = item['id']
                if identity in sections:
                    raise ValueError('A prepared Section has multiple Street owners')
                try:
                    placement = item.get('source_placement') or {}
                    if placement.get('status') == 'UNRESOLVED':
                        raise ValueError(placement.get('reason', 'Unresolved original OSM placement'))
                    points = tuple(tuple(p[:2]) + (float(p[2]) if len(p) > 2 else 0.,) for p in item['points'])
                    sections[identity] = Section(identity, sid, points, float(item['width']),
                        float(item.get('offset') or 0.), rise, interpolation, 0., item['category'])
                except (TypeError, ValueError) as exc:
                    issues.append(dict(kind='Section', id=identity, reason=str(exc)))
        sites, attached = [], set()
        connections = []

        def add_site(identity, refs, kind='Intersection'):
            requested = tuple(End(int(sid), int(end)) for sid, end in refs if categories.get(sid) in self.categories)
            available = tuple(e for e in requested if e.section_id in sections)
            if len(requested) < 2:
                return
            if len(available) != len(requested):
                # The ground-only provider constructs the induced source graph,
                # not an invented connection to an elevated road. Retain every
                # deferred incidence explicitly. Invalid/missing source records
                # are different from intentional grade exclusions and still
                # block the site. All incidences must share their original node.
                deferred = tuple(e for e in requested if e.section_id not in sections)
                missing = [e.section_id for e in deferred]
                original_positions = [originals[e.section_id].get('points', ())[0 if e.end == 0 else -1]
                    for e in requested if originals[e.section_id].get('points')]
                valid_node = (len(original_positions) == len(requested)
                    and all(len(p) >= 2 and all(math.isfinite(float(v)) for v in p[:2]) for p in original_positions)
                    and all(math.dist(original_positions[0][:2], p[:2]) <= 1e-6 for p in original_positions))
                if all(sid in excluded_grade for sid in missing) and valid_node:
                    if len(set(requested)) != len(requested) or attached.intersection(available):
                        raise ValueError('Prepared endpoint incidence is duplicated')
                    if not available:
                        return  # Entirely outside this provider's explicit scope.
                    first = sections[available[0].section_id].points[0 if available[0].end == 0 else -1]
                    deferred_records = [dict(section=e.section_id, end=e.end,
                        reason=excluded_grade[e.section_id],
                        source_way_ids=originals[e.section_id].get('source_way_ids', ())) for e in deferred]
                    if len(available) == 1:
                        endpoint, = available
                        section = sections[endpoint.section_id]
                        sections[endpoint.section_id] = replace(section,
                            scope_ends=tuple(sorted((*section.scope_ends, endpoint.end))))
                        attached.add(endpoint)
                        partner = deferred[0]
                        issues.append(dict(kind='ScopeBoundary', id=identity, source_kind=kind,
                            reason='Ground surface ends at an original grade transition',
                            section=endpoint.section_id, end=endpoint.end,
                            deferred_section=partner.section_id, deferred_end=partner.end,
                            deferred_reason=excluded_grade[partner.section_id],
                            deferred=deferred_records, location=first))
                        return
                    issues.append(dict(kind='ScopeJunction', id=identity, source_kind=kind,
                        reason='Ground-only source junction; elevated approaches are explicitly deferred',
                        available=[dict(section=e.section_id, end=e.end) for e in available],
                        deferred=deferred_records, location=first))
                else:
                    issues.append(dict(kind=kind, id=identity, reason='Required source approaches absent', sections=missing))
                    blocked.update(e.section_id for e in available)
                    return
            if len(set(available)) != len(available) or attached.intersection(available):
                raise ValueError('Prepared endpoint incidence is duplicated')
            attached.update(available)
            if len(available) == 2:
                connections.append(Connection(len(connections), available))
            else:
                if any(site.id == identity for site in sites):
                    raise ValueError('Conflicting original intersection identities')
                radii = tuple(float(self.category_properties[sections[e.section_id].category]['radius']) for e in available)
                sites.append(Intersection(identity, available, radii[0], radii, ('STYLE',) * len(available)))

        for site in snapshot.junction_origins.get('physical', ()):
            refs = [(e['terminal_id'], e['end']) for connector in site['ring'] for e in connector['endpoints']
                    if e['terminal_kind'] == 'Section']
            if refs:
                add_site(site['junction_id'], refs)
        for street in snapshot.streets.values():
            items = street['items']
            for item in items:
                if item['kind'] == 'Intersection' and item.get('surface_treatment') == 'NATIVE_JUNCTION':
                    refs = list(item.get('through', ()))
                    refs.extend((a['section_id'], a['end']) for a in item.get('attachments', ())
                                if a.get('category') in self.categories)
                    add_site(item['id'], refs)
            # A Street's explicit Section sequence provides seam identity.
            # Never connect unrelated coincident controls from different roads.
            ordered = [(i, item) for i, item in enumerate(items) if item['kind'] == 'Section']
            for (ia, a), (ib, b) in zip(ordered, ordered[1:]):
                ea, eb = End(a['id'], 1), End(b['id'], 0)
                if ea in attached or eb in attached:
                    continue
                between = items[ia + 1:ib]
                if any(i.get('surface_treatment') == 'NATIVE_JUNCTION' for i in between):
                    continue
                if a['id'] not in sections or b['id'] not in sections:
                    # Explicit item order also records grade transitions that
                    # do not allocate a separate physical connector ring.
                    add_site(street['id'], ((a['id'], 1), (b['id'], 0)), 'StreetSequence')
                    continue
                connections.append(Connection(len(connections), (ea, eb)))
                attached.update((ea, eb))
        for seam in snapshot.junction_origins.get('seams', ()):
            refs = [(e['terminal_id'], e['end']) for e in seam['endpoints'] if e['terminal_kind'] == 'Section']
            available = [End(*ref) for ref in refs if ref[0] in sections]
            if available and not attached.intersection(available):
                add_site(seam['junction_id'], refs, 'SectionSeam')
        bundles = []
        for bid, source in sorted(snapshot.bundles.items()):
            topology = source.get('topology', {})
            members = topology.get('members', ())
            member_ids = {m['street_id'] for m in members}
            mapped = topology.get('status') == 'MAPPED' and len(members) >= 2
            if not mapped:
                # An unresolved mapper can have no member descriptors at all.
                # Its original ownership still prevents the same ground roads
                # from silently becoming unrelated full-crown surfaces.
                member_ids.update(source.get('owned_streets', ()))
                member_ids.update(source.get('head', ()))
                member_ids.update(source.get('tail', ()))
                member_ids.update(sid for sid, street in snapshot.streets.items()
                                  if street.get('bundle_id') == bid)
            selected = [s for s in sections.values() if s.street_id in member_ids]
            if not selected:
                continue
            if not mapped:
                blocked.update(s.id for s in selected)
                issues.append(dict(kind='Bundle', id=bid, code='UNMAPPED_BUNDLE_MEMBERSHIP',
                    reason=topology.get('reason') or 'Bundle membership is not mapped',
                    sections=sorted(s.id for s in selected)))
                continue
            by_street = {m['street_id']: m for m in members}
            chains = topology.get('chains') or [[m['street_id']] for m in sorted(members, key=lambda m: m['head_order'])]
            lateral_members = []
            try:
                for chain in chains:
                    path = []
                    complete = True
                    for street_id in chain:
                        member = by_street[street_id]
                        source_sections = [i for i in snapshot.streets[street_id]['items'] if i['kind'] == 'Section']
                        forward = member['forward']
                        path.extend(Traversal(i['id'], forward) for i in
                                    (source_sections if forward else reversed(source_sections)))
                        complete &= member['complete']
                    lateral_members.append(BundleMember(tuple(path), complete))
                scope_record = None
                extents = ((tuple(lateral_members), None),)
                if any(t.section_id not in sections for member in lateral_members for t in member.path):
                    from .bundle_scope import ground_bundle_scopes
                    extents = ground_bundle_scopes(
                        bid, tuple(lateral_members), sections, excluded_grade)
                interior_items = [[i for i in snapshot.streets[sid]['items'] if i['kind'] == 'Section']
                                  for sid in topology.get('interior_streets', ())]
                interior_records = []
                if extents[0][1] is not None:
                    from .bundle_scope_interiors import scoped_interior_roles
                    runs = []
                    for items in interior_items:
                        run = []
                        for item in (*items, None):
                            if item is not None and item['id'] in sections:
                                run.append(Traversal(item['id']))
                                continue
                            if run:
                                runs.append(tuple(run))
                                run = []
                            if item is not None:
                                if item['id'] not in excluded_grade:
                                    raise ValueError('Scoped Bundle interior includes an unexplained absent Section')
                                interior_records.append(dict(section_ids=[item['id']], role='OUTSIDE_GROUND_SCOPE',
                                                             reason=excluded_grade[item['id']]))
                    interiors, retained_records = scoped_interior_roles(
                        tuple(members for members, _ in extents), runs, sections, (*sites, *connections))
                    interior_records.extend(retained_records)
                else:
                    interiors = (tuple(i['id'] for items in interior_items for i in items if i['id'] in sections),)
                pending, records = [], []
                for index, ((members, record), interior_ids) in enumerate(zip(extents, interiors)):
                    derived = len(extents) > 1
                    ordinal = index+1 if derived else None
                    identity = ground_extent_id(bid, ordinal) if derived else bid
                    if derived and (identity in snapshot.bundles or any(b.id == identity for b in bundles)):
                        raise ValueError('Derived ground Bundle extent collides with an existing runtime identity')
                    pending.append(Bundle(identity, members, interior_ids,
                                          source_bundle_id=bid if derived else None, scope_ordinal=ordinal))
                    if record is not None:
                        record.update(runtime_id=identity, source_bundle_id=bid,
                                      scope_ordinal=ordinal, interior_source_roles=interior_records)
                        records.append(record)
                bundles.extend(pending)
                issues.extend(records)
            except (KeyError, ValueError) as exc:
                blocked.update(s.id for s in selected)
                issues.append(dict(kind='Bundle', id=bid, reason=str(exc)))
        network = Network(network_uuid, tuple(sections.values()), tuple(sites), tuple(bundles), settings=settings,
                          connections=tuple(connections))
        return Adaptation(network, tuple(issues), frozenset(blocked))


def components(network):
    """Partition by declared topology/Bundle membership, retaining source IDs."""
    adjacency = {section.id: set() for section in network.sections}
    for group in (*network.intersections, *network.connections, *network.bundles):
        ids = tuple(e.section_id for e in group.ends) if hasattr(group, 'ends') else group.sections
        for sid in ids:
            adjacency[sid].update(ids)
    seen = set()
    for sid in sorted(adjacency):
        if sid in seen:
            continue
        found, pending = set(), [sid]
        while pending:
            at = pending.pop()
            if at in found:
                continue
            found.add(at)
            pending.extend(adjacency[at] - found)
        seen.update(found)
        yield replace(network, uuid=f'{network.uuid}/component/{min(found)}',
            sections=tuple(s for s in network.sections if s.id in found),
            intersections=tuple(j for j in network.intersections if j.ends[0].section_id in found),
            connections=tuple(c for c in network.connections if c.ends[0].section_id in found),
            bundles=tuple(b for b in network.bundles if b.sections[0] in found),
            furniture=tuple(p for p in network.furniture if p.section_id in found))
