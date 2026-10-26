"""Transport prepared Blosm records into the numerical backend without bpy.

Original physical incidence, retained Sections and explicit Street item order
remain authoritative. Coordinate proximity never invents a connection.
"""
from dataclasses import dataclass, replace
from .model import Section, End, Intersection, Connection, Bundle, BundleMember, Traversal, Network


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
        categories = {i['id']: i.get('category') for s in snapshot.streets.values()
                      for i in s['items'] if i['kind'] == 'Section'}
        for sid, street in sorted(snapshot.streets.items()):
            for item in street['items']:
                if item['kind'] != 'Section' or item.get('category') not in self.categories:
                    continue
                reason = self.ground_policy(item.get('tags', {}))
                if reason:
                    issues.append(dict(kind='Section', id=item['id'], reason=reason))
                    continue
                identity = item['id']
                if identity in sections:
                    raise ValueError('A prepared Section has multiple Street owners')
                try:
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
                # Excluded transport categories do not supply road width. A
                # missing eligible road is a blocked junction, never a free end.
                missing = [e.section_id for e in requested if e.section_id not in sections]
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
                if ea in attached or eb in attached or a['id'] not in sections or b['id'] not in sections:
                    continue
                between = items[ia + 1:ib]
                if any(i.get('surface_treatment') == 'NATIVE_JUNCTION' for i in between):
                    continue
                connections.append(Connection(len(connections), (ea, eb)))
                attached.update((ea, eb))
        for seam in snapshot.junction_origins.get('seams', ()):
            refs = [(e['terminal_id'], e['end']) for e in seam['endpoints'] if e['terminal_kind'] == 'Section']
            available = [End(*ref) for ref in refs if ref[0] in sections]
            if len(available) == 2 and not attached.intersection(available):
                add_site(seam['junction_id'], refs, 'SectionSeam')
        bundles = []
        for bid, source in sorted(snapshot.bundles.items()):
            topology = source.get('topology', {})
            members = topology.get('members', ())
            selected = [s for s in sections.values() if s.street_id in {m['street_id'] for m in members}]
            if not selected:
                continue
            if topology.get('status') != 'MAPPED' or len(members) < 2:
                blocked.update(s.id for s in selected)
                issues.append(dict(kind='Bundle', id=bid, reason='Bundle membership is not mapped'))
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
                        if any(i['id'] not in sections for i in source_sections):
                            raise ValueError('Bundle member includes an excluded or invalid Section')
                        forward = member['forward']
                        path.extend(Traversal(i['id'], forward) for i in
                                    (source_sections if forward else reversed(source_sections)))
                        complete &= member['complete']
                    lateral_members.append(BundleMember(tuple(path), complete))
                interiors = tuple(i['id'] for sid in topology.get('interior_streets', ())
                                  for i in snapshot.streets[sid]['items'] if i['kind'] == 'Section' and i['id'] in sections)
                bundles.append(Bundle(bid, tuple(lateral_members), interiors))
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
