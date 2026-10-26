"""Original-incidence plans for separators interrupted by interior Streets.

Plans describe existing Section extents only. They do not modify Bundle
membership, substitute profiles, or infer attachment from nearby coordinates.
"""
from dataclasses import dataclass
from .model import BundleMember, End, Traversal


@dataclass(frozen=True)
class InteriorCrossing:
    street_id: int
    path: tuple
    sites: tuple

    @property
    def sections(self):
        return tuple(t.section_id for t in self.path)


def endpoint_sites(network):
    return {e: j for j in (*network.intersections, *network.connections) for e in j.ends}


def street_chain(network, first, *, allowed=None):
    """Walk one unbranched Street run between original physical endpoints.

    A physical junction terminates the walk. Internal Section connections must
    be explicit two-port incidences. With an interior inventory, every Section
    between these ports must be declared; the same Street may continue outside
    either physical junction. Without one, retain the complete-Street rule used
    by a crossing-to-free continuation.
    The terminal endpoint's owner is returned, including None for a true free
    end. Merely stopping at a clipped subset is never treated as a free end.
    """
    sections, sites = network.section_index, endpoint_sites(network)
    street_id = sections[first.section_id].street_id
    inventory = {s.id for s in network.sections if s.street_id == street_id}
    if allowed is not None:
        declared = set(allowed)
        if not declared <= inventory:
            raise ValueError('An interior crossing inventory must belong to its source Street')
        inventory = declared
    path, used, at = [], set(), first
    while True:
        if at.section_id in used:
            raise ValueError('A separator boundary Street has cyclic or repeated Sections')
        used.add(at.section_id)
        path.append(Traversal(at.section_id, at.end == 0))
        last = End(at.section_id, 1 - at.end)
        site = sites.get(last)
        following = [] if site is None or len(site.ends) != 2 else [e for e in site.ends
            if e != last and sections[e.section_id].street_id == street_id]
        if not following:
            break
        at, = following
    if used != inventory:
        raise ValueError('Separator boundary requires every Section of one complete, unbranched original Street run')
    return tuple(path), last, site


def interior_crossings(network, bundle):
    """Validate transverse Streets whose two ports attach to lateral members.

    The inventory covers the run between the two physical member ports, even
    when its Street continues outside the Bundle. Independent islands, branching
    interior roads and incomplete runs remain explicit unsupported cases. Valid
    crossings keep their ordinary carriageway profiles; only the two lateral
    members receive half profiles.
    """
    if not bundle.interior_section_ids:
        return ()
    sections, sites = network.section_index, endpoint_sites(network)
    roles = {t.section_id: role for role, m in enumerate(bundle.members) for t in m.path}
    inventory = {}
    for sid in bundle.interior_section_ids:
        inventory.setdefault(sections[sid].street_id, set()).add(sid)
    result = []
    for street_id, ids in sorted(inventory.items()):
        terminals = []
        for sid in sorted(ids):
            for end in (0, 1):
                endpoint = End(sid, end)
                site = sites.get(endpoint)
                if site is not None and len(site.ends) > 2:
                    member_roles = {roles[e.section_id] for e in site.ends if e.section_id in roles}
                    if member_roles:
                        if len(member_roles) != 1:
                            raise ValueError('An interior crossing port needs one lateral member side')
                        terminals.append((next(iter(member_roles)), endpoint, site))
        terminals.sort(key=lambda row: row[0])
        if not terminals:
            raise ValueError('Bundle interior Streets need two explicit opposite-member crossing ports')
        path, last, site = street_chain(network, terminals[0][1], allowed=ids)
        if len(terminals) != 2 or tuple(row[0] for row in terminals) != (0, 1):
            raise ValueError('Bundle interior Streets need two explicit opposite-member crossing ports')
        if last != terminals[1][1] or site.key != terminals[1][2].key:
            raise ValueError('Interior Street does not connect its declared opposite-member ports')
        result.append(InteriorCrossing(street_id, path, tuple(row[2].key for row in terminals)))
    return tuple(result)


def continuation_sections(network, bundle):
    """Potential complete third Streets at an interior crossing's member end.

    This conservative, local dependency closure precedes geometry evaluation.
    Ambiguous candidates are included here and rejected by the actual planner.
    """
    sites, sections = endpoint_sites(network), network.section_index
    streets = set()
    for member in bundle.members:
        a, b = member.path[0], member.path[-1]
        for endpoint in (End(a.section_id, 0 if a.forward else 1), End(b.section_id, 1 if b.forward else 0)):
            site = sites.get(endpoint)
            if site is not None and any(e.section_id in bundle.interior_section_ids for e in site.ends):
                streets.update(sections[e.section_id].street_id for e in site.ends
                               if e.section_id not in bundle.sections)
    return {s.id for s in network.sections if s.street_id in streets}


@dataclass(frozen=True)
class CrossingFreeContinuation:
    crossing: InteriorCrossing
    side: int
    end: int
    member: BundleMember
    source_street_id: int

    @classmethod
    def build(cls, network, bundle, crossing, side, end):
        """Reuse the GN checkpoint's crossing-to-free incidence policy.

        The finished lateral member has exactly one other, unowned Street at
        the crossing. Its complete source path must reach a genuine free end.
        It bounds a derived region but remains outside the original Bundle.
        """
        sections = network.section_index
        site = next(j for j in network.intersections if j.key == crossing.sites[side])
        own = bundle.members[side].path[0 if end == 0 else -1]
        endpoint = End(own.section_id, (0 if own.forward else 1) if end == 0 else (1 if own.forward else 0))
        required = {endpoint, *(e for e in site.ends if e.section_id in crossing.sections)}
        candidates = [e for e in site.ends if e not in required]
        if endpoint not in site.ends or len(required) != 2 or len(candidates) != 1:
            raise ValueError('Crossing-to-free continuation needs one unambiguous third original Street')
        first, = candidates
        street_id = sections[first.section_id].street_id
        path, last, terminal = street_chain(network, first)
        if terminal is not None:
            raise ValueError('Crossing continuation requires a genuine free original endpoint')
        ids = {t.section_id for t in path}
        if any(ids.intersection(other.sections) for other in network.bundles):
            raise ValueError('Crossing continuation has competing original Bundle ownership')
        if end == 0:
            path = tuple(Traversal(t.section_id, not t.forward) for t in reversed(path))
        return cls(crossing, side, end, BundleMember(path), street_id)

    def record(self):
        return dict(interior_street_id=self.crossing.street_id, boundary_street_id=self.source_street_id,
                    member_side=self.side, bundle_end=self.end,
                    path=[dict(section_id=t.section_id, forward=t.forward) for t in self.member.path],
                    membership_changed=False)
