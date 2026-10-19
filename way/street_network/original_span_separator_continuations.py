"""Explicit candidate regions beyond an original crossing-capped Bundle.

One original member may finish at the crossing while its partner continues to
a genuine free end. At the finished member's original junction, a unique third
road can provide the other boundary of a further region. This opt-in planner
records that derived relationship without changing the original Bundle or
claiming geometry, profile compatibility, or complete coverage.
"""
from dataclasses import dataclass

from .original_span_separator_regions import CrossingSeparatorRegionPlanner
from .physical_ports import SectionEnd


@dataclass(frozen=True)
class CrossingFreeSeparatorRegion:
    bundle_id: int
    original_member_street_ids: tuple[int, int]
    boundary_street_ids: tuple[int, int]
    boundary_span_ids: tuple[tuple[int, int, int], tuple[int, int, int]]
    continuing_member_side: int
    attached_sites: tuple[tuple[int, ...], tuple[int, ...]]
    attached_ends: tuple[SectionEnd, SectionEnd]
    free_ends: tuple[SectionEnd, SectionEnd]
    crossing_span_id: tuple[int, int, int]
    crossing_ends: tuple[SectionEnd, SectionEnd]


@dataclass(frozen=True)
class CrossingFreePreparation:
    regions: tuple[CrossingFreeSeparatorRegion, ...]
    pending: tuple[dict, ...]


class CrossingFreeRegionPlanner(CrossingSeparatorRegionPlanner):
    """Find a unique original continuation; native geometry remains mandatory.

    This is a new, explicit preparation policy, not an automatic change to the
    established region dispatcher. It uses the original incidence ring rather
    than names, coordinate proximity, apparent parallelism, or sample IDs.
    The new road's current profile must not silently be replaced by a Bundle
    half-profile; native profile reconciliation is a separate construction gate.
    """

    def prepare_continuations(self, snapshot, component, *, prefixes):
        """Revalidate supplied prefixes and retain every rejected extent explicitly."""
        prefixes=tuple(prefixes);ids=[row.bundle_id for row in prefixes]
        if len(ids)!=len(set(ids)):
            raise ValueError('Continuation preparation needs unique original prefix declarations')
        regions=[];pending=[]
        for prefix in prefixes:
            bid=prefix.bundle_id
            if type(bid) is not int or bid not in snapshot.bundles:
                raise ValueError('Continuation preparation lost its original Bundle')
            try:
                verified=super().candidate(snapshot.bundles[bid],snapshot,component)
                if verified!=prefix:
                    raise ValueError('Supplied crossing prefix differs from current original topology')
                if not any(prefix.remaining_member_span_ids):
                    continue
                regions.append(self.continuation_candidate(snapshot.bundles[bid],snapshot,component))
            except ValueError as exc:
                pending.append(dict(bundle_id=bid,reason=str(exc)))
        return CrossingFreePreparation(tuple(regions),tuple(pending))

    @staticmethod
    def _identity(route):
        return (route.span.street_id, route.span.section_ids[0], route.span.section_ids[-1])

    def _free_route(self, route, attached_junction):
        if route.kind != 'APPROACH' or not route.span.in_scope:
            raise ValueError('A crossing-to-free boundary needs an in-scope original approach')
        ends = tuple(self.route_end(route, at, allow_free=True) for at in (0, 1))
        fixed = [i for i, end in enumerate(ends) if end.junction_id == attached_junction
                 and end.owner_kind != 'FREE']
        if len(fixed) != 1 or ends[1-fixed[0]].owner_kind != 'FREE':
            raise ValueError('Boundary must attach at the exact crossing and have a genuine free end')
        at = fixed[0]
        return ends[at].endpoint, ends[1-at].endpoint

    def continuation_candidate(self, bundle, snapshot, component):
        # Revalidate the whole Bundle, including all suffix routes and source
        # extent records. A caller cannot pass an unverified prefix declaration.
        prefix = super().candidate(bundle, snapshot, component)
        remaining = prefix.remaining_member_span_ids
        sides = [i for i, spans in enumerate(remaining) if spans]
        if len(sides) != 1 or len(remaining[sides[0]]) != 1:
            raise ValueError('This continuation branch requires exactly one remaining original approach')
        side = sides[0]; other = 1-side
        routes = component.network.routes
        selected = [r for r in routes if self._identity(r) == remaining[side][0]]
        if len(selected) != 1:
            raise ValueError('Remaining original approach is missing or repeated')
        rest = selected[0]
        crossing, = [r for r in routes if self._identity(r) == prefix.crossing_span_id]
        cross_ends = (crossing.span.start, crossing.span.end)
        ordered_cross = [next(end for end in cross_ends if end.endpoint == key)
                         for key in prefix.crossing_ends]
        fixed, free = self._free_route(rest, ordered_cross[side].junction_id)

        jid = ordered_cross[other].junction_id
        sites = [s for s in component.ports.sites if s.junction_id == jid]
        if len(sites) != 1:
            raise ValueError('Continuation needs one recorded original physical junction')
        ports = sites[0].ports
        endpoints = [p.endpoint for p in ports]
        required = {prefix.tail_ends[other], prefix.crossing_ends[other]}
        if (len(endpoints) != len(set(endpoints)) or not required <= set(endpoints)
                or any(p.junction_id != jid for p in ports)):
            raise ValueError('Original junction lost or duplicated the member/crossing endpoints')
        candidates = [p for p in ports if p.endpoint not in required]
        if len(candidates) != 1:
            raise ValueError('Crossing continuation must be unambiguous in the original incidence ring')
        endpoint = candidates[0].endpoint
        sid = endpoint.street_id
        if (sid in (*prefix.member_street_ids, prefix.crossing_span_id[0]) or
                sid not in component.ports.source_streets or sid not in snapshot.streets):
            raise ValueError('Continuation must be a separate original source Street in this component')
        source = snapshot.streets[sid]
        if source.get('bundle_id') is not None:
            raise ValueError('Continuation already belongs to a Bundle and needs a separate ownership policy')
        # Check every Bundle inventory too: a missing Street.bundle_id must not
        # conceal a competing complete or partial member declaration.
        for candidate in snapshot.bundles.values():
            declared = {m['street_id'] for m in candidate.get('topology', {}).get('members', ())}
            if sid in declared or sid in candidate.get('owned_streets', ()):
                raise ValueError('Continuation has competing original Bundle ownership')
        sections = [s for s in source['items'] if s['kind'] == 'Section']
        if not sections or any(s.get('category') not in self.categories for s in sections):
            raise ValueError('Every continuation Section must be in the declared carriageway scope')
        next_route = self.complete_route(dict(street_id=sid), snapshot, component, 'APPROACH')
        new_fixed, new_free = self._free_route(next_route, jid)
        if new_fixed != endpoint:
            raise ValueError('Continuation approach does not match its original junction port')
        ids = list(prefix.member_street_ids); ids[other] = sid
        spans = [None, None]; spans[side] = self._identity(rest); spans[other] = self._identity(next_route)
        fixed_ends = [None, None]; fixed_ends[side] = fixed; fixed_ends[other] = new_fixed
        free_ends = [None, None]; free_ends[side] = free; free_ends[other] = new_free
        return CrossingFreeSeparatorRegion(prefix.bundle_id, prefix.member_street_ids,
            tuple(ids), tuple(spans), side, prefix.tail_sites, tuple(fixed_ends),
            tuple(free_ends), prefix.crossing_span_id, prefix.crossing_ends)
