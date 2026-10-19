"""Discover connected original Section spans without geometric matching.

Each reached junction contributes its complete original ring. Failed rings
remain explicit frontiers; they never turn into free ends or partial rings.
This is import/topology preparation, not part of interactive GN evaluation.
"""
from collections import deque
from dataclasses import dataclass


@dataclass(frozen=True)
class PortConnectionStep:
    owner_kind: str
    junction_id: int
    via_span: tuple


@dataclass(frozen=True)
class PortConnectionIssue:
    owner_kind: str
    junction_id: int
    span_ids: tuple
    reason: str


@dataclass(frozen=True)
class PortConnections:
    plan: object
    physical_ids: tuple
    interior_ids: tuple
    added_sites: tuple[PortConnectionStep, ...]
    frontiers: tuple[PortConnectionIssue, ...]
    components: tuple
    cycle_rank: int
    seam_ids: tuple = ()


def _graph(plan):
    """Cycle rank of selected original incidences, not a mesh-hole certificate."""
    neighbors = {s.junction_id: set() for s in plan.sites}
    ports = {p.endpoint: p.junction_id for s in plan.sites for p in s.ports}
    edges = 0
    for span in plan.spans:
        if not span.in_scope:continue
        for end in (span.start, span.end):
            if end.junction_id in neighbors and (
                ports.get(end.endpoint) != end.junction_id or end.owner_kind not in ('PHYSICAL', 'INTERIOR', 'SEAM')
            ):
                raise ValueError('Selected span owner needs an exact original port incidence')
        a, b = span.start.junction_id, span.end.junction_id
        if a in neighbors and b in neighbors:
            # Count parallel spans and loops separately, even though traversal
            # only needs each neighbor once.
            edges += 1
            neighbors[a].add(b)
            neighbors[b].add(a)
    unseen = set(neighbors)
    components = []
    while unseen:
        pending = [min(unseen)]
        component = set()
        while pending:
            at = pending.pop()
            if at in component:
                continue
            component.add(at)
            pending.extend(neighbors[at] - component)
        unseen.difference_update(component)
        components.append(tuple(sorted(component)))
    return tuple(components), edges - len(neighbors) + len(components)


class PreparedPortConnections:
    def __init__(self, port_reader):
        self.reader = port_reader

    @staticmethod
    def _frontiers(plan):
        selected = {s.junction_id for s in plan.sites}
        found = {}
        for span in plan.spans:
            if not span.in_scope:continue
            for near, far in ((span.start, span.end), (span.end, span.start)):
                if near.junction_id not in selected or far.junction_id is None or far.junction_id in selected:
                    continue
                key = far.owner_kind, far.junction_id
                found.setdefault(key, set()).add(span.identity)
        return {key: tuple(sorted(values)) for key, values in sorted(found.items())}

    def expand(self, snapshot, *, physical_ids=(), interior_ids=(), seam_ids=()):
        """Follow exact owners reachable from the requested seed junctions.

        Validate seeds together; a bad explicit seed is an error. Cache each
        discovered complete ring once, then validate the complete union. A
        ring's local rejection does not suppress independent reachable routes.
        Unselected spans elsewhere on a retained whole Street do not cause a
        jump across an unsupported junction. Managed topology changes require
        a new preparation; this method neither mutates nor samples the sources.
        """
        physical_ids = self.reader._selection(physical_ids, 'Physical junction')
        interior_ids = self.reader._selection(interior_ids, 'Interior junction')
        seam_ids = self.reader._selection(seam_ids, 'Original seam')
        physical = set(physical_ids)
        interior = set(interior_ids)
        seams = set(seam_ids)
        initial = self.reader.prepare(snapshot, physical_ids=tuple(sorted(physical)),
                                      interior_ids=tuple(sorted(interior)), seam_ids=tuple(sorted(seams)))
        selected = ({('PHYSICAL', jid) for jid in physical} | {('INTERIOR', jid) for jid in interior}
                    | {('SEAM', jid) for jid in seams})
        scheduled = set(selected)
        pending = deque()
        failures = {}
        added = []

        def enqueue(plan):
            for key, spans in self._frontiers(plan).items():
                if key not in scheduled:
                    scheduled.add(key)
                    pending.append((key, spans[0]))

        enqueue(initial)
        while pending:
            (kind, jid), via = pending.popleft()
            if kind not in ('PHYSICAL', 'INTERIOR') and not (kind == 'SEAM' and self.reader.seam_connections):
                failures[(kind, jid)] = ('Original two-ended Street seam requires a native connection provider'
                    if kind == 'SEAM' else 'No exact physical ring or native interior junction record')
                continue
            try:
                local = self.reader.prepare(snapshot,
                    physical_ids=(jid,) if kind == 'PHYSICAL' else (),
                    interior_ids=(jid,) if kind == 'INTERIOR' else (),
                    seam_ids=(jid,) if kind == 'SEAM' else ())
            except ValueError as exc:
                failures[(kind, jid)] = str(exc)
                continue
            selected.add((kind, jid))
            {'PHYSICAL':physical, 'INTERIOR':interior, 'SEAM':seams}[kind].add(jid)
            added.append(PortConnectionStep(kind, jid, via))
            enqueue(local)

        # Cross-ring conflicts must reject the union, not be hidden as a local
        # unsupported case. The original reader verifies once-owned Sections,
        # ports, knot ranges and complete source Street consistency here.
        physical_ids, interior_ids = tuple(sorted(physical)), tuple(sorted(interior))
        seam_ids = tuple(sorted(seams))
        plan = self.reader.prepare(snapshot, physical_ids=physical_ids, interior_ids=interior_ids, seam_ids=seam_ids)
        frontiers = []
        for (kind, jid), spans in self._frontiers(plan).items():
            if (kind, jid) not in failures:
                raise ValueError('Connected original span has an unexplained unselected owner')
            frontiers.append(PortConnectionIssue(kind, jid, spans, failures[(kind, jid)]))
        components, cycle_rank = _graph(plan)
        return PortConnections(plan, physical_ids, interior_ids, tuple(added),
                               tuple(frontiers), components, cycle_rank, seam_ids)
