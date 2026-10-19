"""Retain explicit junctions when minor-node concatenation would close a loop."""
from dataclasses import dataclass


@dataclass(frozen=True)
class MinorStreetCycle:
    streets: tuple
    intersections: tuple


class MinorCycleReview:
    """Find directed through-road cycles before any item links are changed.

    Minor classification normally allows a through Street to own its junction
    as an internal item. A closed chain has no first or last Section, however;
    applying the linear-list mutators creates a cyclic list. Keep its original
    Streets and the existing complete connector rings as ordinary junctions.
    This is a topology representation decision, not a claim of GN support.
    """
    def __init__(self, manager):
        self.manager = manager

    def find(self):
        streets = {row[3] for row in self.manager.waymap.iterEdges()}
        following = {}
        for node in self.manager.minorIntersections.values():
            if not node.isMinor:
                continue
            a, b = node.arriving, node.leaving
            if a not in streets or b not in streets:
                continue
            if a.dst != node.location or b.src != node.location:
                continue
            if a in following:
                raise ValueError('Ambiguous minor through-road successor')
            following[a] = (b, node)
        finished, cycles = set(), []
        for seed in following:
            current, path, indices = seed, [], {}
            while current in following and current not in finished:
                if current in indices:
                    members = tuple(path[indices[current]:])
                    cycles.append(MinorStreetCycle(
                        members, tuple(following[s][1] for s in members)))
                    break
                indices[current] = len(path)
                path.append(current)
                current = following[current][0]
            finished.update(path)
        return tuple(cycles)

    def apply(self):
        cycles = self.find()
        # Validate every ring before changing classification. transformToMinor()
        # retains this original ring; its copied side lists are not sufficient.
        for cycle in cycles:
            for node in cycle.intersections:
                connectors, visited, current = [], set(), node.startConnector
                for _ in node.leaveWays:
                    if current is None or id(current) in visited:
                        break
                    visited.add(id(current))
                    connectors.append(current)
                    current = current.succ
                if (len(connectors) != len(node.leaveWays) or not connectors
                        or current is not node.startConnector
                        or {(c.item, c.leaving) for c in connectors}
                        != {(way.street, way.leaving) for way in node.leaveWays}):
                    raise ValueError('A cyclic minor junction has no complete original connector ring')
                if any(c.intersection is not node or c.pred is None
                       or c.succ.pred is not c or c.pred.succ is not c for c in connectors):
                    raise ValueError('A cyclic minor junction has an inconsistent connector ring')
        records = list(getattr(self.manager, 'minorCycleAudit', ()))
        for cycle in cycles:
            for node in cycle.intersections:
                node.isMinor = False
                self.manager.minorIntersections.pop(node.location)
                self.manager.majorIntersections[node.location] = node
            records.append(dict(street_ids=[s.id for s in cycle.streets],
                                intersection_ids=[n.id for n in cycle.intersections],
                                treatment='RETAIN_EXPLICIT_JUNCTIONS'))
        self.manager.minorCycleAudit = records
        return cycles
