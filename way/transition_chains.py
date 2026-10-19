"""Plan and assemble lane-transition Streets without changing the input graph
while discovering its chains. Section identity and graph ownership stay unique.
"""
from dataclasses import dataclass
import re

from .item import Intersection, SideLane, Street, SymLane


@dataclass(frozen=True)
class TransitionJoin:
    intersection: object
    arriving: object
    leaving: object
    has_turns: bool


@dataclass(frozen=True)
class TransitionChain:
    streets: tuple
    joins: tuple


class TransitionStreetAssembler:
    """Join directed two-way nodes; retain cyclic or incompatible topology.

    Discovery uses the original graph and LeavingWay references, never the
    Section.street back references that assembly necessarily changes. All chains
    are planned before any Section list or graph edge is rewritten.
    """
    def __init__(self, waymap, do_pixel_correction):
        self.waymap = waymap
        self.do_pixel_correction = do_pixel_correction
        self.edges = tuple(waymap.edges(data='object', keys=True))
        self.streets = tuple(dict.fromkeys(row[3] for row in self.edges))

    def plan(self):
        available = set(self.streets)
        following, preceding = {}, {}
        for location, node in self.waymap.iterNodes(Intersection):
            if len(node.leaveWays) != 2:
                continue
            a, b = (way.street for way in node.leaveWays)
            if a is b or a not in available or b not in available:
                continue
            if a.dst == location and b.src == location:
                arriving, leaving = a, b
            elif b.dst == location and a.src == location:
                arriving, leaving = b, a
            else:
                continue
            if arriving in following or leaving in preceding:
                raise ValueError('Ambiguous directed lane-transition chain')
            join = TransitionJoin(node, arriving, leaving,
                                  bool(re.search(r'[^N]', leaving.head.lanePatterns[0])))
            following[arriving] = join
            preceding[leaving] = join

        chains, visited = [], set()
        for seed in self.streets:
            if seed in visited:
                continue
            current, backwards = seed, set()
            while current in preceding and current not in backwards:
                backwards.add(current)
                current = preceding[current].arriving
            if current in backwards:
                # No distinguished start/end: keep the original edges/nodes.
                # Never mutate the cycle and then pretend to have undone it.
                visited.update(backwards)
                continue
            members, joins = [current], []
            while current in following:
                join = following[current]
                joins.append(join)
                current = join.leaving
                members.append(current)
            visited.update(members)
            if len(members) > 1:
                chains.append(TransitionChain(tuple(members), tuple(joins)))
        return tuple(chains)

    def assemble(self, symmetric, sided):
        chains = self.plan()
        replacements, result, removed = {}, [], []
        for chain in chains:
            combined = Street(chain.streets[0].src, chain.streets[-1].dst)
            combined.insertStreetEnd(chain.streets[0])
            combined.pred = chain.streets[0].pred
            for join, street in zip(chain.joins, chain.streets[1:]):
                if join.has_turns:
                    transition = SideLane(self.do_pixel_correction, join.intersection.location,
                                          join.arriving.tail, join.leaving.head)
                    sided.append(transition)
                else:
                    transition = SymLane(join.intersection.location,
                                         join.arriving.tail, join.leaving.head)
                    symmetric.append(transition)
                combined.insertEnd(transition)
                combined.insertStreetEnd(street)
                removed.append(join.intersection)
            replacements.update((street, combined) for street in chain.streets)
            result.append(combined)

        # An outer junction can meet several new chains. Remap each incidence
        # to its own replacement, never to the last chain visiting that node.
        for _location, node in self.waymap.iterNodes(Intersection):
            for way in node.leaveWays:
                way.street = replacements.get(way.street, way.street)
        for src, dst, key, street in self.edges:
            if street in replacements:
                self.waymap.remove_edge(src, dst, key)
        for node in removed:
            if self.waymap.degree(node.location):
                raise ValueError('Lane-transition node retains an unexpected graph incidence')
            self.waymap.removeNode(node.location)
        for street in result:
            self.waymap.addEdge(street)
        return tuple(result)
