"""Bounded physical-rim cuts at substantial reflex curb corners.

Each concave run contributes an inward bisector. Endpoints are existing physical
samples and every cut is verified by CompoundPartitioner before chart mapping.
This complements visibility-only splitting when its large positive blocks have
poor quad grading. No boundary/profile positions or source incidences change.
"""
import math
import numpy as np
from .compounds import rim_label
from .mesh import cross
from .partition import BoundaryRef, CutPolicy, JunctionCut


class ReflexCurbCuts(CutPolicy):
    def __init__(self, *, minimum_turn=1.3, maximum_cuts=16):
        if (not math.isfinite(minimum_turn) or not 0. < minimum_turn < math.pi
                or type(maximum_cuts) is not int or maximum_cuts < 1):
            raise ValueError('Invalid reflex-curb cut budgets')
        self.minimum_turn, self.maximum_cuts = minimum_turn, maximum_cuts

    @staticmethod
    def eligible(group):
        if (len(group.loops) != 1 or group.loops[0].hole
                or (group.physical_layout is None and (group.claims or len(group.members) != 1))):
            return False
        points = group.loops[0].xy
        edges = np.roll(points, -1, axis=0)-points
        unit = edges/np.linalg.norm(edges, axis=1)[:, None]
        turns = np.arctan2(cross(np.roll(unit, 1, axis=0), unit),
                           np.sum(np.roll(unit, 1, axis=0)*unit, axis=1))
        # A semicircular return/cap is particularly unsuitable for one radial
        # kernel. Ordinary shallow curved corners keep the normal fast path.
        mass, maximum = 0., 0.
        for angle in np.r_[turns, turns]:
            mass = mass-float(angle) if angle < -1e-6 else 0.
            maximum = max(maximum, mass)
        return maximum > 2.5

    @staticmethod
    def profile_endpoint(points, labels, edge, fraction):
        n = len(points)
        index = (edge+int(fraction > .5)) % n
        label = labels[edge]
        if label.kind != 'mouth':
            return index
        start = edge
        for _ in range(n):
            previous = labels[(start-1) % n]
            if previous.kind != 'mouth' or previous.index != label.index:
                break
            start = (start-1) % n
        else:
            raise ValueError('A physical profile must have distinct boundary endpoints')
        if (index-start) % n % 2:
            target = points[edge]+fraction*(points[(edge+1) % n]-points[edge])
            index = min(((index-1) % n, (index+1) % n),
                        key=lambda vertex: np.linalg.norm(points[vertex]-target))
        return index

    def cuts(self, group, context):
        if (len(group.loops) != 1 or group.loops[0].hole
                or (group.physical_layout is None and (group.claims or len(group.members) != 1))):
            raise ValueError('Reflex cuts require one certified physical exterior')
        loop = group.loops[0]
        points = np.asarray(loop.xy, dtype=float)
        edge = np.roll(points, -1, axis=0)-points
        lengths = np.linalg.norm(edge, axis=1)
        unit = edge/lengths[:, None]
        turns = np.arctan2(cross(np.roll(unit, 1, axis=0), unit),
                           np.sum(np.roll(unit, 1, axis=0)*unit, axis=1))
        sources = {p.region.owner: p for p in group.footprints}
        labels = [rim_label(span, sources) for span in loop.spans]
        clusters = []
        for vertex in np.flatnonzero(turns < -1e-6):
            if not clusters or vertex != clusters[-1][-1]+1:
                clusters.append([])
            clusters[-1].append(int(vertex))
        if len(clusters) > 1 and clusters[0][0] == 0 and clusters[-1][-1] == len(points)-1:
            clusters[0] = clusters.pop()+clusters[0]
        cuts, seen = [], set()
        for indices in clusters:
            mass = -float(turns[indices].sum())
            if mass < self.minimum_turn:
                continue
            index = indices[int(np.searchsorted(np.cumsum(-turns[indices]), mass*.5))]
            # The boundary is CCW. Its left normal points into pavement.
            tangent = unit[index]+unit[index-1]
            direction = np.array((-tangent[1], tangent[0]))
            size = np.linalg.norm(direction)
            if size < 1e-8:
                raise ValueError('A reversing curb has no stable inward bisector')
            direction /= size
            broadcast = np.broadcast_to(direction, edge.shape)
            determinant = cross(broadcast, edge)
            valid = abs(determinant) > 1e-12*lengths
            delta = points-points[index]
            distance = np.divide(cross(delta, edge), determinant,
                                 out=np.full(len(points), np.inf), where=valid)
            fraction = np.divide(cross(delta, broadcast), determinant,
                                 out=np.full(len(points), np.inf), where=valid)
            hits = np.flatnonzero((distance > 1e-7) & (fraction >= 0.) & (fraction <= 1.))
            if not len(hits):
                raise ValueError('A reflex curb has no inward physical boundary contact')
            target = int(hits[np.argmin(distance[hits])])
            last = self.profile_endpoint(points, labels, target, fraction[target])
            token = tuple(sorted((index, last)))
            if token in seen:
                continue
            seen.add(token)
            cuts.append(JunctionCut(None, group.members[0], group.members[0],
                BoundaryRef(0, index, 0.), BoundaryRef(0, last, 0.)))
        if not cuts:
            raise ValueError('Physical exterior has no substantial reflex curb runs')
        if len(cuts) > self.maximum_cuts:
            raise ValueError('Reflex-curb partition exceeds its cut budget')
        return tuple(cuts)
