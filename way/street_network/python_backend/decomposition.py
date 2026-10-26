"""Bounded visibility cuts for an indented simple boundary, without triangles."""
import math
import numpy as np
from .mesh import cross
from .quad_charts import visibility_center
from .regions import _classify, signed_area


class VisibilitySplitter:
    def __init__(self, *, epsilon=1e-7, max_candidates=768, max_cuts=16):
        if (not math.isfinite(epsilon) or epsilon <= 0. or type(max_candidates) is not int
                or max_candidates < 1 or type(max_cuts) is not int or max_cuts < 1):
            raise ValueError('Invalid visibility decomposition budgets')
        self.epsilon, self.max_candidates, self.max_cuts = epsilon, max_candidates, max_cuts

    def visible(self, polygon, first, last):
        a, b = polygon[first], polygon[last]
        d = b - a
        length = np.linalg.norm(d)
        if length <= self.epsilon:
            return False
        edge = np.roll(polygon, -1, axis=0) - polygon
        lengths = np.linalg.norm(edge, axis=1)
        determinant = cross(d, edge)
        valid = abs(determinant) > 1e-12 * length * lengths
        t = np.divide(cross(polygon - a, edge), determinant, out=np.zeros(len(edge)), where=valid)
        u = np.divide(cross(polygon - a, d), determinant, out=np.zeros(len(edge)), where=valid)
        internal = (t * length > self.epsilon) & ((1. - t) * length > self.epsilon)
        if np.any(valid & internal & (u * lengths >= -self.epsilon) & ((1. - u) * lengths >= -self.epsilon)):
            return False
        for at in np.flatnonzero(~valid & (abs(cross(polygon - a, d)) <= self.epsilon * length)):
            ends = np.sort(np.array(((polygon[at] - a) @ d, (polygon[at] + edge[at] - a) @ d)) / (d @ d))
            if min(1., ends[-1]) - max(0., ends[0]) > self.epsilon / length:
                return False
        hit, boundary = _classify(np.array(((a + b) * .5,)), polygon, self.epsilon * .1)
        return bool(hit[0] and boundary[0] < 0)

    @staticmethod
    def cut_angle(polygon, first, last):
        angles = []
        for at, other in ((first, last), (last, first)):
            following = polygon[(at + 1) % len(polygon)] - polygon[at]
            previous = polygon[at - 1] - polygon[at]
            chord = polygon[other] - polygon[at]
            angles.extend((math.atan2(float(cross(following, chord)), float(following @ chord)) % math.tau,
                           math.atan2(float(cross(chord, previous)), float(chord @ previous)) % math.tau))
        return min(angles)

    def split(self, boundary, *, preferred=()):
        """Return noncrossing diagonals between existing boundary vertices.

        A cut separates complete polygonal blocks. Stop as soon as each block
        has a strict visibility kernel; neither triangle faces nor a triangular
        intermediate mesh are created. Every split keeps both boundary chains.
        """
        boundary = np.asarray(boundary, dtype=float)
        if (boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 4
                or not np.all(np.isfinite(boundary)) or signed_area(boundary) <= 0.):
            raise ValueError('Visibility splitting needs a finite CCW boundary')
        cache = {}

        def kernel(indices):
            indices = tuple(indices)
            if indices not in cache:
                try:
                    polygon = boundary[list(indices)]
                    center = visibility_center(polygon)
                    edge = np.roll(polygon, -1, axis=0) - polygon
                    clearance = float(np.min(cross(edge, center - polygon) / np.linalg.norm(edge, axis=1)))
                    cache[indices] = clearance
                except ValueError:
                    cache[indices] = None
            return cache[indices]

        pending = [tuple(range(len(boundary)))]
        cuts, attempts = [], 0
        while pending:
            indices = pending.pop()
            if kernel(indices) is not None:
                continue
            if len(cuts) >= self.max_cuts:
                raise ValueError('Visibility decomposition exceeded its cut budget')
            polygon = boundary[list(indices)]
            edge = np.roll(polygon, -1, axis=0) - polygon
            turning = cross(np.roll(edge, 1, axis=0), edge)
            reflex = np.flatnonzero(turning < -self.epsilon ** 2)
            if not len(reflex):
                raise ValueError('Nonvisible boundary has no usable reflex cut')
            priority = [i for i in reflex if indices[i] in preferred]
            ordered = sorted(reflex, key=lambda i: (turning[i] / (np.linalg.norm(edge[i]) * np.linalg.norm(edge[i - 1])), int(i)))
            starts = list(dict.fromkeys((*priority, *ordered[:16])))
            ends = np.unique(np.r_[np.linspace(0, len(indices) - 1, min(32, len(indices)), dtype=int), starts])
            best, partial, qualified = None, None, 0
            area = signed_area(polygon)
            for i in starts:
                candidates = sorted(ends, key=lambda j: (np.linalg.norm(polygon[i] - polygon[j]), int(j)))
                for j in candidates:
                    # Retain real blocks on both sides, avoiding a succession
                    # of tiny ears beside a densely sampled curved boundary.
                    distance = (int(j) - int(i)) % len(indices)
                    if min(distance, len(indices) - distance) < 3:
                        continue
                    attempts += 1
                    if attempts > self.max_candidates:
                        break
                    if not self.visible(polygon, int(i), int(j)):
                        continue
                    angle = self.cut_angle(polygon, int(i), int(j))
                    if angle < math.radians(5.):
                        continue
                    first = tuple(indices[(int(i) + k) % len(indices)] for k in range(distance + 1))
                    second = tuple(indices[(int(j) + k) % len(indices)] for k in range(len(indices) - distance + 1))
                    a, b = signed_area(boundary[list(first)]), signed_area(boundary[list(second)])
                    if min(a, b) <= self.epsilon ** 2:
                        continue
                    ka, kb = kernel(first), kernel(second)
                    chord = float(np.linalg.norm(polygon[i] - polygon[j]))
                    if ka is not None and kb is not None:
                        score = (min(ka, kb) / math.sqrt(area) + .25 * math.sin(min(math.pi * .5, angle))
                                 + .2 * min(a, b) / area - .01 * chord / math.sqrt(area))
                        option = score, (indices[int(i)], indices[int(j)]), first, second
                        if best is None or option[0] > best[0]:
                            best = option
                        qualified += 1
                        if qualified >= 8:
                            break
                    elif ka is not None or kb is not None:
                        score = (a if ka is not None else b) / area - .01 * chord / math.sqrt(area)
                        if partial is None or score > partial[0]:
                            partial = score, (indices[int(i)], indices[int(j)]), first, second
                if qualified >= 8 or attempts > self.max_candidates:
                    break
            option = best if best is not None else partial
            if option is None:
                raise ValueError('Visibility decomposition needs another compatible cut family')
            _, cut, first, second = option
            cuts.append(cut)
            pending.extend((first, second))
        return tuple(cuts)
