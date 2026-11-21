"""Bounded planar union of sampled physical footprints, without face meshes.

This is preparation for compound junction construction. It returns oriented
boundary loops and their original curve-segment provenance, never triangles or
an accepted pavement surface. Input loops must be simple and counterclockwise.
"""
from dataclasses import dataclass
import math
import numpy as np
from .mesh import cross
from .updates import drain_steps
from .region_pairs import PairBudget, SegmentPairIndex


class BoundaryIncidenceError(ValueError):
    """A retained boundary vertex does not have one incoming and one outgoing edge."""


def signed_area(points):
    local = points - points[0]
    return float(np.sum(cross(local, np.roll(local, -1, axis=0))) * .5)


@dataclass(frozen=True)
class Region:
    owner: object
    xy: np.ndarray
    # A consolidated curve lookup may describe only an exterior while the
    # physical source union also has holes. It is never an area operand.
    boundary_only: bool = False

    def __post_init__(self):
        xy = np.array(self.xy, dtype=float, copy=True)
        if (xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 3
                or not np.all(np.isfinite(xy)) or signed_area(xy) <= 1e-12
                or np.min(np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1)) <= 1e-8):
            raise ValueError('A region requires a finite, simple, CCW boundary with distinct vertices')
        edge = np.roll(xy, -1, axis=0) - xy
        previous = np.roll(edge, 1, axis=0)
        if np.any((np.sum(edge * previous, axis=1) < 0.)
                  & (abs(cross(edge, previous)) <= 1e-12 * np.linalg.norm(edge, axis=1) * np.linalg.norm(previous, axis=1))):
            raise ValueError('A physical footprint boundary backtracks over itself')
        xy.flags.writeable = False
        object.__setattr__(self, 'xy', xy)


@dataclass(frozen=True)
class SourceSpan:
    owner: object
    segment: int
    first: float
    last: float


@dataclass(frozen=True)
class BoundaryLoop:
    xy: np.ndarray
    spans: tuple

    @property
    def area(self):
        return signed_area(self.xy)

    @property
    def hole(self):
        return self.area < 0.


def _classify(points, polygon, epsilon):
    """Even/odd interior plus nearest coincident boundary edge, in batches."""
    return _classify_edges(points, polygon, np.roll(polygon, -1, axis=0), epsilon)


def _classify_edges(points, a, b, epsilon):
    """Ray parity and closest boundary index for a supplied set of original edges."""
    edge = b - a
    length2 = np.sum(edge * edge, axis=1)
    inside, boundary = np.empty(len(points), bool), np.full(len(points), -1, np.int32)
    # Indexing can decline unusually large polygons. Bound the largest
    # (point, edge, coordinate) temporary in that exhaustive path as well.
    batch = max(1, min(128, 32768 // max(1, len(a))))
    for start in range(0, len(points), batch):
        p = points[start:start + batch]
        delta = p[:, None] - a
        fraction = np.clip(np.sum(delta * edge, axis=-1) / length2, 0., 1.)
        distance2 = np.sum((delta - fraction[..., None] * edge) ** 2, axis=-1)
        nearest = np.argmin(distance2, axis=1)
        hit = distance2[np.arange(len(p)), nearest] <= epsilon ** 2
        boundary[start:start + len(p)] = np.where(hit, nearest, -1)
        crosses = (a[:, 1] > p[:, None, 1]) != (b[:, 1] > p[:, None, 1])
        offset = np.divide((p[:, None, 1] - a[:, 1]) * edge[:, 0], edge[:, 1],
                           out=np.zeros_like(distance2), where=crosses)
        intersections = a[:, 0] + offset
        inside[start:start + len(p)] = np.count_nonzero(crosses & (p[:, None, 0] < intersections), axis=1) % 2 == 1
    return inside, boundary


def _coincident_interval(first, last, polygon, epsilon):
    """Certify a whole segment inside the existing boundary welding tubes.

    Each finite polygon edge contributes its exact strip and endpoint disks.
    Their intervals along the source segment must cover it without a gap.
    This also covers unequal sampling across several almost-collinear chords.
    """
    vector=last-first
    length=float(np.linalg.norm(vector))
    along=vector/length
    a,b=polygon,np.roll(polygon,-1,axis=0)
    edge=b-a
    lengths=np.linalg.norm(edge,axis=1)
    direction=edge/lengths[:,None]
    delta=first-a
    slope=cross(along,direction)
    distance=cross(delta,direction)
    projection=np.sum(delta*direction,axis=1)
    speed=direction@along
    intervals=[]
    for i in range(len(polygon)):
        low,high=0.,length
        for value,change,minimum,maximum in ((distance[i],slope[i],-epsilon,epsilon),
                                           (projection[i],speed[i],0.,lengths[i])):
            if change==0.:
                if not minimum<=value<=maximum:
                    high=-1.
                    break
            else:
                p,q=(minimum-value)/change,(maximum-value)/change
                low,high=max(low,min(p,q)),min(high,max(p,q))
        if low<=high:
            intervals.append((low,high))
    # A shared vertex disk is identical for both adjacent finite edges.
    centers=-(delta@along)
    perpendicular=cross(delta,along)
    for center,distance in zip(centers,perpendicular):
        if abs(distance)<=epsilon:
            radius=math.sqrt(max(0.,epsilon**2-distance**2))
            low,high=max(0.,center-radius),min(length,center+radius)
            if low<=high:
                intervals.append((low,high))
    covered=0.
    for low,high in sorted(intervals):
        if low>covered:
            return False
        covered=max(covered,high)
        if covered>=length:
            return True
    return False


class _Nodes:
    def __init__(self, epsilon):
        self.epsilon, self.points, self.buckets = epsilon, [], {}

    def add(self, point):
        cell = tuple(np.floor(point / self.epsilon).astype(np.int64))
        for x in range(cell[0] - 1, cell[0] + 2):
            for y in range(cell[1] - 1, cell[1] + 2):
                for index in self.buckets.get((x, y), ()):
                    if np.linalg.norm(self.points[index] - point) <= self.epsilon:
                        return index
        index = len(self.points)
        self.points.append(point)
        self.buckets.setdefault(cell, []).append(index)
        return index


class RegionUnion:
    """Segment arrangement with bounded NumPy scratch and explicit provenance.

    epsilon is a numerical welding distance, not visual curve tolerance. The
    caller samples the intended curves independently before invoking this step.
    Point-only contacts need a dedicated contact pattern and are diagnosed.
    """
    def __init__(self, *, epsilon=1e-7, max_segments=8192, max_pair_tests=8_000_000):
        if (not math.isfinite(epsilon) or epsilon <= 0. or type(max_segments) is not int or max_segments < 3
                or type(max_pair_tests) is not int or max_pair_tests < 1):
            raise ValueError('Invalid region union budget')
        self.epsilon = epsilon
        self.max_segments, self.max_pair_tests = max_segments, max_pair_tests

    def _classification_candidates(self, points, bounds):
        minimum, maximum = bounds
        # This is only a broad phase. Preserve narrow-phase cancellation at
        # large coordinate ranges by conservatively overselecting its finite
        # arithmetic error; the actual classification epsilon stays unchanged.
        rounding=32.*np.finfo(float).eps*np.maximum(1.,np.maximum(abs(minimum),abs(maximum)))
        return np.all(points >= minimum-self.epsilon-rounding, axis=1) & np.all(points <= maximum+self.epsilon+rounding, axis=1)

    def _classify_batch(self, middle, endpoints, source_edges, vectors, region, regions, owner_bounds, classifier):
        keep=np.ones(len(middle),bool)
        for other in regions:
            if other.owner==region.owner:continue
            possible=np.flatnonzero(keep & self._classification_candidates(middle,owner_bounds[other.owner]))
            if not len(possible):continue
            classify = classifier if len(other.xy) >= 128 else _classify
            inside,boundary=classify(middle[possible],other.xy,self.epsilon)
            # A midpoint can be close to a shallow transverse crossing even
            # though the retained interval leaves the welding neighbourhood.
            # Owner ties describe whole coincident intervals. Use the exact
            # interior result if the interval leaves the boundary tubes;
            # otherwise two crossing exterior pieces can both survive.
            hits=np.flatnonzero(boundary>=0)
            if len(hits):
                first=other.xy[boundary[hits]]
                direction=other.xy[(boundary[hits]+1)%len(other.xy)]-first
                delta=endpoints[possible[hits]]-first[:,None,:]
                fraction=np.clip(np.sum(delta*direction[:,None,:],axis=-1)
                    /np.sum(direction*direction,axis=-1)[:,None],0.,1.)
                distance2=np.sum((delta-fraction[...,None]*direction[:,None,:])**2,axis=-1)
                for index in np.flatnonzero(np.any(distance2>self.epsilon**2,axis=1)):
                    at=hits[index]
                    if not _coincident_interval(*endpoints[possible[at]],other.xy,self.epsilon):
                        boundary[at]=-1
            keep[possible]&=~(inside & (boundary<0))
            for index in np.flatnonzero(boundary>=0):
                edge=boundary[index]
                direction=other.xy[(edge+1)%len(other.xy)]-other.xy[edge]
                at=possible[index]
                # Keep the original two-vector dot arithmetic and lower-owner
                # tie. Batching must not alter coincident span provenance.
                alignment=direction@vectors[source_edges[at]]
                if alignment<0. or (alignment>0. and other.owner<region.owner):keep[at]=False
        return keep

    def _split_pair(self, first, second, offsets, cuts, *, same=False):
        epsilon = self.epsilon
        aa, bb = first.xy, second.xy
        da, db = np.roll(aa, -1, axis=0) - aa, np.roll(bb, -1, axis=0) - bb
        la, lb = np.linalg.norm(da, axis=1), np.linalg.norm(db, axis=1)
        amin, amax = np.minimum(aa, aa + da), np.maximum(aa, aa + da)
        bmin, bmax = np.minimum(bb, bb + db), np.maximum(bb, bb + db)
        batch = max(1, min(128, 32768 // len(bb)))
        for start in range(0, len(aa), batch):
            a, d = aa[start:start + batch], da[start:start + batch]
            nearby = np.all(amax[start:start + len(a), None] + epsilon >= bmin, axis=-1)
            nearby &= np.all(bmax + epsilon >= amin[start:start + len(a), None], axis=-1)
            if same:
                row, column = np.arange(start, start + len(a))[:, None], np.arange(len(bb))[None]
                nearby &= column > row + 1
                nearby &= ~((row == 0) & (column == len(bb) - 1))
            delta = bb[None] - a[:, None]
            determinant = cross(d[:, None], db[None])
            parallel = abs(determinant) <= 1e-12 * la[start:start + len(a), None] * lb
            ta = np.divide(cross(delta, db), determinant, out=np.zeros_like(determinant), where=~parallel)
            tb = np.divide(cross(delta, d[:, None]), determinant, out=np.zeros_like(determinant), where=~parallel)
            valid = nearby & ~parallel
            valid &= (ta >= -epsilon / la[start:start + len(a), None]) & (ta <= 1. + epsilon / la[start:start + len(a), None])
            valid &= (tb >= -epsilon / lb) & (tb <= 1. + epsilon / lb)
            if same and np.any(valid):
                raise ValueError('A physical footprint boundary intersects or touches itself')
            for i, j in zip(*np.nonzero(valid)):
                cuts[offsets[0] + start + i].append(float(np.clip(ta[i, j], 0., 1.)))
                cuts[offsets[1] + j].append(float(np.clip(tb[i, j], 0., 1.)))
            collinear = nearby & parallel & (abs(cross(delta, d[:, None])) <= epsilon * la[start:start + len(a), None])
            if same and np.any(collinear):
                raise ValueError('A physical footprint boundary overlaps or touches itself')
            for i, j in zip(*np.nonzero(collinear)):
                ia = start + i
                for point in (bb[j], bb[j] + db[j]):
                    t = float((point - aa[ia]) @ da[ia] / la[ia] ** 2)
                    if -epsilon / la[ia] <= t <= 1. + epsilon / la[ia]:
                        cuts[offsets[0] + ia].append(float(np.clip(t, 0., 1.)))
                for point in (aa[ia], aa[ia] + da[ia]):
                    t = float((point - bb[j]) @ db[j] / lb[j] ** 2)
                    if -epsilon / lb[j] <= t <= 1. + epsilon / lb[j]:
                        cuts[offsets[1] + j].append(float(np.clip(t, 0., 1.)))

    def _split_pair_indexed(self, first, second, offsets, cuts, budget, *, same=False):
        epsilon = self.epsilon
        aa, bb = first.xy, second.xy
        da, db = np.roll(aa, -1, axis=0) - aa, np.roll(bb, -1, axis=0) - bb
        la, lb = np.linalg.norm(da, axis=1), np.linalg.norm(db, axis=1)
        amin, amax = np.minimum(aa, aa + da), np.maximum(aa, aa + da)
        bmin, bmax = np.minimum(bb, bb + db), np.maximum(bb, bb + db)
        index = SegmentPairIndex(bmin, bmax, epsilon)
        for row, column in index.batches(amin, amax, budget, same=same):
            a, d, b, e = aa[row], da[row], bb[column], db[column]
            delta = b - a
            determinant = cross(d, e)
            parallel = abs(determinant) <= 1e-12 * la[row] * lb[column]
            ta = np.divide(cross(delta, e), determinant, out=np.zeros_like(determinant), where=~parallel)
            tb = np.divide(cross(delta, d), determinant, out=np.zeros_like(determinant), where=~parallel)
            valid = ~parallel
            valid &= (ta >= -epsilon / la[row]) & (ta <= 1. + epsilon / la[row])
            valid &= (tb >= -epsilon / lb[column]) & (tb <= 1. + epsilon / lb[column])
            if same and np.any(valid):
                raise ValueError('A physical footprint boundary intersects or touches itself')
            for hit in np.flatnonzero(valid):
                cuts[offsets[0] + row[hit]].append(float(np.clip(ta[hit], 0., 1.)))
                cuts[offsets[1] + column[hit]].append(float(np.clip(tb[hit], 0., 1.)))
            collinear = parallel & (abs(cross(delta, d)) <= epsilon * la[row])
            if same and np.any(collinear):
                raise ValueError('A physical footprint boundary overlaps or touches itself')
            for hit in np.flatnonzero(collinear):
                i, j = row[hit], column[hit]
                for point in (bb[j], bb[j] + db[j]):
                    t = float((point - aa[i]) @ da[i] / la[i] ** 2)
                    if -epsilon / la[i] <= t <= 1. + epsilon / la[i]:
                        cuts[offsets[0] + i].append(float(np.clip(t, 0., 1.)))
                for point in (aa[i], aa[i] + da[i]):
                    t = float((point - bb[j]) @ db[j] / lb[j] ** 2)
                    if -epsilon / lb[j] <= t <= 1. + epsilon / lb[j]:
                        cuts[offsets[1] + j].append(float(np.clip(t, 0., 1.)))

    def build(self, regions):
        return drain_steps(self.build_steps(regions))

    def build_steps(self, regions):
        regions = tuple(sorted(regions, key=lambda r: r.owner))
        if any(getattr(region, 'boundary_only', False) for region in regions):
            raise ValueError('A boundary-only region cannot be a physical union operand')
        if not regions:
            return ()
        if len({r.owner for r in regions}) != len(regions):
            raise ValueError('Union region owners must be distinct')
        counts = [len(r.xy) for r in regions]
        if sum(counts) > self.max_segments:
            raise ValueError('Compound boundary exceeds its segment budget')
        offsets = np.r_[0, np.cumsum(counts)]
        cuts = [[0., 1.] for _ in range(int(offsets[-1]))]
        tested = sum(count * (count - 1) // 2 for count in counts)
        bounds = [(region.xy.min(axis=0), region.xy.max(axis=0)) for region in regions]
        owner_bounds = {region.owner: bound for region,bound in zip(regions,bounds)}
        estimate = tested
        for i in range(len(regions)):
            if estimate > self.max_pair_tests:
                break
            for j in range(i + 1, len(regions)):
                if not (np.any(bounds[i][1] + self.epsilon < bounds[j][0])
                        or np.any(bounds[j][1] + self.epsilon < bounds[i][0])):
                    estimate += counts[i] * counts[j]
                    if estimate > self.max_pair_tests:
                        break
        # Preserve the original arithmetic and dispatch for every previously
        # accepted input. Only over-budget arrangements use spatial pruning.
        budget = PairBudget(self.max_pair_tests) if estimate > self.max_pair_tests else None
        for i, region in enumerate(regions):
            yield
            if budget is None:
                self._split_pair(region, region, (offsets[i], offsets[i]), cuts, same=True)
            else:
                self._split_pair_indexed(region, region, (offsets[i], offsets[i]), cuts, budget, same=True)
        for i, first in enumerate(regions):
            for j in range(i + 1, len(regions)):
                second = regions[j]
                if (np.any(first.xy.max(axis=0) + self.epsilon < second.xy.min(axis=0))
                        or np.any(second.xy.max(axis=0) + self.epsilon < first.xy.min(axis=0))):
                    continue
                yield
                if budget is None:
                    self._split_pair(first, second, (offsets[i], offsets[j]), cuts)
                else:
                    self._split_pair_indexed(first, second, (offsets[i], offsets[j]), cuts, budget)
        # These immutable polygons can be queried by many other operands.
        # Keep the spatial index local to this one union, including cancellation.
        from .prepared_source_classifier import SourceProofClassifier
        classifier = SourceProofClassifier()
        nodes, retained = _Nodes(self.epsilon), {}
        for ri, region in enumerate(regions):
            yield
            vectors=np.roll(region.xy,-1,axis=0)-region.xy
            pending=[];middles=[];endpoints=[];source_edges=[];pending_count=0
            def retain_batch():
                nonlocal pending_count
                if not pending:return
                keep=self._classify_batch(np.concatenate(middles),np.concatenate(endpoints),np.concatenate(source_edges),
                    vectors,region,regions,owner_bounds,classifier)
                offset=0
                for si,a,b,intervals in pending:
                    selected=keep[offset:offset+len(intervals)];offset+=len(intervals)
                    for first,last in intervals[selected]:
                        u,v=nodes.add(a+first*(b-a)),nodes.add(a+last*(b-a))
                        if u!=v:retained.setdefault((u,v),SourceSpan(region.owner,si,float(first),float(last)))
                pending.clear();middles.clear();endpoints.clear();source_edges.clear();pending_count=0
            for si, (a, b) in enumerate(zip(region.xy, np.roll(region.xy, -1, axis=0))):
                length = np.linalg.norm(b - a)
                fractions = sorted(cuts[offsets[ri] + si])
                unique = [fractions[0]]
                for t in fractions[1:]:
                    if (t - unique[-1]) * length > self.epsilon:
                        unique.append(t)
                intervals = np.column_stack((unique[:-1], unique[1:]))
                # Bound temporary classification scratch independently of the
                # source arrangement size. Retention order stays edge/fraction
                # order, including when one source edge crosses batch limits.
                for first in range(0,len(intervals),8192):
                    part=intervals[first:first+8192]
                    if pending_count+len(part)>8192:
                        yield
                        retain_batch()
                    pending.append((si,a,b,part))
                    middles.append(a+part.mean(axis=1)[:,None]*(b-a))
                    endpoints.append(a+part[...,None]*(b-a))
                    source_edges.append(np.full(len(part),si,np.intp));pending_count+=len(part)
            retain_batch()
        outgoing, incoming = {}, {}
        for u, v in retained:
            outgoing.setdefault(u, []).append(v)
            incoming.setdefault(v, []).append(u)
        if any(len(outgoing.get(i, ())) != 1 or len(incoming.get(i, ())) != 1 for i in outgoing.keys() | incoming.keys()):
            raise BoundaryIncidenceError('Compound boundary has a point-only contact or unresolved numerical incidence')
        pending, loops = set(retained), []
        while pending:
            first = min(pending)
            edge, ids, spans = first, [], []
            while True:
                if edge not in pending:
                    raise ValueError('Compound boundary does not close into disjoint loops')
                pending.remove(edge)
                ids.append(edge[0])
                spans.append(retained[edge])
                edge = (edge[1], outgoing[edge[1]][0])
                if edge == first:
                    break
            xy = np.array([nodes.points[i] for i in ids])
            start = min(range(len(xy)), key=lambda i: tuple(xy[i]))
            xy = np.roll(xy, -start, axis=0)
            xy.flags.writeable = False
            spans = tuple(spans[start:] + spans[:start])
            if abs(signed_area(xy)) <= self.epsilon ** 2:
                raise ValueError('Compound boundary encloses no physical area')
            loops.append(BoundaryLoop(xy, spans))
        return tuple(sorted(loops, key=lambda loop: (loop.hole, tuple(loop.xy[0]))))
