"""Prove that an unresolved internal fillet has no exposed source curbs.

An unrelated physical median hole cannot invalidate an already covered local
corner. Conversely, mere overlap at its node does not establish coverage.
Every complete facing original curb must remain in the pavement closure and
have no positive-length coincidence with a hole rim. Farther portions of the
same original road may still contribute an unrelated exterior curb.
"""
import numpy as np
from .finite_containment import contains_segment, segment_on_polygon_side
from .region_pairs import PairBudget, SegmentPairIndex


class HiddenCornerExposureError(ValueError):
    """Identify the original unresolved pair whose physical curb stays exposed."""
    def __init__(self,message,first,last):
        self.first,self.last=first,last
        super().__init__(message)


def certify_hidden_corners(unresolved, approaches, partitions, loops, consumed, *, epsilon=1e-7):
    outer = [loop for loop in loops if not loop.hole]
    if len(outer) != 1:
        raise ValueError('Covered internal corners require one certified physical exterior')
    bounds = [(loop.xy, np.roll(loop.xy, -1, axis=0)) for loop in loops if loop.hole]
    indices = [SegmentPairIndex(np.minimum(a,b), np.maximum(a,b), epsilon) for a,b in bounds]
    budget = PairBudget(8_000_000)
    records = []
    for first, last, failure in unresolved:
        chains = []
        for endpoint, side in ((first,-1), (last,1)):
            approach = approaches[endpoint]
            for piece in approach.path.pieces:
                sid = piece.section.id
                if sid not in consumed:
                    raise ValueError('Covered internal corner escaped its complete original source interval')
                offset = piece.section.offset + .5*piece.section.width*side*(1 if piece.forward else -1)
                points = piece.path.point(partitions[sid], offset)
                if _has_exposed_interval(points, bounds, indices, budget, epsilon):
                    raise HiddenCornerExposureError('A source-cycle hole retains an unresolved original facing source curb',first,last)
                for a,b in zip(points,points[1:]):
                    if (not contains_segment(a,b,outer[0].xy,tolerance=epsilon) or
                            any(not segment_on_polygon_side(a,b,loop.xy,inside=False,tolerance=epsilon)
                                for loop in loops if loop.hole)):
                        raise HiddenCornerExposureError('An unresolved original internal fillet lacks complete facing-curb coverage',first,last)
                chains.append(dict(section=sid, original_side=1 if offset>piece.section.offset else -1,
                                   intervals=len(points)-1))
        records.append(dict(first=dict(section=first.section_id,end=first.end),
            last=dict(section=last.section_id,end=last.end), original_failure=failure,
            policy='COMPLETE_FACING_SOURCE_CURBS_EXCLUDED_FROM_HOLES', source_curbs=chains,
            exposed_hole_intervals=0, continuous_pavement_closure=True, epsilon=epsilon))
    return records


def partial_corner_candidate(first, last, pair, internal, consumed):
    """One complete internal arm and one finite original external route only."""
    if (first in internal)+(last in internal) != 1:
        return False
    for endpoint, approach in zip((first, last), pair):
        route = approach.path
        pieces = getattr(route, 'pieces', ())
        if (not pieces or not callable(getattr(route, 'trims', None))
                or hasattr(approach, 'original_approaches')
                or getattr(route, 'absorbed_sites', ())
                or getattr(route, 'mouth_piece', None) is not None
                or getattr(route, 'deferred_profile_reserve', None) is not None):
            return False
        if endpoint in internal:
            if any(piece.section.id not in consumed for piece in pieces):return False
        elif any(piece.section.id in consumed for piece in pieces):return False
    return True


def _owned_prefix(approach, mouths, consumed):
    """Use only the real complete mouth's original global source intervals."""
    result = []
    for piece in approach.path.pieces:
        sid = piece.section.id
        if sid in consumed:
            raise ValueError('A hidden partial corner cannot relabel a consumed route as an external prefix')
        key = sid, int(not piece.forward)
        mouth = mouths.get(key)
        if mouth is None:break
        first, last = ((0., mouth.station) if piece.forward else (mouth.station, piece.path.length))
        if not 0. <= first <= last <= piece.path.length:
            raise ValueError('A hidden partial corner cannot extend its real original mouth')
        if first == last:break
        result.append((piece, float(first), float(last)))
        if last-first < piece.path.length:break
    if not result:
        raise ValueError('A hidden partial corner needs an actual positive owned source prefix')
    return result


def certify_partially_hidden_corners(unresolved, approaches, partitions, loops, consumed, mouths, *, epsilon=1e-7):
    """Prove complete internal curbs and actual owned external-prefix curbs hidden.

    Every tested chord stays continuously inside the pavement closure and off
    every positive-length physical boundary interval, including its exterior.
    The former two-consumed policy is deliberately unchanged above.
    """
    outer = [loop for loop in loops if not loop.hole]
    if len(outer) != 1:
        raise ValueError('Covered partial corners require one certified physical exterior')
    bounds = [(loop.xy, np.roll(loop.xy, -1, axis=0)) for loop in loops]
    indices = [SegmentPairIndex(np.minimum(a,b), np.maximum(a,b), epsilon) for a,b in bounds]
    budget = PairBudget(8_000_000)
    records = []
    for first, last, failure in unresolved:
        internal = {end for end in (first, last) if all(
            piece.section.id in consumed for piece in approaches[end].path.pieces)}
        if not partial_corner_candidate(first, last, (approaches[first], approaches[last]), internal, consumed):
            raise ValueError('Covered partial corners need exactly one complete consumed original arm')
        chains = []
        for endpoint, side in ((first,-1), (last,1)):
            approach = approaches[endpoint]
            ranges = ([(piece, 0., piece.path.length) for piece in approach.path.pieces]
                      if endpoint in internal else _owned_prefix(approach, mouths, consumed))
            for piece, low, high in ranges:
                sid = piece.section.id
                canonical = partitions[sid]
                stations = np.r_[low, canonical[(canonical > low+1e-9) & (canonical < high-1e-9)], high]
                offset = piece.section.offset+.5*piece.section.width*side*(1 if piece.forward else -1)
                points = piece.path.point(stations, offset)
                if _has_exposed_interval(points, bounds, indices, budget, epsilon, conservative=True):
                    raise HiddenCornerExposureError(
                        'A source-cycle physical boundary retains an unresolved partial facing source curb',first,last)
                for a,b in zip(points,points[1:]):
                    if (not contains_segment(a,b,outer[0].xy,tolerance=epsilon) or
                            any(not segment_on_polygon_side(a,b,loop.xy,inside=False,tolerance=epsilon)
                                for loop in loops if loop.hole)):
                        raise HiddenCornerExposureError(
                            'An unresolved partial fillet lacks continuous owned facing-curb coverage',first,last)
                chains.append(dict(section=sid, original_side=1 if offset>piece.section.offset else -1,
                    original_interval=[low,high], intervals=len(points)-1,
                    scope='complete_consumed_source' if endpoint in internal else 'actual_owned_source_prefix'))
        records.append(dict(first=dict(section=first.section_id,end=first.end),
            last=dict(section=last.section_id,end=last.end), original_failure=failure,
            policy='COMPLETE_CONSUMED_AND_OWNED_PREFIX_CURBS_EXCLUDED_FROM_ALL_BOUNDARIES',
            source_curbs=chains, exposed_physical_boundary_intervals=0,
            continuous_pavement_closure=True, external_remainder_unchanged=True, epsilon=epsilon))
    return records


def _has_exposed_interval(points, boundaries, indices, budget, epsilon, *, conservative=False):
    first,last = points[:-1],points[1:]
    delta = last-first
    length = np.linalg.norm(delta,axis=1)
    for (a,b),index in zip(boundaries,indices):
        edge = b-a
        for rows,columns in index.batches(np.minimum(first,last),np.maximum(first,last),budget):
            d,e = delta[rows],edge[columns]
            offset = a[columns]-first[rows]
            cross = lambda x,y:x[:,0]*y[:,1]-x[:,1]*y[:,0]
            if conservative:
                # Both ends in the existing boundary band conservatively
                # exclude near-coincident short clipped arrangement edges.
                aligned = (abs(cross(offset,d)) <= epsilon*length[rows]) & (
                    abs(cross(b[columns]-first[rows],d)) <= epsilon*length[rows])
            else:
                aligned = abs(cross(d,e)) <= 1e-12*length[rows]*np.linalg.norm(e,axis=1)
                aligned &= abs(cross(offset,d)) <= epsilon*length[rows]
            if not np.any(aligned):continue
            r,c = rows[aligned],columns[aligned]
            squared = np.sum(delta[r]*delta[r],axis=1)
            t = np.sum((a[c]-first[r])*delta[r],axis=1)/squared
            u = np.sum((b[c]-first[r])*delta[r],axis=1)/squared
            overlap = np.minimum(1.,np.maximum(t,u))-np.maximum(0.,np.minimum(t,u))
            if np.any(overlap*length[r] > epsilon):return True
    return False
