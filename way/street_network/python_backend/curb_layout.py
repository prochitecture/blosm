"""Finite source leads for a through curb whose carriageway width changes.

Circular fillets retain their contacts. Only an inconclusive/failed width
continuation may request longer leads; no source point, profile or radius moves.
The ownership planner still has to reconcile any newly overlapping footprints.
"""
from .continuations import ContinuationReachError


def finite_boundaries(junction, approaches, corners, stations, settings, *, maximum_steps=16):
    from .intersections import BoundaryCurve, IntersectionProvider, ApproachExtentError
    stations = list(stations)
    count = len(stations)
    for step in range(maximum_steps + 1):
        boundaries, failed = [], []
        for i, corner in enumerate(corners):
            try:
                boundaries.append(BoundaryCurve(corner,
                    (stations[i], stations[(i + 1) % count]), settings.curve_tolerance))
            except ContinuationReachError:
                if not corner.continuation:
                    raise
                failed.append(i)
        if not failed:
            return stations, boundaries
        if step == maximum_steps:
            raise ContinuationReachError('Width continuation exhausts its bounded source-lead budget')
        increments = {}
        for i in failed:
            corner = corners[i]
            reach = max(corner.radius, .5*max(corner.first.width, corner.second.width))
            increment = max(.125, reach/8.)
            for index in (i, (i + 1) % count):
                increments[index] = max(increments.get(index, 0.), increment)
        exhausted = []
        for index, increment in increments.items():
            approach = approaches[index]
            candidate = IntersectionProvider.mouth_station(approach, stations[index] + increment, settings)
            if not IntersectionProvider.mouth_fits(approach, candidate):
                exhausted.append((approach, candidate))
            else:
                stations[index] = candidate
        if exhausted:
            # The caller may form a source-owned compound, but may never
            # extrapolate a real neighboring site or known scope boundary.
            raise ApproachExtentError(junction, exhausted)
