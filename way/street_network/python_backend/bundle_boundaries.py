"""Source-owned Bundle curbs, including certified dimensional Section seams.

The collector restricts the original metric paths and inserts the transition's
actual quintic curbs. It never interpolates new source controls, substitutes
virtual Sections, or infers a topological connection from coincident positions.
"""
import numpy as np
from .model import End, Key


def internal_sites(network, bundle):
    """Consecutive member incidences, with their explicitly declared seam."""
    sites = {e: j for j in (*network.intersections, *network.connections) for e in j.ends}
    result = []
    for member in bundle.members:
        for previous, following in zip(member.path, member.path[1:]):
            exit = End(previous.section_id, 1 if previous.forward else 0)
            entry = End(following.section_id, 0 if following.forward else 1)
            site = sites.get(exit)
            if site is not None and entry in site.ends and len(site.ends) == 2:
                result.append(site)
    return tuple(result)


def dimensional_sites(network, bundle):
    from .transitions import required
    return tuple(site for site in internal_sites(network, bundle) if required(network, site))


def candidate_owner(network, site):
    """Potential composite owner; its evaluated plan must also claim the seam."""
    return next((bundle for bundle in network.bundles
                 if any(j.key == site.key for j in dimensional_sites(network, bundle))), None)


def _certify_increasing(coefficients, *, coordinate='axial'):
    """Positive Bernstein derivative hull, with a bounded subdivision proof."""
    from .bundle_axis import AxialChartUnavailable
    pending, checks, lower = [(np.asarray(coefficients), 0)], 0, float('inf')
    while pending:
        values, depth = pending.pop()
        checks += 1
        if checks > 4096:
            raise AxialChartUnavailable(f'Bundle transition {coordinate} proof exceeded its budget')
        if values.min() > 1e-8:
            lower = min(lower, float(values.min()))
            continue
        if values.max() <= 1e-8 or depth >= 16:
            raise AxialChartUnavailable(f'Bundle transition curb is not continuously {coordinate}-monotone')
        work, left, right = values, [values[0]], [values[-1]]
        while len(work) > 1:
            work = .5 * (work[:-1] + work[1:])
            left.append(work[0]); right.append(work[-1])
        pending.extend(((np.asarray(left), depth + 1), (np.asarray(right[::-1]), depth + 1)))
    return lower


class MemberBoundary:
    def __init__(self, context, member, axis):
        self.context, self.member, self.axis = context, member, axis
        self.exact_corners = []
        self.maximum_join_error = 0.
        self.sites = {e: j for j in (*context.network.intersections, *context.network.connections) for e in j.ends}
        self.transitions = {}
        from .transitions import required
        from .bundle_axis import AxialChartUnavailable
        for previous, following in zip(member.path, member.path[1:]):
            exit = End(previous.section_id, 1 if previous.forward else 0)
            entry = End(following.section_id, 0 if following.forward else 1)
            site = self.sites.get(exit)
            if site is not None and entry in site.ends and required(context.network, site):
                transition = context[site.layout_key].transition
                if transition is not None:
                    self.transitions[(previous.section_id, following.section_id)] = transition
        if self.transitions:
            for previous, following in zip(member.path, member.path[1:]):
                exit = End(previous.section_id, 1 if previous.forward else 0)
                entry = End(following.section_id, 0 if following.forward else 1)
                site = self.sites.get(exit)
                if site is None or entry not in site.ends or len(site.ends) != 2:
                    raise AxialChartUnavailable('Bundle transition chains need explicit two-ended Section seams')

    @property
    def owners(self):
        return tuple(sorted({plan.owner for plan in self.transitions.values()}))

    def _certify_path(self, path, traversal):
        from .bundle_axis import AxialChartUnavailable
        if not path.exact:
            # Same exact quadratic derivative minimum used by the original
            # constant-width chart. Offset regularity is certified by Path.
            coefficients = path.coefficients @ self.axis * (1 if traversal.forward else -1)
            a, b, c = coefficients[:, 1], coefficients[:, 2], coefficients[:, 3]
            critical = np.divide(-b, 3. * c, out=np.zeros_like(b), where=abs(c) > 1e-15)
            middle = a + critical * (2. * b + 3. * critical * c)
            minimum = np.minimum(a, a + 2. * b + 3. * c)
            minimum = np.minimum(minimum, np.where((critical > 0.) & (critical < 1.), middle, np.inf))
            if np.min(minimum) <= 1e-8:
                raise AxialChartUnavailable('Cubic axial derivative is not continuously positive')

    def _path_boundary(self, path, traversal, side):
        from .bundle_axis import _simplify
        settings = self.context.network.settings
        source = path.section
        q = side * (1 if traversal.forward else -1)
        self._certify_path(path, traversal)
        stations = path.samples(0., path.length, settings.spacing, tolerance=settings.curve_tolerance * .1)
        controls = np.interp(path.knots, path.parameters, path.stations)
        stations = np.unique(np.r_[stations, controls])
        xy = path.point(stations, source.offset + q * .5 * source.width)
        protected = np.searchsorted(stations, controls)
        if not traversal.forward:
            xy = xy[::-1]
            protected = len(xy) - 1 - protected
        tolerance = 1e-9 if path.exact else settings.curve_tolerance * .1
        points = _simplify(xy, tolerance, protected)
        if path.exact and len(points) > 2:
            from .mesh import cross
            delta = np.diff(points, axis=0)
            bends = abs(cross(delta[:-1], delta[1:])) > 1e-10 * np.linalg.norm(delta[:-1], axis=1) * np.linalg.norm(delta[1:], axis=1)
            self.exact_corners.extend(points[1:-1][bends])
        return points

    def _restricted_path(self, traversal, index):
        from .bundle_ends import PathRange
        from .bundle_axis import AxialChartUnavailable
        view = self.context[Key('path', traversal.section_id)]
        path = view.path if isinstance(view, PathRange) else view
        first = view.first if isinstance(view, PathRange) else 0.
        last = first + view.length
        neighbors = self.member.path
        pairs = ((neighbors[index-1], traversal),) if index else ()
        pairs += ((traversal, neighbors[index+1]),) if index+1 < len(neighbors) else ()
        for a, b in pairs:
            transition = self.transitions.get((a.section_id, b.section_id))
            if transition is None:
                continue
            own = transition.source_index(traversal.section_id)
            approach, distance = transition.approaches[own], transition.distances[own]
            if approach.end == 0:
                first = max(first, distance)
            else:
                last = min(last, path.length - distance)
        if last <= first + 1e-7:
            raise AxialChartUnavailable('Bundle dimensional transitions consume a complete Section; shared transition ownership is required')
        return PathRange(path, first, last) if first > 1e-8 or last < path.length-1e-8 else path

    def _certify_transition(self, curve, direction):
        if hasattr(curve, 'segments'):
            for i,piece in enumerate(curve.segments):
                span = curve.knots[i+1]-curve.knots[i]
                _certify_increasing(piece.first_derivative @ self.axis * direction/span)
        else:
            _certify_increasing(curve.first_derivative @ self.axis * direction)

    def _transition_boundary(self, transition, previous, side):
        from .bundle_axis import AxialChartUnavailable
        own = transition.source_index(previous.section_id)
        # q of the first approach is reversed by SectionTransition.surface.
        q = side if own == 0 else -side
        curve = transition.curbs[int(q == 1)]
        direction = 1 if own == 0 else -1
        settings = self.context.network.settings
        self._certify_transition(curve, direction)
        if hasattr(curve, 'segments'):
            t = curve.sampling_parameters(settings.spacing,settings.curve_tolerance*.1,max_intervals=200000)
        else:
            count = curve.intervals(settings.spacing, settings.curve_tolerance * .1)
            if count > 200_000:
                raise AxialChartUnavailable('Bundle transition boundary exceeds its sampling budget')
            t = np.linspace(0., 1., count+1)
        return curve(t if direction == 1 else t[::-1])

    def points(self, side):
        from .bundle_axis import AxialChartUnavailable
        pieces = []

        def append(xy):
            if pieces:
                error = float(np.linalg.norm(pieces[-1][-1] - xy[0]))
                self.maximum_join_error = max(self.maximum_join_error,error)
                if error > getattr(self,'join_tolerance',1e-7):
                    raise AxialChartUnavailable('Axial Bundle member needs compatible explicit Section seams')
                xy = xy[1:]
            pieces.append(xy)

        for index, traversal in enumerate(self.member.path):
            path = self._restricted_path(traversal, index) if self.transitions else self.context[Key('path', traversal.section_id)]
            append(self._path_boundary(path, traversal, side))
            if index+1 < len(self.member.path):
                following = self.member.path[index+1]
                transition = self.transitions.get((traversal.section_id, following.section_id))
                if transition is not None:
                    append(self._transition_boundary(transition, traversal, side))
        return np.concatenate(pieces)
