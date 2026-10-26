"""Actual dimensional transition curbs in a certified angular Bundle chart."""
import numpy as np
from .bundle_boundaries import MemberBoundary


class AngularMemberBoundary(MemberBoundary):
    def __init__(self, context, member, frame, *, join_tolerance=1e-7):
        if not 0.<join_tolerance<=1e-5:
            raise ValueError('Angular curb reconciliation exceeds its ten-micrometre limit')
        self.join_tolerance = join_tolerance
        self.frame = frame
        self.transition_certificates = []
        super().__init__(context, member, None)

    def _certify_path(self, path, traversal):
        # A restriction retains its original source polynomial. Certifying the
        # whole source is conservative and never certifies a replacement fit.
        original = getattr(path, 'path', path)
        self.frame.certify(original, forward=traversal.forward)

    def _path_boundary(self, path, traversal, side):
        """Retain the source sampling phase when a dimensional seam trims it.

        A changed trim adds its actual endpoint without redistributing rows on
        the rest of the source. Exact controls stay compulsory; Smooth offsets
        retain their continuous curvature-based error bound. The whole-source
        allocation is checked before each sampling array is created.
        """
        from .mesh import cross
        original = getattr(path, 'path', path)
        first = getattr(path, 'first', 0.)
        last = first+path.length
        self._certify_path(path, traversal)
        settings = self.context.network.settings
        stations = original.samples(0., original.length, settings.spacing,
            tolerance=.1*settings.curve_tolerance, max_samples=200000)
        stations = np.r_[first, stations[(stations>first+1e-10)&(stations<last-1e-10)], last]
        q = side*(1 if traversal.forward else -1)
        points = original.point(stations, original.section.offset+q*.5*original.section.width)
        if not traversal.forward:
            points = points[::-1]
        if path.exact and len(points)>2:
            delta = np.diff(points, axis=0)
            bends = abs(cross(delta[:-1], delta[1:]))>1e-10*np.linalg.norm(delta[:-1], axis=1)*np.linalg.norm(delta[1:], axis=1)
            self.exact_corners.extend(points[1:-1][bends])
        return points

    def _certify_transition(self, curve, direction):
        lower = self.frame.certify_curve(curve, forward=direction==1)
        self.transition_certificates.append(dict(lower=lower, direction=direction,
            pieces=len(curve.segments) if hasattr(curve, 'segments') else 1))
