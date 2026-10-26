"""Certified angular coordinates and physical curb queries for curved Bundles.

Coordinates choose no road shape: original Cartesian offset chords, authored
Section IDs and complete physical end cuts remain authoritative. No bpy import
or fixture data is required by the numerical construction.
"""
import math
import numpy as np
from .mesh import cross
from .bundle_axis import AxialChartUnavailable, AxialPlan
from .slope_projection import corner_matrices


class AngularCoordinates:
    def __init__(self, paths):
        points = np.concatenate([path.points for path in paths])
        origin = points.mean(axis=0)
        local = points - origin
        solution, _, rank, _ = np.linalg.lstsq(
            np.c_[2. * local, np.ones(len(local))], np.sum(local * local, axis=1), rcond=None)
        if rank != 3:
            raise AxialChartUnavailable('The Bundle has no determinate angular coordinate centre')
        self.center = origin + solution[:2]
        self.scale = float(np.mean(np.linalg.norm(points - self.center, axis=1)))
        middle = paths[0].frame(.5 * paths[0].length)[0] - self.center
        self.angle = math.atan2(middle[1], middle[0])
        angles = np.unwrap(np.arctan2(paths[0].points[:, 1] - self.center[1],
                                      paths[0].points[:, 0] - self.center[0]))
        self.direction = float(np.sign(angles[-1] - angles[0]))
        reach = max(getattr(path, 'maximum_reach', abs(path.section.offset) + .5 * path.section.width)
                    for path in paths)
        if (self.direction == 0. or not math.isfinite(self.scale) or self.scale <= 1e-8
                or np.min(np.linalg.norm(points - self.center, axis=1)) <= 2. * reach):
            raise AxialChartUnavailable('Angular coordinates have insufficient radial clearance')

    def local(self, points):
        delta = np.asarray(points) - self.center
        angle = np.arctan2(delta[..., 1], delta[..., 0]) - self.angle
        angle = np.arctan2(np.sin(angle), np.cos(angle))
        return np.stack((self.direction * self.scale * angle,
            -self.direction * (np.linalg.norm(delta, axis=-1) - self.scale)), axis=-1)

    def world(self, coordinates):
        coordinates = np.asarray(coordinates)
        angle = self.angle + self.direction * coordinates[..., 0] / self.scale
        radius = self.scale - self.direction * coordinates[..., 1]
        if np.any(radius <= 1e-8):
            raise ValueError('An angular quad chart reaches its coordinate origin')
        return self.center + radius[..., None] * np.stack((np.cos(angle), np.sin(angle)), axis=-1)

    def certify(self, path, *, forward=True):
        """Whole-source positive angular progress, including physical offsets."""
        direction = self.direction if forward else -self.direction
        if path.exact:
            bounds = []
            for side in (-1., 1.):
                xy = path.point(path.knots, path.section.offset + side * .5 * path.section.width)
                bounds.extend(direction * cross(xy[:-1] - self.center, np.diff(xy, axis=0)))
            lower = float(min(bounds))
        else:
            from .bundles import _beziers
            lower = math.inf
            reach = abs(path.section.offset) + .5 * path.section.width
            for controls in _beziers(path):
                position, velocity = controls - self.center, 3. * np.diff(controls, axis=0)
                coefficients = np.zeros(6)
                for i, p in enumerate(position):
                    for j, v in enumerate(velocity):
                        coefficients[i + j] += (math.comb(3, i) * math.comb(2, j)
                            / math.comb(5, i + j) * cross(p, v))
                margin = 128. * np.finfo(float).eps * max(1.,
                    float(abs(position).max()) * float(abs(velocity).max()))
                lower = min(lower, float(np.min(direction * coefficients))
                    - reach * float(np.linalg.norm(velocity, axis=1).max()) - margin)
        if lower <= 1e-8:
            raise AxialChartUnavailable('No continuous positive angular-progress certificate')
        return lower

    def longitudinal_bound(self, mesh):
        """Bound all bilinear gradients over each quad's angular direction cone."""
        angle = self.local(mesh.xyz[mesh.quads, :2])[..., 0] / self.scale
        first, last = angle.min(axis=1), angle.max(axis=1)
        middle, half = .5 * (last + first), .5 * (last - first)
        if half.max() >= .5 * math.pi:
            raise ValueError('A quad has no acute angular tangent cone')
        theta = self.angle + self.direction * middle
        tangent = np.c_[-np.sin(theta), np.cos(theta)]
        normal = np.c_[np.cos(theta), np.sin(theta)]
        # Subtracting the datum avoids cancellation for translated profiles.
        heights = mesh.xyz[mesh.quads, 2] - float(mesh.xyz[:, 2].min())
        gradients = np.einsum('fcij,fj->fci', corner_matrices(mesh), heights)
        along = abs(np.sum(gradients * tangent[:, None], axis=-1))
        across = abs(np.sum(gradients * normal[:, None], axis=-1))
        cosine, sine = np.cos(half)[:, None], np.sin(half)[:, None]
        bound = np.where(across * cosine <= along * sine,
            np.linalg.norm(gradients, axis=-1), along * cosine + across * sine)
        return float(bound.max())

    def certify_curve(self, curve, *, forward=True):
        """Prove positive angular progress of each complete polynomial curb."""
        from .bundle_boundaries import _certify_increasing
        if hasattr(curve, 'segments'):
            return min(self.certify_curve(piece, forward=forward)/(curve.knots[i+1]-curve.knots[i])
                       for i, piece in enumerate(curve.segments))
        position = curve.controls-self.center
        velocity = curve.first_derivative
        degree, derivative_degree = len(position)-1, len(velocity)-1
        coefficients = np.zeros(degree+derivative_degree+1)
        for i, point in enumerate(position):
            for j, tangent in enumerate(velocity):
                coefficients[i+j] += (math.comb(degree,i)*math.comb(derivative_degree,j)
                    /math.comb(degree+derivative_degree,i+j))*cross(point,tangent)
        direction = self.direction if forward else -self.direction
        margin = 256.*np.finfo(float).eps*max(1., float(abs(position).max())*float(abs(velocity).max()),
                                             float(abs(curve.controls).max())*float(abs(velocity).max()))
        return _certify_increasing(direction*coefficients-margin, coordinate='angular')


class RayCurb:
    """Unique ray intersection with an angularly monotone physical polyline."""
    def __init__(self, frame, points):
        self.frame = frame
        self.points = np.asarray(points, dtype=float)
        self.stations = frame.local(self.points)[:, 0]
        if np.any(np.diff(self.stations) <= 1e-9):
            raise AxialChartUnavailable('A physical curb is not angularly monotone')
        if np.ptp(self.stations) / frame.scale >= 2. * math.pi - .01:
            raise AxialChartUnavailable('The Bundle needs another angular branch cut')

    def __call__(self, station):
        station = np.asarray(station)
        index = np.clip(np.searchsorted(self.stations, station, side='right') - 1, 0, len(self.points) - 2)
        angle = self.frame.angle + self.frame.direction * station / self.frame.scale
        direction = np.stack((np.cos(angle), np.sin(angle)), axis=-1)
        a, b = self.points[index], self.points[index + 1]
        delta = b - a
        determinant = cross(delta, direction)
        if np.any(abs(determinant) < 1e-10):
            raise AxialChartUnavailable('A radial fibre is parallel to its curb chord')
        fraction = -cross(a - self.frame.center, direction) / determinant
        result = a + fraction[..., None] * delta
        result = np.where((abs(station - self.stations[0]) < 1e-10)[..., None], self.points[0], result)
        result = np.where((abs(station - self.stations[-1]) < 1e-10)[..., None], self.points[-1], result)
        if np.any(np.linalg.norm(result - self.frame.center, axis=-1) <= 1e-8):
            raise AxialChartUnavailable('A physical curb reaches the angular origin')
        return result

    def jet(self, station):
        """Local transverse coordinate and dy/ds on the actual curb chord."""
        station = np.asarray(station)
        index = np.clip(np.searchsorted(self.stations, station, side='right')-1, 0, len(self.points)-2)
        edge = self.points[index+1]-self.points[index]
        xy = self(station)
        delta = xy-self.frame.center
        radius = np.linalg.norm(delta, axis=-1)
        radial = delta/radius[..., None]
        denominator = cross(edge, radial)
        if np.any(abs(denominator) < 1e-10):
            raise ValueError('No finite angular curb jet')
        slope = radius/self.frame.scale*np.sum(edge*radial, axis=-1)/denominator
        return self.frame.local(xy)[..., 1], slope


class AngularCrown:
    """Existing dynamic peak policy evaluated on actual Cartesian curb rays."""
    fill_micro_gap = False
    peaks = AxialPlan.peaks

    def __init__(self, frame, curves, *, datum, rise):
        self.frame = frame
        self.curves = tuple(RayCurb(frame, points) for points in curves)
        self.datum, self.rise = float(datum), float(rise)
        self.ceiling = self.datum + self.rise

    def edges(self, station):
        return np.asarray([self.frame.local(curve(station))[..., 1] for curve in self.curves])

    def values(self, points):
        local = self.frame.local(points)
        ao, ai, bi, bo = self.edges(local[..., 0])
        high, low = np.maximum(ao, bi), np.minimum(ai, bo)
        first, last = self.peaks(local[..., 0])
        y = local[..., 1]
        distance = np.where(y >= first, (y - first) / np.maximum(high - first, 1e-12),
            np.where(y <= last, (y - last) / np.maximum(last - low, 1e-12), 0.))
        return self.datum + self.rise * np.clip(1. - distance * distance, 0., 1.)

    def gradient(self, points):
        """World-XY gradient of the physical crown, including moving peaks.

        Curb samples define actual Cartesian chords, so their ray derivatives
        are analytic. At a curb return the incoming one-sided polynomial
        gradient; differentiating clipped values would incorrectly halve it.
        At an Exact chord knot the selected chord gives the corresponding
        one-sided derivative. Existing value/sampling arithmetic is unchanged.
        """
        points = np.asarray(points)
        local = self.frame.local(points)
        jets = [curve.jet(local[..., 0]) for curve in self.curves]
        (ao, da), (ai, di), (bi, db), (bo, do) = jets
        high, dh = np.maximum(ao, bi), np.where(ao >= bi, da, db)
        low, dl = np.minimum(ai, bo), np.where(ai <= bo, di, do)
        inner_high, dih = np.maximum(ai, bo), np.where(ai >= bo, di, do)
        inner_low, dil = np.minimum(ao, bi), np.where(ao <= bi, da, db)
        gap, dg = inner_high-inner_low, dih-dil
        wa, wb = ao-ai, bi-bo
        width = np.maximum(1e-10, np.minimum(wa, wb))
        dw = np.where(np.minimum(wa, wb) > 1e-10, np.where(wa <= wb, da-di, db-do), 0.)
        if self.fill_micro_gap:
            gap, dg = np.zeros_like(gap), np.zeros_like(dg)
        fraction = np.clip(-gap/width, 0., 1.)
        df = np.where((fraction > 0.) & (fraction < 1.), -dg/width+gap*dw/width**2, 0.)
        blend = fraction**3*(10.+fraction*(-15.+6.*fraction))
        derivative = 30.*fraction**2*(1.-fraction)**2*df
        inner, outer = .5*(inner_high+inner_low), .5*(high+low)
        peak = (1.-blend)*inner+blend*outer
        dp = ((1.-blend)*.5*(dih+dil)+blend*.5*(dh+dl)
              +derivative*(outer-inner))
        first, dfirst = np.where(gap > 0., inner_high, peak), np.where(gap > 0., dih, dp)
        last, dlast = np.where(gap > 0., inner_low, peak), np.where(gap > 0., dil, dp)
        y = local[..., 1]
        upper, lower = y >= first, y <= last
        peak_value = np.where(upper, first, last)
        peak_slope = np.where(upper, dfirst, dlast)
        side_width = np.where(upper, high-first, last-low)
        side_slope = np.where(upper, dh-dfirst, dlast-dl)
        active = upper | lower
        valid = side_width > 1e-12
        side_slope = np.where(valid, side_slope, 0.)
        side_width = np.maximum(side_width, 1e-12)
        distance = np.where(active, (y-peak_value)/side_width, 0.)
        ds = (-peak_slope-distance*side_slope)/side_width
        z_s = np.where(active, -2.*self.rise*distance*ds, 0.)
        z_y = np.where(active, -2.*self.rise*distance/side_width, 0.)
        delta = points-self.frame.center
        radius = np.linalg.norm(delta, axis=-1)
        if np.any(radius <= 1e-8):
            raise ValueError('An angular profile gradient reaches its coordinate origin')
        radial = delta/radius[..., None]
        angular = np.stack((-radial[..., 1], radial[..., 0]), axis=-1)
        return self.frame.direction*(z_s[..., None]*(self.frame.scale/radius)[..., None]*angular
                                     -z_y[..., None]*radial)


def polyline_distance(points, curve):
    """Bounded scratch distance used only to audit declared physical rims."""
    answer = np.full(len(points), np.inf)
    a, delta = curve[:-1], np.diff(curve, axis=0)
    length2 = np.sum(delta * delta, axis=1)
    for start in range(0, len(points), 64):
        difference = points[start:start + 64, None] - a
        fraction = np.clip(np.sum(difference * delta, axis=-1) / length2, 0., 1.)
        answer[start:start + 64] = np.sqrt(np.min(np.sum(
            (difference - fraction[..., None] * delta) ** 2, axis=-1), axis=1))
    return answer
