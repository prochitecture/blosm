"""Vectorized Exact/Smooth centerlines, physical stations and offset guards.

Smooth uses a bounded source envelope: gentle interpolating cubics are retained,
and excessive span-wide bows are replaced by local tangent cubic bends with
straight source runs. The continuous offset-curvature certificate is unchanged.
Arc length uses bounded
adaptive intervals and Gauss integration, rather than a fixed dense subdivision
for every segment. No SciPy or Blender API is required.
"""
import math
from dataclasses import replace
import numpy as np
from .mesh import cross

_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(8)


def right(vector):
    return np.stack((vector[..., 1], -vector[..., 0]), axis=-1)


def unit(vector):
    length = np.linalg.norm(vector, axis=-1, keepdims=True)
    if np.any(length < 1e-10):
        raise ValueError('A path has an undefined tangent')
    return vector / length


class PathSampleBudgetExceeded(ValueError):
    """The requested sampling would exceed the supplied point allocation."""


class Path:
    def __init__(self, section, tolerance=.005, endpoint_directions=None, *, deviation_limit=None,
                 endpoint_bends=None, endpoint_alignment_limits=None, returning=False, _coalescing=True,
                 _radius_factor=1.5):
        self.coalescing = None
        self.radius_reserve_fallback = None
        try:
            self._initialize(section, tolerance, endpoint_directions, deviation_limit=deviation_limit,
                endpoint_bends=endpoint_bends, endpoint_alignment_limits=endpoint_alignment_limits,
                returning=returning, minimum_radius_factor=_radius_factor)
        except ValueError as exc:
            if (not _coalescing or section.interpolation != 3 or len(section.points) < 3
                    or not any(message in str(exc) for message in
                        ('Local Smooth', 'Smooth offset curvature'))):
                raise
            candidate = self._coalesced(section, tolerance, endpoint_directions,
                deviation_limit, endpoint_bends, endpoint_alignment_limits, returning)
            if candidate is None:
                # Fit the largest reserve first. This preference is separate
                # from mandatory continuous curvature and envelope guards.
                for factor in (1.45, 1.4, 1.35, 1.3, 1.25):
                    try:
                        candidate = Path(section, tolerance, endpoint_directions,
                            deviation_limit=deviation_limit, endpoint_bends=endpoint_bends,
                            endpoint_alignment_limits=endpoint_alignment_limits, returning=returning,
                            _coalescing=False, _radius_factor=factor)
                    except ValueError:
                        continue
                    candidate.radius_reserve_fallback = dict(minimum_radius_factor=factor,
                        original_failure=str(exc))
                    break
                if candidate is None:
                    raise exc
            else:
                candidate.coalescing['original_failure'] = str(exc)
            self.__dict__.update(candidate.__dict__)

    def _coalesced(self, section, tolerance, endpoint_directions, deviation_limit,
                   endpoint_bends, endpoint_alignment_limits, returning):
        from .local_smoothing import coalesced_source_indices, within_polyline_envelope
        original = np.asarray(section.points, dtype=float)
        bound = max(4.*tolerance, section.width if deviation_limit is None else deviation_limit)
        directions = dict(endpoint_directions or {})
        # Coalescing must not independently change an authored endpoint frame.
        directions.setdefault(0, unit(original[1, :2]-original[0, :2]))
        directions.setdefault(1, unit(original[-1, :2]-original[-2, :2]))
        protected = tuple(1 if end == 0 else len(original)-2 for end in (endpoint_bends or {}))
        seen = set()
        for fraction in (.01, .025, .05, .1, .2, .25, .3, .4):
            budget = min(fraction*section.width, .4*bound)
            indices = coalesced_source_indices(original[:, :2], budget, protected)
            if len(indices) == len(original) or indices in seen:
                continue
            seen.add(indices)
            source = replace(section, points=tuple(map(tuple, original[list(indices)])))
            try:
                candidate = Path(source, tolerance, directions, deviation_limit=bound-budget,
                    endpoint_bends=endpoint_bends, endpoint_alignment_limits=endpoint_alignment_limits,
                    returning=returning, _coalescing=False)
            except ValueError:
                continue
            if not within_polyline_envelope(candidate.coefficients, original[:, :2], bound):
                continue
            # Preserve the authored source identity and complete controls in
            # every downstream consumer and report; only evaluation data differ.
            candidate.section = section
            candidate.source_points = original[:, :2].copy()
            candidate.deviation_limit = bound
            candidate.coalescing = dict(policy='CERTIFIED_DENSE_CONTROLS_V1',
                retained_indices=indices, source_control_count=len(original),
                coalescing_limit=budget, source_deviation_limit=bound)
            return candidate
        return None

    def _initialize(self, section, tolerance=.005, endpoint_directions=None, *, deviation_limit=None,
                    endpoint_bends=None, endpoint_alignment_limits=None, returning=False,
                    minimum_radius_factor=1.5):
        self.section = section
        self.tolerance = tolerance
        self.points = np.asarray(section.points, dtype=float)[:, :2]
        self.source_points = self.points.copy()
        self.local_bends = ()
        self.localized = False
        self.alignment_fallback = None
        self.endpoint_bends = endpoint_bends or {}
        if deviation_limit is not None and (not math.isfinite(deviation_limit) or deviation_limit <= 0.):
            raise ValueError('Source-deviation limit must be positive and finite')
        self.deviation_limit = max(4. * tolerance, section.width if deviation_limit is None else deviation_limit)
        # Keep retained span-wide interpolation in the central half of the
        # source corridor. Physical local rounds may need the larger envelope;
        # using that envelope to accept an interpolating cubic allowed an
        # entire long straight run to bow sideways by a road width.
        self.alignment_limit = max(4. * tolerance, .25 * self.deviation_limit)
        if np.max(abs(np.asarray(section.points)[:, 2])) > 1e-7:
            raise ValueError('Source centerlines must remain in the construction XY plane; terrain placement derives elevation separately')
        self.returning = bool(returning)
        if np.linalg.norm(self.points[0] - self.points[-1]) < 1e-7 and not self.returning:
            raise ValueError('Closed Section needs explicit start/end incidence at its closure')
        delta = np.diff(self.points, axis=0)
        self.chords = np.linalg.norm(delta, axis=1)
        self.alignment_limits = np.full(len(self.chords), self.alignment_limit)
        for end, limit in (endpoint_alignment_limits or {}).items():
            self.alignment_limits[0 if end == 0 else -1] = min(self.deviation_limit, max(self.alignment_limit, limit))
        if np.min(self.chords) < 1e-7:
            raise ValueError('A path has coincident XY controls')
        self.knots = np.r_[0., np.cumsum(self.chords)]
        directions = delta / self.chords[:, None]
        self.exact = section.interpolation == 0
        if self.exact:
            dot = np.sum(directions[:-1] * directions[1:], axis=1)
            if np.any(dot < -.95):
                raise ValueError('Exact polyline has a near-reversal miter')
            normals = np.vstack((right(directions[0]),
                (right(directions[:-1]) + right(directions[1:])) / (1. + dot)[:, None],
                right(directions[-1])))
            if np.max(np.linalg.norm(normals, axis=1)) > 5.:
                raise ValueError('Exact miter exceeds the supported reach')
            self.normals = normals
            self.directions = directions
            self.length = float(self.knots[-1])
            self.stations = self.knots.copy()
            self.parameters = self.knots.copy()
            self.certificate_intervals = 0
            return
        tangent = unit(np.vstack((directions[0], directions[:-1] + directions[1:], directions[-1])))
        for end, value in (endpoint_directions or {}).items():
            tangent[0 if end == 0 else -1] = unit(np.asarray(value))
        self.coefficients = self._coefficients(delta, tangent)
        certified = self._certified()
        if not certified:
            # Established spacing-weighted retry; the interpolation still goes
            # through every authored point and retains endpoint directions.
            tangent = unit(np.vstack((directions[0],
                directions[:-1] / self.chords[:-1, None] + directions[1:] / self.chords[1:, None], directions[-1])))
            for end, value in (endpoint_directions or {}).items():
                tangent[0 if end == 0 else -1] = unit(np.asarray(value))
            self.coefficients = self._coefficients(delta, tangent)
            certified = self._certified()
        from .local_smoothing import LocalCubicSmoothing, within_source_envelope
        if endpoint_bends or not certified or not within_source_envelope(
                self.coefficients, self.source_points, self.alignment_limits):
            reach = abs(section.offset) + .5 * section.width
            original = self.points, self.chords, self.coefficients, self.knots
            original_certified = certified
            try:
                self.points, self.chords, self.coefficients, self.local_bends = LocalCubicSmoothing(
                    self.source_points, reach, self.deviation_limit, endpoint_directions, endpoint_bends,
                    minimum_radius_factor=minimum_radius_factor).build()
                self.knots = np.r_[0., np.cumsum(self.chords)]
                self.localized = True
                certified = self._certified()
                if not certified:
                    raise ValueError('Smooth offset curvature could fold this carriageway')
            except ValueError as exc:
                # Alignment is a preference; the original physical guards are
                # mandatory. A sharp/dimensional transition may need a wider,
                # continuously certified curve than a local circular-style bend.
                # Fixed shared halves must never fall back independently.
                if (endpoint_bends or not original_certified or not within_source_envelope(
                        original[2], self.source_points, self.deviation_limit)):
                    raise
                self.points, self.chords, self.coefficients, self.knots = original
                self.local_bends, self.localized = (), False
                self.alignment_fallback = str(exc)
                certified = self._certified()
        if not certified:
            raise ValueError('Smooth offset curvature could fold this carriageway')
        parameters = [0.]
        # Midpoint-to-chord error controls the interval seed. Continuous offset
        # safety is supplied separately by the Bernstein curvature certificate.
        for i in range(len(self.chords)):
            stack = [(self.knots[i], self.knots[i + 1], 0)]
            while stack:
                a, b, depth = stack.pop()
                t = np.linspace(a, b, 5)
                p = self.jet(t)[0]
                linear = p[0] + np.linspace(0., 1., 5)[:, None] * (p[-1] - p[0])
                if np.max(np.linalg.norm(p - linear, axis=1)) > min(tolerance * .25, .001):
                    if depth >= 14:
                        raise ValueError('Arc-length sampling exceeded its refinement budget')
                    middle = .5 * (a + b)
                    stack.extend(((middle, b, depth + 1), (a, middle, depth + 1)))
                else:
                    parameters.append(b)
        self.parameters = np.asarray(parameters)
        lengths = self._integral(self.parameters[:-1], self.parameters[1:])
        self.stations = np.r_[0., np.cumsum(lengths)]
        self.length = float(self.stations[-1])
        if self.localized:
            self.piece_stations = self.station(self.knots)

    def _coefficients(self, delta, tangent):
        a, b = tangent[:-1] * self.chords[:, None], tangent[1:] * self.chords[:, None]
        return np.stack((self.points[:-1], a, 3. * delta - 2. * a - b,
                         -2. * delta + a + b), axis=1)

    def _certified(self, reach=None):
        record = reach is None
        a, b, c = [self.coefficients[:, None, k, :] for k in (1, 2, 3)]
        reach = abs(self.section.offset) + self.section.width * .5 if reach is None else reach
        for count in ((32, 64) if self.localized else (8, 16, 32, 64)):
            t0 = (np.arange(count) / count)[None, :, None]
            t1 = t0 + 1. / count
            tm = t0 + .5 / count
            velocity = lambda t: a + t * (2. * b + 3. * t * c)
            acceleration = lambda t: 2. * b + 6. * t * c
            d0, d2 = velocity(t0), velocity(t1)
            acc0, acc1 = acceleration(t0), acceleration(t1)
            d1 = d0 + acc0 * (.5 / count)
            axis = unit(velocity(tm))
            lower = np.minimum.reduce([np.sum(d * axis, axis=-1) for d in (d0, d1, d2)])
            upper = np.maximum.reduce([abs(v) for v in (cross(d0, acc0),
                (2. * cross(d1, acc0) + cross(d0, acc1)) / 3.,
                (cross(d2, acc0) + 2. * cross(d1, acc1)) / 3., cross(d2, acc1))])
            minimum_speed = 1e-5 * self.chords[:, None] if self.localized else 1e-5
            denominator = np.maximum(lower, 1e-30 if self.localized else 1e-8)
            if np.all(lower > minimum_speed) and np.all(reach * upper / denominator ** 3 < .98):
                if record:
                    self.certificate_intervals = count
                    self.piece_curvature_bounds = np.max(upper / lower ** 3, axis=1)
                    self.curvature_bound = float(self.piece_curvature_bounds.max())
                return True
        return False

    def jet(self, parameter):
        parameter = np.asarray(parameter)
        at = np.clip(np.searchsorted(self.knots, parameter, side='right') - 1, 0, len(self.chords) - 1)
        length = self.chords[at]
        u = ((parameter - self.knots[at]) / length)[..., None]
        if self.exact:
            p = self.points[at] + u * (self.points[at + 1] - self.points[at])
            return p, self.directions[at], np.zeros_like(p)
        p, a, b, c = np.moveaxis(self.coefficients[at], -2, 0)
        return (p + u * (a + u * (b + u * c)),
                (a + u * (2. * b + 3. * u * c)) / length[..., None],
                (2. * b + 6. * u * c) / length[..., None] ** 2)

    def _integral(self, first, last):
        middle, half = (first + last) * .5, (last - first) * .5
        t = middle[..., None] + half[..., None] * _GAUSS_X
        speed = np.linalg.norm(self._velocity(t), axis=-1)
        return half * np.sum(speed * _GAUSS_W, axis=-1)

    def _velocity(self, parameter):
        """The exact jet derivative without unused positions/accelerations.

        Arc-length quadrature and Newton steps require only speed. Preserve
        the jet's interval lookup, dtype behavior and arithmetic order. A
        subclass overriding the full jet keeps that extension point.
        """
        if type(self).jet is not Path.jet:
            return self.jet(parameter)[1]
        parameter = np.asarray(parameter)
        at = np.clip(np.searchsorted(self.knots, parameter, side='right') - 1, 0, len(self.chords) - 1)
        if self.exact:
            return self.directions[at]
        length = self.chords[at]
        u = ((parameter - self.knots[at]) / length)[..., None]
        _, a, b, c = np.moveaxis(self.coefficients[at], -2, 0)
        return (a + u * (2. * b + 3. * u * c)) / length[..., None]

    def parameter(self, station):
        station = np.asarray(station)
        if np.any(station < -1e-8) or np.any(station > self.length + 1e-8):
            raise ValueError('Station falls outside the source path')
        if self.exact:
            return np.clip(station, 0., self.length)
        at = np.clip(np.searchsorted(self.stations, station, side='right') - 1, 0, len(self.stations) - 2)
        lo, hi = self.parameters[at], self.parameters[at + 1]
        t = lo + (hi - lo) * (station - self.stations[at]) / (self.stations[at + 1] - self.stations[at])
        for _ in range(4):
            error = self._integral(lo, t) - (station - self.stations[at])
            speed = np.linalg.norm(self._velocity(t), axis=-1)
            t = np.clip(t - error / speed, lo, hi)
        return t

    def station(self, parameter):
        """Metric inverse of parameter(), including derived polynomial joins."""
        parameter = np.asarray(parameter)
        if np.any(parameter < -1e-8) or np.any(parameter > self.knots[-1] + 1e-8):
            raise ValueError('Parameter falls outside the source path')
        parameter = np.clip(parameter, 0., self.knots[-1])
        if self.exact:
            return parameter
        at = np.clip(np.searchsorted(self.parameters, parameter, side='right') - 1,
                     0, len(self.parameters) - 2)
        return self.stations[at] + self._integral(self.parameters[at], parameter)

    def frame(self, station):
        parameter = self.parameter(station)
        p, d, _ = self.jet(parameter)
        if self.exact:
            at = np.clip(np.searchsorted(self.knots, parameter, side='right') - 1, 0, len(self.chords) - 1)
            t = ((parameter - self.knots[at]) / self.chords[at])[..., None]
            normal = (1. - t) * self.normals[at] + t * self.normals[at + 1]
        else:
            normal = right(unit(d))
        return p, unit(d), normal

    def visible_curvature_bound(self, first, last):
        """Bernstein proof restricted to the actual retained Smooth interval.

        A removed junction/seam span can contain the source's largest bend.
        Its global bound remains valid but can force excessive visible rows.
        Eight bounded intervals per original cubic tighten that proof without
        changing the centerline, offsets or sampling error budget.
        """
        if self.exact:
            return 0.
        p0, p1 = self.parameter(np.array((first, last)))
        margin = max(1e-8, 128. * np.finfo(float).eps * max(1., self.knots[-1]))
        p0, p1 = max(0., p0-margin), min(self.knots[-1], p1+margin)
        maximum, eps = 0., np.finfo(float).eps
        for i, length in enumerate(self.chords):
            lo, hi = max(0., (p0-self.knots[i])/length), min(1., (p1-self.knots[i])/length)
            if hi <= lo:
                continue
            a, b, c = self.coefficients[i, 1:]
            t0 = np.linspace(lo, hi, 9)[:-1, None]
            t1 = t0+(hi-lo)/8.
            d0, d2 = a+t0*(2.*b+3.*t0*c), a+t1*(2.*b+3.*t1*c)
            acc0, acc1 = 2.*b+6.*t0*c, 2.*b+6.*t1*c
            d1 = d0+.5*(t1-t0)*acc0
            middle = .5*(t0+t1)
            axis = unit(a+middle*(2.*b+3.*middle*c))
            lower = np.minimum.reduce([np.sum(v*axis, axis=1) for v in (d0, d1, d2)])
            numerator = np.stack((cross(d0, acc0), (2.*cross(d1, acc0)+cross(d0, acc1))/3.,
                (cross(d2, acc0)+2.*cross(d1, acc1))/3., cross(d2, acc1)))
            speed_scale = max(1., float(abs(np.r_[d0, d1, d2]).max()))
            determinant_scale = max(1., speed_scale*float(abs(np.r_[acc0, acc1]).max()))
            lower -= 128.*eps*speed_scale
            if lower.min() <= 0.:
                return self.curvature_bound
            upper = np.max(abs(numerator), axis=0)+128.*eps*determinant_scale
            maximum = max(maximum, float(np.max(upper/lower**3)))
        return min(self.curvature_bound, maximum)

    def samples(self, first, last, spacing, *, tolerance=None, local_curvature=False, piecewise=False,
                max_samples=None):
        if max_samples is not None and (type(max_samples) is not int or max_samples<3):
            raise ValueError('A path sample budget must be an integer of at least three')
        if last <= first + 1e-7:
            raise ValueError('Intersection trims consume the complete carriageway')
        if self.localized or (piecewise and not self.exact):
            return self._localized_samples(first, last, spacing,
                self.tolerance if tolerance is None else tolerance, max_samples, local_curvature)
        if not self.exact:
            # A continuous curvature bound controls physical offset-curve chord
            # error. It permits regular metric rows instead of copying the fine
            # arc-length lookup table into the final mesh. This also avoids
            # tiny rows when a moving mouth almost coincides with a lookup knot.
            reach = abs(self.section.offset) + .5 * self.section.width
            curvature = self.visible_curvature_bound(first, last) if local_curvature else self.curvature_bound
            if curvature > 1e-12:
                offset_curvature = curvature / (1. - reach * curvature)
                offset_speed = 1. + reach * curvature
                allowed = min(math.sqrt(8. * (self.tolerance if tolerance is None else tolerance) / offset_curvature),
                              math.pi / (2. * offset_curvature)) / offset_speed
                spacing = min(spacing, allowed)
            count = max(2, int(np.ceil((last - first) / spacing)))
            count += count % 2
            if max_samples is not None and count+1>max_samples:
                raise PathSampleBudgetExceeded('Path sampling exceeds its allocation budget')
            return np.linspace(first, last, count + 1)
        # Exact control vertices are compulsory samples, even at sharp corners.
        fixed = np.r_[first, self.stations[(self.stations > first) & (self.stations < last)], last]
        # Subdivide each compulsory interval independently. Unioning a uniform
        # grid with compulsory samples creates arbitrarily close neighboring
        # rows even on an otherwise regular road.
        pieces, total = [], 1
        for a, b in zip(fixed, fixed[1:]):
            count = max(2, int(np.ceil((b - a) / spacing)))
            count += count % 2
            total += count
            if max_samples is not None and total>max_samples:
                raise PathSampleBudgetExceeded('Path sampling exceeds its allocation budget')
            pieces.append(np.linspace(a, b, count + 1)[:-1])
        return np.concatenate(pieces + [np.array([last])])

    def _localized_samples(self, first, last, spacing, tolerance, max_samples, local_curvature):
        """Refine only nominal row intervals meeting a certified local bend.

        Nominal intervals have equal length. Subdividing each independently
        avoids tiny leftovers at arc/line joins, while bounding every interval
        by the maximum curvature of all polynomial pieces it intersects.
        """
        if self.radius_reserve_fallback is not None:
            return self._offset_samples(first, last, spacing, tolerance, max_samples)
        count = max(2, int(np.ceil((last - first) / spacing)))
        count += count % 2
        if max_samples is not None and count + 1 > max_samples:
            raise PathSampleBudgetExceeded('Path sampling exceeds its allocation budget')
        base = np.linspace(first, last, count + 1)
        curvature = np.zeros(count)
        stations = self.piece_stations if self.localized else self.station(self.knots)
        for a, b, bound in zip(stations, stations[1:], self.piece_curvature_bounds):
            if b < first or a > last:
                continue
            lo = max(0, int(np.searchsorted(base, a, side='left')) - 1)
            hi = min(count, int(np.searchsorted(base, b, side='right')))
            curvature[lo:hi] = np.maximum(curvature[lo:hi], bound)
        if local_curvature:
            curvature = np.minimum(curvature, self.visible_curvature_bound(first, last))
        reach = abs(self.section.offset) + .5 * self.section.width
        allowed = np.full(count, spacing)
        curved = curvature > 1e-12
        k = curvature[curved]
        offset_curvature = k / (1. - reach * k)
        allowed[curved] = np.minimum(allowed[curved],
            np.minimum(np.sqrt(8. * tolerance / offset_curvature),
                       math.pi / (2. * offset_curvature)) / (1. + reach * k))
        parts = np.maximum(1, np.ceil(np.diff(base) / allowed).astype(np.int64))
        # Numerical roundoff must not split a nominal straight interval twice.
        parts[~curved] = 1
        if int(parts.sum()) % 2:
            parts[int(np.argmax(parts))] += 1
        total = int(parts.sum())
        if max_samples is not None and total + 1 > max_samples:
            raise PathSampleBudgetExceeded('Path sampling exceeds its allocation budget')
        starts = np.r_[0, np.cumsum(parts)[:-1]]
        offsets = np.arange(total) - np.repeat(starts, parts)
        rows = np.repeat(base[:-1], parts) + offsets * np.repeat(np.diff(base) / parts, parts)
        return np.r_[rows, last]

    def _signed_piece_bounds(self):
        """Signed curvature enclosure from the existing absolute certificate.

        cross(P', P'') is quadratic for a cubic. Its Bernstein coefficients
        establish each piece's possible signs without numerical sampling.
        """
        a, b, c = (self.coefficients[:, k, :] for k in (1, 2, 3))
        p0, p1, p2 = 2.*cross(a, b), 6.*cross(a, c), 6.*cross(b, c)
        bernstein = np.stack((p0, p0+.5*p1, p0+p1+p2), axis=1)
        scale = np.maximum(1., np.max(abs(bernstein), axis=1))
        margin = 128.*np.finfo(float).eps*scale
        return (np.where(bernstein.max(axis=1)+margin > 0., self.piece_curvature_bounds, 0.),
                np.where(bernstein.min(axis=1)-margin < 0., self.piece_curvature_bounds, 0.))

    def _offset_chord_bound(self, first, last, signed=None):
        """Continuous error bound for every lateral offset in the pavement.

        For rightward offset r, ds_r=(1+r*kappa)ds and its exact arc length is
        ds+r*dtheta. The offset-curvature bound times this length squared / 8
        bounds chord error. Each signed branch is convex in r, so checking the
        two outermost lateral values also covers every intermediate profile row.
        """
        positive, negative = self._signed_piece_bounds() if signed is None else signed
        lo = max(0, int(np.searchsorted(self.piece_stations, first, side='right'))-1)
        hi = min(len(positive), int(np.searchsorted(self.piece_stations, last, side='left'))+1)
        pos = float(positive[lo:hi].max(initial=0.))
        neg = float(negative[lo:hi].max(initial=0.))
        ds = last-first
        if max(pos, neg)*ds >= math.pi*.5:
            return math.inf  # Cannot infer the unwrapped tangent turn yet.
        tangents = self.frame(np.array((first, last)))[1]
        turn = math.atan2(float(cross(tangents[0], tangents[1])), float(tangents[0] @ tangents[1]))
        margin = 256.*np.finfo(float).eps*max(1., ds, self.length)
        result = 0.
        for offset in (self.section.offset-.5*self.section.width, self.section.offset+.5*self.section.width):
            denominators = 1.+offset*pos, 1.-offset*neg
            if min(denominators) <= 0.:
                return math.inf
            curvature = max(pos/denominators[0], neg/denominators[1])
            length = ds+offset*turn
            if length <= 0.:
                return math.inf
            result = max(result, curvature*(length+margin)**2/8.+margin)
        return result

    def _offset_samples(self, first, last, spacing, tolerance, max_samples):
        """Avoid combining an inner curb's curvature with an outer curb's speed.

        Used only for a newly recovered tight bend. Previously accepted paths
        retain their complete row layout. Every new interval is certified;
        there is no sampled-error relaxation or larger tolerance.
        """
        signed = self._signed_piece_bounds()
        rows = [float(first)]
        target = .99*tolerance
        while rows[-1] < last-1e-10:
            start, finish = rows[-1], min(last, rows[-1]+spacing)
            if self._offset_chord_bound(start, finish, signed) > target:
                low, high = start, finish
                for unused in range(18):
                    mid = .5*(low+high)
                    if self._offset_chord_bound(start, mid, signed) <= target:
                        low = mid
                    else:
                        high = mid
                finish = low
            if finish-start <= 1e-8:
                raise ValueError('Certified offset sampling cannot advance within its finite budget')
            rows.append(float(finish))
            if max_samples is not None and len(rows) > max_samples:
                raise PathSampleBudgetExceeded('Path sampling exceeds its allocation budget')
        rows[-1] = float(last)
        # Distribute a short terminal remainder across a bounded tail. Splitting
        # only the final pair can produce two unnecessarily short inner rows.
        if len(rows) >= 3 and rows[-1]-rows[-2] < .75*(rows[-2]-rows[-3]):
            best, minimum = None, rows[-1]-rows[-2]
            for count in range(2, min(12, len(rows)-1)+1):
                tail = np.linspace(rows[-count-1], rows[-1], count+1)
                step = float(tail[1]-tail[0])
                if step > minimum and all(self._offset_chord_bound(a, b, signed) <= tolerance
                                          for a, b in zip(tail, tail[1:])):
                    best, minimum = tail, step
            if best is not None:
                rows[-len(best):] = best.tolist()
        return np.asarray(rows)

    def point(self, station, offset=0.):
        point, _, normal = self.frame(station)
        return point + np.asarray(offset)[..., None] * normal
