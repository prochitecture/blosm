"""Bounded, batched boundary-jet fields for compound surface construction."""
import math
import numpy as np
from .updates import drain_steps


def _compile_rule(rule, fit_tolerance):
    """The existing residual-qualified polynomial branch, unchanged."""
    polynomial = getattr(rule, 'polynomial', None)
    supplied = polynomial() if polynomial is not None else None
    for degree in (5, 10, 20):
        test = np.linspace(0., 1., 2 * degree + 7)
        h, g = rule.sample(test)
        actual = np.column_stack((h, g))
        if supplied is not None:
            c = np.asarray(supplied, dtype=float)
            if (c.ndim == 2 and c.shape[1] == 3 and 1 <= len(c) <= degree + 1
                    and np.isfinite(c).all()):
                error = float(np.max(abs(np.polynomial.polynomial.polyval(test, c).T - actual)))
                if np.isfinite(error) and error <= fit_tolerance:
                    break
            supplied = None
        nodes = .5 * (1. - np.cos(np.linspace(0., np.pi, degree + 1)))
        height, gradient = rule.sample(nodes)
        c = np.polynomial.polynomial.polyfit(nodes, np.column_stack((height, gradient)), degree)
        approximate = np.polynomial.polynomial.polyval(test, c).T
        error = float(np.max(abs(approximate - actual)))
        if np.isfinite(error) and error <= fit_tolerance:
            break
    else:
        raise ValueError('Boundary jet polynomial needs a different approximation branch')
    return c, error, float(np.linalg.norm(g, axis=1).max())


class _RulePiece:
    def __init__(self, rule, first, last):
        self.rule, self.first, self.last = rule, first, last

    def sample(self, t):
        sample_piece = getattr(self.rule, 'sample_piece', None)
        if sample_piece is not None:
            return sample_piece(t, self.first, self.last)
        return self.rule.sample(self.first+np.asarray(t)*(self.last-self.first))


class BoundaryJetField:
    """Blend bounded local Taylor predictions, preserving boundary values/jets.

    The field is defined over its supplied pavement region. Closest-edge
    projection and bounded Taylor branches can introduce internal derivative
    seams: this is not a global C1 certificate. Publication must independently
    check the sampled mesh, incoming joins, elevation and physical slopes.
    """
    def __init__(self, boundaries, *, fit_tolerance=1e-8, scratch_bytes=8 * 1024 * 1024):
        if (not boundaries or not math.isfinite(fit_tolerance) or fit_tolerance <= 0.
                or type(scratch_bytes) is not int or scratch_bytes < 4096):
            raise ValueError('Invalid boundary field settings')
        self.datum, self.ceiling = boundaries[0].datum, boundaries[0].ceiling
        if not all(map(math.isfinite, (self.datum, self.ceiling))) or self.ceiling < self.datum:
            raise ValueError('Boundary field has an invalid elevation envelope')
        if any(abs(b.datum - self.datum) > 1e-8 or abs(b.ceiling - self.ceiling) > 1e-8 for b in boundaries):
            raise ValueError('Boundary field needs one elevation envelope')
        self.origins = np.concatenate([b.xy for b in boundaries])
        self.edges = np.concatenate([np.roll(b.xy, -1, axis=0) - b.xy for b in boundaries])
        self.length2 = np.sum(self.edges ** 2, axis=1)
        if sum(len(b.rules) for b in boundaries) != len(self.edges):
            raise ValueError('Every boundary field edge needs one jet rule')
        if (len(self.edges) > 20000 or not np.all(np.isfinite(self.origins))
                or np.min(self.length2) <= 1e-16):
            raise ValueError('Boundary field has excessive, invalid or collapsed edges')
        self.batch_rows = max(1, min(128, scratch_bytes // (len(self.edges) * 160)))
        rules = tuple(rule for boundary in boundaries for rule in boundary.rules)
        intervals = []
        branches = 0
        for rule in rules:
            declared = getattr(rule, 'polynomial_breaks', None)
            if declared is not None and not callable(declared):
                raise ValueError('Boundary jet breaks require an optional explicit declaration')
            breaks = () if declared is None else declared()
            if (not isinstance(breaks, tuple) or len(breaks) > 19999
                    or any(not np.isfinite(value) or not 0. < value < 1. for value in breaks)
                    or any(first >= last for first, last in zip(breaks, breaks[1:]))):
                raise ValueError('Boundary jet pieces need finite ordered declared interior breaks')
            branches += len(breaks)+1
            if branches > 20000:
                raise ValueError('Boundary jet pieces exceed the existing prediction branch budget')
            intervals.append((0., *breaks, 1.))
        coefficients, pieces, maximum_error, maximum_slope = [], [], 0., 0.
        for edge, (rule, bounds) in enumerate(zip(rules, intervals)):
            if len(bounds) == 2:
                c, error, slope = _compile_rule(rule, fit_tolerance)
                maximum_error = max(maximum_error, error)
                maximum_slope = max(maximum_slope, slope)
                coefficients.append(c)
                continue
            local = []
            for first, last in zip(bounds, bounds[1:]):
                c, error, slope = _compile_rule(_RulePiece(rule, first, last), fit_tolerance)
                maximum_error = max(maximum_error, error)
                maximum_slope = max(maximum_slope, slope)
                local.append(c)
            array = np.zeros((max(map(len, local)), len(local), 3))
            for index, c in enumerate(local):
                array[:len(c), index] = c
            bound_array = np.asarray(bounds)
            array.flags.writeable = bound_array.flags.writeable = False
            pieces.append((edge, bound_array, array))
            coefficients.append(np.zeros((1, 3)))
        count = max(map(len, coefficients))
        self.coefficients = np.zeros((count, len(coefficients), 3))
        for edge, c in enumerate(coefficients):
            self.coefficients[:len(c), edge] = c
        self.maximum_boundary_slope = maximum_slope
        self._pieces = tuple(pieces)
        self.metadata = dict(maximum_sampled_fit_error=maximum_error,
                             maximum_boundary_slope=maximum_slope, batch_rows=self.batch_rows)
        if pieces:
            self.metadata['declared_polynomial_pieces'] = dict(
                physical_edges=len(self.edges), prediction_branches=branches,
                split_edges=len(pieces), maximum_prediction_branches=20000,
                physical_projection_and_weights_unchanged=True)
        for values in (self.origins, self.edges, self.length2, self.coefficients):
            values.flags.writeable = False

    def __call__(self, xy):
        return drain_steps(self.evaluate_steps(xy))

    def evaluate_steps(self, xy):
        xy = np.asarray(xy, dtype=float)
        if xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy)):
            raise ValueError('Boundary field queries need finite XY rows')
        result = np.empty(len(xy))
        for start in range(0, len(xy), self.batch_rows):
            if start:
                yield
            query = xy[start:start + self.batch_rows]
            delta = query[:, None] - self.origins
            t = np.clip(np.sum(delta * self.edges, axis=-1) / self.length2, 0., 1.)
            offset_xy = delta - t[..., None] * self.edges
            distance2 = np.sum(offset_xy ** 2, axis=-1)
            jet = np.zeros(t.shape + (3,))
            for c in self.coefficients[::-1]:
                jet = jet * t[..., None] + c
            for edge, bounds, coefficients in self._pieces:
                at = t[:, edge]
                branch = np.clip(np.searchsorted(bounds, at, side='right')-1, 0, len(bounds)-2)
                local = (at-bounds[branch])/(bounds[branch+1]-bounds[branch])
                value = np.zeros((len(query), 3))
                for c in coefficients[::-1]:
                    value = value*local[:, None]+c[branch]
                jet[:, edge] = value
            height = np.clip(jet[..., 0], self.datum, self.ceiling)
            offset_z = np.sum(jet[..., 1:] * offset_xy, axis=-1)
            room = np.where(offset_z >= 0., self.ceiling - height, height - self.datum)
            ratio = np.divide(offset_z, room, out=np.zeros_like(offset_z), where=room > 1e-14)
            prediction = height + room * np.tanh(ratio)
            # Scale inverse-fourth-distance weights before division to avoid
            # overflow at boundary points; at least one weight stays nonzero.
            scale = np.maximum(distance2.min(axis=1), 1e-20)
            weights = (scale[:, None] / np.maximum(distance2, 1e-20)) ** 2 * np.sqrt(self.length2)
            result[start:start + len(query)] = np.sum(prediction * weights, axis=1) / weights.sum(axis=1)
        return np.clip(result, self.datum, self.ceiling)
