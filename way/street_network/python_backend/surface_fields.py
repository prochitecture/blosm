"""Bounded, batched boundary-jet fields for compound surface construction."""
import math
import numpy as np


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
        coefficients, maximum_error, maximum_slope = [], 0., 0.
        for boundary in boundaries:
            for rule in boundary.rules:
                for degree in (5, 10, 20):
                    nodes = .5 * (1. - np.cos(np.linspace(0., np.pi, degree + 1)))
                    height, gradient = rule.sample(nodes)
                    c = np.polynomial.polynomial.polyfit(nodes, np.column_stack((height, gradient)), degree)
                    test = np.linspace(0., 1., 2 * degree + 7)
                    h, g = rule.sample(test)
                    actual = np.column_stack((h, g))
                    approximate = np.polynomial.polynomial.polyval(test, c).T
                    error = float(np.max(abs(approximate - actual)))
                    if np.isfinite(error) and error <= fit_tolerance:
                        break
                else:
                    raise ValueError('Boundary jet polynomial needs a different approximation branch')
                maximum_error = max(maximum_error, error)
                maximum_slope = max(maximum_slope, float(np.linalg.norm(g, axis=1).max()))
                coefficients.append(c)
        count = max(map(len, coefficients))
        self.coefficients = np.zeros((count, len(coefficients), 3))
        for edge, c in enumerate(coefficients):
            self.coefficients[:len(c), edge] = c
        self.maximum_boundary_slope = maximum_slope
        self.metadata = dict(maximum_sampled_fit_error=maximum_error,
                             maximum_boundary_slope=maximum_slope, batch_rows=self.batch_rows)
        for values in (self.origins, self.edges, self.length2, self.coefficients):
            values.flags.writeable = False

    def __call__(self, xy):
        xy = np.asarray(xy, dtype=float)
        if xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy)):
            raise ValueError('Boundary field queries need finite XY rows')
        result = np.empty(len(xy))
        for start in range(0, len(xy), self.batch_rows):
            query = xy[start:start + self.batch_rows]
            delta = query[:, None] - self.origins
            t = np.clip(np.sum(delta * self.edges, axis=-1) / self.length2, 0., 1.)
            offset_xy = delta - t[..., None] * self.edges
            distance2 = np.sum(offset_xy ** 2, axis=-1)
            jet = np.zeros(t.shape + (3,))
            for c in self.coefficients[::-1]:
                jet = jet * t[..., None] + c
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
