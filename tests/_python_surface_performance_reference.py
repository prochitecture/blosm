"""Frozen dev78 construction oracles for optimization equivalence tests.

Copied from b9517de; deliberately retain the original scans and clipping.
Not imported by the runtime backend.
"""
import numpy as np
import math
from tested_python_streets.patterns import Grid, PatternLibrary
from tested_python_streets.mesh import cross
from tested_python_streets.regions import signed_area


class ReferencePatternLibrary(PatternLibrary):
    def quad(self, bottom, right, top, left):
        """A four-sided direct-quad atlas retaining every requested rim sample.

        Boundary partitions remain explicit (possibly nonuniform) parameters.
        Callers map them to their existing source rim vertices, so refinements
        do not create geometric cracks or force neighboring meshes to resample.
        """
        counts = (bottom, right, top, left)
        if (any(type(v) is not int or v < 1 for v in counts)
                or (bottom - top) % 2 or (right - left) % 2):
            raise ValueError('Quad rim counts require compatible opposite-edge parity')
        if ((min(bottom, top) + 1) * (min(left, right) + 1)
                + 2 * (abs(bottom - top) + abs(right - left)) > 2_000_000):
            raise ValueError('Requested quad atlas exceeds the per-part vertex budget')
        if bottom == top and left == right:
            return self.grid(bottom, left)
        key = ('quad', *counts)
        saved = self._lookup(key)
        if saved is not None:
            return saved
        base = self.grid(min(bottom, top), min(left, right))
        points, faces = list(base.uv.copy()), list(map(tuple, base.quads))
        owners = {}

        def side_of(a, b):
            if a[1] == b[1] == 0.:
                return 0
            if a[0] == b[0] == 1.:
                return 1
            if a[1] == b[1] == 1.:
                return 2
            if a[0] == b[0] == 0.:
                return 3
            return None

        def record(fi):
            face = faces[fi]
            for a, b in zip(face, (*face[1:], face[0])):
                side = side_of(points[a], points[b])
                if side is not None:
                    owners[tuple(sorted((a, b)))] = (side, fi)

        for fi in range(len(faces)):
            record(fi)
        for side, target in enumerate(counts):
            current = sum(s == side for s, fi in owners.values())
            for _ in range((target - current) // 2):
                edges = [edge for edge, (s, fi) in owners.items() if s == side]
                edge = max(edges, key=lambda e: np.linalg.norm(points[e[1]] - points[e[0]]))
                owner = owners[edge][1]
                face = faces[owner]
                for x, y in zip(face, (*face[1:], face[0])):
                    owners.pop(tuple(sorted((x, y))), None)
                at = next(i for i in range(4) if tuple(sorted((face[i], face[(i + 1) % 4]))) == edge)
                a0, a3, b1, b0 = face[at:] + face[:at]
                a1, a2 = len(points), len(points) + 1
                points.extend(((2. * points[a0] + points[a3]) / 3., (points[a0] + 2. * points[a3]) / 3.))
                p, q = len(points), len(points) + 1
                points.extend((.5 * (points[a1] + (2. * points[b0] + points[b1]) / 3.),
                               .5 * (points[a2] + (points[b0] + 2. * points[b1]) / 3.)))
                added = ((a0, a1, p, b0), (a1, a2, q, p),
                         (a2, a3, b1, q), (p, q, b1, b0))
                indices = (owner, len(faces), len(faces) + 1, len(faces) + 2)
                faces[owner] = added[0]
                faces.extend(added[1:])
                for fi in indices:
                    record(fi)
        result = Grid(np.asarray(points), np.asarray(faces, dtype=np.int32), max(bottom, top), max(left, right))
        self._positive(result)
        return self._remember(key, result)


def visibility_center(polygon, epsilon=1e-10):
    """Centroid of the strict visibility kernel, from bounded half-plane clips."""
    origin = polygon[0]
    local = polygon - origin
    lo, hi = local.min(axis=0), local.max(axis=0)
    result = np.array((lo, (hi[0], lo[1]), hi, (lo[0], hi[1])))
    for a, b in zip(local, np.roll(local, -1, axis=0)):
        if len(result) < 3:
            raise ValueError('Block has no visibility kernel; decompose it into local charts')
        edge = b - a
        distance = cross(edge, result - a) / np.linalg.norm(edge)
        points = []
        for i, (p, q) in enumerate(zip(result, np.roll(result, -1, axis=0))):
            first, last = distance[i], distance[(i + 1) % len(result)]
            inside, next_inside = first >= -epsilon, last >= -epsilon
            if inside:
                points.append(p)
            if inside != next_inside:
                fraction = np.clip(first / (first - last), 0., 1.)
                points.append(p + fraction * (q - p))
        if not points:
            raise ValueError('Block has no visibility kernel; decompose it into local charts')
        result = np.asarray(points)
        keep = np.linalg.norm(result - np.roll(result, 1, axis=0), axis=1) > epsilon
        result = result[keep]
    if len(result) < 3 or abs(signed_area(result)) <= epsilon ** 2:
        raise ValueError('Block has no strict visibility kernel')
    weights = cross(result, np.roll(result, -1, axis=0))
    center = np.sum((result + np.roll(result, -1, axis=0)) * weights[:, None], axis=0) / (3. * np.sum(weights))
    edge = np.roll(local, -1, axis=0) - local
    if np.min(cross(edge, center - local) / np.linalg.norm(edge, axis=1)) <= epsilon:
        raise ValueError('Block visibility kernel has insufficient numerical clearance')
    return center + origin


from tested_python_streets.slope_projection import SlopeProjection, corner_matrices, face_colors
from tested_python_streets.mesh import Mesh


class ReferenceSlopeProjection(SlopeProjection):
    def solve(self, mesh, fixed, *, datum, ceiling, maximum_slope, minimum_heights=None,
              feasibility_slope=None):
        """Project to a slope ball, optionally stopping at a certified feasible mesh.

        By default the displacement must converge as well as the constraints.
        A caller needing feasibility rather than the nearest height vector may
        supply a larger acceptance slope and project towards a stricter inner
        target. Fixed vertices and the height box still hold at every return;
        the caller must independently certify its physical surface constraints.
        """
        mesh.validate()
        fixed = np.asarray(fixed)
        if fixed.dtype != np.bool_ or fixed.shape != (len(mesh.xyz),):
            raise ValueError('Slope projection needs a Boolean fixed-vertex mask')
        if (not all(map(math.isfinite, (datum, ceiling, maximum_slope)))
                or ceiling < datum or maximum_slope <= 0.):
            raise ValueError('Invalid slope-projection envelope')
        if feasibility_slope is not None and (not math.isfinite(feasibility_slope)
                or feasibility_slope < maximum_slope):
            raise ValueError('Feasibility slope must be finite and at least the projection target')
        lower = datum
        if minimum_heights is not None:
            lower = np.asarray(minimum_heights, dtype=float)
            if (lower.shape != (len(mesh.xyz),) or not np.all(np.isfinite(lower))
                    or np.any(lower < datum) or np.any(lower > ceiling)):
                raise ValueError('Invalid per-vertex minimum heights')
            if np.any(mesh.xyz[fixed, 2] < lower[fixed] - 1e-10):
                raise ValueError('Fixed boundary heights violate minimum heights')
            lower = lower.copy()
            lower[fixed] = mesh.xyz[fixed, 2]
        estimated = len(mesh.quads) * 2048 + len(mesh.xyz) * 96
        if len(mesh.quads) > self.max_cells or estimated > self.memory_bytes:
            raise ValueError('Slope projection exceeds its mesh/memory budget')
        initial = mesh.xyz[:, 2]
        if np.any(initial[fixed] < datum - 1e-10) or np.any(initial[fixed] > ceiling + 1e-10):
            raise ValueError('Fixed boundary heights exceed the crown envelope')
        coefficient = corner_matrices(mesh)
        faces = mesh.quads
        def gradients(z):
            return np.linalg.norm(np.einsum('fcij,fj->fci', coefficient, z[faces]), axis=-1)
        before = float(gradients(initial).max())
        if before <= maximum_slope and np.all((initial >= lower) & (initial <= ceiling)):
            return mesh, dict(iterations=0, maximum=before, maximum_height_change=0., colors=0,
                              estimated_bytes=estimated, construction='fixed-chart height projection')
        free = ~fixed[faces]
        movable = coefficient * free[:, None, None, :]
        gram = np.einsum('fcij,fckj->fcik', movable, movable)
        values, vectors = np.linalg.eigh(gram)
        values = np.maximum(values, 0.)
        groups = face_colors(faces, fixed, self.max_colors)
        corrections = np.zeros((len(faces), 4, 4))
        box_correction = np.zeros(len(initial))
        heights = initial.copy()
        aim = maximum_slope
        maximum, movement = before, float('inf')
        for iteration in range(self.max_iterations):
            previous = heights.copy()
            for group in groups:
                ids = faces[group]
                for corner in range(4):
                    x = heights[ids] + corrections[group, corner]
                    projected = self.project(x, coefficient[group, corner], movable[group, corner],
                        values[group, corner], vectors[group, corner], aim)
                    projected[~free[group]] = initial[ids][~free[group]]
                    corrections[group, corner] = (x - projected) * free[group]
                    heights[ids] = projected
            x = heights + box_correction
            heights = np.clip(x, lower, ceiling)
            heights[fixed] = initial[fixed]
            box_correction = x - heights
            box_correction[fixed] = 0.
            maximum = float(gradients(heights).max())
            movement = float(np.max(abs(heights - previous)))
            converged = (maximum <= maximum_slope + self.slope_tolerance and movement <= self.height_tolerance
                         if feasibility_slope is None else maximum <= feasibility_slope)
            if converged:
                xyz = mesh.xyz.copy()
                xyz[:, 2] = heights
                return Mesh(xyz, faces, mesh.material, dict(mesh.seams)), dict(iterations=iteration + 1,
                    initial_maximum=before, maximum=maximum,
                    maximum_height_change=float(np.max(abs(heights - initial))),
                    squared_displacement=float(np.sum((heights - initial) ** 2)),
                    colors=len(groups), estimated_bytes=estimated,
                    **({} if feasibility_slope is None else dict(convergence='feasible_surface')),
                    construction='fixed-chart height projection')
        raise ValueError(f'Height projection did not converge within its iteration budget: slope {maximum:.6f}, step {movement:.3g}')


class ReferenceQuadGradient:
    def __init__(self, faces, vertices):
        self.targets = np.stack([
            np.stack((faces[:, i], faces[:, (i + 1) % 4], faces[:, (i + 2) % 4]))
            for i in range(4)]).reshape(-1).astype(np.intp, copy=False)
        self.contributions = np.empty((4, 3, len(faces), 2))
        self.right_a = np.empty((len(faces), 4, 2))
        self.right_b = np.empty_like(self.right_a)
        self.gradient = np.zeros((vertices, 2))

    def evaluate(self, energy):
        a, b = energy.edge, energy.next_edge
        b[:, :-1] = a[:, 1:]
        b[:, -1] = a[:, 0]
        self.right_a[..., 0], self.right_a[..., 1] = a[..., 1], -a[..., 0]
        self.right_b[..., 0], self.right_b[..., 1] = b[..., 1], -b[..., 0]
        denominator = energy.jac[..., None]
        numerator = energy.numerator[..., None]
        da = 2. * a / denominator - numerator * self.right_b / denominator ** 2
        db = 2. * b / denominator + numerator * self.right_a / denominator ** 2
        if energy.power == 2:
            factor = numerator / (2. * denominator)
            da, db = da * factor, db * factor
        self.contributions[:, 0] = -da.transpose(1, 0, 2)
        self.contributions[:, 1] = (da - db).transpose(1, 0, 2)
        self.contributions[:, 2] = db.transpose(1, 0, 2)
        # Both coordinate streams retain the original corner/contribution/face
        # accumulation order. minlength also clears unreferenced vertices.
        values = self.contributions.reshape(-1, 2)
        self.gradient[:, 0] = np.bincount(self.targets, weights=values[:, 0], minlength=len(self.gradient))
        self.gradient[:, 1] = np.bincount(self.targets, weights=values[:, 1], minlength=len(self.gradient))
        return self.gradient


class ReferenceBoundaryJetField:
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


class ReferenceQuadEnergy:
    def __init__(self, faces, power):
        self.faces, self.power = faces, power
        self.corners = np.empty((len(faces), 4, 2))
        self.edge = np.empty_like(self.corners)
        self.next_edge = np.empty_like(self.corners)
        self.jac = np.empty((len(faces), 4))
        self.numerator = np.empty_like(self.jac)
        self.value = float('inf')
        self.rejected_corner = None

    def evaluate(self, points):
        # A corner rejected on a previous trial is a cheap exact witness. Read
        # its current coordinates/connectivity: it is never a cached geometric
        # value or an approximate acceptance test. Only finite evaluations have
        # usable scratch buffers; an accepted proposal still evaluates all faces.
        if self.rejected_corner is not None:
            face, corner = self.rejected_corner
            ids = self.faces[face]
            a = points[ids[(corner + 1) % 4]] - points[ids[corner]]
            b = points[ids[(corner + 2) % 4]] - points[ids[(corner + 1) % 4]]
            if a[0] * b[1] - a[1] * b[0] <= 1e-10:
                self.value = float('inf')
                return self.value
        np.take(points, self.faces, axis=0, out=self.corners)
        np.subtract(self.corners[:, 1:], self.corners[:, :-1], out=self.edge[:, :-1])
        np.subtract(self.corners[:, 0], self.corners[:, -1], out=self.edge[:, -1])
        self.next_edge[:, :-1] = self.edge[:, 1:]
        self.next_edge[:, -1] = self.edge[:, 0]
        self.jac[:] = cross(self.edge, self.next_edge)
        if self.jac.min() <= 1e-10:
            self.rejected_corner = divmod(int(np.argmin(self.jac)), 4)
            self.value = float('inf')
            return self.value
        # Preserve the original squared-edge addition and reduction order.
        # next_edge is scratch after evaluation; QuadGradient restores it.
        np.square(self.edge, out=self.corners)
        np.square(self.next_edge, out=self.next_edge)
        np.add(self.corners, self.next_edge, out=self.corners)
        np.sum(self.corners, axis=-1, out=self.numerator)
        shape = self.numerator / self.jac
        self.value = float(np.sum(shape) if self.power == 1 else np.sum(shape ** 2) / 4.)
        return self.value
