"""Direct quad layers for visibility-star-shaped planar blocks.

This is a mesh chart building block, not a compound surface provider. A fixed
physical rim feeds positive radial cells and cached 3:1 transition topology.
Bounded local quality improvement moves only internal vertices. Heights and
compound decomposition are separate, mandatory construction stages.
"""
from dataclasses import dataclass
import math
import numpy as np
from .mesh import Mesh, cross
from .model import Key
from .quad_workspace import QuadEnergy, QuadGradient
from .regions import Region, RegionUnion, signed_area


def _quality(points):
    edge = np.roll(points, -1, axis=-2) - points
    jac = cross(edge, np.roll(edge, -1, axis=-2))
    lengths = np.linalg.norm(edge, axis=-1)
    return jac, lengths.max(axis=-1) / lengths.min(axis=-1)


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


class RadialTransition:
    """Place the existing eight-vertex, four-quad 3:1 template in physical XY.

    The two internal points lie on the two intermediate rays and on a line
    parallel to the coarse inner edge. This gives three positive radial quads
    and one positive trapezoid, including along concave exterior curb samples.
    """
    def __init__(self, patterns):
        self.pattern = patterns.strip(1, 3, 1)
        self.fine = patterns.boundary(self.pattern, 3)
        self.coarse = patterns.boundary(self.pattern, 1)
        inside = np.flatnonzero((self.pattern.uv[:, 0] > 0.) & (self.pattern.uv[:, 0] < 1.))
        self.inside = inside[np.argsort(self.pattern.uv[inside, 1])]
        self.quads = self.pattern.quads[:, ::-1]

    def candidates(self, outer, inner):
        """Coordinates are relative to the chart's visibility centre."""
        edge = inner[1] - inner[0]
        length = np.linalg.norm(edge)
        if length <= 1e-10:
            return ()
        normal = np.array((edge[1], -edge[0])) / length
        lower, upper = float(normal @ inner[0]), float(np.min(outer[1:3] @ normal))
        if lower <= 0. or upper <= lower + 1e-8:
            return ()
        direction = outer[1:3] / np.linalg.norm(outer[1:3], axis=1)[:, None]
        denominator = direction @ normal
        if np.min(denominator) <= 1e-12:
            return ()
        result = []
        for fraction in (.25, .5, .75):
            points = np.empty((8, 2))
            points[self.fine], points[self.coarse] = outer, inner
            points[self.inside] = (lower + fraction * (upper - lower)) * direction / denominator[:, None]
            jac, ratio = _quality(points[self.quads])
            if np.min(jac) > 1e-10:
                result.append((float(ratio.max()), points))
        return tuple(result)


class FixedRimQuality:
    """Bounded local shape-energy descent; every accepted step stays positive."""
    def __init__(self, *, iterations=300, target_ratio=7.5, power=1):
        if (type(iterations) is not int or iterations < 0 or not math.isfinite(target_ratio)
                or target_ratio < 1. or type(power) is not int or power not in (1, 2)):
            raise ValueError('Invalid local quad-quality settings')
        self.iterations, self.target_ratio = iterations, target_ratio
        self.power = power

    def improve(self, mesh):
        xy, faces = mesh.xyz[:, :2].copy(), mesh.quads
        edge_ids = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edge_ids, axis=0, return_counts=True)
        fixed = np.zeros(len(xy), bool)
        fixed[unique[counts == 1].ravel()] = True
        src, dst = np.r_[unique[:, 0], unique[:, 1]], np.r_[unique[:, 1], unique[:, 0]]
        degree = np.bincount(dst, minlength=len(xy))

        # A rejected proposal must not overwrite the accepted state's buffers.
        # Workspaces are bounded by this chart and released after this call.
        current, trial = QuadEnergy(faces, self.power), QuadEnergy(faces, self.power)
        scatter = QuadGradient(faces, len(xy))
        value = current.evaluate(xy)
        if not np.isfinite(value):
            raise ValueError('Quality improvement requires a positive initial mesh')
        accepted_steps = 0
        for iteration in range(self.iterations):
            edge, jac = current.edge, current.jac
            lengths = np.linalg.norm(edge, axis=-1)
            sine = jac / (lengths * np.roll(lengths, -1, axis=1))
            if np.max(lengths.max(axis=1) / lengths.min(axis=1)) <= self.target_ratio and sine.min() >= .1:
                break
            gradient = scatter.evaluate(current)
            length2 = np.sum((xy[src] - xy[dst]) ** 2, axis=1)
            mean2 = np.bincount(dst, weights=length2, minlength=len(xy)) / degree
            move = -gradient * mean2[:, None] / degree[:, None]
            move[fixed] = 0.
            norm = np.linalg.norm(move, axis=1)
            cap = .2 * np.sqrt(mean2)
            move *= np.minimum(1., np.divide(cap, norm, out=np.ones_like(norm), where=norm > 0.))[:, None]
            for power in range(16):
                proposal = xy + (2. ** -power) * move
                proposed = trial.evaluate(proposal)
                if proposed < value * (1. - 1e-12):
                    xy, value = proposal, proposed
                    current, trial = trial, current
                    accepted_steps += 1
                    break
            else:
                break
        result = Mesh(np.column_stack((xy, mesh.xyz[:, 2])), faces, mesh.material, dict(mesh.seams))
        return result, dict(steps=accepted_steps, energy=value,
            boundary_error=float(np.max(np.linalg.norm(xy[fixed] - mesh.xyz[fixed, :2], axis=1))))


@dataclass
class ChartResult:
    mesh: Mesh
    metadata: dict


class QuadChartQualityError(ValueError):
    """A positive, bounded chart that failed its final shape-quality guard."""
    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


class LayeredQuadChart:
    def __init__(self, patterns, *, spacing=.7, max_cells=50000, max_layers=256,
                 quality_iterations=300, maximum_ratio=8.):
        if (not math.isfinite(spacing) or spacing <= 0. or type(max_cells) is not int or max_cells < 4
                or type(max_layers) is not int or max_layers < 1 or type(quality_iterations) is not int
                or quality_iterations < 0 or not math.isfinite(maximum_ratio) or maximum_ratio < 1.):
            raise ValueError('Invalid direct-quad chart settings')
        self.patterns, self.spacing = patterns, spacing
        self.max_cells, self.max_layers = max_cells, max_layers
        self.quality = FixedRimQuality(iterations=quality_iterations, target_ratio=min(7.5, maximum_ratio))
        self.maximum_ratio = maximum_ratio

    def build(self, boundary, *, require_quality=True):
        boundary = np.asarray(boundary, dtype=float)
        RegionUnion().build((Region(Key('chart', 0), boundary),))
        if len(boundary) % 2 or len(boundary) < 4:
            raise ValueError('Direct quad chart needs an even physical rim; refine an owned curb edge')
        if len(boundary) // 2 - 1 > self.max_cells:
            raise ValueError('Direct quad chart rim exceeds its cell budget')
        center = visibility_center(boundary)
        relative = boundary - center
        points, faces = list(relative.copy()), []
        front = list(range(len(points)))
        transition = RadialTransition(self.patterns)
        coarsened, layers = 0, 0
        for _ in range(self.max_layers):
            n = len(front)
            a = np.asarray([points[i] for i in front])
            radii = np.linalg.norm(a, axis=1)
            direction = a / radii[:, None]
            gap = np.arctan2(cross(direction, np.roll(direction, -1, axis=0)),
                            np.sum(direction * np.roll(direction, -1, axis=0), axis=1))
            if np.min(gap) <= 1e-10:
                raise ValueError('Block rays cross or touch the physical boundary more than once')
            edge_lengths = np.linalg.norm(np.roll(a, -1, axis=0) - a, axis=1)
            fraction = max(1e-4, min(.30, self.spacing / radii.max(),
                float(np.min(3. * np.minimum(edge_lengths, np.roll(edge_lengths, 1)) / radii))))
            b = (1. - fraction) * a
            final = n <= 6 and radii.max() <= self.spacing
            if final:
                b = .45 * radii.min() * direction
            candidates = []
            if n > 6:
                for i in range(n):
                    indices = [(i + k) % n for k in range(4)]
                    if sum(gap[indices[:3]]) >= math.pi - .01:
                        continue
                    inner = b[[indices[0], indices[3]]]
                    length = float(np.linalg.norm(inner[1] - inner[0]))
                    if length > 2. * self.spacing or np.max(edge_lengths[indices[:3]]) >= 1.05 * self.spacing:
                        continue
                    variants = transition.candidates(a[indices], inner)
                    if variants:
                        ratio, nodes = min(variants, key=lambda variant: round(variant[0], 7))
                        score = (length / self.spacing - 1.) ** 2 + .2 * ratio
                        candidates.append((score, i, nodes))
            selected, occupied = {}, set()
            for score, i, nodes in sorted(candidates, key=lambda candidate: (round(candidate[0], 7), candidate[1])):
                edges = {(i + k) % n for k in range(3)}
                if not edges & occupied and n - 2 * (len(selected) + 1) >= 6:
                    selected[i] = nodes
                    occupied.update(edges)
            start = min(selected) if selected else 0
            groups, consumed, index = [], 0, start
            while consumed < n:
                count = 3 if index in selected else 1
                groups.append((index, count))
                consumed += count
                index = (index + count) % n
            if consumed != n:
                raise ValueError('Cyclic transition selection does not partition the front')
            if len(faces) + sum(4 if count == 3 else 1 for i, count in groups) + 2 > self.max_cells:
                raise ValueError('Direct quad chart exceeds its cell budget')
            new_ids = {}
            for i in sorted(i for i, count in groups):
                new_ids[i] = len(points)
                points.append(b[i])
            for i, count in groups:
                end = (i + count) % n
                if count == 1:
                    faces.append((front[i], front[end], new_ids[end], new_ids[i]))
                else:
                    nodes = selected[i]
                    ids = np.empty(8, dtype=np.int32)
                    ids[transition.fine] = [front[(i + k) % n] for k in range(4)]
                    ids[transition.coarse] = (new_ids[i], new_ids[end])
                    ids[transition.inside] = (len(points), len(points) + 1)
                    points.extend(nodes[transition.inside])
                    faces.extend(ids[transition.quads])
                    coarsened += 1
            front = [new_ids[i] for i, count in groups]
            layers += 1
            if len(faces) > self.max_cells:
                raise ValueError('Direct quad chart exceeds its cell budget')
            if final:
                break
        else:
            raise ValueError('Direct quad chart exceeds its layer budget')
        if len(front) == 4:
            faces.append(tuple(front))
        elif len(front) == 6:
            variants = []
            for i in range(6):
                ring = front[i:] + front[:i]
                pair = np.array((ring[:4], (ring[0], ring[3], ring[4], ring[5])))
                jac, ratio = _quality(np.asarray(points)[pair])
                if np.min(jac) > 1e-10:
                    variants.append((float(ratio.max()), i, pair))
            if not variants:
                raise ValueError('No positive six-vertex direct quad core')
            faces.extend(min(variants, key=lambda variant: (round(variant[0], 7), variant[1]))[2])
        else:
            raise ValueError('Direct quad core has unsupported parity')
        mesh = Mesh(np.column_stack((np.asarray(points) + center, np.zeros(len(points)))), np.asarray(faces))
        # The local numerical frame must not introduce even a rounding change
        # at a shared source rim. These coordinates remain authoritative.
        mesh.xyz[:len(boundary), :2] = boundary
        mesh.seams = {0: np.arange(len(boundary))}
        initial = mesh.validate()
        mesh, improvement = self.quality.improve(mesh)
        quality = mesh.validate()
        rotations = None
        if require_quality and self.maximum_ratio < quality['maximum_edge_ratio'] <= 2. * self.maximum_ratio:
            from .quad_rotations import QuadEdgeRotations
            mesh, rotations = QuadEdgeRotations(maximum_ratio=self.maximum_ratio,
                max_cells=self.max_cells).improve(mesh,
                    FixedRimQuality(iterations=min(100, self.quality.iterations)))
            quality = mesh.validate()
        points = mesh.xyz[mesh.quads, :2]
        edges = np.roll(points, -1, axis=1) - points
        lengths = np.linalg.norm(edges, axis=-1)
        quality['minimum_corner_sine'] = float(np.min(cross(edges, np.roll(edges, -1, axis=1))
            / (lengths * np.roll(lengths, -1, axis=1))))
        result = ChartResult(mesh, dict(layers=layers, coarsened=coarsened, initial=initial,
            quality=quality, improvement=improvement, edge_rotations=rotations,
            construction='direct radial layers with cached 3:1 quad blocks'))
        if require_quality and quality['maximum_edge_ratio'] > self.maximum_ratio:
            raise QuadChartQualityError(
                f'Direct quad chart needs better grading: edge ratio {quality["maximum_edge_ratio"]:.3f}', result)
        return result
