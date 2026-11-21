"""Grade sparse defects in a positive planar pavement with its physical rims fixed."""
import numpy as np

from .four_side_charts import quality
from .mesh import Mesh, cross, publication_coordinates


def _bad(xy, faces, *, ratio=8., sine=.03):
    edges = np.roll(xy[faces], -1, axis=1) - xy[faces]
    lengths = np.linalg.norm(edges, axis=2)
    turns = cross(edges, np.roll(edges, -1, axis=1))
    return ((lengths.max(1) > ratio * lengths.min(1)) |
            (turns <= 1e-10).any(1) |
            (turns < sine * lengths * np.roll(lengths, -1, axis=1)).any(1))


def _objective(base, faces, free, scale, offsets, target_sine):
    xy = base.copy()
    xy[free] += offsets.reshape(-1, 2) * scale[:, None]
    edges = np.roll(xy[faces], -1, axis=1) - xy[faces]
    lengths = np.sqrt(np.sum(edges * edges, axis=2) + 1e-24)
    unit = edges / lengths[:, :, None]
    following = np.roll(unit, -1, axis=1)
    # Penalizing every violating edge pair gives a continuous first derivative
    # when two edges exchange their order as the longest or shortest edge.
    logs = np.log(lengths)
    excess = np.maximum(0., logs[:, :, None] - logs[:, None, :] - np.log(7.75))
    grad = 2. * (excess.sum(2) - excess.sum(1))[:, :, None] * unit / lengths[:, :, None]
    sine = cross(unit, following)
    angle = np.minimum(0., sine - target_sine)
    factor = 150. * angle[:, :, None]
    right = np.stack((unit[:, :, 1], -unit[:, :, 0]), axis=2)
    grad += factor * (np.roll(right, -1, axis=1) - sine[:, :, None]*unit) / lengths[:, :, None]
    grad += np.roll(factor * (-right - sine[:, :, None]*following), 1, axis=1) / lengths[:, :, None]
    vertex_grad = np.roll(grad, 1, axis=1) - grad
    gradient = np.column_stack([np.bincount(faces.ravel(), weights=g.ravel(), minlength=len(base))
                                for g in (vertex_grad[:, :, 0], vertex_grad[:, :, 1])])
    return float(np.sum(excess*excess) + 75.*np.sum(angle*angle)), (gradient[free]*scale[:, None]).ravel()


def _minimize(evaluate, size):
    """Limited-memory BFGS with bounded descent search, using NumPy only."""
    x = np.zeros(size)
    value, gradient = evaluate(x)
    history = []
    for iteration in range(400):
        if value < 1e-14 or np.max(abs(gradient)) < 1e-10:
            break
        direction = gradient.copy()
        weights = []
        for step, change, inverse in reversed(history):
            alpha = inverse * (step @ direction)
            direction -= alpha * change
            weights.append(alpha)
        if history:
            step, change, inverse = history[-1]
            direction /= inverse * (change @ change)
        for (step, change, inverse), alpha in zip(history, reversed(weights)):
            direction += step * (alpha - inverse * (change @ direction))
        descent = float(gradient @ direction)
        if not np.isfinite(descent) or descent <= 0.:
            history.clear()
            direction = gradient.copy()
            descent = float(gradient @ gradient)
        distance = 1.
        for _ in range(24):
            delta = -distance*direction
            next_value, next_gradient = evaluate(x+delta)
            if next_value <= value - 1e-4*distance*descent:
                break
            distance *= .5
        else:
            break
        change = next_gradient-gradient
        curvature = float(change @ delta)
        if curvature > max(1e-20, 1e-12*np.linalg.norm(delta)*np.linalg.norm(change)):
            history.append((delta, change, 1./curvature))
            if len(history) > 8:
                history.pop(0)
        x += delta
        value, gradient = next_value, next_gradient
    return x, iteration+1


def _adjust(xy, faces, fixed, bad, rings, target_sine):
    affected = bad.copy()
    for _ in range(rings):
        movable = np.zeros(len(xy), dtype=bool)
        movable[faces[affected].ravel()] = True
        movable[fixed] = False
        affected = movable[faces].any(1)
    if not affected.any() or affected.sum() > 8192:
        return 0
    ids, inverse = np.unique(faces[affected], return_inverse=True)
    local = inverse.reshape(-1, 4)
    base = xy[ids] - xy[ids[0]]
    free = np.flatnonzero(movable[ids])
    lengths = np.linalg.norm(np.roll(base[local], -1, axis=1) - base[local], axis=2)
    scale = np.full(len(base), np.inf)
    np.minimum.at(scale, local, lengths)
    np.minimum.at(scale, np.roll(local, -1, axis=1), lengths)
    step, iterations = _minimize(lambda x: _objective(base, local, free, scale[free], x, target_sine), 2*len(free))
    xy[ids[free]] += step.reshape(-1, 2)*scale[free, None]
    return iterations


def _collars(xy, faces, marked):
    # One quad per existing edge and a central quad. All four original edges
    # remain intact, so neighbors and any attached profiles still share them.
    boundary = faces[marked]
    corners = xy[boundary]
    inner = corners + .32*(corners.mean(1)[:, None, :] - corners)
    ids = np.arange(len(xy), len(xy)+4*len(marked), dtype=np.int32).reshape(-1, 4)
    rings = np.stack((boundary, np.roll(boundary, -1, axis=1), np.roll(ids, -1, axis=1), ids), axis=2)
    faces = faces.copy()
    faces[marked] = ids
    return np.vstack((xy, inner.reshape(-1, 2))), np.vstack((faces, rings.reshape(-1, 4)))


def _fixed_rim(mesh, max_cells):
    mesh.validate()
    if (mesh.uv_layers or np.ptp(mesh.xyz[:, 2]) != 0. or len(mesh.quads) > max_cells
            or not mesh.seams):
        raise ValueError('Local grading requires a bounded planar chart before UV construction')
    faces = mesh.quads
    edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=2).reshape(-1, 2), axis=1)
    edges, counts = np.unique(edges, axis=0, return_counts=True)
    boundary = np.unique(edges[counts == 1])
    seams = np.unique(np.concatenate(list(mesh.seams.values())))
    if counts.max() > 2 or not np.array_equal(boundary, seams):
        raise ValueError('Local grading requires every physical rim and a manifold pavement')
    return boundary


def _independent_vertices(faces, fixed):
    """Vertices of one color have disjoint incident faces, even at diagonals."""
    incident = [[] for _ in fixed]
    for fi, face in enumerate(faces):
        for corner, vertex in enumerate(face):
            incident[vertex].append((fi, corner))
    if any(not row for row in incident):
        raise ValueError('Independent vertex relaxation requires every vertex to belong to a quad')
    colors = np.full(len(fixed), -1, dtype=np.int32)
    for vertex in np.flatnonzero(~fixed):
        occupied = {colors[other] for fi, _ in incident[vertex] for other in faces[fi]}
        colors[vertex] = next(i for i in range(len(occupied)+1) if i not in occupied)
    for color in range(int(colors.max())+1):
        vertices = np.flatnonzero(colors == color)
        fi, corner, owner = np.asarray([(fi, corner, i) for i, vertex in enumerate(vertices)
            for fi, corner in incident[vertex]]).T
        yield vertices, faces[fi], corner, owner


def grade_dense_steps(mesh, max_cells):
    """Relax a positive pavement with independent step lengths, then repair sparse defects."""
    from .quality_correction import QuadShapeCorrection
    from .updates import numerical_time_slices
    boundary = _fixed_rim(mesh, max_cells)
    fixed = np.zeros(len(mesh.xyz), dtype=bool)
    fixed[boundary] = True
    xy, faces = mesh.xyz[:, :2].copy(), mesh.quads
    initial = int(_bad(xy, faces).sum())
    if mesh.material != 'pavement' or initial <= 256:
        raise ValueError('Dense relaxation requires pavement beyond the sparse repair limit')
    families = list(_independent_vertices(faces, fixed))
    quality_guard = QuadShapeCorrection(max_cells=max_cells)
    def relax():
        for sweep in range(64):
            changed = False
            for vertices, cells, corner, owner in families:
                yield
                _, _, gradient = quality_guard.evaluate(xy, cells, gradient=True)
                if gradient is None:
                    raise ValueError('Independent vertex relaxation requires positive noncollapsed cells')
                points = xy[cells]
                rows = np.arange(len(cells))
                pivot = points[rows, corner]
                before = points[rows, (corner-1)%4]-pivot
                after = points[rows, (corner+1)%4]-pivot
                degree = np.bincount(owner)
                scale2 = np.bincount(owner, weights=.5*np.sum(before**2+after**2, axis=1))/degree
                move = -gradient[vertices]*(scale2/degree)[:, None]
                magnitude = np.linalg.norm(move, axis=1)
                move *= np.minimum(1., .1*np.sqrt(scale2)/np.maximum(magnitude, 1e-100))[:, None]
                previous = np.bincount(owner, weights=quality_guard.face_losses(points))
                active = (magnitude>1e-12) & (previous>1e-16)
                for _ in range(18):
                    if not active.any():
                        break
                    yield
                    touched = np.flatnonzero(active[owner])
                    trial = points[touched].copy()
                    trial[np.arange(len(touched)),corner[touched]] += move[owner[touched]]
                    losses = np.bincount(owner[touched], weights=quality_guard.face_losses(trial),
                                         minlength=len(vertices))
                    accepted = active & (losses < previous*(1.-1e-12))
                    xy[vertices[accepted]] += move[accepted]
                    changed |= bool(accepted.any())
                    active &= ~accepted
                    move[active] *= .5
            if not changed:
                break
        return sweep+1
    sweeps = yield from numerical_time_slices(relax())
    relaxed = mesh.with_positions(np.column_stack((xy,mesh.xyz[:,2])))
    result, details = grade(relaxed, max_cells)
    if not np.array_equal(result.xyz[boundary], mesh.xyz[boundary]):
        raise ValueError('Independent vertex relaxation moved the physical rim')
    return result, dict(details, independent_vertex_sweeps=sweeps, initial_defects=initial,
                        relaxed_defects=int(_bad(xy,faces).sum()))


def grade(mesh, max_cells, *, narrow_grass_disk=False):
    boundary = _fixed_rim(mesh, max_cells)
    xy = mesh.xyz[:, :2].copy()
    faces = mesh.quads
    disk = mesh.material == 'grass' and len(mesh.seams) == 1 and len(faces) <= 8192
    if narrow_grass_disk and not disk:
        raise ValueError('Narrow-tip grading requires a bounded grass disk')
    # The disk contract requires positive stored corners and the same 8:1
    # edge bound. A small positive angle margin accommodates a tapered tip.
    minimum_sine, target_sine = (.001, .008) if narrow_grass_disk else (.03, .06)
    fixed = np.zeros(len(xy), dtype=bool)
    fixed[boundary] = True
    attempts = []
    added = 0
    for stage in range(2):
        for rings in (3, 6, 10):
            stored = publication_coordinates(xy, 1., physical=True)
            bad = _bad(xy, faces, sine=minimum_sine) | _bad(stored, faces, sine=minimum_sine)
            if not bad.any():
                break
            if not disk and bad.sum() > 256:
                raise ValueError('Local grading exceeds the sparse defect budget')
            iterations = _adjust(xy, faces, fixed, bad, rings, target_sine)
            attempts.append(dict(rings=rings, failed_cells=int(bad.sum()), iterations=iterations))
        stored = publication_coordinates(xy, 1., physical=True)
        marked = np.flatnonzero(_bad(xy, faces, sine=minimum_sine) | _bad(stored, faces, sine=minimum_sine))
        if not len(marked):
            xyz = np.column_stack((xy, np.full(len(xy), mesh.xyz[0, 2])))
            result = Mesh(xyz, faces, mesh.material, dict(mesh.seams))
            physical = quality(result)
            published = quality(result.with_positions(publication_coordinates(xyz, 1., physical=True)))
            if not np.array_equal(xyz[boundary], mesh.xyz[boundary]):
                raise ValueError('Local grading moved a physical rim vertex')
            return result, dict(qualified=True, local_adjustments=attempts, inserted_quads=added,
                                physical_quality=physical, stored_quality=published, exact_fixed_rim=True)
        added = 4*len(marked)
        if stage or len(marked) > 64 or len(faces)+added > max_cells:
            break
        xy, faces = _collars(xy, faces, marked)
        fixed = np.r_[fixed, np.zeros(added, dtype=bool)]
    raise ValueError('Local quad adjustments did not qualify the assembled pavement')
