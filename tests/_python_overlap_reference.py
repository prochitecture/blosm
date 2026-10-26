"""Frozen pre-batching validator for same-process performance comparisons only."""
import numpy as np
from blosm.way.street_network.python_backend.validation import _interiors_overlap


def validate_reference(meshes, *, changed=None, cell_size=2., tolerance=1e-6):
    if not meshes:
        return dict(quads=0, candidate_pairs=0)
    keys = list(meshes)
    changed = set(keys) if changed is None else set(changed)
    arrays = [meshes[key].xyz[meshes[key].quads, :2] for key in keys]
    points = np.concatenate(arrays)
    owners = np.concatenate([np.full(len(a), i, dtype=np.int32) for i, a in enumerate(arrays)])
    selected = np.array([keys[i] in changed for i in owners])
    minimum, maximum = points.min(axis=1), points.max(axis=1)
    cells, pair_count = {}, 0
    for index, (lo, hi) in enumerate(zip(minimum, maximum)):
        first = np.floor(lo / cell_size).astype(np.int64)
        last = np.floor(hi / cell_size).astype(np.int64)
        if np.prod(last - first + 1) > 4096:
            raise ValueError('Surface quad exceeds the overlap-index spatial budget')
        slots = [(x, y) for x in range(int(first[0]), int(last[0]) + 1)
                 for y in range(int(first[1]), int(last[1]) + 1)]
        candidates = set()
        for slot in slots:
            candidates.update(cells.get(slot, ()))
        if candidates:
            partner = np.fromiter(candidates, dtype=np.int32)
            active = selected[partner] | selected[index]
            active &= np.all(maximum[partner] - minimum[index] > tolerance, axis=1)
            active &= np.all(maximum[index] - minimum[partner] > tolerance, axis=1)
            partner = partner[active]
            pair_count += len(partner)
            for start in range(0, len(partner), 1024):
                batch = partner[start:start + 1024]
                hit = _interiors_overlap(np.broadcast_to(points[index], (len(batch), 4, 2)), points[batch], tolerance)
                if np.any(hit):
                    other = int(batch[np.flatnonzero(hit)[0]])
                    raise ValueError(f'Pavement interiors overlap: {keys[owners[index]].text()} '
                                     f'and {keys[owners[other]].text()}')
        for slot in slots:
            cells.setdefault(slot, []).append(index)
    return dict(quads=len(points), candidate_pairs=pair_count)
