"""Spatial broad phase and batched separating-axis tests for pavement quads."""
import numpy as np


def bounds_overlap(a, b):
    return all(a[0][k] <= b[1][k] and b[0][k] <= a[1][k] for k in (0, 1))


def _interiors_overlap(first, second, tolerance):
    # Both inputs have already passed convex CCW quad validation. Subtract a
    # local origin before projection to avoid cancellation on large map frames.
    origin = first[:, :1, :]
    a, b = first - origin, second - origin
    edges = np.concatenate((np.roll(a, -1, axis=1) - a, np.roll(b, -1, axis=1) - b), axis=1)
    axes = np.stack((edges[..., 1], -edges[..., 0]), axis=-1)
    axes /= np.linalg.norm(axes, axis=-1, keepdims=True)
    pa = np.einsum('nac,nvc->nav', axes, a)
    pb = np.einsum('nac,nvc->nav', axes, b)
    depth = np.minimum(pa.max(axis=2), pb.max(axis=2)) - np.maximum(pa.min(axis=2), pb.min(axis=2))
    return np.min(depth, axis=1) > tolerance


_CELL_DTYPE = np.dtype([('x', np.int64), ('y', np.int64)])
_ENTRY_BATCH = 65536
_PAIR_BATCH = 65536
_SAT_BATCH = 1024


def _keys(cells):
    return np.ascontiguousarray(cells, dtype=np.int64).view(_CELL_DTYPE).reshape(-1)


def _entries(first, spans, indices):
    """Expand quad/grid incidences in bounded NumPy batches."""
    counts = spans[indices, 0] * spans[indices, 1]
    ends = np.cumsum(counts, dtype=np.int64)
    at = 0
    while at < len(indices):
        previous = int(ends[at - 1]) if at else 0
        stop = max(at + 1, int(np.searchsorted(ends, previous + _ENTRY_BATCH, side='right')))
        sizes = counts[at:stop]
        owners = np.repeat(indices[at:stop], sizes)
        offset = np.arange(int(sizes.sum()), dtype=np.int64)
        offset -= np.repeat(np.r_[0, np.cumsum(sizes[:-1], dtype=np.int64)], sizes)
        cells = np.column_stack((first[owners, 0] + offset // spans[owners, 1],
                                 first[owners, 1] + offset % spans[owners, 1]))
        yield owners, cells
        at = stop


class _CellGrid:
    """CSR grid containing changed faces only; no unchanged/unchanged work."""
    def __init__(self, first, spans, indices, max_entries):
        count = int(np.sum(spans[indices, 0] * spans[indices, 1]))
        if count > max_entries:
            raise ValueError('Overlap index exceeds its entry budget; increase cell size or partition the network')
        coordinates = np.empty((count, 2), dtype=np.int64)
        owners = np.empty(count, dtype=np.int32)
        at = 0
        for faces, cells in _entries(first, spans, indices):
            stop = at + len(faces)
            coordinates[at:stop], owners[at:stop] = cells, faces
            at = stop
        order = np.lexsort((owners, coordinates[:, 1], coordinates[:, 0]))
        coordinates, owners = coordinates[order], owners[order]
        starts = np.flatnonzero(np.r_[True, np.any(coordinates[1:] != coordinates[:-1], axis=1)])
        self.coordinates = coordinates[starts]
        self.keys = _keys(self.coordinates)
        self.offsets = np.r_[starts, len(owners)]
        self.owners = owners
        self.entries = count

    def pairs(self, first, spans, indices):
        for faces, cells in _entries(first, spans, indices):
            query = _keys(cells)
            positions = np.searchsorted(self.keys, query)
            inside = positions < len(self.keys)
            inside &= self.keys[np.minimum(positions, len(self.keys) - 1)] == query
            faces, positions = faces[inside], positions[inside]
            sizes = self.offsets[positions + 1] - self.offsets[positions]
            # A dense cell can itself exceed the pair workspace. Stream that
            # row rather than allocating its Cartesian product at once.
            ends = np.cumsum(sizes, dtype=np.int64)
            at = 0
            while at < len(faces):
                if sizes[at] > _PAIR_BATCH:
                    start, finish = self.offsets[positions[at]:positions[at] + 2]
                    for lo in range(int(start), int(finish), _PAIR_BATCH):
                        partners = self.owners[lo:min(lo + _PAIR_BATCH, finish)]
                        yield np.full(len(partners), faces[at], dtype=np.int32), partners, positions[at]
                    at += 1
                    continue
                previous = int(ends[at - 1]) if at else 0
                stop = max(at + 1, int(np.searchsorted(ends, previous + _PAIR_BATCH, side='right')))
                counts = sizes[at:stop]
                selected = np.repeat(np.arange(at, stop), counts)
                offset = np.arange(int(counts.sum()), dtype=np.int64)
                offset -= np.repeat(np.r_[0, np.cumsum(counts[:-1], dtype=np.int64)], counts)
                slots = positions[selected]
                partners = self.owners[self.offsets[slots] + offset]
                yield faces[selected], partners, slots
                at = stop


def validate_overlaps(meshes, *, changed=None, cell_size=2., tolerance=1e-6,
                      max_index_entries=8_000_000):
    """Test changed parts against each other, themselves and nearby old parts.

    The index is rebuilt from supplied geometry, so a newly introduced spatial
    interaction cannot be missed because it was absent from an old dependency
    graph. Persistent topology still decides which feature owns each surface.
    """
    if (not np.isfinite(cell_size) or cell_size <= 0. or not np.isfinite(tolerance) or tolerance < 0.
            or type(max_index_entries) is not int or max_index_entries < 1):
        raise ValueError('Invalid overlap-index limits')
    if not meshes:
        return dict(quads=0, candidate_pairs=0)
    keys = list(meshes)
    changed = set(keys) if changed is None else set(changed)
    arrays = [meshes[key].xyz[meshes[key].quads, :2] for key in keys]
    points = np.concatenate(arrays)
    owners = np.concatenate([np.full(len(a), i, dtype=np.int32) for i, a in enumerate(arrays)])
    selected = np.array([key in changed for key in keys], dtype=bool)[owners]
    if not np.any(selected):
        return dict(quads=len(points), candidate_pairs=0)
    minimum, maximum = points.min(axis=1), points.max(axis=1)
    bins = np.floor(np.stack((minimum, maximum)) / cell_size)
    if not np.all(np.isfinite(bins)) or np.max(abs(bins)) > 2 ** 61:
        raise ValueError('Surface coordinates exceed the overlap-index integer range')
    first, last = bins.astype(np.int64)
    spans = last - first + 1
    if np.any(spans > 4096) or np.any(spans[:, 0] * spans[:, 1] > 4096):
        raise ValueError('Surface quad exceeds the overlap-index spatial budget')
    grid = _CellGrid(first, spans, np.flatnonzero(selected), max_index_entries)
    pair_count, batches = 0, 0
    for target, partner, slots in grid.pairs(first, spans, np.arange(len(points), dtype=np.int32)):
        # Choose the lexicographically first common cell, so a face pair is
        # tested once even when its boxes share many cells. Both changed faces
        # use their original order; an unchanged face can query any changed one.
        active = (target > partner) | ~selected[target]
        active &= np.all(np.maximum(first[target], first[partner]) == grid.coordinates[slots], axis=1)
        active &= np.all(maximum[partner] - minimum[target] > tolerance, axis=1)
        active &= np.all(maximum[target] - minimum[partner] > tolerance, axis=1)
        left, right = np.maximum(target[active], partner[active]), np.minimum(target[active], partner[active])
        pair_count += len(left)
        for start in range(0, len(left), _SAT_BATCH):
            a, b = left[start:start + _SAT_BATCH], right[start:start + _SAT_BATCH]
            # Preserve the old SAT argument order and local origin, including
            # comparisons close to the configured overlap tolerance.
            hit = _interiors_overlap(points[a], points[b], tolerance)
            batches += 1
            if np.any(hit):
                at = int(np.flatnonzero(hit)[0])
                raise ValueError(f'Pavement interiors overlap: {keys[owners[a[at]]].text()} '
                                 f'and {keys[owners[b[at]]].text()}')
    return dict(quads=len(points), candidate_pairs=pair_count,
                index_entries=grid.entries, sat_batches=batches)
