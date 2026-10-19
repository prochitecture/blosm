"""Bounded temporary storage for complete generated-surface validation."""
import numpy as np


def minimum_positive_quad_corner(xyz, faces, *, batch_size=16384):
    """Check every corner using the original XY cross product, in bounded batches.

    Batching does not sample or change the surface. Each face uses exactly the
    same arithmetic as the whole-array check; only the lifetime of temporary
    edge and corner arrays changes. Callers retain the original coordinate and
    connectivity arrays for any subsequent ownership checks.
    """
    if type(batch_size) is not int or batch_size <= 0:
        raise ValueError('Quad validation batch size must be a positive integer')
    if faces.ndim != 2 or faces.shape[1] != 4 or not len(faces):
        raise ValueError('Quad validation requires nonempty four-corner faces')
    if not np.all(np.isfinite(xyz)):
        raise ValueError('Original-span surface has non-finite coordinates')
    minimum = float('inf')
    for start in range(0, len(faces), batch_size):
        points = xyz[faces[start:start+batch_size], :2]
        edges = np.roll(points, -1, axis=1)-points
        following = np.roll(edges, -1, axis=1)
        cross = edges[:, :, 0]*following[:, :, 1]-edges[:, :, 1]*following[:, :, 0]
        local = float(np.min(cross))
        if not np.isfinite(local) or local <= 0:
            raise ValueError('Original-span surface has a folded or collapsed quad')
        minimum = min(minimum, local)
    return minimum
