"""Bounded fallback seeds for circular contacts on finite curved approaches.

Sampled center-offset curves supply guesses only. Acceptance always requires
the ordinary exact contact solve, prescribed radius and circular-sector guards.
Previously successful local contacts never enter this search.
"""
import numpy as np
from .mesh import cross


def contact_seeds(first, second, radius, *, samples=129, maximum=24, side_override=None):
    gap = np.arctan2(cross(first.direction, second.direction), first.direction @ second.direction) % (2. * np.pi)
    side = (1. if gap < np.pi else -1.) if side_override is None else side_override
    distances = (first.offset - .5 * first.width - side * radius,
                 second.offset + .5 * second.width + side * radius)
    lower = -4. * max(first.width + radius, second.width + radius)
    axes = [np.unique(np.r_[np.linspace(lower, 0., 9), np.linspace(0., approach.length, samples)])
            for approach in (first, second)]
    curves = [approach.point(axis, distance) for approach, axis, distance in zip((first, second), axes, distances)]
    a, b = curves
    da, db = np.diff(a, axis=0), np.diff(b, axis=0)
    offset = b[None, :-1] - a[:-1, None]
    denominator = cross(da[:, None], db[None])
    valid = abs(denominator) > 1e-12 * np.maximum(1.,
        np.linalg.norm(da, axis=1)[:, None] * np.linalg.norm(db, axis=1)[None])
    u = np.divide(cross(offset, db[None]), denominator, out=np.zeros_like(denominator), where=valid)
    v = np.divide(cross(offset, da[:, None]), denominator, out=np.zeros_like(denominator), where=valid)
    hits = np.argwhere(valid & (u >= 0.) & (u <= 1.) & (v >= 0.) & (v <= 1.))
    candidates = [(axes[0][i] + u[i, j] * (axes[0][i+1] - axes[0][i]),
                   axes[1][j] + v[i, j] * (axes[1][j+1] - axes[1][j])) for i, j in hits]
    # Tangent contacts can lack a transverse sampled intersection. The nearest
    # bounded pairs provide a few additional guesses, never an acceptance test.
    squared = np.sum((a[:, None] - b[None]) ** 2, axis=2)
    nearest = np.argsort(squared.ravel(), kind='stable')[:min(4, maximum)]
    candidates.extend((axes[0][i], axes[1][j]) for i, j in zip(*np.unravel_index(nearest, squared.shape)))
    candidates.sort(key=lambda pair: (max(0., pair[0]) + max(0., pair[1]),
                                       max(pair), pair[0], pair[1]))
    unique = []
    for pair in candidates:
        if not any(np.linalg.norm(np.asarray(pair) - old) < 1e-5 for old in unique):
            unique.append(np.asarray(pair, dtype=float))
        if len(unique) == maximum:
            break
    return tuple(unique)
