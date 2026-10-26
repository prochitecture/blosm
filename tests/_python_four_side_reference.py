"""Scalar corner-search oracle retained from dev48 (source a0543f7d...)."""
import itertools
import math
import numpy as np
from tested_python_streets.mesh import cross


def corners(boundary):
    n = len(boundary)
    if n == 4:
        return [(0., (0, 1, 2, 3))]
    edge = np.roll(boundary, -1, axis=0) - boundary
    lengths = np.linalg.norm(edge, axis=1)
    unit = edge / lengths[:, None]
    previous = np.roll(unit, 1, axis=0)
    turns = np.arctan2(cross(previous, unit), np.sum(previous * unit, axis=1))
    chosen = sorted(np.flatnonzero(turns > .05),
                    key=lambda i: (-round(float(turns[i]), 8), int(i)))[:16]
    choices = sorted({(int(i) + shift) % n for i in chosen for shift in (-1, 0, 1)})
    if len(choices) > 36:
        choices = sorted(sorted(choices, key=lambda i: (-round(float(turns[i]), 8), i))[:36])
    options = []
    for selected in itertools.combinations(choices, 4):
        counts = np.diff((*selected, selected[0] + n))
        if min(counts) < 2 or (counts[0] - counts[2]) % 2 or (counts[1] - counts[3]) % 2:
            continue
        p = boundary[list(selected)]
        chord = np.roll(p, -1, axis=0) - p
        size = np.linalg.norm(chord, axis=1)
        sine = cross(chord, np.roll(chord, -1, axis=0)) / (size * np.roll(size, -1))
        if sine.min() < .15:
            continue
        angles = [math.atan2(float(-cross(-unit[i - 1], unit[i])),
                             float(-unit[i - 1] @ unit[i])) % math.tau for i in selected]
        if min(angles) < .20 or max(angles) > math.pi + .01:
            continue
        balance = np.sum(np.log(size[:2] / size[2:]) ** 2)
        score = np.sum((np.asarray(angles) - math.pi / 2.) ** 2) + 1.5 * np.sum(1. - sine ** 2) + .3 * balance
        options.append((round(float(score), 8), selected))
    return sorted(options)
