"""Bounded graph flooding to find the minimum fill needed for boundary drainage."""
import heapq
import numpy as np


class BoundaryDrainage:
    """Raise interior basins to their lowest spill path; retain fixed outlets.

    The heap visits vertices in increasing spill elevation. Each newly visited
    vertex has a parent at an equal or lower filled height. Consequently every
    reachable vertex has a nonascending edge path to a fixed boundary vertex.
    This is a diagnostic/floor construction, not a replacement for slope checks.
    """
    def __init__(self, *, max_cells=50000, memory_bytes=96 * 1024 * 1024):
        if (type(max_cells) is not int or max_cells < 1
                or type(memory_bytes) is not int or memory_bytes < 4096):
            raise ValueError('Invalid drainage budgets')
        self.max_cells, self.memory_bytes = max_cells, memory_bytes

    def fill(self, mesh, fixed):
        fixed = np.asarray(fixed)
        if fixed.dtype != np.bool_ or fixed.shape != (len(mesh.xyz),) or not np.any(fixed):
            raise ValueError('Drainage needs a Boolean mask with fixed boundary outlets')
        estimated = len(mesh.quads) * 256 + len(mesh.xyz) * 256
        if len(mesh.quads) > self.max_cells or estimated > self.memory_bytes:
            raise ValueError('Drainage exceeds its mesh/memory budget')
        mesh.validate()
        faces = mesh.quads
        edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        edges = np.unique(edges, axis=0)
        source, target = np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]]
        order = np.lexsort((target, source))
        target = target[order]
        offsets = np.r_[0, np.cumsum(np.bincount(source, minlength=len(mesh.xyz)))]
        heights = mesh.xyz[:, 2].copy()
        visited = fixed.copy()
        heap = [(float(heights[i]), int(i)) for i in np.flatnonzero(fixed)]
        heapq.heapify(heap)
        while heap:
            height, first = heapq.heappop(heap)
            for last in target[offsets[first]:offsets[first + 1]]:
                if visited[last]:
                    continue
                visited[last] = True
                heights[last] = max(float(heights[last]), height)
                heapq.heappush(heap, (float(heights[last]), int(last)))
        if not np.all(visited):
            raise ValueError('A pavement component has no fixed drainage outlet')
        change = heights - mesh.xyz[:, 2]
        return heights, dict(raised_vertices=int(np.sum(change > 1e-10)),
            maximum_fill=float(change.max()), estimated_bytes=estimated,
            construction='minimum graph spill elevation to fixed boundary outlets')
