"""Bounded graph flooding to find the minimum fill needed for boundary drainage."""
import heapq
import numpy as np


class ProtectedSpillUnavailable(ValueError):
    """No bounded outlet could be constructed with the prescribed heights."""


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

    def _graph(self, mesh, fixed):
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
        return offsets, target, estimated

    @staticmethod
    def _spill(original, fixed, offsets, target):
        heights = original.copy()
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
        return heights

    def fill(self, mesh, fixed):
        offsets, target, estimated = self._graph(mesh, fixed)
        heights = self._spill(mesh.xyz[:, 2], np.asarray(fixed), offsets, target)
        change = heights - mesh.xyz[:, 2]
        return heights, dict(raised_vertices=int(np.sum(change > 1e-10)),
            maximum_fill=float(change.max()), estimated_bytes=estimated,
            construction='minimum graph spill elevation to fixed boundary outlets')

    def fill_preserving(self, mesh, fixed, prescribed):
        """Open a path for trapped prescribed heights before filling the other basins."""
        offsets, neighbors, _ = self._graph(mesh, fixed)
        heights = mesh.xyz[:, 2].copy()
        protected = np.asarray(prescribed)
        if protected.dtype != np.bool_ or protected.shape != heights.shape or not np.all(protected[fixed]):
            raise ValueError('Prescribed drainage heights must include all boundary outlets')
        filled = self._spill(heights, fixed, offsets, neighbors)
        conflicts = np.flatnonzero(protected & (filled > heights))
        routes = []
        visits = 0
        for start in conflicts[np.argsort(heights[conflicts], kind='stable')]:
            level = heights[start]
            spill = filled[start]
            if spill <= level:
                continue
            # Prefer the lowest pass. If its prescribed heights seal it,
            # a slightly higher, editable vertex can provide a cheaper cut.
            for upper in (spill, np.inf):
                queue = [(0., 0, int(start))]
                costs = {start: (0., 0)}
                parent = {}
                while queue:
                    cost, length, node = heapq.heappop(queue)
                    if (cost, length) != costs[node]:
                        continue
                    visits += 1
                    if visits > 8*len(heights) or len(routes) >= 64:
                        raise ProtectedSpillUnavailable('Prescribed drainage exceeds its local search budget')
                    if filled[node] <= level:
                        break
                    for neighbor in neighbors[offsets[node]:offsets[node+1]]:
                        z = heights[neighbor]
                        if z > upper or (protected[neighbor] and z != level and filled[neighbor] > level):
                            continue
                        score = (cost + max(0., z-level), length+1)
                        if score < costs.get(neighbor, (np.inf, 0)):
                            costs[neighbor] = score
                            parent[neighbor] = node
                            heapq.heappush(queue, (*score, int(neighbor)))
                else:
                    continue
                break
            else:
                raise ProtectedSpillUnavailable('The lowest spill is sealed by prescribed heights')
            end = node
            while node in parent:
                node = parent[node]
                if not protected[node]:
                    heights[node] = level
            filled = self._spill(heights, fixed, offsets, neighbors)
            if filled[start] > level:
                raise ProtectedSpillUnavailable('The proposed route did not drain to a physical outlet')
            routes.append(dict(start=int(start), outlet_path_vertex=end, edges=length,
                               previous_fill=float(spill-level)))
        if not routes or not np.array_equal(filled[protected], mesh.xyz[protected, 2]):
            raise ProtectedSpillUnavailable('Drainage did not retain every prescribed height')
        return filled, dict(routes=routes, visited=visits,
                            maximum_cut=float(np.max(mesh.xyz[:, 2]-filled)))
