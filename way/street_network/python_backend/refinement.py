"""Conforming local quad refinement using cached 1:3 boundary templates."""
import numpy as np
from .mesh import Assembler


class QuadRefiner:
    def __init__(self, patterns, *, max_cells=50000):
        if type(max_cells) is not int or max_cells < 1:
            raise ValueError('Quad refinement needs a positive cell budget')
        self.patterns, self.max_cells = patterns, max_cells

    def refine(self, mesh, selected, field):
        selected = np.asarray(selected)
        if (selected.ndim != 1 or selected.dtype.kind not in 'iu'
                or np.any(selected < 0) or np.any(selected >= len(mesh.quads))):
            raise ValueError('Quad refinement needs valid face indices')
        edges = np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1)
        keys = np.sort(edges, axis=-1)
        unique, counts = np.unique(keys.reshape(-1, 2), axis=0, return_counts=True)
        if not set(counts).issubset({1, 2}):
            raise ValueError('Quad refinement requires a manifold input mesh')
        if len(np.unique(mesh.quads)) != len(mesh.xyz):
            raise ValueError('Quad refinement cannot retain unused input vertices')
        boundary = {tuple(edge) for edge, count in zip(unique, counts) if count == 1}
        split = {tuple(edge) for edge in keys[selected].reshape(-1, 2) if tuple(edge) not in boundary}
        atlases, total = [], 0
        for row in keys:
            pattern = self.patterns.quad(*(3 if tuple(edge) in split else 1 for edge in row))
            total += len(pattern.quads)
            if total > self.max_cells:
                raise ValueError('Adaptive quad refinement exceeds its cell budget')
            atlases.append(pattern)
        assembly = Assembler()
        for fi, (q, pattern) in enumerate(zip(mesh.quads, atlases)):
            u, v = pattern.uv.T
            p = mesh.xyz[q, :2]
            xy = (((1. - u) * (1. - v))[:, None] * p[0] + (u * (1. - v))[:, None] * p[1]
                  + (u * v)[:, None] * p[2] + ((1. - u) * v)[:, None] * p[3])
            addresses = []
            for x, y in pattern.uv:
                if x in (0., 1.) and y in (0., 1.):
                    corner = {(0., 0.): 0, (1., 0.): 1, (1., 1.): 2, (0., 1.): 3}[x, y]
                    addresses.append(('vertex', int(q[corner])))
                    # Avoid roundoff moving original shared/boundary vertices.
                    xy[len(addresses) - 1] = mesh.xyz[q[corner], :2]
                    continue
                side, t = ((0, x) if y == 0. else (1, y) if x == 1. else
                           (2, 1. - x) if y == 1. else (3, 1. - y) if x == 0. else (None, None))
                if side is None:
                    addresses.append(('face', fi, len(addresses)))
                else:
                    a, b = map(int, (q[side], q[(side + 1) % 4]))
                    if a > b:
                        t = 1. - t
                    addresses.append(('edge', min(a, b), max(a, b), int(round(t * 3.))))
            height = field(xy)
            for at, address in enumerate(addresses):
                if address[0] == 'vertex':
                    height[at] = mesh.xyz[address[1], 2]
            assembly.add(np.column_stack((xy, height)), pattern.quads, addresses)
        result = assembly.mesh(mesh.material)
        vertex_map = np.array([assembly.lookup['vertex', i] for i in range(len(mesh.xyz))], dtype=np.int32)
        result.seams = {key: vertex_map[ids] for key, ids in mesh.seams.items()}
        result.validate()
        return result, vertex_map
