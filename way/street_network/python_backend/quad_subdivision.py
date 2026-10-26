"""Uniform direct-quad refinement of an already qualified physical patch.

Refinement does not reselect its footprint, atlas corners or cap decomposition.
Canonical edge addresses keep neighbouring faces conforming. The caller still
qualifies world shape, profiles and slopes after sampling its physical field.
"""
import numpy as np
from .mesh import Assembler


class UniformQuadSubdivision:
    def __init__(self, patterns, *, max_cells=100000):
        if type(max_cells) is not int or max_cells < 1:
            raise ValueError('Quad subdivision needs a positive cell budget')
        self.patterns, self.max_cells = patterns, max_cells

    def build(self, mesh, factor, field, *, datum, ceiling):
        mesh.validate()
        if type(factor) is not int or not 1 <= factor <= 4:
            raise ValueError('Quad subdivision factor must be between one and four')
        if len(mesh.quads)*factor*factor > self.max_cells:
            raise ValueError('Complete quad subdivision exceeds its cell budget')
        if factor == 1:
            return mesh
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if counts.max()>2 or len(np.unique(mesh.quads)) != len(mesh.xyz):
            raise ValueError('Quad subdivision needs a conforming mesh without unused vertices')
        boundary = {tuple(edge) for edge, count in zip(unique, counts) if count==1}
        pattern = self.patterns.grid(factor, factor)
        u, v = pattern.uv.T
        weights = np.c_[(1.-u)*(1.-v), u*(1.-v), u*v, (1.-u)*v]
        assembly = Assembler()
        for face, q in enumerate(mesh.quads):
            xy = weights@mesh.xyz[q, :2]
            addresses, constraints = [], {}
            for i, (x, y) in enumerate(pattern.uv):
                if x in (0., 1.) and y in (0., 1.):
                    corner = {(0., 0.):0, (1., 0.):1, (1., 1.):2, (0., 1.):3}[x, y]
                    vertex = int(q[corner])
                    xy[i] = mesh.xyz[vertex, :2]
                    addresses.append(('vertex', vertex))
                    constraints[i] = mesh.xyz[vertex, 2]
                    continue
                side, t = ((0, x) if y==0. else (1, y) if x==1. else
                    (2, 1.-x) if y==1. else (3, 1.-y) if x==0. else (None, None))
                if side is None:
                    addresses.append(('face', face, i))
                    continue
                a, b = map(int, (q[side], q[(side+1)%4]))
                if a>b:
                    a, b, t = b, a, 1.-t
                step = int(round(t*factor))
                t = step/factor
                # Canonical arithmetic, independent of face traversal direction.
                xy[i] = (1.-t)*mesh.xyz[a, :2]+t*mesh.xyz[b, :2]
                addresses.append(('edge', a, b, step))
                za, zb = mesh.xyz[[a, b], 2]
                if (a, b) in boundary and za==datum and zb==datum:
                    constraints[i] = datum
                elif za==ceiling and zb==ceiling:
                    constraints[i] = ceiling
            height = np.asarray(field(xy), dtype=float)
            if height.shape != (len(xy),):
                raise ValueError('A refinement height field must address every generated vertex')
            for i, z in constraints.items():
                height[i] = z
            assembly.add(np.c_[xy, height], pattern.quads, addresses)
        result = assembly.mesh(mesh.material)
        result.seams = {'crown_guides':np.flatnonzero(result.xyz[:, 2] == ceiling)}
        result.validate()
        return result
