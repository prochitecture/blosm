"""Bounded affine height constraints for an optional fixed-chart projection.

Each row constrains |a.z| <= limit for four existing vertices. Dykstra state
is local to one solve; immutable source profiles cannot move to satisfy a row.
"""
import numpy as np
from .slope_projection import face_colors


class HeightSlabs:
    def __init__(self, vertices, coefficients, limits):
        vertices=np.asarray(vertices)
        coefficients=np.asarray(coefficients,dtype=float)
        limits=np.asarray(limits,dtype=float)
        if (vertices.ndim != 2 or vertices.shape[1] != 4 or not len(vertices)
                or vertices.dtype.kind not in 'iu' or np.any(vertices < 0) or np.any(vertices > np.iinfo(np.int32).max)
                or coefficients.shape != vertices.shape or not np.isfinite(coefficients).all()
                or limits.shape != (len(vertices),) or not np.isfinite(limits).all()
                or np.any(limits < 0.) or any(len(set(row)) != 4 for row in vertices)):
            raise ValueError('Height slabs need finite four-vertex rows and nonnegative limits')
        self.vertices=vertices.astype(np.int32,copy=True)
        self.coefficients,self.limits=coefficients.copy(),limits.copy()

    @property
    def estimated_bytes(self):return 512*len(self.vertices)

    def violation(self, heights):
        value=abs(np.einsum('ij,ij->i',self.coefficients,np.asarray(heights)[self.vertices]))
        return float(np.maximum(0.,value-self.limits).max())

    def prepare(self, mesh, fixed, max_colors):
        if (self.vertices.max() >= len(mesh.xyz) or len(self.vertices) > 8*len(mesh.quads)):
            raise ValueError('Height slabs exceed the original chart index or row budget')
        return _PreparedSlabs(self,mesh,fixed,max_colors)


class _PreparedSlabs:
    def __init__(self, source, mesh, fixed, max_colors):
        self.source=source
        self.free=~fixed[source.vertices]
        self.movable=source.coefficients*self.free
        self.norm2=np.sum(self.movable**2,axis=1)
        self.constant=mesh.xyz[source.vertices,2].copy()
        fixed_gradient=np.einsum('ij,ij->i',source.coefficients,self.constant)
        if np.any((self.norm2 <= 1e-24) & (abs(fixed_gradient) > source.limits+1e-12)):
            raise ValueError('Fixed profile heights make a directional slope constraint infeasible')
        self.groups=face_colors(source.vertices,fixed,max_colors)
        self.corrections=np.zeros_like(source.coefficients)

    def sweep(self, heights):
        source=self.source
        for group in self.groups:
            ids=source.vertices[group]
            free=self.free[group]
            trial=heights[ids]+self.corrections[group]
            value=np.einsum('ij,ij->i',source.coefficients[group],trial)
            delta=np.clip(value,-source.limits[group],source.limits[group])-value
            step=np.divide(delta,self.norm2[group],out=np.zeros_like(delta),where=self.norm2[group] > 1e-24)
            projected=trial+step[:,None]*self.movable[group]
            projected[~free]=self.constant[group][~free]
            self.corrections[group]=(trial-projected)*free
            heights[ids]=projected
