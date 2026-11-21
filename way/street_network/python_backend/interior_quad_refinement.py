"""Bounded anisotropic direct-quad refinement with an unchanged physical rim."""
import numpy as np
from .mesh import Assembler
from .quality_correction import QuadShapeCorrection
from .updates import drain_steps


class InteriorQuadRefinementError(ValueError):
    """A bounded alternative atlas did not meet the unchanged quality limits."""
    def __init__(self, message, result, vertex_map, metadata):
        self.result,self.vertex_map,self.metadata=result,vertex_map,metadata
        super().__init__(message)


class InteriorQuadRefinement:
    def __init__(self, patterns, *, max_cells=50000, max_passes=3, maximum_candidates=512):
        if (type(max_cells) is not int or max_cells<1 or type(max_passes) is not int or
                not 1<=max_passes<=3 or type(maximum_candidates) is not int or
                not 1<=maximum_candidates<=512):
            raise ValueError('Invalid bounded interior quad-refinement allocation')
        self.patterns,self.max_cells=patterns,max_cells
        self.max_passes,self.maximum_candidates=max_passes,maximum_candidates

    @staticmethod
    def _edges(mesh):
        edges=np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1),axis=-1)
        unique,counts=np.unique(edges.reshape(-1,2),axis=0,return_counts=True)
        if not set(counts).issubset({1,2}):
            raise ValueError('Interior quad refinement requires manifold edge incidence')
        return edges,{tuple(edge) for edge in unique[counts==1]}

    def _refine(self, mesh, selected):
        edges,boundary=self._edges(mesh)
        points=mesh.xyz[mesh.quads,:2]
        lengths=np.linalg.norm(np.roll(points,-1,axis=1)-points,axis=2)
        split={tuple(edge) for i in selected for edge,length in zip(edges[i],lengths[i])
               if length>7.5*lengths[i].min() and tuple(edge) not in boundary}
        if not split:
            raise ValueError('Interior refinement cannot split an immutable physical rim edge')
        atlases,total=[],0
        for row in edges:
            pattern=self.patterns.quad(*(3 if tuple(edge) in split else 1 for edge in row))
            total+=len(pattern.quads)
            if total>self.max_cells:
                raise ValueError('Interior quad refinement exceeds its finite cell budget')
            atlases.append(pattern)
        assembly=Assembler()
        for fi,(q,pattern) in enumerate(zip(mesh.quads,atlases)):
            u,v=pattern.uv.T;p=mesh.xyz[q,:2]
            xy=(((1.-u)*(1.-v))[:,None]*p[0]+(u*(1.-v))[:,None]*p[1]
                +(u*v)[:,None]*p[2]+((1.-u)*v)[:,None]*p[3])
            keys=[]
            for at,(x,y) in enumerate(pattern.uv):
                if x in (0.,1.) and y in (0.,1.):
                    corner={(0.,0.):0,(1.,0.):1,(1.,1.):2,(0.,1.):3}[x,y]
                    keys.append(('vertex',int(q[corner])));xy[at]=mesh.xyz[q[corner],:2]
                    continue
                side,t=((0,x) if y==0. else (1,y) if x==1. else (2,1.-x) if y==1.
                        else (3,1.-y) if x==0. else (None,None))
                if side is None:
                    keys.append(('face',fi,at))
                else:
                    a,b=map(int,(q[side],q[(side+1)%4]))
                    if a>b:t=1.-t
                    keys.append(('edge',min(a,b),max(a,b),int(round(t*3.))))
            xyz=np.column_stack((xy,np.full(len(xy),mesh.xyz[0,2])))
            assembly.add(xyz,pattern.quads,keys)
        result=assembly.mesh(mesh.material)
        mapping=np.array([assembly.lookup['vertex',i] for i in range(len(mesh.xyz))],dtype=np.int32)
        result.seams={key:mapping[ids] for key,ids in mesh.seams.items()}
        result.validate()
        return result,mapping,len(split)

    def improve(self, mesh):
        return drain_steps(self.improve_steps(mesh))

    def improve_steps(self, mesh):
        mesh.validate()
        if len(mesh.quads)>self.max_cells or np.ptp(mesh.xyz[:,2])>1e-10:
            raise ValueError('Interior quad refinement requires a bounded planar input mesh')
        if len(np.unique(mesh.quads))!=len(mesh.xyz):
            raise ValueError('Interior quad refinement cannot retain unused vertices')
        _,boundary=self._edges(mesh)
        fixed=np.array(sorted({node for edge in boundary for node in edge}),dtype=np.int32)
        original_rim=mesh.xyz[fixed].copy()
        mapping=np.arange(len(mesh.xyz),dtype=np.int32)
        correction=QuadShapeCorrection(max_cells=self.max_cells)
        details=dict(construction='anisotropic direct 1:3 interior quad templates',
                     original_physical_rim_exact=True,physical_boundary_edges_split=0,
                     maximum_passes=self.max_passes,maximum_cells=self.max_cells,passes=[])
        current=mesh
        for attempt in range(self.max_passes+1):
            _,quality,_=correction.evaluate(current.xyz[:,:2],current.quads)
            if correction.qualified(quality):
                details.update(qualified=True,quality=quality,quads=len(current.quads))
                return current,mapping,details
            if attempt==self.max_passes:break
            points=current.xyz[current.quads,:2]
            lengths=np.linalg.norm(np.roll(points,-1,axis=1)-points,axis=2)
            selected=np.flatnonzero(lengths.max(axis=1)/lengths.min(axis=1)>8.)
            if not 1<=len(selected)<=self.maximum_candidates:
                break
            try:
                refined,remap,split=self._refine(current,selected)
            except ValueError as failure:
                details['refinement_failure']=str(failure)
                break
            mapping=remap[mapping];fixed=remap[fixed]
            if not np.array_equal(refined.xyz[fixed],original_rim):
                raise ValueError('Interior refinement moved an original physical boundary')
            current,report=yield from correction.improve_steps(refined)
            if not np.array_equal(current.xyz[fixed],original_rim):
                raise ValueError('Interior correction moved an original physical boundary')
            details['passes'].append(dict(selected_cells=len(selected),split_interior_edges=split,
                quads=len(current.quads),correction=report))
        details.update(qualified=False,quality=quality,quads=len(current.quads))
        raise InteriorQuadRefinementError('Bounded interior quad refinement did not meet the physical quality contract',
                                         current,mapping,details)
