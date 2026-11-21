"""Bounded fixed-rim collar atlases for a certified composed scope leaf."""
import numpy as np
from .conforming_boundary_refinement import _graded_face_steps,_parent_parameters
from .four_side_charts import FourSidedQuadChart,quality
from .mesh import Assembler,Mesh
from .quad_charts import ChartResult
from .quality_correction import QuadShapeCorrection
from .updates import drain_steps


def grade_fixed_cells(mesh,patterns,*,max_cells=50000):
    return drain_steps(grade_fixed_cells_steps(mesh,patterns,max_cells=max_cells))


def grade_fixed_cells_steps(mesh,patterns,*,max_cells=50000):
    """Repair only unqualified positive cells, preserving every old vertex."""
    mesh.validate()
    if type(max_cells) is not int or max_cells<len(mesh.quads) or np.ptp(mesh.xyz[:,2])>1e-10:
        raise ValueError('Fixed-cell grading needs a bounded planar mesh')
    correction=QuadShapeCorrection(max_cells=max_cells);assembly=Assembler();total=0;reports=[]
    layers={name:[] for name in mesh.uv_layers}
    for fi,quad in enumerate(mesh.quads):
        if fi and fi%64==0:
            yield
        p=mesh.xyz[quad];faces=np.array(((0,1,2,3),),dtype=np.int32)
        _,double,_=correction.evaluate(p[:,:2],faces)
        _,stored,_=correction.evaluate(p[:,:2].astype(np.float32).astype(float),faces)
        if correction.qualified(double) and correction.qualified(stored):
            assembly.add(p,faces,[('vertex',int(v)) for v in quad]);total+=1
            for name,values in mesh.uv_layers.items():layers[name].append(values[fi])
            continue
        local,details=yield from _graded_face_steps(p[:,:2],patterns,
            max_cells-total-(len(mesh.quads)-fi-1),np.mean(p[:,:2],axis=0))
        local.xyz[:,2]=p[0,2];local.xyz[:4]=p
        assembly.add(local.xyz,local.quads,[('vertex',int(v)) for v in quad]+
            [('local_face',fi,i) for i in range(4,len(local.xyz))])
        total+=len(local.quads);reports.append(dict(face=fi,quads=len(local.quads),**details))
        if total>max_cells:raise ValueError('Fixed-cell collar grading exceeds its finite cell budget')
        if layers:
            uv=_parent_parameters(local.xyz[:,:2],p);u,v=uv.T
            for name,values in mesh.uv_layers.items():
                a=values[fi];mapped=(((1.-u)*(1.-v))[:,None]*a[0]+(u*(1.-v))[:,None]*a[1]
                    +(u*v)[:,None]*a[2]+((1.-u)*v)[:,None]*a[3])
                layers[name].extend(mapped[local.quads])
    output=assembly.mesh(mesh.material)
    mapping=np.asarray([assembly.lookup['vertex',i] for i in range(len(mesh.xyz))],dtype=np.int32)
    output.seams={key:mapping[ids] for key,ids in mesh.seams.items()}
    for name,values in layers.items():output=output.with_uv_layer(name,np.asarray(values))
    if not np.array_equal(output.xyz[mapping],mesh.xyz):raise ValueError('Fixed-cell grading moved an original source vertex')
    report=quality(output)
    _,stored,_=correction.evaluate(output.xyz[:,:2].astype(np.float32).astype(float),output.quads)
    if not correction.qualified(report) or not correction.qualified(stored):
        raise ValueError('Fixed-cell grading failed its unchanged double/stored shape limits')
    return output,mapping,dict(construction='bounded positive-cell collars with all original vertices fixed',
        quality=report,stored_quality=stored,local_collars=reports,original_vertices_exact=True,
        physical_rim_exact=True,all_direct_quads=True)


def fixed_rim_chart(boundary,patterns,*,max_cells=50000,max_candidates=6):
    return drain_steps(fixed_rim_chart_steps(boundary,patterns,max_cells=max_cells,max_candidates=max_candidates))


def fixed_rim_chart_steps(boundary,patterns,*,max_cells=50000,max_candidates=6):
    """Return a positive, qualified direct-quad atlas, with its rim first.

    Try only existing four-sided corner choices. A positive mapped face whose
    immutable edge lengths prevent grading gets a bounded local collar; its
    complete original face rim remains fixed. Actual world float32 coordinates
    must meet the ordinary 8/.03 shape limits too. This is a computational
    atlas for a separately certified physical footprint, never new ownership.
    """
    if (type(max_cells) is not int or max_cells<4 or type(max_candidates) is not int
            or not 1<=max_candidates<=6):raise ValueError('Invalid finite scope atlas allocation')
    boundary=np.asarray(boundary,dtype=float)
    if boundary.ndim!=2 or boundary.shape[1]!=2 or len(boundary)<4 or len(boundary)%2 or not np.isfinite(boundary).all():
        raise ValueError('Finite scope atlas needs a finite even physical rim')
    charts=FourSidedQuadChart(patterns,max_cells=max_cells,max_candidates=max_candidates)
    correction=QuadShapeCorrection(max_cells=max_cells)
    failures=[]
    for _,corners in charts.corners(boundary)[:max_candidates]:
        try:
            seed,counts=charts.map(boundary,corners);initial=quality(seed)
            assembly=Assembler();graded=[];total=0
            for fi,quad in enumerate(seed.quads):
                if fi % 64 == 0:
                    yield
                p=seed.xyz[quad];faces=np.array(((0,1,2,3),),dtype=np.int32)
                _,double,_=correction.evaluate(p[:,:2],faces)
                _,stored,_=correction.evaluate(p[:,:2].astype(np.float32).astype(float),faces)
                if correction.qualified(double) and correction.qualified(stored):
                    total+=1
                    assembly.add(p,faces,[('vertex',int(v)) for v in quad])
                else:
                    local,details=yield from _graded_face_steps(p[:,:2],patterns,max_cells-total,np.mean(p[:,:2],axis=0))
                    total+=len(local.quads)
                    assembly.add(local.xyz,local.quads,[('vertex',int(v)) for v in quad]+
                        [('local_face',fi,i) for i in range(4,len(local.xyz))])
                    graded.append(dict(face=fi,quads=len(local.quads),**details))
                if total>max_cells:raise ValueError('Finite scope collar atlas exceeds its finite cell budget')
            candidate=assembly.mesh()
            source_rim=np.asarray([assembly.lookup['vertex',int(v)] for v in seed.seams[0]],dtype=np.int32)
            interior=np.flatnonzero(~np.isin(np.arange(len(candidate.xyz)),source_rim))
            order=np.r_[source_rim,interior];remap=np.empty(len(order),dtype=np.int32)
            remap[order]=np.arange(len(order),dtype=np.int32)
            result=Mesh(candidate.xyz[order],remap[candidate.quads],seams={0:np.arange(len(boundary),dtype=np.int32)})
            if not np.array_equal(result.xyz[:len(boundary),:2],boundary):
                raise ValueError('Finite scope atlas changed an original physical boundary point')
            report=quality(result)
            _,stored,_=correction.evaluate(result.xyz[:,:2].astype(np.float32).astype(float),result.quads)
            if not charts.accepted(report) or not correction.qualified(stored):
                raise ValueError('Finite scope atlas failed its unchanged double/stored shape limits')
            return ChartResult(result,dict(construction='fixed-rim direct quad atlas with bounded local collars',
                initial=initial,quality=report,stored_quality=stored,local_collars=graded,
                atlas_corners=corners,rim_intervals=counts,boundary_first=True,
                exact_physical_rim=True,all_direct_quads=True))
        except ValueError as failure:failures.append(str(failure))
    raise ValueError('No qualified finite fixed-rim collar atlas: '+'; '.join(failures))
