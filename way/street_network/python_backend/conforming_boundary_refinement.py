"""Insert prescribed separator-rim points with conforming direct quads.

All original vertices and outer/profile chords are immutable. An odd number
of added vertices on a face propagates through original *internal* dual edges,
so every refined face has an even boundary. Cached positive templates cover
both opposite and adjacent even side counts. No triangles are constructed.
"""
from collections import deque
import numpy as np
from .compounds import rim_label
from .mesh import Assembler,Mesh
from .partition import BoundaryRef
from .patterns import Grid
from .regions import signed_area
from .updates import drain_steps


def _template(patterns,counts):
    counts=tuple(map(int,counts))
    if min(counts)<1 or sum(counts)%2:
        raise ValueError('Conforming quad template needs a positive even rim')
    if sum(counts)>8192:raise ValueError('Conforming template exceeds its finite rim budget')
    if (counts[0]-counts[2])%2==0 and (counts[1]-counts[3])%2==0:
        return patterns.quad(*counts)
    key=('conforming_boundary_quad',*counts)
    saved=patterns._lookup(key)
    if saved is not None:return saved
    # Exactly two adjacent even sides remain. Their midpoints and opposite
    # corner join one interior degree-three vertex in three positive quads.
    even={i for i,n in enumerate(counts) if n%2==0}
    rotation=next(i for i in range(4) if {i,(i+1)%4}==even)
    uv=np.array(((0.,0.),(.5,0.),(1.,0.),(1.,.5),(1.,1.),(0.,1.),(.5,.5)))
    for _ in range(rotation):uv=np.column_stack((1.-uv[:,1],uv[:,0]))
    points=list(uv)
    faces=[(6,1,2,3),(6,3,4,5),(6,5,0,1)]

    def side_of(a,b):
        if a[1]==b[1]==0.:return 0
        if a[0]==b[0]==1.:return 1
        if a[1]==b[1]==1.:return 2
        if a[0]==b[0]==0.:return 3
        return None

    for side,target in enumerate(counts):
        current=2 if side in even else 1
        for _ in range((target-current)//2):
            candidates=[]
            for fi,face in enumerate(faces):
                for at,(a,b) in enumerate(zip(face,(*face[1:],face[0]))):
                    if side_of(points[a],points[b])==side:
                        candidates.append((float(np.linalg.norm(points[a]-points[b])),fi,at))
            _,owner,at=max(candidates,key=lambda row:row[0])
            face=faces[owner];a0,a3,b1,b0=face[at:]+face[:at]
            a1,a2=len(points),len(points)+1
            points.extend(((2.*points[a0]+points[a3])/3.,(points[a0]+2.*points[a3])/3.))
            p,q=len(points),len(points)+1
            points.extend((.5*(points[a1]+(2.*points[b0]+points[b1])/3.),
                           .5*(points[a2]+(points[b0]+2.*points[b1])/3.)))
            added=((a0,a1,p,b0),(a1,a2,q,p),(a2,a3,b1,q),(p,q,b1,b0))
            faces[owner]=added[0];faces.extend(added[1:])
    result=Grid(np.asarray(points),np.asarray(faces,dtype=np.int32),max(counts[::2]),max(counts[1::2]))
    patterns._positive(result)
    return patterns._remember(key,result)


def _edge(a,b):return (min(int(a),int(b)),max(int(a),int(b)))


def _hexagon_chart(boundary, max_cells):
    """Tile a six-vertex rim with one diagonal or three quads around a center."""
    from .quality_correction import QuadShapeCorrection
    quality = QuadShapeCorrection()
    indices = np.arange(6)
    for rotation in range(3):
        ring = np.roll(indices, -rotation)
        layouts = [(boundary, (ring[:4], ring[[3, 4, 5, 0]]))]
        if rotation < 2:
            # Joining alternating rim vertices to their centroid splits their
            # triangle into three pieces. Each absorbs one outer corner.
            center = boundary[ring[::2]].mean(axis=0)
            layouts.append((np.vstack((boundary, center)),
                [(*np.roll(ring, -offset)[:3], 6) for offset in (0, 2, 4)]))
        for xy, cells in layouts:
            if len(cells) > max_cells:
                continue
            quads = np.asarray(cells, dtype=np.int32)
            _, physical, _ = quality.evaluate(xy, quads)
            if not quality.qualified(physical):
                continue
            _, stored, _ = quality.evaluate(xy.astype(np.float32).astype(float), quads)
            if quality.qualified(stored):
                mesh = Mesh(np.column_stack((xy, np.zeros(len(xy)))), quads, seams={0:indices})
                mesh.validate()
                return mesh, dict(construction='direct six-vertex rim', quality=physical, stored_quality=stored)
    return None


def _graded_face(boundary,patterns,max_cells,center):
    return drain_steps(_graded_face_steps(boundary,patterns,max_cells,center))


def _graded_face_steps(boundary,patterns,max_cells,center):
    """Separate immutable unequal corner edges before local shape repair.

    A quad containing two fixed boundary edges with ratio above eight cannot
    be repaired by moving its other vertices. Positive radial collars give
    each prescribed edge its own face. Their graded inner core then uses the
    ordinary direct-quad charts, rotations and anisotropic refinement.
    """
    from .quad_charts import LayeredQuadChart,FixedRimQuality
    from .quad_rotations import QuadEdgeRotations
    from .quality_correction import QuadShapeCorrection
    from .interior_quad_refinement import InteriorQuadRefinement,InteriorQuadRefinementError
    if len(boundary) == 6:
        direct = _hexagon_chart(np.asarray(boundary), max_cells)
        if direct is not None:
            return direct
    original=np.asarray(boundary).copy();origin=original[0].copy()
    boundary=original-origin;center=np.asarray(center)-origin
    limit=min(2000,max_cells)
    if limit<4:raise ValueError('Prescribed local grading exceeds its finite cell budget')
    lengths=np.linalg.norm(np.roll(boundary,-1,axis=0)-boundary,axis=1)
    count=len(boundary);best=None
    for factor in (1.,1.5,2.,.75,.5):
        for rings in (1,2,3):
            yield
            p=boundary.copy();points=list(p);faces=[]
            try:
                for layer in range(rings):
                    radii=np.linalg.norm(p-center,axis=1)
                    local=np.linalg.norm(np.roll(p,-1,axis=0)-p,axis=1)
                    step=np.minimum(.3*radii,factor*np.sqrt(local*np.roll(local,1)))
                    inner=p+(center-p)*(step/radii)[:,None]
                    base,new=layer*count,(layer+1)*count
                    faces.extend((base+i,base+(i+1)%count,new+(i+1)%count,new+i) for i in range(count))
                    points.extend(inner);p=inner
                core=LayeredQuadChart(patterns,spacing=float(np.median(lengths)),max_cells=limit,
                    quality_iterations=10).build(p,require_quality=False).mesh
                remap=np.r_[np.arange(rings*count,(rings+1)*count),
                    np.arange(len(points),len(points)+len(core.xyz)-count)]
                points.extend(core.xyz[count:,:2]);faces.extend(remap[core.quads])
                if len(faces)>limit:continue
                candidate=Mesh(np.column_stack((points,np.zeros(len(points)))),np.asarray(faces),
                    seams={0:np.arange(count,dtype=np.int32)})
                candidate.validate()
                correction=QuadShapeCorrection(max_cells=limit,max_iterations=500)
                candidate,details=yield from correction.improve_steps(candidate)
                if not details['qualified']:
                    candidate,_=yield from QuadEdgeRotations(max_cells=limit).improve_steps(candidate,FixedRimQuality(iterations=100))
                    candidate,details=yield from QuadShapeCorrection(max_cells=limit,max_iterations=1000).improve_steps(candidate)
                if not details['qualified']:
                    try:candidate,_,_=yield from InteriorQuadRefinement(patterns,max_cells=limit).improve_steps(candidate)
                    except InteriorQuadRefinementError as failure:candidate=failure.result
                _,report,_=correction.evaluate(candidate.xyz[:,:2],candidate.quads)
                if correction.qualified(report):
                    if not np.array_equal(candidate.xyz[:count,:2],boundary):
                        raise ValueError('Local prescribed grading changed its original face rim')
                    world=candidate.xyz.copy();world[:,:2]+=origin;world[:count,:2]=original
                    _,stored,_=correction.evaluate(world[:,:2].astype(np.float32).astype(float),candidate.quads)
                    if not correction.qualified(stored):
                        # Reserve numerical room only when the actual stored
                        # coordinates need it. The final limits remain 8/.03.
                        candidate,_=yield from QuadShapeCorrection(max_cells=limit,max_iterations=1500,
                            maximum_ratio=7.6,target_ratio=7.4).improve_steps(candidate)
                        world=candidate.xyz.copy();world[:,:2]+=origin;world[:count,:2]=original
                        _,stored,_=correction.evaluate(world[:,:2].astype(np.float32).astype(float),candidate.quads)
                    if correction.qualified(stored):
                        _,report,_=correction.evaluate(candidate.xyz[:,:2],candidate.quads)
                        candidate.xyz=world
                        return candidate,dict(radial_collars=rings,collar_factor=factor,quality=report,
                            stored_quality=stored)
                if best is None or report['maximum_edge_ratio']<best['maximum_edge_ratio']:best=report
            except ValueError:continue
    raise ValueError('Prescribed face needs another finite grading atlas: '+str(best))


def _parent_parameters(points,quad):
    """Bilinear source-face coordinates for retained heights and UV seams."""
    p=np.asarray(quad)[:,:2];query=np.asarray(points)
    uv=np.full((len(query),2),.5)
    for _ in range(16):
        u,v=uv.T
        current=((1.-u)*(1.-v))[:,None]*p[0]+(u*(1.-v))[:,None]*p[1]+(u*v)[:,None]*p[2]+((1.-u)*v)[:,None]*p[3]
        a=(1.-v)[:,None]*(p[1]-p[0])+v[:,None]*(p[2]-p[3])
        b=(1.-u)[:,None]*(p[3]-p[0])+u[:,None]*(p[2]-p[1])
        delta=query-current;det=a[:,0]*b[:,1]-a[:,1]*b[:,0]
        uv[:,0]+=(delta[:,0]*b[:,1]-delta[:,1]*b[:,0])/det
        uv[:,1]+=(a[:,0]*delta[:,1]-a[:,1]*delta[:,0])/det
        if float(np.max(abs(delta)))<1e-10:break
    if np.any(uv< -1e-8) or np.any(uv>1.+1e-8):
        raise ValueError('Prescribed local grading left its original finite source face')
    return np.clip(uv,0.,1.)


def _requests(mesh,refs,group,requested,tolerance):
    sources={p.region.owner:p for p in group.footprints}
    insertions={};new_points={}
    for li,raw in sorted(requested.items()):
        if type(li) is not int or not 0<=li<len(group.loops) or not group.loops[li].hole or li not in mesh.seams:
            raise ValueError('Prescribed refinement needs an existing hole-boundary loop')
        points=np.asarray(raw,dtype=float)
        if points.ndim!=2 or points.shape[1]!=2 or len(points)<4 or not np.isfinite(points).all() or signed_area(points)<=0.:
            raise ValueError('Prescribed separator boundary must be finite and counterclockwise')
        if np.min(np.linalg.norm(np.roll(points,-1,axis=0)-points,axis=1))<=1e-8:
            raise ValueError('Prescribed separator rim contains a collapsed chord')
        ids=np.asarray(mesh.seams[li],dtype=np.int32)
        old=mesh.xyz[ids,:2]
        if signed_area(old)>=0.:raise ValueError('Pavement hole seam must be clockwise')
        # Start both clockwise cycles at exactly the same old point. The
        # requested cycle must retain every original point, in source order.
        matches=np.flatnonzero(np.all(points==old[0],axis=1))
        if len(matches)!=1:raise ValueError('Prescribed rim lost an original boundary vertex')
        clockwise=points[::-1]
        start=int(np.flatnonzero(np.all(clockwise==old[0],axis=1))[0])
        clockwise=np.roll(clockwise,-start,axis=0)
        cursor=0
        for i,(a,b) in enumerate(zip(ids,np.roll(ids,-1))):
            a,b=int(a),int(b);first,last=mesh.xyz[a,:2],mesh.xyz[b,:2]
            if not np.array_equal(clockwise[cursor%len(clockwise)],first):
                raise ValueError('Prescribed rim changed original vertex order or coordinates')
            cursor+=1;fractions=[];edge=last-first;length2=float(edge@edge)
            while cursor<len(clockwise) and not np.array_equal(clockwise[cursor],last):
                point=clockwise[cursor];fraction=float((point-first)@edge/length2)
                if (not 0.<fraction<1. or fractions and fraction<=fractions[-1]
                        or np.linalg.norm(point-first-fraction*edge)>tolerance):
                    raise ValueError('Prescribed separator point leaves its original boundary chord')
                fractions.append(fraction)
                canonical=fraction if a<b else 1.-fraction
                new_points[_edge(a,b),canonical]=point.copy()
                cursor+=1
            if i<len(ids)-1 and cursor>=len(clockwise):
                raise ValueError('Prescribed rim lost an original boundary vertex')
            if fractions:
                first_ref=refs[a]
                if first_ref.loop!=li or rim_label(group.loops[li].spans[first_ref.edge],sources).kind!='curb':
                    raise ValueError('Prescribed refinement cannot split an incoming profile')
                insertions[_edge(a,b)]=np.asarray(sorted(f if a<b else 1.-f for f in fractions))
        if cursor!=len(clockwise):raise ValueError('Prescribed rim contains unmatched source vertices')
    return insertions,new_points


def _propagate(mesh,incidence,insertions):
    odd=np.zeros(len(mesh.quads),dtype=bool)
    adjacency=[[] for _ in mesh.quads]
    for edge,faces in incidence.items():
        if edge in insertions and len(insertions[edge])%2:
            for face in faces:odd[face]^=True
        if len(faces)==2:
            a,b=faces;adjacency[a].append((b,edge));adjacency[b].append((a,edge))
    parents={};order=[];roots=[]
    for first in range(len(mesh.quads)):
        if first in parents:continue
        roots.append(first);parents[first]=(None,None);queue=deque((first,))
        while queue:
            face=queue.popleft();order.append(face)
            for other,edge in sorted(adjacency[face]):
                if other not in parents:parents[other]=(face,edge);queue.append(other)
    added=0
    for face in reversed(order):
        parent,edge=parents[face]
        if parent is None:
            if odd[face]:raise ValueError('Prescribed hole refinement has incompatible component boundary parity')
        elif odd[face]:
            if edge in insertions:raise ValueError('Parity propagation repeated an interior split')
            insertions[edge]=np.array((.5,));odd[parent]^=True;added+=1
    return added,len(roots)


def refine(mesh,refs,group,patterns,requested_CCW_hole_boundaries,*,max_cells=50000,max_insertions=20000,tolerance=1e-7):
    """Return (mesh, boundary_refs, old_vertex_map, construction_metadata).

    Requested CCW hole rims include every old CW pavement rim point byte-exact.
    Extra points lie on those same finite curb chords. The returned mesh is
    positive and conforming; final physical ratio/slope checks remain mandatory.
    """
    if (type(max_cells) is not int or max_cells<1 or type(max_insertions) is not int or max_insertions<1
            or not np.isfinite(tolerance) or tolerance<=0. or not isinstance(requested_CCW_hole_boundaries,dict)):
        raise ValueError('Invalid bounded prescribed-boundary refinement allocation')
    mesh.validate()
    if any(type(key) is not int for key in requested_CCW_hole_boundaries):
        raise ValueError('Prescribed separator loops need integer source indices')
    if sum(len(points) for points in requested_CCW_hole_boundaries.values())>max_insertions+sum(len(ids) for ids in mesh.seams.values()):
        raise ValueError('Prescribed refinement exceeds its finite edge-insertion budget')
    if len(mesh.quads)>max_cells or len(np.unique(mesh.quads))!=len(mesh.xyz):
        raise ValueError('Prescribed refinement needs a bounded mesh with no unused vertices')
    incidence={}
    for face,quad in enumerate(mesh.quads):
        for a,b in zip(quad,np.roll(quad,-1)):incidence.setdefault(_edge(a,b),[]).append(face)
    if any(len(faces) not in (1,2) for faces in incidence.values()):
        raise ValueError('Prescribed refinement requires manifold edge incidence')
    declared={_edge(a,b) for ids in mesh.seams.values() for a,b in zip(ids,np.roll(ids,-1))}
    if declared!={edge for edge,faces in incidence.items() if len(faces)==1} or set(refs)!={int(v) for ids in mesh.seams.values() for v in ids}:
        raise ValueError('Prescribed refinement requires complete authoritative boundary seams and references')
    insertions,new_points=_requests(mesh,refs,group,requested_CCW_hole_boundaries,tolerance)
    requested_count=sum(map(len,insertions.values()))
    propagated,components=_propagate(mesh,incidence,insertions)
    if sum(map(len,insertions.values()))>max_insertions:
        raise ValueError('Prescribed refinement exceeds its finite edge-insertion budget')
    atlases=[];total=0;refined_faces=[]
    for fi,quad in enumerate(mesh.quads):
        counts=tuple(1+len(insertions.get(_edge(a,b),())) for a,b in zip(quad,np.roll(quad,-1)))
        if any(count>1 for count in counts):refined_faces.append(fi)
        pattern=_template(patterns,counts);atlases.append(pattern);total+=len(pattern.quads)
        if total>max_cells:raise ValueError('Prescribed refinement exceeds its finite cell budget')
    assembly=Assembler();uv_layers={name:[] for name in mesh.uv_layers};local_grading=[];refined_set=set(refined_faces)
    for fi,(quad,pattern) in enumerate(zip(mesh.quads,atlases)):
        u,v=pattern.uv.T;p=mesh.xyz[quad]
        curves=[];fractions=[];side_fractions=[]
        for side,(a,b) in enumerate(zip(quad,np.roll(quad,-1))):
            a,b=int(a),int(b);canonical=np.r_[0.,insertions.get(_edge(a,b),()),1.]
            local=canonical if a<b else 1.-canonical[::-1]
            ids=patterns.boundary(pattern,side)
            raw=pattern.uv[ids,0 if side%2==0 else 1]
            target=local if side in (0,1) else 1.-local[::-1]
            f=np.interp(u if side%2==0 else v,raw,target)
            if side in (2,3):f=1.-f
            curves.append((1.-f)[:,None]*p[side]+f[:,None]*p[(side+1)%4])
            side_fractions.append(f)
            fractions.append(local)
        bilinear=(((1.-u)*(1.-v))[:,None]*p[0]+(u*(1.-v))[:,None]*p[1]
            +(u*v)[:,None]*p[2]+((1.-u)*v)[:,None]*p[3])
        xyz=(1.-v)[:,None]*curves[0]+u[:,None]*curves[1]+v[:,None]*curves[2]+(1.-u)[:,None]*curves[3]-bilinear
        addresses=[]
        for at,(x,y) in enumerate(pattern.uv):
            if x in (0.,1.) and y in (0.,1.):
                corner={(0.,0.):0,(1.,0.):1,(1.,1.):2,(0.,1.):3}[x,y]
                node=int(quad[corner]);addresses.append(('vertex',node));xyz[at]=mesh.xyz[node]
            else:
                side,t=((0,x) if y==0. else (1,y) if x==1. else (2,1.-x) if y==1. else (3,1.-y) if x==0. else (None,None))
                if side is None:addresses.append(('face',fi,at));continue
                a,b=int(quad[side]),int(quad[(side+1)%4]);edge=_edge(a,b)
                ids=patterns.boundary(pattern,side)
                raw=pattern.uv[ids,0 if side%2==0 else 1]
                raw=raw if side in (0,1) else 1.-raw[::-1]
                index=int(np.argmin(abs(raw-t)));fraction=float(fractions[side][index])
                canonical=fraction if a<b else 1.-fraction
                position=int(np.argmin(abs(insertions[edge]-canonical)))
                canonical=float(insertions[edge][position])
                addresses.append(('edge',*edge,position))
                xyz[at]=(1.-canonical)*mesh.xyz[edge[0]]+canonical*mesh.xyz[edge[1]]
                if (edge,canonical) in new_points:xyz[at,:2]=new_points[edge,canonical]
        faces=pattern.quads;parent=None
        if fi in refined_set:
            from .quality_correction import QuadShapeCorrection
            _,shape,_=QuadShapeCorrection().evaluate(xyz[:,:2],faces)
            _,stored_shape,_=QuadShapeCorrection().evaluate(xyz[:,:2].astype(np.float32).astype(float),faces)
            if not QuadShapeCorrection().qualified(shape) or not QuadShapeCorrection().qualified(stored_shape):
                ids=[patterns.boundary(pattern,s) for s in range(4)]
                trace=np.r_[ids[0][:-1],ids[1][:-1],ids[2][::-1][:-1],ids[3][::-1][:-1]]
                rim=xyz[trace].copy()
                try:graded,details=_graded_face(rim[:,:2],patterns,max_cells-total+len(faces),np.mean(p[:,:2],axis=0))
                except ValueError as failure:raise ValueError(f'Prescribed source face {fi}: {failure}') from failure
                total+=len(graded.quads)-len(faces)
                if total>max_cells:raise ValueError('Prescribed refinement exceeds its finite cell budget')
                parent=_parent_parameters(graded.xyz[:,:2],p)
                gu,gv=parent.T
                z=((1.-gu)*(1.-gv))*p[0,2]+(gu*(1.-gv))*p[1,2]+(gu*gv)*p[2,2]+((1.-gu)*gv)*p[3,2]
                xyz=np.column_stack((graded.xyz[:,:2],z));xyz[:len(trace)]=rim
                addresses=[addresses[int(i)] for i in trace]+[('graded_face',fi,i) for i in range(len(trace),len(xyz))]
                faces=graded.quads;local_grading.append(dict(face=fi,quads=len(faces),**details))
        assembly.add(xyz,faces,addresses)
        for name,values in mesh.uv_layers.items():
            if parent is not None:
                gu,gv=parent.T;a=values[fi]
                uv=(((1.-gu)*(1.-gv))[:,None]*a[0]+(gu*(1.-gv))[:,None]*a[1]
                    +(gu*gv)[:,None]*a[2]+((1.-gu)*gv)[:,None]*a[3])
                uv_layers[name].extend(uv[faces]);continue
            a=values[fi];base=(((1.-u)*(1.-v))[:,None]*a[0]+(u*(1.-v))[:,None]*a[1]
                +(u*v)[:,None]*a[2]+((1.-u)*v)[:,None]*a[3])
            curves_uv=[(1.-f)[:,None]*a[s]+f[:,None]*a[(s+1)%4] for s,f in enumerate(side_fractions)]
            uv=(1.-v)[:,None]*curves_uv[0]+u[:,None]*curves_uv[1]+v[:,None]*curves_uv[2]+(1.-u)[:,None]*curves_uv[3]-base
            uv_layers[name].extend(uv[pattern.quads])
    output=assembly.mesh(mesh.material)
    mapping=np.asarray([assembly.lookup['vertex',i] for i in range(len(mesh.xyz))],dtype=np.int32)
    updated={int(mapping[i]):ref for i,ref in refs.items()};seams={}
    for li,ids in mesh.seams.items():
        values=[]
        for a,b in zip(ids,np.roll(ids,-1)):
            a,b=int(a),int(b);values.append(int(mapping[a]));edge=_edge(a,b)
            fractions=insertions.get(edge,())
            indexes=range(len(fractions)) if a<b else reversed(range(len(fractions)))
            for index in indexes:
                canonical=float(fractions[index]);fraction=canonical if a<b else 1.-canonical
                node=assembly.lookup[('edge',*edge,index)];values.append(node)
                first,last=refs[a],refs[b]
                if first.loop!=last.loop:raise ValueError('Boundary refinement crosses physical loops')
                end=last.fraction if first.edge==last.edge else 1.
                updated[node]=BoundaryRef(first.loop,first.edge,first.fraction+fraction*(end-first.fraction))
        seams[li]=np.asarray(values,dtype=np.int32)
    output.seams=seams
    for name,values in uv_layers.items():output=output.with_uv_layer(name,np.asarray(values))
    quality=output.validate()
    from .quality_correction import QuadShapeCorrection
    _,stored_quality,_=QuadShapeCorrection().evaluate(output.xyz[:,:2].astype(np.float32).astype(float),output.quads)
    if not np.array_equal(output.xyz[mapping],mesh.xyz):raise ValueError('Prescribed refinement moved an original source vertex')
    affected={int(v) for fi in refined_faces for v in mesh.quads[fi]}
    local={int(v) for quad in mesh.quads if affected.intersection(map(int,quad)) for v in quad}
    return output,updated,mapping,dict(construction='prescribed finite-hole chords with conforming direct quad templates',
        requested_boundary_vertices=requested_count,propagated_interior_edges=propagated,
        dual_components=components,quads=len(output.quads),original_vertices_exact=True,
        outer_profile_chords_unchanged=True,all_direct_quads=True,quality=quality,stored_quality=stored_quality,
        refined_original_faces=refined_faces,relaxation_old_vertex_ids=sorted(local),local_grading=local_grading)
