"""Transport a finite source entrance to an independent certified run chart."""
import math
import numpy as np
from .mesh import Mesh
from .entrance_slopes import EntranceSlopeGuard
from .slope_projection import corner_matrices
from .drainage import BoundaryDrainage
from .four_side_charts import quality


def _blend(t):return t**3*(10.+t*(-15.+6.*t))


def _ranges(field,points):
    values=field.frame.local(points)[...,0]
    return all(np.min(values)>=ray.stations[0]-1e-8 and np.max(values)<=ray.stations[-1]+1e-8
               for ray in field.curves)


def _cut(field,station,sign,fractions,sides):
    edges=field.edges(station);ao,ai,bi,bo=map(float,edges)
    if ao<=bi or ai<=bo or max(ai,bo)>=min(ao,bi)-1e-7:
        raise ValueError('Chart transport needs an overlapping, consistently ordered source pair')
    high,low=max(ao,bi),min(ai,bo)
    first,last=map(float,field.peaks(station))
    if abs(first-last)>1e-9:raise ValueError('Chart transport cannot merge a real separator gap')
    peak=-sign*(2.*first-high-low)/(high-low)
    q=peak+fractions*np.where(sides,1.+peak,1.-peak)
    y=.5*(high+low)-sign*.5*(high-low)*q
    return field.frame.world(np.c_[np.full(len(q),station),y])


def _trace(ray,first,last,frame,sign):
    a,b=map(float,frame.local(np.array((first,last)))[:,0])
    if sign*(b-a)<=1e-7:raise ValueError('Chart transport does not advance on a real original curb')
    if max(np.linalg.norm(ray(a)-first),np.linalg.norm(ray(b)-last))>1e-7:
        raise ValueError('Run and finite entrance fields do not share the same original curb chords')
    low,high=sorted((a,b));middle=ray.points[(ray.stations>low+1e-9)&(ray.stations<high-1e-9)]
    if sign<0.:middle=middle[::-1]
    points=np.vstack((first,middle,last));lengths=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
    if lengths[-1]<=1e-8:raise ValueError('Chart transport has no positive original curb lead')
    return points,lengths


def _at(points,lengths,t):
    s=np.asarray(t)*lengths[-1]
    return np.stack([np.interp(s,lengths,points[:,i]) for i in (0,1)],axis=-1)


def _rows(traced,count,exact,tolerance):
    """Coalesce incidental dual-curb knots within the source curve budget."""
    regular=np.linspace(0.,1.,count+1)
    values=np.unique(np.concatenate((regular,*(lengths/lengths[-1] for points,lengths in traced))))
    coordinates=np.column_stack([_at(points,lengths,values) for points,lengths in traced])
    protected=set(int(np.searchsorted(values,t)) for t in regular)
    for (points,lengths),is_exact in zip(traced,exact):
        if not is_exact or len(points)<3:continue
        edge=np.diff(points,axis=0);edge/=np.linalg.norm(edge,axis=1)[:,None]
        for t in (lengths[1:-1]/lengths[-1])[np.linalg.norm(np.diff(edge,axis=0),axis=1)>1e-9]:
            protected.add(int(np.searchsorted(values,t)))
    keep=set(protected);ordered=sorted(protected);pending=list(zip(ordered,ordered[1:]))
    while pending:
        a,b=pending.pop()
        if b<=a+1:continue
        fraction=(values[a+1:b]-values[a])/(values[b]-values[a])
        error=np.linalg.norm(coordinates[a+1:b]-(coordinates[a]+fraction[:,None]*(coordinates[b]-coordinates[a])),axis=1)
        i=int(np.argmax(error))
        if error[i]>tolerance:
            index=a+1+i;keep.add(index);pending.extend(((a,index),(index,b)))
    result=values[sorted(keep)];maximum=0.
    for a,b in zip(sorted(keep),sorted(keep)[1:]):
        fraction=(values[a:b+1]-values[a])/(values[b]-values[a])
        maximum=max(maximum,float(np.linalg.norm(coordinates[a:b+1]
            -(coordinates[a]+fraction[:,None]*(coordinates[b]-coordinates[a])),axis=1).max()))
    return result,maximum


def _contains(mesh,approach,projection):
    """Every complete quad remains in the actual finite original strip union."""
    from .regions import Region,RegionUnion,signed_area,_classify
    from .model import Key
    from .finite_containment import contains_segment,segment_on_polygon_side
    from .source_strip_containment import _distance
    operands=[]
    for role in (0,1):
        first,last=approach.field.curves[2*role:2*role+2]
        xy=np.vstack((first.points,last.points[::-1]))
        operands.append(Region(Key('bridge_source_strip',role),xy if signed_area(xy)>0. else xy[::-1]))
    loops=RegionUnion(max_segments=20000).build(operands)
    outer=[l.xy for l in loops if not l.hole];holes=[l.xy for l in loops if l.hole]
    if len(outer)!=1:raise ValueError('Chart transport needs one connected finite original strip union')
    band=projection+1e-7
    for face in mesh.quads:
        polygon=mesh.xyz[face,:2]
        for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
            if not contains_segment(a,b,outer[0],tolerance=band):
                raise ValueError('Chart transport leaves its finite original source strips')
            for hole in holes:
                if not segment_on_polygon_side(a,b,hole,inside=False,tolerance=band):
                    raise ValueError('Chart transport fills an original finite source gap')
        sides=np.stack((polygon,np.roll(polygon,-1,axis=0)),axis=1)
        for boundary in (outer[0],*holes):
            selected=np.flatnonzero(np.all(boundary>=polygon.min(axis=0)-band,axis=1)
                &np.all(boundary<=polygon.max(axis=0)+band,axis=1))
            if not len(selected):continue
            inside,near=_classify(boundary[selected],polygon,1e-9)
            if np.any((inside | (near>=0))&(_distance(boundary[selected],sides)>band)):
                raise ValueError('Chart transport hides a source boundary inside a quad')
    return dict(proof='continuous_complete_quad_edges_inside_finite_original_strip_union',
        source_chord_error=approach.source_curve_tolerance,derived_polyline_projection_error=projection,
        boundary_band=band,maximum_combined_error=approach.source_curve_tolerance+band,
        holes=len(holes),complete_quads=len(mesh.quads))


def _same_source_curbs(approach,target,first,last):
    """Continuous equality of all four prescribed source-chord chains."""
    frame=approach.frame_chart;low,high=sorted((float(first),float(last)));maximum=0.;checks=0
    for original,actual in zip(approach.field.curves,target.curves):
        values=frame.local(actual.points)[:,0];covered=low
        for p,q,a,b in zip(actual.points,actual.points[1:],values,values[1:]):
            if max(a,b)<low-1e-9 or min(a,b)>high+1e-9:continue
            if b<=a+1e-12:raise ValueError('Independent run field reverses an original finite source curb')
            start,stop=max(low,a),min(high,b)
            if stop<start:continue
            if start>covered+1e-8:raise ValueError('Independent run field lost a finite original curb interval')
            covered=max(covered,stop)
            nodes=np.r_[start,original.stations[(original.stations>start)&(original.stations<stop)],stop]
            xy=p+(nodes-a)[:,None]/(b-a)*(q-p)
            error=float(np.linalg.norm(xy-original(nodes),axis=1).max())
            maximum=max(maximum,error);checks+=len(nodes)
            if checks>20000:raise ValueError('Independent run source-curb identity proof exceeds its bounded allocation')
            if error>1e-7:raise ValueError('Independent run field changed an original finite source curb chain')
        if covered<high-1e-8:raise ValueError('Independent run field does not cover all finite original curb intervals')
    return dict(proof='continuous_shared_original_piecewise_linear_curbs',maximum_error=maximum,checks=checks)


def _terminal_piece(context, approach, inventory):
    """A finite terminal connector can serve its exact incident run piece.

    The caller's whole assembly retains all other pieces separately. This
    admission never changes the full original run used for progress proofs.
    """
    from .model import End
    from .profiles import memberships
    whole = tuple(tuple(member.visits) for member in approach.run.members)
    full = {visit.section_id for chain in whole for visit in chain}
    declared = set(getattr(inventory, 'source_sections', full))
    if declared == full: return None
    piece = getattr(inventory, 'source_visits', None)
    if (getattr(inventory, 'entrance', None) is not approach or piece is None
            or len(piece) != 2 or not all(piece) or approach.end_index not in (0, 1)):
        raise ValueError('Chart transport run inventory changed the complete original source sections')
    piece = tuple(tuple(chain) for chain in piece)
    used = [visit.section_id for chain in piece for visit in chain]
    if set(used) != declared or len(used) != len(set(used)) or not declared < full:
        raise ValueError('Terminal piece inventory repeats or invents original source sections')
    by_end = {end: site for site in (*context.network.intersections, *context.network.connections) for end in site.ends}
    membership = memberships(context.network)
    bundles = set()
    for role, (original, selected) in enumerate(zip(whole, piece)):
        expected = original[:len(selected)] if approach.end_index == 0 else original[-len(selected):]
        if expected != selected:
            raise ValueError('Terminal source piece must be its exact incident original prefix or suffix')
        root = selected[0 if approach.end_index == 0 else -1]
        endpoint = End(root.section_id, int(not root.forward) if approach.end_index == 0 else int(root.forward))
        if endpoint != approach.original_ends[role]:
            raise ValueError('Terminal source piece changed its original physical root incidence')
        for visit in selected:
            record = membership.get(visit.section_id)
            if record is None or record[1:] != (role, visit.forward):
                raise ValueError('Terminal source piece changed original Bundle membership or direction')
            bundles.add(record[0].id)
        for a, b in zip(selected, selected[1:]):
            exit, entry = End(a.section_id, int(a.forward)), End(b.section_id, int(not b.forward))
            site = by_end.get(exit)
            if site is None or entry not in site.ends:
                raise ValueError('Terminal source piece skipped an original source incidence')
    if len(bundles) != 1:
        raise ValueError('Terminal source piece combines unrelated original Bundles')
    return dict(proof='exact_terminal_incident_original_run_piece',
        original_run_sections=[v.section_id for chain in whole for v in chain],
        piece_sections=used, terminal_end=approach.end_index,
        whole_original_run_progress_certified=True, requires_remaining_run_ownership_audit=True)


def _run_progress(context,approach,inventory):
    """Certify explicit original run visits without inventing interior Ends."""
    from .model import Key,End
    if hasattr(approach,'run'):
        visits=tuple(tuple(m.visits) for m in approach.run.members)
        if hasattr(inventory,'source_sections') and set(inventory.source_sections)!=set(approach.source_run_sections):
            _terminal_piece(context, approach, inventory)
    else:
        flags=getattr(approach,'original_incidence',())
        near=getattr(approach,'near_role',None)
        if (len(flags)!=2 or flags.count(True)!=1 or near not in (0,1) or not flags[near]
                or getattr(inventory,'entrance',None) is not approach):
            raise ValueError('Interior chart transport needs its explicit original/derived source incidence contract')
        run=getattr(inventory,'run',None)
        visits=(tuple(tuple(m.visits) for m in run.members) if run is not None
                else getattr(inventory,'source_visits',None))
        if visits is None or len(visits)!=2 or not all(visits):
            raise ValueError('Interior chart transport requires explicit original visits for both certified piece roles')
        declared=set(getattr(inventory,'source_sections',()))
        used=[v.section_id for role in visits for v in role]
        if not declared or set(used)!=declared or len(used)!=len(set(used)):
            raise ValueError('Interior chart transport piece inventory lost or repeats original source visits')
        if not set(approach.source_run_sections)<=declared:
            raise ValueError('Interior chart transport local original sources leave its certified piece inventory')
        from .profiles import memberships
        membership=memberships(context.network)
        bundles=set()
        for role,local in enumerate(approach.visits):
            if local not in visits[role]:raise ValueError('Interior chart transport changed a local original visit or role')
            record=membership.get(local.section_id)
            if record is None or record[1:]!=(role,local.forward):
                raise ValueError('Interior chart transport local visits lost their original Bundle provenance')
            bundles.add(record[0].id)
        if len(bundles)!=1:raise ValueError('Interior chart transport combines unrelated original Bundles')
        by_end={e:j for j in (*context.network.intersections,*context.network.connections) for e in j.ends}
        near_end=approach.original_ends[near]
        if near_end not in by_end:raise ValueError('Interior chart transport has no real original near incidence')
        for role,chain in enumerate(visits):
            for visit in chain:
                record=membership.get(visit.section_id)
                if record is None or record[0].id not in bundles or record[1:]!=(role,visit.forward):
                    raise ValueError('Interior chart transport piece changed an original Bundle role or orientation')
            for a,b in zip(chain,chain[1:]):
                exit=End(a.section_id,int(a.forward));entry=End(b.section_id,int(not b.forward))
                site=by_end.get(exit)
                if site is None or entry not in site.ends:
                    raise ValueError('Interior chart transport visits lack their original source incidence')
    return [inventory.frame.certify(context[Key('path',v.section_id)],forward=v.forward)
        for member in visits for v in member]


def _overlap(field,first,last):
    """Linear source-curb events certify the whole local overlap interval."""
    a,b=sorted((float(first),float(last)))
    nodes=np.unique(np.concatenate((np.array((a,b)),
        *(r.stations[(r.stations>a)&(r.stations<b)] for r in field.curves))))
    edges=field.edges(nodes);ao,ai,bi,bo=edges
    if np.any(ao<=bi) or np.any(ai<=bo) or np.any(np.maximum(ai,bo)>=np.minimum(ao,bi)-1e-7):
        raise ValueError('Chart transport crosses a separated or role-swapped original source interval')
    return len(nodes)


def build(context,approach,station,mouth,target_inventory,*,max_cells=2048):
    """Return a bounded direct-quad connector and exact front/back interfaces.

    Both fields use the same original sampled curb chains over the connector.
    The local field is queried only inside its certified finite range. A
    quintic field blend has zero first/second endpoint jets, retaining the
    physical front and independent body field at their exact complete cuts.
    """
    target=target_inventory.field
    if target_inventory.frame is not target.frame or target_inventory.frame is approach.frame_chart:
        raise ValueError('Finite entrance transport needs an independent, consistent certified run frame')
    if (not hasattr(target.frame,'certify') or target_inventory.field.datum!=approach.field.datum
            or abs(target_inventory.field.rise-approach.field.rise)>1e-9):
        raise ValueError('Chart transport needs the run certificate and original shared crown')
    if hasattr(target_inventory,'entrance') and target_inventory.entrance is not approach:
        raise ValueError('Chart transport run inventory belongs to another prepared physical entrance')
    progress=_run_progress(context,approach,target_inventory)
    mouth=np.asarray(mouth,float)
    if mouth.ndim!=2 or mouth.shape[1]!=3 or len(mouth)<3 or not np.isfinite(mouth).all():
        raise ValueError('Chart transport needs one actual finite whole mouth')
    initial=approach.cross_section(station);delta=mouth[-1,:2]-mouth[0,:2]
    q=2.*((mouth[:,:2]-mouth[0,:2])@delta/(delta@delta))-1.;q[[0,-1]]=-1.,1.
    if np.any(np.diff(q)<=0.) or not np.any(abs(q-initial.peak)<1e-8):
        raise ValueError('Chart transport needs an ordered whole mouth including its shared crest')
    expected=approach.point(station,.5*approach.width*q)
    if max(abs(expected-mouth[:,:2]).max(),abs(approach.field.values(mouth[:,:2])-mouth[:,2]).max())>1e-8:
        raise ValueError('Chart transport disagrees with its authoritative finite source mouth')
    sides=q<=initial.peak
    fractions=(q-initial.peak)/np.where(sides,1.+initial.peak,1.-initial.peak)
    across=len(q)-1;sign=approach.sign
    projected=target.frame.local(mouth[:,:2])[:,0]
    first=float(projected.max() if sign>0. else projected.min())
    failures=[]
    source_budget=.1*context.network.settings.curve_tolerance
    projection_budget=min(approach.source_curve_tolerance,source_budget-approach.source_curve_tolerance-1e-7)
    if projection_budget<=0. or source_budget>.002:
        raise ValueError('Chart transport has no remaining original source-curve error allocation')
    for attempt,factor in enumerate((2.,3.,4.,6.,8.)):
        at=first+sign*factor*context.network.settings.spacing
        try:
            back_xy=_cut(target,at,sign,fractions,sides)
            if not _ranges(target,back_xy) or not _ranges(approach.field,back_xy):
                raise ValueError('Independent body cut exceeds the certified finite original entrance range')
            rays=(approach.field.curves[0],approach.field.curves[3]) if sign>0. else (approach.field.curves[3],approach.field.curves[0])
            traced=[_trace(ray,mouth[i,:2],back_xy[i],approach.frame_chart,sign) for ray,i in zip(rays,(0,-1))]
            length=max(t[-1] for points,t in traced)
            count=max(2,int(math.ceil(length/context.network.settings.spacing)))
            paths=approach.paths if sign>0. else approach.paths[::-1]
            rows,projection=_rows(traced,count,[p.exact for p in paths],projection_budget)
            if (len(rows)-1)*across>max_cells:raise ValueError('Chart transport exceeds its bounded direct-quad allocation')
            left,right=(_at(points,lengths,rows) for points,lengths in traced)
            weight=_blend(rows);center=.5*(left+right);edge=right-left;squared=np.sum(edge*edge,axis=1)
            if np.min(squared)<=1e-12:raise ValueError('Chart transport collapses its complete profile')
            crest=[]
            for field in (approach.field,target):
                x=field.frame.local(center)[:,0]
                if not _ranges(field,np.vstack((left,right,center))):
                    raise ValueError('Chart transport would extrapolate one of its prescribed fields')
                p1,p2=field.peaks(x)
                if np.max(abs(p1-p2))>1e-9:raise ValueError('Chart transport reaches a true split crown')
                peak=field.frame.world(np.c_[x,p1])
                crest.append(np.sum((peak-left)*edge,axis=1)/squared)
            peak=(1.-weight)*crest[0]+weight*crest[1]
            if np.any(peak<=1e-8) or np.any(peak>=1.-1e-8):raise ValueError('Chart transport loses a positive crown half')
            u=peak[:,None]+fractions[None,:]*np.where(sides[None,:],peak[:,None],1.-peak[:,None])
            xy=(left[:,None]+u[...,None]*edge[:,None]).reshape(-1,2)
            xy[:across+1]=mouth[:,:2];xy[-across-1:]=back_xy
            if not _ranges(approach.field,xy) or not _ranges(target,xy):
                raise ValueError('Chart transport leaves the common finite source-field support')
            x=approach.frame_chart.local(xy)[:,0]
            events=_overlap(approach.field,x.min(),x.max())
            identity=_same_source_curbs(approach,target,x.min(),x.max())
            source_z=approach.field.values(xy);target_z=target.values(xy)
            w=np.repeat(weight,across+1)
            z=(1.-w)*source_z+w*target_z
            xyz=np.c_[xy,z];front=np.arange(across+1);back=np.arange(len(xyz)-across-1,len(xyz))
            xyz[front]=mouth;xyz[back,2]=target.values(back_xy)
            side_ids=np.r_[np.arange(len(rows))*(across+1),np.arange(len(rows))*(across+1)+across]
            if np.max(abs(xyz[side_ids,2]-target.datum))>1e-8:
                raise ValueError('Chart transport changed an original outer curb datum')
            interior_sides=side_ids[~np.isin(side_ids,np.r_[front,back])]
            xyz[interior_sides,2]=target.datum
            faces=context.patterns.grid(len(rows)-1,across).quads[:,::-1]
            mesh=Mesh(xyz,faces);report=quality(mesh)
            stored=Mesh(xyz.astype(np.float32).astype(float),faces);published=quality(stored)
            mesh.validate();stored.validate()
            if max(report['maximum_edge_ratio'],published['maximum_edge_ratio'])>8. or min(report['minimum_corner_sine'],published['minimum_corner_sine'])<.03:
                raise ValueError('Chart transport fails unchanged double/stored quad quality '
                    f"({published['maximum_edge_ratio']:.6g}/{published['minimum_corner_sine']:.6g})")
            gradient=np.einsum('fcij,fj->fci',corner_matrices(stored),stored.xyz[stored.quads,2]-target.datum)
            slope=float(np.linalg.norm(gradient,axis=-1).max())
            along=max(approach.frame_chart.longitudinal_bound(stored),target.frame.longitudinal_bound(stored))
            if slope>.1+1e-6 or along>.03+1e-8:raise ValueError('Chart transport exceeds physical slope limits')
            velocity=np.broadcast_to(approach.direction,(len(xyz),2)).copy()
            velocity[front]=approach.jet(station,.5*approach.width*q)[1]
            epsilon=1e-4
            velocity[back]=(_cut(target,at+epsilon,sign,fractions,sides)-_cut(target,at-epsilon,sign,fractions,sides))/(2.*epsilon)
            joins=[EntranceSlopeGuard().validate(stored,ids,velocity,stored.xyz[:,2]-target.datum) for ids in (front,back)]
            fixed=np.zeros(len(xyz),bool);fixed[np.r_[front,back,side_ids]]=True
            if BoundaryDrainage(max_cells=max_cells).fill(stored,fixed)[1]['maximum_fill']>1e-6:
                raise ValueError('Chart transport creates an interior drainage minimum')
            if not np.array_equal(mesh.xyz[front],mouth):raise ValueError('Chart transport moved the actual junction mouth')
            containment=_contains(mesh,approach,projection)
            if containment['maximum_combined_error']>min(.002,source_budget)+1e-12:
                raise ValueError('Chart transport exceeds its composed original source-curve error budget')
            report = dict(construction='finite source-curb transport between independent charts',
                quads=len(faces),quality=report,published_quality=published,exact_mouth=True,
                physical_slope=slope,chart_longitudinal=along,entrance_bounds=joins,
                target_station=at,start=float(station),attempts=attempt+1,continuous_overlap_events=events,
                local_field_extrapolation=False,original_outer_curb_chords_exact=projection<=1e-12,
                maximum_source_polyline_projection_error=projection,source_polyline_projection_allowance=projection_budget,
                source_containment=containment,original_source_curve_budget=source_budget,
                source_curve_identity=identity,whole_original_run_progress_lower=min(progress),
                source_sections=[p.section.id for p in approach.paths],shared_crest_column=int(np.argmin(abs(q-initial.peak))),
                all_direct_quads=True,previous_attempt_failures=failures)
            if hasattr(approach, 'run'):
                piece = _terminal_piece(context, approach, target_inventory)
                if piece is not None: report['original_terminal_run_piece'] = piece
            return mesh,front,back,report
        except ValueError as failure:failures.append(str(failure))
    raise ValueError('No qualified finite chart-aware merged profile bridge: '+'; '.join(failures))
