"""Extend an interior computational diagonal to its first physical curb.

Only computational lines change. Every returned reference addresses an
unchanged physical boundary chord; profiles, source curves and holes remain
authoritative. A line cannot pass a physical boundary to reach a later bank.
"""
import numpy as np
from .mesh import cross
from .partition import BoundaryRef
from .finite_containment import segment_on_polygon_side


def extend(group, first, last, *, curb_edges, first_ref=None, last_ref=None, max_segments=20000):
    """Return two ordinary BoundaryRefs for one visible provisional diagonal.

    Existing physical endpoints retain their exact references. An interior
    endpoint extends along the same line to the closest physical boundary in
    that direction, and that first boundary must be an eligible true curb.
    """
    first,last=np.asarray(first,float),np.asarray(last,float)
    if (first.shape!=(2,) or last.shape!=(2,) or not np.isfinite((first,last)).all()
            or type(max_segments) is not int or not 1<=max_segments<=20000):
        raise ValueError('Computational cut extension needs finite endpoints and a bounded rim')
    delta=last-first;length=float(np.linalg.norm(delta));epsilon=1e-7
    if length<=epsilon:raise ValueError('Computational diagonal has no resolved length')
    if sum(len(loop.xy) for loop in group.loops)>max_segments:
        raise ValueError('Computational cut extension exceeds its physical rim segment budget')
    if sum(not loop.hole for loop in group.loops)!=1:
        raise ValueError('Computational cut extension needs one physical exterior')
    curb_edges=frozenset(curb_edges);hits=[];parallel=[];dt=epsilon/length
    for li,loop in enumerate(group.loops):
        a=loop.xy;b=np.roll(a,-1,axis=0);edge=b-a
        denominator=cross(delta,edge);offset=a-first
        valid=abs(denominator)>1e-12*length*np.linalg.norm(edge,axis=1)
        t=np.divide(cross(offset,edge),denominator,out=np.zeros(len(a)),where=valid)
        u=np.divide(cross(offset,delta),denominator,out=np.zeros(len(a)),where=valid)
        for ei in np.flatnonzero(valid & (u>=-1e-10) & (u<=1.+1e-10)):
            hits.append((float(t[ei]),BoundaryRef(li,int(ei),float(np.clip(u[ei],0.,1.)))))
        for ei in np.flatnonzero(~valid & (abs(cross(offset,delta))<=epsilon*length)):
            values=np.array(((a[ei]-first)@delta,(b[ei]-first)@delta))/(length*length)
            parallel.append((float(values.min()),float(values.max())))
    if not hits:raise ValueError('Computational diagonal has no finite physical boundary contacts')
    def point(ref):
        if (not isinstance(ref,BoundaryRef) or not 0<=ref.loop<len(group.loops)
                or not 0<=ref.edge<len(group.loops[ref.loop].xy) or not 0.<=ref.fraction<=1.):
            raise ValueError('Computational cut lost an authoritative physical boundary reference')
        polygon=group.loops[ref.loop].xy
        return polygon[ref.edge]+ref.fraction*(polygon[(ref.edge+1)%len(polygon)]-polygon[ref.edge])
    def endpoint(ref,origin,lower):
        if ref is not None:
            if np.linalg.norm(point(ref)-origin)>epsilon:
                raise ValueError('Computational diagonal moved an existing physical endpoint')
            if (ref.loop,ref.edge) not in curb_edges:
                raise ValueError('Computational cut cannot split an incoming physical profile')
            return (0. if lower else 1.),ref
        choices=[hit for hit in hits if hit[0]<-dt] if lower else [hit for hit in hits if hit[0]>1.+dt]
        if not choices:raise ValueError('Interior computational endpoint has no first finite physical curb')
        station=(max if lower else min)(value for value,ref in choices)
        tied=[ref for value,ref in choices if abs(value-station)<=dt]
        # Boundary vertices can have two incident edge records. No profile
        # edge may be replaced by a nearby curb or a farther intersection.
        if any((ref.loop,ref.edge) not in curb_edges for ref in tied):
            raise ValueError('First physical cut contact reaches an incoming profile')
        return station,min(tied,key=lambda ref:(ref.loop,ref.edge,ref.fraction))
    lo,a=endpoint(first_ref,first,True);hi,b=endpoint(last_ref,last,False)
    if hi-lo<=dt:raise ValueError('Extended computational cut has no positive pavement interval')
    if any(lo+dt<t<hi-dt for t,ref in hits):
        raise ValueError('Computational extension crosses a physical hole or boundary before its first curb')
    if any(min(hi,q)-max(lo,p)>dt for p,q in parallel):
        raise ValueError('Computational extension overlaps a physical boundary interval')
    left,right=point(a),point(b)
    for loop in group.loops:
        if not segment_on_polygon_side(left,right,loop.xy,inside=not loop.hole,tolerance=epsilon):
            raise ValueError('Extended computational cut leaves the complete physical pavement')
    return (a,b),dict(original_length=length,extended_length=float(np.linalg.norm(right-left)),
        first_extension_m=max(0.,-lo*length),last_extension_m=max(0.,(hi-1.)*length),
        first_ref_preserved=first_ref is not None,last_ref_preserved=last_ref is not None,
        first_physical_curb_only=True,physical_rim_changed=False,profiles_changed=False,
        physical_segments_tested=sum(len(loop.xy) for loop in group.loops))
