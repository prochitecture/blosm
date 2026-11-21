"""Remove a proven redundant internal lobe from a finite source fillet patch.

No physical curb is removed. A single crossing of construction centerline
seams may enclose a lobe already covered by one of those actual source strips.
Only that continuously contained internal lobe can be omitted. The resulting
curb-bearing patch must still pass normal region and complete-union checks.
"""
import math
import numpy as np
from .mesh import cross
from .regions import Region, RegionUnion, signed_area
from .finite_containment import contains_segment, _capsule_intervals
from .model import Key


def _source_intervals(corner):
    for walk in corner.source_contact_walks:
        for piece in walk.prefix:
            yield piece.path, 0., piece.path.length
    for walk, distance in zip(corner.source_contact_walks, corner.stations):
        approach=walk.approach
        if not hasattr(approach.path,'pieces'):
            yield approach.path, *((0.,distance) if approach.end==0 else
                                   (approach.path.length-distance,approach.path.length))
            continue
        for index,piece in enumerate(approach.path.pieces):
            extent=min(distance,float(approach.path.stations[index+1]))-float(approach.path.stations[index])
            if extent<=1e-9:break
            yield piece.path, *((0.,extent) if piece.forward else (piece.path.length-extent,piece.path.length))


def _covered(lobe,corner,spacing,tolerance):
    maximum=4096
    epsilon=min(.00025,tolerance)
    used=0
    for path,first,last in _source_intervals(corner):
        if not 0.<=first<last<=path.length:
            raise ValueError('Internal closure proof cannot extrapolate an original source interval')
        section=path.section
        reach=abs(section.offset)+.5*section.width
        if path.exact:
            knots=np.asarray(path.knots)
            at=np.r_[first,knots[(knots>first)&(knots<last)],last]
            chord_error=0.
        else:
            k=float(path.curvature_bound)
            if reach*k>=1.:
                raise ValueError('Internal closure proof needs unfolded original offset curves')
            curvature=k/(1.-reach*k);speed=1.+reach*k
            step=min(.3,spacing)
            if curvature>1e-12:
                step=min(step,math.sqrt(8.*epsilon/curvature)/speed,math.pi/(16.*curvature*speed))
            if (last-first)/step>maximum-used:
                raise ValueError('Internal closure source proof exceeds its finite cell budget')
            at=path.samples(first,last,step,tolerance=epsilon,max_samples=maximum-used+1)
            chord_error=epsilon
        used+=len(at)-1
        if used>maximum:
            raise ValueError('Internal closure source proof exceeds its finite cell budget')
        high=path.point(at,section.offset+.5*section.width)
        low=path.point(at,section.offset-.5*section.width)
        polygon=np.vstack((high,low[::-1]))
        region=Region(Key('internal_closure_source',0),polygon if signed_area(polygon)>0. else polygon[::-1])
        RegionUnion().build((region,))
        edges=tuple(zip(lobe,np.roll(lobe,-1,axis=0)))
        if not all(contains_segment(a,b,region.xy,tolerance=1e-10) for a,b in edges):
            continue
        # The transverse endpoints are exact source cuts. Longitudinal curbs
        # require enough interior clearance for their analytical chord error.
        banks=(*zip(high[:-1],high[1:]),*zip(low[:-1],low[1:]))
        if any(_capsule_intervals(a,b-a,np.array((x,y)),chord_error+1e-7)
               for a,b in edges for x,y in banks):
            continue
        return dict(section=section.id,source_interval=[first,last],source_cells=used,
            longitudinal_chord_allowance=chord_error,continuous_internal_containment=True,
            no_physical_curb_removed=True,original_source_dimensions=True)
    raise ValueError('A redundant internal closure lobe is not covered by an actual finite source strip')


def normalize(points,kinds,corner,spacing,tolerance,*,maximum_pairs=65536):
    points=np.asarray(points,float)
    if len(points)!=len(kinds) or not getattr(corner,'source_contact_walks',None):
        raise ValueError('Internal closure needs its explicit actual source contact walks')
    internal=np.flatnonzero(np.asarray(kinds)=='internal')
    hits=[];tested=0
    for at,i in enumerate(internal):
        js=internal[at+1:]
        js=js[(js!=i+1)&~((i==0)&(js==len(points)-1))]
        tested+=len(js)
        if tested>maximum_pairs:
            raise ValueError('Internal closure intersections exceed their finite pair budget')
        if not len(js):continue
        a=points[i];d=points[(i+1)%len(points)]-a
        b=points[js];e=points[(js+1)%len(points)]-b
        determinant=cross(d,e)
        proper=abs(determinant)>1e-12*np.linalg.norm(d)*np.linalg.norm(e,axis=1)
        ta=np.divide(cross(b-a,e),determinant,out=np.zeros(len(js)),where=proper)
        tb=np.divide(cross(b-a,d),determinant,out=np.zeros(len(js)),where=proper)
        valid=proper&(ta>1e-8)&(ta<1.-1e-8)&(tb>1e-8)&(tb<1.-1e-8)
        hits.extend((int(i),int(js[k]),a+ta[k]*d) for k in np.flatnonzero(valid))
        if len(hits)>1:
            raise ValueError('Internal closure needs a different policy for multiple crossing seams')
    if len(hits)!=1:
        raise ValueError('Internal closure has no single proper construction-seam crossing')
    i,j,point=hits[0]
    loops=((np.vstack((point,points[i+1:j+1])),[kinds[i],*kinds[i+1:j],kinds[j]]),
           (np.vstack((point,points[j+1:],points[:i+1])),[kinds[j],*kinds[j+1:],*kinds[:i],kinds[i]]))
    discarded=[index for index,(xy,labels) in enumerate(loops) if all(label=='internal' for label in labels)]
    if len(discarded)!=1:
        raise ValueError('An internal closure crossing would change a physical source curb')
    drop=discarded[0]
    lobe=loops[drop][0]
    proof=_covered(lobe,corner,spacing,tolerance)
    result,labels=loops[1-drop]
    if signed_area(result)<0.:
        result=result[::-1];labels=list(np.roll(np.asarray(labels)[::-1],-1))
    RegionUnion().build((Region(Key('internal_closure_result',0),result),))
    proof.update(discarded_area=abs(signed_area(lobe)),source_seam_pair=[i,j],
                 intersections_tested=tested,derived_internal_cut=point.tolist())
    return result,labels,proof
