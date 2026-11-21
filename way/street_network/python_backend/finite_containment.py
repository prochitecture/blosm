"""Continuous containment of finite linear source boundaries.

The physical exterior must already be certified simple and without holes.
Every boundary intersection partitions a source segment into intervals whose
interior classification is constant. Testing those intervals also catches a
source exit shorter than any fixed diagnostic sampling spacing.
"""
import numpy as np
from .mesh import cross
from .regions import _classify


def contains_segment(first,last,polygon,*,tolerance=1e-7,classifier=None):
    """Include every point of a finite segment in a simple physical exterior.

    Coincident endpoints and boundary-aligned intervals are permitted. Source
    endpoints are never displaced or extrapolated to obtain containment.
    """
    return segment_on_polygon_side(first,last,polygon,inside=True,tolerance=tolerance,classifier=classifier)


def separated_segment(first,last,polygon,*,clearance):
    """Keep a complete segment outside a polygon and its boundary capsules."""
    if not np.isfinite(clearance) or clearance<=0.:
        raise ValueError('Finite segment separation needs a positive numerical clearance')
    if not segment_on_polygon_side(first,last,polygon,inside=False,tolerance=1e-10):return False
    first,last=np.asarray(first,float),np.asarray(last,float)
    return not _capsule_intervals(first,last-first,np.asarray(polygon,float),clearance)


def _capsule_intervals(first,delta,polygon,radius):
    """Exact line intervals inside the union of finite boundary capsules."""
    result=[]
    speed2=float(delta@delta)
    next_points=np.roll(polygon,-1,axis=0)
    low=np.minimum(first,first+delta)-radius
    high=np.maximum(first,first+delta)+radius
    nearby=np.all(np.maximum(polygon,next_points)>=low,axis=1)&np.all(np.minimum(polygon,next_points)<=high,axis=1)
    for a,b in zip(polygon[nearby],next_points[nearby]):
        edge=b-a;length=float(np.linalg.norm(edge))
        if length<=1e-15:continue
        axis=edge/length;normal=np.array((-axis[1],axis[0]))
        low,high=0.,1.
        for value,speed,minimum,maximum in (
                (float((first-a)@axis),float(delta@axis),0.,length),
                (float((first-a)@normal),float(delta@normal),-radius,radius)):
            if abs(speed)<=1e-15:
                if value<minimum or value>maximum:low,high=1.,0.;break
            else:
                p,q=(minimum-value)/speed,(maximum-value)/speed
                low=max(low,min(p,q));high=min(high,max(p,q))
        if low<=high:result.append((low,high))
        if speed2<=1e-30:continue
        for center in (a,b):
            offset=first-center
            center_at=-float(offset@delta)/speed2
            nearest=offset+center_at*delta
            remaining=radius*radius-float(nearest@nearest)
            if remaining<0.:continue
            reach=np.sqrt(remaining/speed2)
            low=max(0.,center_at-reach);high=min(1.,center_at+reach)
            if low<=high:result.append((low,high))
    return sorted(result)


def _intersection_stations(first, delta, polygon, edge, lengths):
    offset = polygon-first
    determinant = cross(delta, edge)
    valid = abs(determinant) > 1e-12*np.linalg.norm(delta)*lengths
    t = np.divide(cross(offset, edge), determinant, out=np.zeros(len(edge)), where=valid)
    u = np.divide(cross(offset, delta), determinant, out=np.zeros(len(edge)), where=valid)
    hits = t[valid & (t>0.) & (t<1.) & (u>=-1e-10) & (u<=1.+1e-10)]
    return np.r_[0., np.sort(hits), 1.]


def segment_on_polygon_side(first,last,polygon,*,inside=True,tolerance=1e-7,classifier=None):
    """Continuously stay on one polygon side, with a bounded boundary band.

    For each interval on the wrong side, finite segment/disk capsules cover
    the *entire* interval. A single near-boundary midpoint cannot excuse a
    longer excursion elsewhere in that interval.
    """
    first,last,polygon=(np.asarray(value,float) for value in (first,last,polygon))
    if (first.shape!=(2,) or last.shape!=(2,) or polygon.ndim!=2 or polygon.shape[1]!=2
            or len(polygon)<3 or not all(np.isfinite(value).all() for value in (first,last,polygon))
            or not np.isfinite(tolerance) or tolerance<=0.):
        raise ValueError('Finite source containment needs valid planar inputs')
    delta=last-first
    edge=np.roll(polygon,-1,axis=0)-polygon
    stations=_intersection_stations(first,delta,polygon,edge,np.linalg.norm(edge,axis=1))
    middle=.5*(stations[:-1]+stations[1:])
    classify=_classify if classifier is None else classifier
    hit,boundary=classify(first+middle[:,None]*delta,polygon,1e-10)
    wrong=np.flatnonzero((hit!=inside)&(boundary<0))
    endpoints,near=classify(np.array((first,last)),polygon,tolerance)
    if np.any((endpoints!=inside)&(near<0)):return False
    return _within_capsules(first,delta,polygon,tolerance,stations,wrong)


def _within_capsules(first, delta, polygon, tolerance, stations, wrong):
    if not len(wrong):return True
    capsules=_capsule_intervals(first,delta,polygon,tolerance)
    for index in wrong:
        covered=float(stations[index]);end=float(stations[index+1])
        for low,high in capsules:
            if high<covered-1e-12:continue
            if low>covered+1e-12:break
            covered=max(covered,high)
            if covered>=end-1e-12:break
        if covered<end-1e-12:return False
    return True


def segments_on_polygon_side(first, last, polygon, *, inside, tolerance, classifier):
    """Apply the scalar proof to edges sharing one polygon.

    Each intersection partitions only its own segment. Classify endpoints
    and interval midpoints together, then use the original finite capsules
    on each wrong-side interval. The offsets map witnesses back to segments.
    """
    first, last, polygon = (np.asarray(value, float) for value in (first, last, polygon))
    if (first.ndim != 2 or first.shape[1] != 2 or last.shape != first.shape
            or polygon.ndim != 2 or polygon.shape[1] != 2 or len(polygon) < 3
            or not all(np.isfinite(value).all() for value in (first, last, polygon))
            or not np.isfinite(tolerance) or tolerance <= 0.):
        raise ValueError('Finite source containment needs valid planar inputs')
    edge = np.roll(polygon, -1, axis=0)-polygon
    lengths = np.linalg.norm(edge, axis=1)
    result = np.empty(len(first), bool)
    # Bound the number of interval witnesses, even for a highly concave rim.
    batch = max(1, min(64, 32768 // len(polygon)))
    for start in range(0, len(first), batch):
        a, b = first[start:start+batch], last[start:start+batch]
        points, indices = np.unique(np.concatenate((a, b)), axis=0, return_inverse=True)
        hit, near = classifier(points, polygon, tolerance)
        valid = ((hit == inside) | (near >= 0))[indices.reshape(2, len(a))].all(axis=0)
        delta = b-a
        stations = [_intersection_stations(p, d, polygon, edge, lengths) for p, d in zip(a, delta)]
        mid = [.5*(t[:-1]+t[1:]) for t in stations]
        counts = [len(values) for values in mid]
        offsets = np.r_[0, np.cumsum(counts)]
        segment = np.repeat(np.arange(len(a)), counts)
        hit, near = classifier(a[segment]+np.concatenate(mid)[:,None]*delta[segment], polygon, 1e-10)
        wrong = (hit != inside) & (near < 0)
        for i in np.unique(segment[wrong & valid[segment]]):
            intervals = np.flatnonzero(wrong[offsets[i]:offsets[i+1]])
            valid[i] = _within_capsules(a[i], delta[i], polygon, tolerance, stations[i], intervals)
        result[start:start+len(a)] = valid
    return result


def contains_polygon(source,exterior,*,tolerance=1e-7):
    """Contain a finite source polygon inside a certified hole-free exterior."""
    source=np.asarray(source,float)
    if source.ndim!=2 or source.shape[1]!=2 or len(source)<3:
        raise ValueError('Finite source footprint needs a planar polygon')
    return all(contains_segment(first,last,exterior,tolerance=tolerance)
        for first,last in zip(source,np.roll(source,-1,axis=0)))
