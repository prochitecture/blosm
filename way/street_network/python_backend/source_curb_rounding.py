"""Round an exposed intersection of two explicitly owned original curbs.

The small added wedge closes along the actual source curbs to their physical
intersection. It contains one certified circular arc and no invented node or
centreline bridge. Source support and all original profile boundaries remain
the caller's mandatory ownership checks.
"""
from types import SimpleNamespace
from dataclasses import dataclass
import math
import numpy as np
from .curb_contact import CircularCurbContact
from .compounds import RimSegment,SampledFootprint
from .finite_free_caps import PreparedCapCurb
from .regions import Region,RegionUnion,signed_area


@dataclass(frozen=True)
class FiniteSourceBank:
    """An original bank and its declared finite native support intervals.

    These are geometric parameters, not new junction incidences. Ownership of
    the intervals and exposure of the bank switch remain caller proof gates.
    Disjoint owned ranges must stay disjoint; a fillet cannot bridge the gap.
    """
    path: object
    side: int
    intervals: tuple

    def __post_init__(self):
        from .paths import Path
        if type(self.path) is not Path or type(self.side) is not int or self.side not in (-1,1):
            raise ValueError('Finite source bank needs an original Path and declared native side')
        intervals=tuple(tuple(row) for row in self.intervals)
        if (not 1<=len(intervals)<=20000 or any(len(row)!=2 for row in intervals)
                or any(not all(math.isfinite(v) for v in row)
                       or not 0.<=row[0]<row[1]<=self.path.length for row in intervals)
                or any(a[1]>b[0] for a,b in zip(intervals,intervals[1:]))):
            raise ValueError('Finite source bank needs ordered disjoint original station ranges')
        object.__setattr__(self,'intervals',intervals)

    def contains(self,first,last):
        low,high=sorted((float(first),float(last)))
        return any(a<=low and high<=b for a,b in self.intervals)


@dataclass(frozen=True)
class SourceBankFillet:
    corner: object
    intersection: tuple
    banks: tuple


def fit_banks(context,first,last,radius,mode,guess):
    """Fit the exposed CCW bank transition with the existing circle solver.

    The arriving bank is traversed toward its right-bank Approach; the leaving
    bank uses its left-bank Approach. Native sides determine those parameter
    orientations, independently of distant original junctions. No Route may
    silently acquire additional Sections beyond either finite declared range.
    """
    from .model import Key
    from .intersections import Approach
    from .corner_radii import CornerRadiusPolicy,CornerContactError
    if type(first) is not FiniteSourceBank or type(last) is not FiniteSourceBank:
        raise ValueError('Source bank fillet needs two explicit finite native banks')
    if (not math.isfinite(radius) or radius<=0. or mode not in ('STRICT','STYLE')
            or np.asarray(guess).shape!=(2,) or not np.isfinite(guess).all()):
        raise ValueError('Source bank fillet needs a finite crossing and positive prescribed radius')
    banks=(first,last);approaches=[]
    for bank,end in zip(banks,(int(first.side==1),int(last.side==-1))):
        if context[Key('path',bank.path.section.id)] is not bank.path:
            raise ValueError('Source bank fillet changed its original Path identity')
        approaches.append(Approach(bank.path,end,context[Key('profile',bank.path.section.id)]))
    def native(arm,station):return station if arm.end==0 else arm.length-station
    accepted={}
    def accept(corner):
        if corner.center is None:return False
        try:crossing=intersection_stations(corner,np.asarray(guess))
        except ValueError:return False
        if not all(bank.contains(native(arm,at),native(arm,cut))
                   for bank,arm,at,cut in zip(banks,approaches,corner.stations,crossing)):
            return False
        accepted['intersection']=crossing
        return True
    corner=CornerRadiusPolicy.solve(*approaches,radius,mode,0.,accept=accept)
    if 'intersection' not in accepted:
        raise CornerContactError('Source bank fillet has no crossing in the declared original ranges')
    return SourceBankFillet(corner,accepted['intersection'],banks)


def intersection_stations(corner,guess):
    if corner.center is None or not math.isfinite(corner.radius) or corner.radius<=0.:
        raise ValueError('Source curb rounding needs an actual positive circle')
    approaches=(corner.first,corner.second)
    lower=np.maximum(0.,np.asarray(corner.stations)-4.*corner.radius)
    upper=np.minimum([a.length for a in approaches],np.asarray(corner.stations)+4.*corner.radius)
    initial=[]
    for i,a in enumerate(approaches):
        at=np.linspace(lower[i],upper[i],65)
        points=a.point(at,a.offset+(-.5 if i==0 else .5)*a.width)
        initial.append(float(at[np.argmin(np.linalg.norm(points-guess,axis=1))]))
    solver=CircularCurbContact(*approaches,0.,1.)
    solver.tolerance=max(1e-11,8.*np.finfo(float).eps*max(1.,float(np.max(abs(guess)))))
    values,point,contacts=solver.solve(initial,lower,upper)
    if np.linalg.norm(contacts[0]-contacts[1])>1e-8:
        raise ValueError('Source curb intersection lacks one finite original contact')
    return tuple(map(float,values))


def patch(owner,corner,intersection,spacing,tolerance,*,partitions=None,max_samples=8192):
    if (not all(math.isfinite(v) and v>0. for v in (spacing,tolerance))
            or type(max_samples) is not int or max_samples<4):raise ValueError('Invalid bounded curb rounding samples')
    approaches=(corner.first,corner.second)
    for a,contact,cut in zip(approaches,corner.stations,intersection):
        if not (0.<=contact<=a.length and 0.<=cut<=a.length):
            raise ValueError('Source curb rounding cannot extrapolate a source endpoint')
    from .finite_scope_spurs import _source_stations
    chains=[]
    for i,(a,first,last) in enumerate(zip(approaches,intersection,corner.stations)):
        at=(_source_stations(a,max(first,last),min(.2,spacing),tolerance)
            if partitions is None else np.asarray(partitions[i]))
        values=np.asarray(sorted({first,last,*(float(v) for v in at if min(first,last)<v<max(first,last))}))
        if last<first:values=values[::-1]
        chains.append(a.point(values,a.offset+(-.5 if i==0 else .5)*a.width))
    if np.linalg.norm(chains[0][0]-chains[1][0])>1e-8:
        raise ValueError('Source rounding needs the same original curb intersection')
    points=[chains[0][0]];kinds=[]
    def append(chain,kind):
        for point in chain:
            if np.linalg.norm(point-points[-1])<=1e-8:continue
            points.append(np.asarray(point));kinds.append(kind)
            if len(points)>max_samples:raise ValueError('Source curb rounding exceeds its finite sample budget')
    append(chains[0][1:],'internal')
    count=max(2,int(math.ceil(abs(corner.sweep)*corner.radius/min(spacing,math.sqrt(8.*tolerance*corner.radius)))))
    if count+len(points)>max_samples:raise ValueError('Source curb rounding exceeds its finite sample budget')
    theta=np.linspace(corner.start_angle,corner.start_angle+corner.sweep,count+1)
    arc=corner.center+corner.radius*np.column_stack((np.cos(theta),np.sin(theta)))
    contacts=np.asarray((chains[0][-1],chains[1][-1]))
    scale=max(1.,corner.radius,float(np.max(abs(corner.center))),float(np.max(abs(contacts))))
    precision=max(1e-11,32.*np.spacing(scale))
    if np.max(np.linalg.norm(arc[[0,-1]]-contacts,axis=1))>precision:
        raise ValueError('Source-bank circle endpoints lack representable native contact agreement')
    # These are the explicitly solved same-source contacts, not nearby points
    # found by an XY merge. Use their canonical native coordinates at both
    # ends; a trig conversion can otherwise lose one ulp and clip a nearly
    # tangent arc's declared endpoint. Interior circle samples stay untouched.
    arc[0],arc[-1]=contacts
    append(arc[1:],'curb');append(chains[1][::-1][1:],'internal')
    if np.linalg.norm(points[-1]-points[0])>1e-8:
        raise ValueError('Source curb rounding needs an invented closing bridge')
    points=np.asarray(points[:-1])
    if signed_area(points)<0.:
        points=points[::-1];kinds=list(np.roll(np.asarray(kinds)[::-1],-1))
    start=next((i for i in range(len(points)) if kinds[i]!=kinds[i-1]),0)
    points=np.roll(points,-start,axis=0);kinds=kinds[start:]+kinds[:start]
    region=Region(owner,points);RegionUnion().build((region,))
    curves=[];labels=[];first=0
    while first<len(points):
        last=first+1
        while last<len(points) and kinds[last]==kinds[first]:last+=1
        curve=PreparedCapCurb(np.vstack((points[first:last],points[last%len(points)])))
        index=len(curves);curves.append(curve)
        labels.extend(RimSegment(kinds[first],index,float(a/curve.length),float(b/curve.length))
            for a,b in zip(curve.stations,curve.stations[1:]))
        first=last
    return SampledFootprint(region,tuple(labels),SimpleNamespace(approaches=(),boundaries=tuple(curves),
        rows=len(points),mouths={},corners=(corner,)))
