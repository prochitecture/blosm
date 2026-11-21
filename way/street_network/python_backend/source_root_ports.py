"""Actual original source ports used by failed cyclic-root construction.

One source bank is represented once. A second polygon approximation of the
same native bank is not a valid independent physical operand.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import math
import numpy as np

from .model import End,Key
from .compounds import RimSegment,SampledFootprint
from .regions import Region,RegionUnion
from .source_curb_coordinates import CurbCoordinates
from .source_identity import matches_current_section


class OriginalRootPortUnavailable(ValueError):
    """No complete original native-port reference is available."""


@dataclass(frozen=True)
class OriginalRootPort:
    footprint: object
    approach: object
    station: float

    @classmethod
    def build(cls,context,source,approach):
        if getattr(getattr(approach, 'profile', None), 'half', False):
            raise OriginalRootPortUnavailable('Root reference needs its actual centred full original source port')
        return cls._from_native_geometry(context, source, approach)

    @classmethod
    def _from_native_geometry(cls,context,source,approach):
        from .compound_cycles import PolylineCurb
        from .intersections import Approach
        layout=source.layout
        if (type(approach) is not Approach or tuple(getattr(layout,'approaches',()))!=(approach,)
                or len(getattr(layout,'outward_stations',()))!=1 or len(getattr(layout,'boundaries',()))!=3
                or not getattr(layout,'original_root_caps',None)
                or not matches_current_section(context,approach.source)
                or approach.width!=approach.source.width or approach.offset!=0.
                or approach.profile!=context[Key('profile',approach.source.id)].oriented(approach.end)):
            raise OriginalRootPortUnavailable('Root reference needs its actual centred full original source port')
        if any(getattr(getattr(approach,name),'__func__',None) is not getattr(Approach,name)
                for name in ('frame','point','jet','physical_gradient')):
            raise OriginalRootPortUnavailable('Root reference cannot substitute a derived original-source frame')
        station=layout.outward_stations[0]
        if not math.isfinite(station) or not 0.<station<approach.length:
            raise OriginalRootPortUnavailable('Root reference lost its finite complete mouth')
        curves=layout.boundaries
        for curve in curves:
            coordinates=getattr(curve,'source_coordinates',None)
            if (type(curve) is not PolylineCurb or type(coordinates) is not CurbCoordinates
                    or coordinates.path is not approach.path or coordinates.approach is not approach):
                raise OriginalRootPortUnavailable('Root reference lost its exact native point ledger')
            coordinates.validate(curve.points,context)
        at=curves[0].source_coordinates.stations
        if (at[0]!=0. or at[-1]!=station or np.any(np.diff(at)<=0.)
                or not np.array_equal(curves[1].source_coordinates.stations,at[::-1])
                or not np.array_equal(curves[2].source_coordinates.stations,np.zeros(2))):
            raise OriginalRootPortUnavailable('Root reference changed its original bank intervals')
        radial=layout.radial;u=np.linspace(1.,0.,2*radial+1)
        mouth=approach.point(station,approach.offset+approach.width*(u-.5))
        points=np.vstack((curves[0].points[:-1],mouth[:-1],curves[1].points[:-1],curves[2].points[:-1]))
        if not np.array_equal(points,source.region.xy):
            raise OriginalRootPortUnavailable('Root reference is not the declared native port polygon')
        expected=[]
        for index,curve in enumerate(curves):
            if index==1:expected.extend(RimSegment('mouth',0,float(x),float(y)) for x,y in zip(u,u[1:]))
            expected.extend(RimSegment('internal' if index==2 else 'curb',index,float(x/curve.length),float(y/curve.length))
                for x,y in zip(curve.stations,curve.stations[1:]))
        if tuple(expected)!=source.segments:
            raise OriginalRootPortUnavailable('Root reference changed its actual segment provenance')
        return cls(source,approach,float(station))

    @property
    def endpoint(self):return End(self.approach.source.id,self.approach.end)


class NativeRolePort(OriginalRootPort):
    """Original port with its actual full or Bundle half profile at this End."""

    @classmethod
    def build(cls, context, source, approach):
        from .ordinary_mouth_clearance import _authenticate_native_half_profiles
        if not _authenticate_native_half_profiles(context, SimpleNamespace(approaches=(approach,))):
            raise OriginalRootPortUnavailable('Native port needs its original oriented carriageway profile')
        return cls._from_native_geometry(context, source, approach)


def native_bank_curvature(approach,station,side,*,reverse=False):
    """Actual Smooth offset-bank curvature, with an explicit Route orientation.

    Exact banks are affine on each declared source interval. A true miter or
    Route seam needs its own one-sided provider and is not guessed here.
    """
    from .routes import Route
    from .intersections import Approach
    from .mesh import cross
    if (type(approach) is not Approach or type(approach.path) is not Route or side not in (-1,1)
            or not math.isfinite(station) or not 0.<station<approach.length):
        raise OriginalRootPortUnavailable('Native root contact needs a finite original Route bank')
    route=approach.path
    if np.any(route.stations[1:-1]==station):
        raise OriginalRootPortUnavailable('Native root curvature cannot invent a Route-seam jet')
    index=int(np.searchsorted(route.stations,station,side='right')-1)
    piece=route.pieces[index];local=station-float(route.stations[index]);path=piece.path
    native=local if piece.forward else path.length-local
    parameter=path.parameter(native)
    if np.any(path.knots[1:-1]==parameter):
        raise OriginalRootPortUnavailable('Native root curvature needs a declared one-sided control-knot jet')
    if path.exact:return 0.
    _,velocity,acceleration=path.jet(parameter)
    k=float(cross(velocity,acceleration)/np.linalg.norm(velocity)**3)*(1. if piece.forward else -1.)
    distance=.5*side*approach.width
    denominator=1.+distance*k
    if not math.isfinite(denominator) or denominator<=0.:
        raise OriginalRootPortUnavailable('Native root bank has no regular offset curvature')
    return k/denominator*(-1. if reverse else 1.)


def refine_contact_port(context,source,contact,*,maximum_samples=4096,root=None):
    """Refine only the owned side before a true contact, inside the old budget.

    These proof knots lie inside the accepted root/circle patches. The exposed
    contact-to-mouth banks retain their existing approximation and mesh cuts.
    """
    a=source.layout.approaches[0]
    reference=OriginalRootPort.build(context,source,a)
    return refine_native_contact(context,reference,contact,maximum_samples=maximum_samples,root=root)


def refine_native_contact(context,reference,contact,*,maximum_samples=4096,root=None):
    """Subdivide a previously authenticated port at its physical contact."""
    from .compound_cycles import PolylineCurb
    source,a=reference.footprint,reference.approach
    if (isinstance(contact,(bool,np.bool_)) or not isinstance(contact,(int,float,np.integer,np.floating))
            or not math.isfinite(contact) or not 0.<contact<=reference.station
            or type(maximum_samples) is not int or not 3<=maximum_samples<=4096):
        raise OriginalRootPortUnavailable('Root refinement needs its actual finite original contact')
    old=source.layout.boundaries[0].source_coordinates.stations
    low,high=max(0.,contact-.5),contact
    tolerance=.1*RegionUnion().epsilon
    additions=[];bounds=[];allocated=0
    for i,piece in enumerate(a.path.pieces):
        begin,finish=float(a.path.stations[i]),float(a.path.stations[i+1])
        start,end=max(low,begin),min(high,finish)
        if end<=start:continue
        first,last=((start-begin,end-begin) if piece.forward else
                    (piece.path.length-(end-begin),piece.path.length-(start-begin)))
        remaining=maximum_samples-allocated
        if remaining<3:raise OriginalRootPortUnavailable('Root contact sampling exhausts its existing finite allocation')
        values=piece.path.samples(first,last,min(.5,context.network.settings.spacing),
            tolerance=.99*tolerance,max_samples=remaining,local_curvature=True)
        allocated+=len(values)
        additions.extend(begin+(values if piece.forward else piece.path.length-values))
        bounds.append((piece,begin,first,last))
    if not additions:raise OriginalRootPortUnavailable('Root contact has no positive original refinement interval')
    at=np.unique(np.r_[old,additions,contact])
    replaced=()
    if root is not None:
        at,replaced=root.contact_grid(context,reference,at,float(contact))
    if len(at)>20001:raise OriginalRootPortUnavailable('Root native point ledger exceeds its original sample allocation')
    right=a.point(at,a.offset+.5*a.width);left=a.point(at,a.offset-.5*a.width)[::-1]
    radial=source.layout.radial;u=np.linspace(1.,0.,2*radial+1)
    mouth=a.point(reference.station,a.offset+a.width*(u-.5));cap=source.layout.boundaries[2].points
    curves=tuple(PolylineCurb(points) for points in (right,left,cap))
    declarations=(
        CurbCoordinates(a.path,at,np.full(len(at),a.offset+.5*a.width),np.ones(len(at)),approach=a,
            maximum_chord_error=source.layout.boundaries[0].source_coordinates.maximum_chord_error,evaluation_grid=at),
        CurbCoordinates(a.path,at[::-1],np.full(len(at),a.offset-.5*a.width),-np.ones(len(at)),approach=a,
            maximum_chord_error=source.layout.boundaries[1].source_coordinates.maximum_chord_error,evaluation_grid=at),
        source.layout.boundaries[2].source_coordinates)
    for curve,coordinates in zip(curves,declarations):
        curve.source_coordinates=coordinates;coordinates.validate(curve.points,context)
    points=np.vstack((right[:-1],mouth[:-1],left[:-1],cap[:-1]));labels=[]
    for index,curve in enumerate(curves):
        if index==1:labels.extend(RimSegment('mouth',0,float(x),float(y)) for x,y in zip(u,u[1:]))
        labels.extend(RimSegment('internal' if index==2 else 'curb',index,float(x/curve.length),float(y/curve.length))
            for x,y in zip(curve.stations,curve.stations[1:]))
    layout=SimpleNamespace(approaches=(a,),outward_stations=(reference.station,),boundaries=curves,
        rows=len(points),radial=radial,mouths={},corners=(),
        original_root_caps={len(labels)-1:reference.endpoint})
    result=SampledFootprint(Region(source.region.owner,points),tuple(labels),layout)
    if replaced and not set(map(tuple,source.region.xy)).issubset(map(tuple,points)):
        raise OriginalRootPortUnavailable('Contact sampling lost an original port coordinate')
    maximum=0.
    for piece,begin,first,last in bounds:
        values=at[(at>=begin)&(at<=begin+piece.path.length)]-begin
        values=np.sort(values if piece.forward else piece.path.length-values)
        values=values[(values>=first)&(values<=last)]
        k=piece.path.visible_curvature_bound(first,last) if not piece.path.exact else 0.
        reach=abs(piece.section.offset)+.5*piece.section.width
        if reach*k>=1.:raise OriginalRootPortUnavailable('Refined native bank lacks a regular strip certificate')
        bound=(k/(1.-reach*k))*((1.+reach*k)*np.diff(values).max())**2/8.
        bound+=256.*np.finfo(float).eps*max(1.,float(np.abs(points).max()))
        maximum=max(maximum,float(bound))
    if maximum>tolerance:
        raise OriginalRootPortUnavailable('Root contact refinement failed its unchanged analytical error target')
    return result,dict(section=a.source.id,end=a.end,contact=float(contact),owned_interval=(low,high),
        old_knots=len(old),new_knots=len(at),added_knots=len(at)-len(old),
        maximum_local_continuous_chord_error=maximum,maximum_local_requested_error=tolerance,
        intended_original_curves_unchanged=True,original_mouth_plane_unchanged=True,
        replaced_generated_stations=replaced)
