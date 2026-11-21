"""Exact source parameters retained beside prepared curb polylines.

The arrays describe the existing vertices; they never resample or move them.
Route seams are explicit parameter slices of an existing chord, not invented
source nodes or nearest-source matches.
"""
from dataclasses import dataclass
import numpy as np
from .model import Key


SOURCE_CURB_COORDINATE_POLICY='6-explicit-owned-prefix-physical-native-ledgers'


class CurbCoordinatesUnavailable(ValueError):
    """A derived approach needs its own original-source parameter contract."""


@dataclass(frozen=True)
class OriginalCurbPiece:
    edge: int
    first_fraction: float
    last_fraction: float
    path: object
    first_station: float
    last_station: float
    first_offset: float
    last_offset: float
    authored_side: int | None

    @property
    def section(self):return self.path.section


@dataclass(frozen=True)
class _VertexEvaluation:
    path: object
    approach: object
    stations: object
    offsets: object
    sides: object
    points: object


def _readonly(value):
    result=np.array(value,copy=True);result.flags.writeable=False
    return result


class CurbCoordinates:
    """One exact station/offset pair for each existing physical curve vertex."""
    def __init__(self,path,stations,offsets,sides,*,approach=None,maximum_chord_error=None,evaluation_grid=None):
        from .paths import Path
        from .routes import Route
        from .intersections import Approach
        from .owned_prefixes import OwnedPrefixRoute
        route=type(path) in (Route,OwnedPrefixRoute)
        if (type(path) not in (Path,Route,OwnedPrefixRoute) or route and approach is None
                or approach is not None and (type(approach) is not Approach
                    or approach.path is not path or approach.end not in (0,1))):
            raise CurbCoordinatesUnavailable('Prepared source curb needs an explicit original Path or Route')
        if type(path) is OwnedPrefixRoute:
            from .owned_prefix_source_curbs import validate_declaration
            validate_declaration(path,approach)
        stations,offsets,sides=(np.asarray(v) for v in (stations,offsets,sides))
        if stations.ndim==1 and len(stations)>20001:
            raise ValueError('Curb source parameter ledger exceeds its finite vertex budget')
        if (stations.ndim!=1 or len(stations)<2
                or offsets.shape!=stations.shape or sides.shape!=stations.shape
                or not np.isfinite(stations).all() or not np.isfinite(offsets).all()
                or not np.isin(sides,(-1,1)).all()
                or np.any(stations<0.) or np.any(stations>path.length)):
            raise ValueError('Prepared source curb needs finite original vertex parameters')
        center=path.section.offset if approach is None else approach.offset
        width=path.section.width if approach is None else approach.width
        if not np.array_equal(offsets,center+.5*width*sides):
            raise ValueError('Curb bank side does not match its declared source offset')
        if (maximum_chord_error is not None and
                (isinstance(maximum_chord_error,(bool,np.bool_))
                 or not isinstance(maximum_chord_error,(int,float,np.integer,np.floating))
                 or not np.isfinite(maximum_chord_error) or maximum_chord_error<=0.)):
            raise ValueError('Curb source sampling needs its positive finite chord-error bound')
        self.path,self.approach=path,approach
        self._maximum_chord_error=(None if maximum_chord_error is None else float(maximum_chord_error))
        self.stations=stations.astype(float,copy=True)
        self.offsets=offsets.astype(float,copy=True)
        self.sides=sides.astype(np.int8,copy=True)
        for array in (self.stations,self.offsets,self.sides):array.flags.writeable=False
        if evaluation_grid is not None:
            grid=np.asarray(evaluation_grid,dtype=float)
            if (grid.ndim!=1 or not 2<=len(grid)<=20001 or not np.isfinite(grid).all()
                    or grid[0]<0. or grid[-1]>path.length or np.any(np.diff(grid)<=0.)):
                raise ValueError('Curb source evaluation needs its finite original monotone batch')
            indices=np.searchsorted(grid,self.stations)
            if np.any(indices>=len(grid)) or not np.array_equal(grid[indices],self.stations):
                raise ValueError('Curb source parameters do not belong to their declared evaluation batch')
            evaluate=path if approach is None else approach
            # The builder evaluates both entire banks in the same ascending
            # batch, then takes/reverses existing vertices. A shorter SIMD
            # batch may differ by one ulp; retain that exact declared source
            # evaluation instead of reevaluating a subset with another shape.
            banks=np.stack((evaluate.point(grid,center+.5*width),
                            evaluate.point(grid,center-.5*width)))
            points=banks[(self.sides<0).astype(int),indices]
            self._vertex_evaluation=_VertexEvaluation(path,approach,_readonly(self.stations),
                _readonly(self.offsets),_readonly(self.sides),_readonly(points))

    @property
    def maximum_chord_error(self):
        return getattr(self,'_maximum_chord_error',None)

    def subset(self,indices):
        indices=np.asarray(indices)
        if (indices.ndim!=1 or indices.dtype.kind not in 'iu' or len(indices)<2
                or np.any(indices<0) or np.any(indices>=len(self.stations))):
            raise ValueError('Curb source subset needs explicit retained vertex indices')
        result=type(self)(self.path,self.stations[indices],self.offsets[indices],
                          self.sides[indices],approach=self.approach,
                          maximum_chord_error=self.maximum_chord_error)
        reference=getattr(self,'_vertex_evaluation',None)
        if reference is not None:
            result._vertex_evaluation=_VertexEvaluation(reference.path,reference.approach,
                _readonly(reference.stations[indices]),_readonly(reference.offsets[indices]),
                _readonly(reference.sides[indices]),_readonly(reference.points[indices]))
        return result

    def validate(self,points,context=None):
        """Check the declaration directly; no bank search or inverse fitting."""
        from .routes import Route
        from .owned_prefixes import OwnedPrefixRoute
        if type(self.path) is OwnedPrefixRoute:
            from .owned_prefix_source_curbs import validate_declaration
            validate_declaration(self.path,self.approach,context)
        points=np.asarray(points)
        if points.shape!=(len(self.stations),2):
            raise ValueError('Curb source parameters do not match its physical vertices')
        reference=getattr(self,'_vertex_evaluation',None)
        if reference is not None:
            if (self.path is not reference.path or self.approach is not reference.approach
                    or not np.array_equal(self.stations,reference.stations)
                    or not np.array_equal(self.offsets,reference.offsets)
                    or not np.array_equal(self.sides,reference.sides)):
                raise ValueError('Curb source parameters changed their declared source evaluation')
            expected=reference.points
        else:
            expected=(self.path.point(self.stations,self.offsets) if self.approach is None
                      else self.approach.point(self.stations,self.offsets))
        if not np.array_equal(points,expected):
            raise ValueError('Curb source parameters changed their original physical vertices')
        paths=(tuple(p.path for p in self.path.pieces) if type(self.path) in (Route,OwnedPrefixRoute) else (self.path,))
        if context is not None and any(context[Key('path',p.section.id)] is not p for p in paths):
            raise ValueError('Curb source parameters changed original Path ownership')
        return True

    def edge_sources(self,edge):
        """Split parameter ancestry at real Route seams without adding vertices."""
        from .routes import Route
        from .owned_prefixes import OwnedPrefixRoute
        if (isinstance(edge,(bool,np.bool_)) or not isinstance(edge,(int,np.integer))
                or not 0<=edge<len(self.stations)-1):
            raise ValueError('Curb source edge is outside its explicit vertex map')
        first,last=map(float,self.stations[edge:edge+2])
        offsets=self.offsets[edge:edge+2]
        side=int(self.sides[edge]) if self.sides[edge]==self.sides[edge+1] else None
        if type(self.path) not in (Route,OwnedPrefixRoute):
            direction=-1 if self.approach is not None and self.approach.end==1 else 1
            native=np.array((first,last)) if direction==1 else self.path.length-np.array((first,last))
            return (OriginalCurbPiece(int(edge),0.,1.,self.path,float(native[0]),float(native[1]),
                                     direction*float(offsets[0]),direction*float(offsets[1]),
                                     None if side is None else side*direction),)
        low,high=sorted((first,last))
        begin=np.searchsorted(self.path.stations,low,side='right')
        end=np.searchsorted(self.path.stations,high,side='left')
        seams=self.path.stations[begin:end]
        fractions=(np.r_[0.,np.sort((seams-first)/(last-first)),1.]
                   if first!=last else np.array((0.,1.)))
        result=[]
        for a,b in zip(fractions,fractions[1:]):
            at=first+.5*(a+b)*(last-first)
            index=int(np.clip(np.searchsorted(self.path.stations,at,side='right')-1,
                              0,len(self.path.pieces)-1))
            piece=self.path.pieces[index];base=self.path.stations[index]
            local=np.clip(np.array((first+a*(last-first),first+b*(last-first)))-base,
                          0.,piece.path.length)
            direction=1 if piece.forward else -1
            native=local if piece.forward else piece.path.length-local
            lateral=(offsets[0]+np.array((a,b))*(offsets[1]-offsets[0]))*direction
            result.append(OriginalCurbPiece(int(edge),float(a),float(b),piece.path,
                float(native[0]),float(native[1]),float(lateral[0]),float(lateral[1]),
                None if side is None else side*direction))
        return tuple(result)


class NativeBankCurb:
    """One prepared bank evaluated on its declared original analytical source.

    This is used only when a finite exposed contact identifies that very bank
    through another analytical wrapper. Every old physical vertex and its
    parameter remains compulsory. New chords keep an equal or smaller certified
    error bound; ordinary prepared polylines are otherwise left untouched.
    """
    continuation=None

    def __init__(self,curve):
        from .compound_cycles import PolylineCurb
        from .sampling import CertifiedCurbSampler
        if type(curve) is not PolylineCurb:
            raise ValueError('Native bank parameterization needs an original prepared polyline')
        coordinates=getattr(curve,'source_coordinates',None)
        if type(coordinates) is not CurbCoordinates or coordinates.approach is None:
            raise CurbCoordinatesUnavailable('Native bank needs its explicit original Approach ledger')
        if coordinates.maximum_chord_error is None:
            raise CurbCoordinatesUnavailable('Native bank needs its original declared source sampling bound')
        coordinates.validate(curve.points)
        delta=np.diff(coordinates.stations)
        if (len(set(coordinates.sides))!=1 or len(set(coordinates.offsets))!=1
                or not (np.all(delta>0.) or np.all(delta<0.))):
            raise CurbCoordinatesUnavailable('Native bank needs one monotone original bank interval')
        self.original=curve
        self.source_coordinates=coordinates
        self.approach=coordinates.approach
        self.offset=float(coordinates.offsets[0])
        self.declared_chord_error_bound=coordinates.maximum_chord_error
        self.parameters=coordinates.stations
        self.cumulative=curve.stations
        self.points,self.stations,self.length=curve.points,curve.stations,curve.length
        self.corner,self.contact_reference_points=curve.corner,curve.points
        metrics=CertifiedCurbSampler(1.,1.,20000).lead_metrics(self.approach,
            self.parameters,self.cumulative,self.offset,0.,self.length)
        self.curvature_bound,self.maximum_speed_bound,_=metrics
        self.original_chord_error_bound=(self.curvature_bound*
            (self.maximum_speed_bound*float(np.max(np.diff(self.cumulative))))**2/8.)
        # Preserve exact original vertex parameters after normalization. Source
        # Route seams and Exact knots are semantic cuts too; dividing a curve
        # distance and multiplying it back must not move such a cut by an ulp.
        knots={float(t):float(s) for t,s in zip(self.stations/self.length,self.parameters)}
        route=self.approach.path
        from .paths import Path
        if type(route) is Path:
            pieces=((route,self.approach.end==0,0.),)
        else:
            pieces=tuple((piece.path,piece.forward,float(route.stations[i]))
                         for i,piece in enumerate(route.pieces))
        xp,fp=(self.parameters,self.stations/self.length)
        if xp[-1]<xp[0]:xp,fp=xp[::-1],fp[::-1]
        for path,forward,base in pieces:
            native=path.stations if path.exact else np.array((0.,path.length))
            for station in native:
                outward=base+(station if forward else path.length-station)
                if not xp[0]<=outward<=xp[-1]:continue
                t=float(np.interp(outward,xp,fp))
                if t in knots and knots[t]!=outward:
                    raise CurbCoordinatesUnavailable('Distinct native semantic knots collapse in the prepared parameter')
                knots[t]=float(outward)
        self.evaluation_knots=np.array(sorted(knots))
        self.evaluation_stations=np.array([knots[t] for t in self.evaluation_knots])
        self.evaluation_knots.flags.writeable=False
        self.evaluation_stations.flags.writeable=False

    def outward_station(self,t):
        return np.interp(np.asarray(t),self.evaluation_knots,self.evaluation_stations)

    def __call__(self,t):
        return self.approach.point(self.outward_station(t),self.offset)

    def certified_parameters(self,sampler,first,last):
        from .compound_curbs import OffsetCurb
        from .sampling import CertifiedCurbSampler
        tolerance=min(sampler.tolerance,self.declared_chord_error_bound)
        if self.original_chord_error_bound>0.:
            tolerance=min(tolerance,self.original_chord_error_bound)
        bounded=CertifiedCurbSampler(tolerance,sampler.spacing,sampler.max_samples)
        sampled=OffsetCurb.sample(self,bounded,first,last)
        knots=self.evaluation_knots
        values=np.unique(np.r_[sampled,knots[(knots>first)&(knots<last)]])
        if len(values)>sampler.max_samples:
            raise ValueError('Native bank sampling exceeds its unchanged sample budget')
        return values


def original_bank_intervals(curve,lower,upper,initial,*,prepared=False):
    """Finite native parameter maps for explicit banks and composite offsets.

    Prepared chord maps are available only to the eligibility check. They must
    not make a chord and a curved bank into an alleged exact curve intersection.
    The caller first replaces that eligible chord bank with NativeBankCurb.
    """
    from .compound_curbs import CompositeCurb,OffsetCurb,SamplingBudgetCurb
    from .compound_cycles import PolylineCurb
    from .intersections import Approach
    from .paths import Path
    from .routes import Route
    if not (np.isfinite((lower,upper,initial)).all() and 0.<=lower<=upper<=1.):return ()
    depth=0
    while isinstance(curve,SamplingBudgetCurb):
        depth+=1
        if depth>16:raise ValueError('Native bank exceeds its finite wrapper depth')
        curve=curve.curve
    result=[]

    def append(path,offset,s0,s1,t0,t1):
        if t1<t0:s0,s1,t0,t1=s1,s0,t1,t0
        low,high=max(lower,t0),min(upper,t1)
        if high<=low or s0==s1:return
        a=s0+(low-t0)/(t1-t0)*(s1-s0)
        b=s0+(high-t0)/(t1-t0)*(s1-s0)
        seed=s0+(float(np.clip(initial,low,high))-t0)/(t1-t0)*(s1-s0)
        parameters=np.array((s0,s1));cumulative=np.array((t0,t1))*curve.length
        if parameters[1]<parameters[0]:parameters,cumulative=parameters[::-1],cumulative[::-1]
        result.append(dict(path=path,offset=float(offset),low=min(a,b),high=max(a,b),
            seed=float(seed),forward=True,base=0.,length=path.length,
            parameters=parameters,cumulative=cumulative,curve_base=0.,curve_length=curve.length,
            explicit_bank=True))

    if isinstance(curve,NativeBankCurb) or prepared and type(curve) is PolylineCurb:
        coordinates=getattr(curve,'source_coordinates',None)
        if type(coordinates) is not CurbCoordinates or coordinates.approach is None:return ()
        if (len(set(coordinates.sides))!=1 or len(set(coordinates.offsets))!=1
                or not (np.all(np.diff(coordinates.stations)>0.)
                        or np.all(np.diff(coordinates.stations)<0.))):return ()
        coordinates.validate(curve.points)
        for edge in range(len(coordinates.stations)-1):
            a,b=curve.stations[edge:edge+2]/curve.length
            if b<lower or a>upper:continue
            for piece in coordinates.edge_sources(edge):
                if piece.authored_side is not None and piece.first_offset==piece.last_offset:
                    append(piece.path,piece.first_offset,piece.first_station,piece.last_station,
                        a+piece.first_fraction*(b-a),a+piece.last_fraction*(b-a))
    elif isinstance(curve,CompositeCurb):
        for i,piece in enumerate(curve.pieces):
            if not isinstance(piece,OffsetCurb):continue
            approach=piece.approach
            if (type(approach) is not Approach or type(approach.path) not in (Path,Route)
                    or piece.offset not in (approach.offset-.5*approach.width,
                                             approach.offset+.5*approach.width)):continue
            route=approach.path
            if type(route) is Route:
                paths=tuple((p.path,bool(p.forward),float(route.stations[j]))
                            for j,p in enumerate(route.pieces))
                seams=route.stations
            else:
                paths=((route,approach.end==0,0.),);seams=np.array((0.,route.length))
            for j,(a,b) in enumerate(zip(piece.parameters,piece.parameters[1:])):
                if a==b:continue
                ca,cb=piece.cumulative[j:j+2]
                ta,tb=(curve.lengths[i]+np.array((ca,cb)))/curve.length
                if tb<lower or ta>upper:continue
                fractions=sorted({0.,1.,*((float(s)-a)/(b-a) for s in seams if min(a,b)<s<max(a,b))})
                for u,v in zip(fractions,fractions[1:]):
                    first,last=a+u*(b-a),a+v*(b-a)
                    at=int(np.clip(np.searchsorted(seams,.5*(first+last),side='right')-1,0,len(paths)-1))
                    path,forward,base=paths[at]
                    s0,s1=first-base,last-base
                    if not forward:s0,s1=path.length-s0,path.length-s1
                    append(path,(1 if forward else -1)*piece.offset,s0,s1,
                        ta+u*(tb-ta),ta+v*(tb-ta))
    return tuple(result)


def canonical_native_banks(group):
    """Upgrade only an exposed contact with a proven coincident native bank."""
    from copy import copy
    from dataclasses import replace
    from .compound_cycles import PolylineCurb
    from .compounds import rim_label
    sources={f.region.owner:f for f in group.footprints};eligible={};records=[]
    for loop in group.loops:
        for index,last in enumerate(loop.spans):
            first=loop.spans[index-1]
            if first.owner==last.owner:continue
            labels=tuple(rim_label(span,sources) for span in (first,last))
            if any(label.kind!='curb' for label in labels):continue
            owners=(first.owner,last.owner)
            curves=tuple(sources[owner].layout.boundaries[label.index] for owner,label in zip(owners,labels))
            initial=(labels[0].last,labels[1].first)
            reach=tuple(2./sources[owner].layout.rows for owner in owners)
            bounds=tuple((max(0.,seed-r),min(1.,seed+r)) for seed,r in zip(initial,reach))
            for side in (0,1):
                curve=curves[side]
                if type(curve) is not PolylineCurb or (owners[side],labels[side].index) in eligible:continue
                coordinates=getattr(curve,'source_coordinates',None)
                if (type(coordinates) is not CurbCoordinates
                        or coordinates.maximum_chord_error is None):continue
                own=original_bank_intervals(curve,*bounds[side],initial[side],prepared=True)
                other=original_bank_intervals(curves[1-side],*bounds[1-side],initial[1-side])
                matches=[]
                for a in own:
                    for b in other:
                        if (a['path'] is b['path'] and a['offset']==b['offset']
                                and min(a['high'],b['high'])>max(a['low'],b['low'])):
                            matches.append(dict(section=a['path'].section.id,offset=a['offset'],
                                interval=(max(a['low'],b['low']),min(a['high'],b['high']))))
                if matches:
                    eligible[owners[side],labels[side].index]=NativeBankCurb(curve)
                    records.append(dict(owner=owners[side].text(),index=labels[side].index,
                        original_sources=matches,
                        declared_chord_error_bound=coordinates.maximum_chord_error,
                        proof='EXPOSED_CONTACT_IDENTICAL_NATIVE_PATH_BANK_AND_FINITE_DOMAIN'))
    if not eligible:return group,()
    footprints=[]
    for source in group.footprints:
        curves=tuple(eligible.get((source.region.owner,i),curve)
                     for i,curve in enumerate(source.layout.boundaries))
        if any(a is not b for a,b in zip(curves,source.layout.boundaries)):
            layout=copy(source.layout);layout.boundaries=curves
            source=replace(source,layout=layout)
        footprints.append(source)
    physical=group.physical_layout
    if physical is not None:
        physical=copy(physical);physical.native_bank_parameterization=tuple(records)
    return replace(group,footprints=tuple(footprints),physical_layout=physical),tuple(records)
