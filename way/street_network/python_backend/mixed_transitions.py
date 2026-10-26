"""Bounded polynomial transitions retaining Exact and Smooth centerlines.

Curved dimensional seams use the existing orientation, sampling, profile and
physical mesh guards. Interior transverse frames are polynomial; complete outer
profiles and both source centerline curves remain authoritative.
"""
import math
import numpy as np
from .exact_transitions import PiecewiseCurb,PolynomialCurb,restricted
from .paths import right
from .transitions import SectionTransition,TransitionClearance,_certify_orientation
from .transition_fitting import ExactRowSampler,ExactTransitionFitter
from .intersections import half_intervals
from .entrance_slopes import EntranceSlopeGuard
from .mixed_heights import MixedHeightCorrection


def elevate(values,degree):
    values=np.array(values,dtype=float,copy=True)
    while len(values)<=degree:
        n=len(values)
        result=np.empty((n+1,*values.shape[1:]))
        result[0],result[-1]=values[0],values[-1]
        t=np.arange(1,n)/n
        t=t.reshape((n-1,)+(1,)*(values.ndim-1))
        result[1:-1]=t*values[:-1]+(1.-t)*values[1:]
        values=result
    return values


def product(vector,scalar):
    a,b=len(vector)-1,len(scalar)-1
    result=np.zeros((a+b+1,2))
    for i,v in enumerate(vector):
        for j,s in enumerate(scalar):
            result[i+j]+=math.comb(a,i)*math.comb(b,j)/math.comb(a+b,i+j)*v*s
    return result


class MixedTransitionGeometry:
    max_pieces=4096

    @staticmethod
    def leg(approach,distance,side):
        path=approach.path
        if getattr(path,'outward',False) or not hasattr(path,'knots'):
            raise ValueError('Mixed transition requires local Section paths')
        lo,hi=(0.,distance) if approach.end==0 else (path.length-distance,path.length)
        p0,p1=map(float,path.parameter(np.array((lo,hi))))
        cuts=np.r_[p0,path.knots[(path.knots>p0+1e-10)&(path.knots<p1-1e-10)],p1]
        pieces=[]
        for first,last in zip(cuts,cuts[1:]):
            index=int(np.clip(np.searchsorted(path.knots,.5*(first+last),side='right')-1,0,len(path.chords)-1))
            chord=float(path.chords[index])
            a,b=(first-path.knots[index])/chord,(last-path.knots[index])/chord
            stations=path.station(np.array((first,last)))
            normal=path.frame(stations)[2]
            if path.exact:
                controls=np.array((path.points[index],path.points[index+1]))
                normals=normal.copy()
            else:
                p,v,c,d=path.coefficients[index]
                controls=np.array((p,p+v/3.,p+(2.*v+c)/3.,p+v+c+d))
                _,velocity,acceleration=path.jet(np.array((first,last)))
                speed=np.linalg.norm(velocity,axis=1)
                derivative=right(acceleration/speed[:,None]-velocity*
                    (np.sum(velocity*acceleration,axis=1)/speed**3)[:,None])
                normals=np.array((normal[0],normal[0]+(last-first)*derivative[0]/3.,
                                  normal[1]-(last-first)*derivative[1]/3.,normal[1]))
            controls=restricted(controls,float(a),float(b))
            pieces.append([controls,normals,float(stations[1]-stations[0]),
                           (int(approach.source.id),float(first),float(last))])
        forward=(approach.end==0)==(side==1)
        if not forward:
            pieces=[[c[::-1].copy(),-n[::-1].copy(),length,(identity,last,first)]
                    for c,n,length,(identity,first,last) in pieces[::-1]]
        if not pieces or len(pieces)>MixedTransitionGeometry.max_pieces:
            raise ValueError('Mixed transition exceeds its source-piece budget')
        return pieces

    @classmethod
    def build(cls,approaches,distances):
        first,last=approaches
        weld=float(np.linalg.norm(first.origin-last.origin))
        if weld>1e-5:raise ValueError('Mixed seam endpoint copies disagree')
        directions=(-first.frame(0.)[1],last.frame(0.)[1])
        dot=float(directions[0]@directions[1])
        if dot<-.95:raise ValueError('Mixed transition has a near-reversal miter')
        shared=(right(directions[0])+right(directions[1]))/(1.+dot)
        legs=[cls.leg(a,d,i) for i,(a,d) in enumerate(zip(approaches,distances))]
        center=.5*(first.origin+last.origin)
        legs[0][-1][0][-1]=center;legs[1][0][0][0]=center
        for controls,at in ((legs[0][-1][1],-1),(legs[1][0][1],0)):
            delta=shared-controls[at]
            controls[at]=shared
            if len(controls)>2:controls[-2 if at==-1 else 1]+=delta
        pieces=legs[0]+legs[1]
        if len(pieces)>cls.max_pieces:raise ValueError('Mixed transition exceeds its source-piece budget')
        reach=max(float(np.linalg.norm(n,axis=1).max()) for _,n,_,_ in pieces)
        if reach>5.:raise ValueError('Mixed transition exceeds the supported miter reach')
        lengths=np.array([p[2] for p in pieces])
        knots=np.r_[0.,np.cumsum(lengths)]/float(lengths.sum())
        knots[-1]=1.
        curbs=[]
        for q in (-1.,1.):
            a,b=-first.offset+q*.5*first.width,last.offset+q*.5*last.width
            widths=np.array((a,a,a,b,b,b))
            segments=[]
            for (centerline,normals,_,_),lo,hi in zip(pieces,knots,knots[1:]):
                offsets=restricted(widths,float(lo),float(hi))
                transverse=product(normals,offsets)
                controls=elevate(centerline,len(transverse)-1)+transverse
                segments.append(PolynomialCurb(controls))
            curbs.append(PiecewiseCurb(knots,segments))
        return tuple(curbs),dict(mixed_miter=True,source_centerlines_preserved=True,
            polynomial_pieces=len(pieces),mandatory_knots=knots.tolist(),maximum_miter_reach=reach,
            source_parameter_pieces=[p[3] for p in pieces],
            endpoint_reconciliation=weld,interpolation_modes=[a.source.interpolation for a in approaches])

    @staticmethod
    def parameter_at_station(plan, side, station):
        """Transport a metric source station into its actual polynomial piece.

        The chart allocates intervals by each piece's arc length, but the
        polynomial within an interval uses the original native parameter.
        Linear station fractions therefore do not locate a Smooth source.
        """
        approach=plan.approaches[side]
        native=float(approach.path.parameter(np.asarray(np.clip(station,0.,approach.length))))
        pieces=plan.certification['source_parameter_pieces']
        for index,(identity,first,last) in enumerate(pieces):
            if identity==approach.source.id and min(first,last)-1e-9<=native<=max(first,last)+1e-9:
                fraction=float(np.clip((native-first)/(last-first),0.,1.))
                lo,hi=plan.curbs[0].knots[index:index+2]
                return float(lo+(hi-lo)*fraction)
        raise ValueError('Furniture station has no owned mixed-transition source piece')


class MixedTransitionFitter:
    max_candidates=64

    @staticmethod
    def certify(context,result,clearance):
        result=SectionTransition.certify(context,result,clearance)
        mesh=result.mesh(context.patterns,exaggeration=1.)
        maximum=0.
        for side,approach in enumerate(result.approaches):
            ids=mesh.seams[side]
            q=np.linspace(-1.,1.,result.across+1)*(-1. if side==0 else 1.)
            velocities=np.zeros((len(mesh.xyz),2))
            velocities[ids]=approach.jet(result.distances[side],approach.offset+.5*approach.width*q)[1]
            maximum=max(maximum,EntranceSlopeGuard().validate(
                mesh,ids,velocities,mesh.xyz[:,2]-result.datum))
        result.certification['maximum_source_entrance_longitudinal_bound']=maximum
        return result

    @classmethod
    def solve(cls,context,site,approaches,clearance=None):
        first,last=approaches
        if abs(first.profile.datum-last.profile.datum)>1e-8:
            raise ValueError('Mixed transition requires a common minimum height')
        differences=np.array([first.profile.value(-q)-last.profile.value(q) for q in (-1.,0.,1.)])
        a=.5*(differences[0]+differences[2])-differences[1]
        b=.5*(differences[2]-differences[0])
        delta=float(abs(differences).max())
        if abs(a)>1e-15:
            vertex=-b/(2.*a)
            if abs(vertex)<1.:delta=max(delta,abs(float(first.profile.value(-vertex)-last.profile.value(vertex))))
        desired=max(2.,4.*abs(first.width-last.width)+2.*abs(first.offset+last.offset),1.875*delta/.03)
        clearance=clearance if clearance is not None else TransitionClearance.build(context,approaches)
        extended,free=ExactTransitionFitter.free_limits(context.network,approaches,clearance)
        across=2*half_intervals(context.network.settings)
        seen=set();failures=[];attempts=0;height_candidates=[]
        for cap in (clearance.limits,extended):
            effective=TransitionClearance(cap,clearance.constraints)
            for factor in (1.,1.5,2.,3.,4.,6.,8.):
                distances=np.minimum(.5*desired*factor,cap)
                identity=tuple(distances)
                if identity in seen:continue
                seen.add(identity)
                try:
                    curbs,metadata=MixedTransitionGeometry.build(approaches,distances)
                    orientation=_certify_orientation(*curbs)
                    sampler=ExactRowSampler(curbs,context.network.settings.curve_tolerance,50000//across)
                    for parameters,bound in sampler.candidates(context.network.settings.resolution):
                        attempts+=1
                        if attempts>cls.max_candidates:raise ValueError('Mixed transition exhausted its candidate budget')
                        certificate=dict(orientation,**metadata,certified_curve_error=bound,
                            sampling_policy='mixed_source_polynomials',
                            extended_free_ends=[r for r in free if distances[r['approach']]>clearance.limits[r['approach']]])
                        result=SectionTransition(site.key,tuple(approaches),tuple(map(float,distances)),curbs,
                            len(parameters)-1,across,parameters,certificate)
                        try:return cls.certify(context,result,effective)
                        except ValueError as exc:
                            failures.append(str(exc))
                            score=MixedHeightCorrection.score(context,result)
                            if score is not None:
                                height_candidates.append((score,result,effective))
                                height_candidates.sort(key=lambda row:(row[0],row[1].rows*row[1].across))
                                del height_candidates[MixedHeightCorrection.max_candidates:]
                except ValueError as exc:failures.append(str(exc))
                if attempts>=cls.max_candidates:break
            if attempts>=cls.max_candidates:break
        # Preserve every accepted primary mesh. Only after its bounded scan
        # fails may a matching crown redistribute interior heights. Full source
        # mouths/neighbor rows, zero curbs, crest, XY and direct quads stay fixed.
        for _,result,effective in height_candidates:
            try:return cls.certify(context,MixedHeightCorrection.solve(context,result),effective)
            except ValueError as exc:failures.append(str(exc))
        raise ValueError('No bounded mixed transition: '+'; '.join(dict.fromkeys(failures[-8:])))
