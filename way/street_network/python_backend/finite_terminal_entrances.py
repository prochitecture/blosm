"""A transverse common fibre on two finite original terminal prefixes."""
import math
import numpy as np
from .bundle_merged_entrances import PairedAngularEntrance, MemberInventory
from .bundle_angular_curves import AngularCrown
from .bundle_axis import AxialChartUnavailable
from .model import End,Key
from .paths import unit
from .routes import Piece
from .source_mouths import SourceMouth


class LinearFrame:
    def __init__(self,origin,axis):
        self.origin=np.asarray(origin,float);self.axis=unit(axis)
        self.left=np.array((-self.axis[1],self.axis[0]))

    def local(self,points):
        delta=np.asarray(points)-self.origin
        return np.stack((delta@self.axis,delta@self.left),axis=-1)

    def world(self,coordinates):
        c=np.asarray(coordinates)
        return self.origin+c[...,0,None]*self.axis+c[...,1,None]*self.left

    def gradient_world(self,points,along,across):
        del points
        return along[...,None]*self.axis+across[...,None]*self.left

    def longitudinal_bound(self,mesh):
        from .slope_projection import corner_matrices
        height=mesh.xyz[mesh.quads,2]-float(mesh.xyz[:,2].min())
        gradients=np.einsum('fcij,fj->fci',corner_matrices(mesh),height)
        return float(abs(gradients@self.axis).max())


class LinearRay:
    def __init__(self,frame,points):
        self.frame=frame;self.points=np.asarray(points,float)
        self.local=frame.local(self.points);self.stations=self.local[:,0]
        if np.any(np.diff(self.stations)<=1e-9):
            raise AxialChartUnavailable('A finite terminal curb is not axially monotone')

    def __call__(self,station):
        station=np.asarray(station)
        if np.any(station<self.stations[0]-1e-8) or np.any(station>self.stations[-1]+1e-8):
            raise AxialChartUnavailable('A finite terminal fibre would extrapolate an original curb')
        y=np.interp(station,self.stations,self.local[:,1])
        return self.frame.world(np.stack(np.broadcast_arrays(station,y),axis=-1))

    def jet(self,station,*,interval=None):
        station=np.asarray(station)
        if interval is None:
            i=np.clip(np.searchsorted(self.stations,station,side='right')-1,0,len(self.points)-2)
        else:
            first,last=sorted(interval)
            i=int(np.clip(np.searchsorted(self.stations,first,side='right')-1,0,len(self.points)-2))
            margin=32.*np.spacing(max(1.,abs(first),abs(last),float(np.max(abs(self.frame.origin)))))
            if (not first<last or first<self.stations[i] or last>self.stations[i+1]
                    or not np.isfinite(station).all() or np.any(station<first-margin) or np.any(station>last+margin)):
                raise ValueError('A one-sided linear curb jet must stay on its declared original chord')
        edge=self.local[i+1]-self.local[i]
        return self.frame.local(self(station))[...,1],edge[...,1]/edge[...,0]


class LinearCrown(AngularCrown):
    def __init__(self,frame,curves,*,datum,rise):
        self.frame=frame;self.curves=tuple(LinearRay(frame,p) for p in curves)
        self.datum,self.rise=float(datum),float(rise);self.ceiling=self.datum+self.rise


class FinitePairedEntrance(PairedAngularEntrance):
    """Two original half-profile aliases on a bounded transverse source chart.

    The chart covers at most eight source widths of each terminal Section.
    Complete run identity survives separately; every mouth is inside all four
    finite, monotonically advancing original curb ranges.
    """
    def __init__(self,context,run,end,*,source_tolerance,source_parameters=None):
        from .bundles import CenterlineRange
        from .bundle_boundaries import _certify_increasing
        if len(run.members)!=2 or end not in (0,1):
            raise ValueError('A finite terminal descriptor requires two original member incidences')
        if not np.isfinite(source_tolerance) or not 0.<source_tolerance<=.00025:
            raise ValueError('Finite terminal source chord error must be positive and at most a quarter millimetre')
        self.run=run;self.source_run_sections=run.source_sections;self.end_index=end
        self.original_ends=tuple(m.ends[end] for m in run.members)
        self.visits=tuple(m.visits[0 if end==0 else -1] for m in run.members)
        if any(e!=End(v.section_id,int(not v.forward) if end==0 else int(v.forward))
                for v,e in zip(self.visits,self.original_ends)):
            raise ValueError('Finite terminal descriptor lost its original source incidence')
        self.paths=tuple(context[Key('path',v.section_id)] for v in self.visits)
        self.path=MemberInventory(tuple(Piece(p,True) for p in self.paths))
        self.sign=1. if end==0 else -1.
        origins=[];directions=[];profiles=[]
        for role,(p,v,e) in enumerate(zip(self.paths,self.visits,self.original_ends)):
            point,tangent,_=p.frame(0. if e.end==0 else p.length)
            origins.append(point);directions.append(tangent*(1. if v.forward else -1.))
            profile=context[Key('profile',v.section_id)].oriented(int(not v.forward))
            if profile.half!=(1 if role==0 else -1):
                raise ValueError('Finite terminal lost an original Bundle half-profile')
            profiles.append(profile)
        if float(unit(directions[0])@unit(directions[1]))<np.cos(np.deg2rad(30.)):
            raise AxialChartUnavailable('Finite paired terminal directions differ by more than thirty degrees')
        self.frame_chart=LinearFrame(np.mean(origins,axis=0),np.sum(directions,axis=0))
        if float((origins[0]-origins[1])@self.frame_chart.left)<=1e-6:
            raise AxialChartUnavailable('Finite terminal roles disagree with their original lateral order')
        self.source_ranges=[];curves=[];source_stations=[];proofs=[]
        reach=8.*max(p.section.width for p in self.paths)
        for p,v,e in zip(self.paths,self.visits,self.original_ends):
            if not p.exact and (abs(p.section.offset)+.5*p.section.width)*p.curvature_bound>=1.:
                raise AxialChartUnavailable('Finite terminal source offsets lack a regular curvature certificate')
            extent=min(p.length,reach)
            first,last=(0.,extent) if e.end==0 else (p.length-extent,p.length)
            self.source_ranges.append((first,last))
            lower=None
            spacing=min(.5,context.network.settings.spacing)
            if not p.exact:
                offset_reach=abs(p.section.offset)+.5*p.section.width
                curvature=p.curvature_bound/(1.-offset_reach*p.curvature_bound)
                speed=1.+offset_reach*p.curvature_bound
                if curvature>1e-12:
                    spacing=min(spacing,math.sqrt(8.*source_tolerance/curvature)/speed,
                        math.pi/(16.*curvature*speed))
                lower=min(_certify_increasing(3.*np.diff(c,axis=0)@self.frame_chart.axis
                    *(1. if v.forward else -1.),coordinate='finite terminal axial')
                    for c in CenterlineRange.build(p,first,last).curves)
            extra=(() if source_parameters is None else source_parameters.get(p.section.id,()))
            if any(not np.isfinite(s) or not 0.<=s<=p.length for s in extra):
                raise ValueError('Finite terminal contact stations must lie on their original finite Section')
            if source_parameters is None or not len(extra):
                samples=p.samples(first,last,spacing,tolerance=source_tolerance)
            else:
                from .paired_source_sampling import sample_source_parameters
                samples=sample_source_parameters(p,first,last,spacing,source_tolerance,extra,
                    maximum_samples=self.maximum_samples)
                first,last=float(samples[0]),float(samples[-1])
                self.source_ranges[-1]=(first,last)
            if len(samples)>self.maximum_samples:raise ValueError('Finite terminal exceeds its source sampling budget')
            for side in (-1.,1.):
                xy=p.point(samples,p.section.offset+side*(1. if v.forward else -1.)*.5*p.section.width)
                curves.append(xy if v.forward else xy[::-1])
            source_stations.append(samples if v.forward else samples[::-1])
            proofs.append(dict(section=p.section.id,interval=(first,last),axial_derivative_lower=lower))
        if max(abs(p.rise-profiles[0].rise)+abs(p.datum-profiles[0].datum) for p in profiles)>1e-9:
            raise ValueError('Finite paired terminal requires its original shared crown')
        self.field=LinearCrown(self.frame_chart,curves,datum=profiles[0].datum,rise=profiles[0].rise)
        self.source_curve_stations=tuple(source_stations)
        self.source_curve_tolerance=source_tolerance;self.source_prefix_proof=proofs
        self.origin_station=0.
        self.length=(min(r.stations[-1] for r in self.field.curves) if self.sign>0.
                     else -max(r.stations[0] for r in self.field.curves))
        minimum=(max(r.stations[0] for r in self.field.curves) if self.sign>0.
                 else -min(r.stations[-1] for r in self.field.curves))+1e-5
        self.minimum_station=max(0.,float(minimum))
        self.source,self.end=self.paths[0].section,self.original_ends[0].end
        high,low=self.bounds(self.minimum_station)
        self.width,self.offset=float(high-low),0.
        self.origin=self.frame_chart.world(np.array((0.,.5*(high+low))))
        self.direction=self.sign*self.frame_chart.axis
        self.profile=self.cross_section(self.minimum_station)

    def cut_plane(self,station):
        return self.frame_chart.world(np.array((self.chart_station(station),0.))),self.sign*self.frame_chart.axis

    def source_mouths(self,station):
        target=float(self.chart_station(station));peak_u=.5*(self.cross_section(station).peak+1.)
        aliases=[]
        for role,(p,v,e,interval) in enumerate(zip(self.paths,self.visits,self.original_ends,self.source_ranges)):
            low,high=interval;first,last=self.frame_chart.local(p.frame(np.array((low,high)))[0])[:,0]
            if not min(first,last)-1e-8<=target<=max(first,last)+1e-8:
                raise ValueError('Finite terminal cut lies outside an original centreline prefix')
            for _ in range(56):
                middle=.5*(low+high);value=float(self.frame_chart.local(p.frame(middle)[0])[0])
                if (value<target)==(last>first):low=middle
                else:high=middle
            at=.5*(low+high)
            if abs(float(self.frame_chart.local(p.frame(at)[0])[0])-target)>1e-7:
                raise ValueError('Finite terminal source-mouth station did not converge')
            outer=float(role if self.sign>0. else 1-role)
            factor=(1 if e.end==0 else -1)*(1 if role==0 else -1)*(1 if v.forward else -1)
            aliases.append(SourceMouth(e,at,outer,peak_u,.5*(1-factor),.5*(1+factor)))
        return tuple(aliases)

    def mouth_parameters(self,radial,spacing):
        del spacing
        left=self.profile.coordinate(-1.,np.linspace(0.,1.,self.profile.intervals(-1.,radial)+1))
        right=self.profile.coordinate(1.,np.linspace(0.,1.,self.profile.intervals(1.,radial)+1))
        return .5*(np.r_[left,right[::-1][1:]]+1.)[::-1]
