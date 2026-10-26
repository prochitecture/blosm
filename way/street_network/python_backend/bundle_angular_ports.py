"""Original-mouth profiles and certified angular coordinates for attached Bundles."""
from dataclasses import dataclass
import math
import numpy as np
from . import model
from .bundle_axis import AxialChartUnavailable
from .bundle_angular_curves import AngularCoordinates,AngularCrown,RayCurb
from .mesh import Mesh,cross
from .entrance_slopes import EntranceSlopeGuard


class TraversedPath:
    def __init__(self,path,forward):
        self.original,self.forward=path,forward
        self.section,self.length=path.section,path.length
        self.points=path.points if forward else path.points[::-1]
    def frame(self,station):
        p,t,n=self.original.frame(station if self.forward else self.length-np.asarray(station))
        return p,t if self.forward else -t,n if self.forward else -n


@dataclass
class AngularPort:
    role:int
    end:int
    path:object
    forward:bool
    station:float
    intervals:int
    profile:object

    def coordinates(self,v):
        q=(2.*np.asarray(v)-1.)*(1 if self.role==0 else -1)
        original=q*(1 if self.forward else -1)
        point,tangent,normal=self.path.frame(self.station)
        xy=point+(self.path.section.offset+.5*self.path.section.width*original)[...,None]*normal
        return xy,np.broadcast_to(tangent,xy.shape),self.profile.datum+self.profile.value(q)


class AngularAttachedDomain:
    """Use original polynomials and actual mouth stations/partitions."""
    def make_frame(self, paths):
        return AngularCoordinates(paths)

    def make_port(self, context, sites, role, end, visit, path, profile):
        source_end=model.End(visit.section_id,end if visit.forward else 1-end)
        site=sites.get(source_end)
        if site is None or len(site.ends)<3:
            raise AxialChartUnavailable('Attached angular member ends need actual multi-arm mouths')
        mouth=context[site.layout_key].mouths.get((source_end.section_id,source_end.end))
        if mouth is None:
            raise AxialChartUnavailable('Attached angular mouth is consumed by another junction')
        return AngularPort(role,end,path,visit.forward,mouth.station,mouth.transverse_intervals,profile)

    def patterns_for(self, port):
        return self.context.patterns

    def __init__(self,context,bundle):
        self.context,self.bundle=context,bundle
        if len(bundle.members)!=2 or any(len(member.path)!=1 for member in bundle.members):
            raise AxialChartUnavailable('Attached angular charts require two complete single-Section members')
        paths=[];self.ports={};curves=[];self.curve_proofs=[]
        sites={end:junction for junction in (*context.network.intersections,*context.network.connections) for end in junction.ends}
        for role,member in enumerate(bundle.members):
            visit=member.path[0];path=context[model.Key('path',visit.section_id)]
            paths.append(TraversedPath(path,visit.forward))
            profile=context[model.Key('profile',visit.section_id)].oriented(0 if visit.forward else 1)
            if profile.half!=(1 if role==0 else -1):
                raise AxialChartUnavailable('Attached angular members require inward-rising shared half crowns')
            for end in (0,1):
                self.ports[role,end]=self.make_port(context,sites,role,end,visit,path,profile)
            samples=path.samples(0.,path.length,context.network.settings.spacing,
                                 tolerance=.1*context.network.settings.curve_tolerance)
            samples=np.unique(np.r_[samples,[self.ports[role,end].station for end in (0,1)]])
            for side in (-1.,1.):
                original=side*(1 if visit.forward else -1)
                xy=path.point(samples,path.section.offset+.5*path.section.width*original)
                curves.append(xy if visit.forward else xy[::-1])
        from .bundles import CenterlineRange,certify_disjoint
        visible=[]
        for role,view in enumerate(paths):
            stations=sorted(self.ports[role,end].station for end in (0,1))
            visible.append((CenterlineRange.build(view.original,*stations),))
        self.proof_checks=certify_disjoint(*visible)
        self.frame=self.make_frame(paths)
        for path in paths:self.curve_proofs.append(self.frame.certify(path.original,forward=path.forward))
        first=self.ports[0,0].profile
        if any(abs(port.profile.rise-first.rise)>1e-8 or abs(port.profile.datum-first.datum)>1e-8 for port in self.ports.values()):
            raise AxialChartUnavailable('Attached angular members require a common crown datum and rise')
        self.field=AngularCrown(self.frame,curves,datum=first.datum,rise=first.rise)
        self.rays=self.field.curves
        self.head=max(float(self.frame.local(port.coordinates(np.array((0.,1.)))[0])[:,0].max())
                      for (role,end),port in self.ports.items() if end==0)
        self.tail=min(float(self.frame.local(port.coordinates(np.array((0.,1.)))[0])[:,0].min())
                      for (role,end),port in self.ports.items() if end==1)

    def sides(self,role):
        return (self.rays[0],self.rays[1]) if role==0 else (self.rays[3],self.rays[2])


def ray_jet(ray,station):
    """Exact dy/ds of each Cartesian chord in the angular coordinate chart."""
    frame=ray.frame;station=np.asarray(station)
    index=np.clip(np.searchsorted(ray.stations,station,side='right')-1,0,len(ray.points)-2)
    edge=ray.points[index+1]-ray.points[index]
    xy=ray(station);delta=xy-frame.center;radius=np.linalg.norm(delta,axis=-1)
    radial=delta/radius[...,None]
    denominator=cross(edge,radial)
    if np.any(abs(denominator)<1e-10):raise ValueError('No finite angular curb jet')
    slope=radius/frame.scale*np.sum(edge*radial,axis=-1)/denominator
    return frame.local(xy)[...,1],slope


class AngularProfileBridge:
    def __init__(self,domain):self.domain=domain

    def build(self,role,end,margin,*,physical_output=False,extra_rows=0):
        domain=self.domain;context=domain.context;frame=domain.frame;port=domain.ports[role,end]
        plane=domain.head+margin if end==0 else domain.tail-margin
        ao,ai,bi,bo=domain.field.edges(plane)
        if max(ai,bo)-min(ao,bi)<=1e-7:raise ValueError('The connector cut needs separated whole profiles')
        outer,inner=domain.sides(role)
        across=port.intervals
        front,_,_=port.coordinates(np.linspace(0.,1.,across+1))
        local=frame.local(front)
        length=np.max(abs(plane-local[:,0]))
        if np.any((plane-local[:,0])*(1 if end==0 else -1)<=0.):raise ValueError('Connector cut precedes its mouth')
        if extra_rows not in (0,1):raise ValueError('Invalid coordinated entrance subdivision')
        rows=max(2,int(math.ceil(length/context.network.settings.spacing)))+extra_rows
        if rows*across>4096:raise ValueError('Angular entrance connector exceeds its cell budget')
        pattern=domain.patterns_for(port).grid(rows,across);u,v=pattern.uv.T
        mouth,tangent,physical=port.coordinates(v);mouth_uv=frame.local(mouth)
        first=mouth_uv[:,0];station=(1.-u)*first+u*plane
        y0,_=ray_jet(outer,station);y1,_=ray_jet(inner,station)
        a,da=ray_jet(outer,first);b,db=ray_jet(inner,first)
        shift=mouth_uv[:,1]-((1.-v)*a+v*b)
        delta=mouth-frame.center;radius=np.linalg.norm(delta,axis=1)
        radial=delta/radius[:,None];angular=np.c_[-radial[:,1],radial[:,0]]
        source_slope=-radius/frame.scale*np.sum(radial*tangent,axis=1)/np.sum(angular*tangent,axis=1)
        derivative=(source_slope-((1.-v)*da+v*db))*(plane-first)
        shift[(v==0.)|(v==1.)]=0.;derivative[(v==0.)|(v==1.)]=0.
        y=(1.-v)*y0+v*y1+(1.-3.*u*u+2.*u*u*u)*shift+(u-2.*u*u+u*u*u)*derivative
        xy=frame.world(np.c_[station,y])
        front_ids=context.patterns.boundary(pattern,3);back_ids=context.patterns.boundary(pattern,1)
        xy[front_ids]=mouth[front_ids]
        for side,ray in ((0,outer),(2,inner)):
            ids=context.patterns.boundary(pattern,side);xy[ids]=ray(station[ids])
        # Reapply the prescribed whole mouth after the sampled lateral chains.
        xy[front_ids]=mouth[front_ids]
        physical_heights=physical
        exaggeration=1. if physical_output else context.network.settings.exaggeration
        xyz=np.c_[xy,port.profile.datum+exaggeration*(physical_heights-port.profile.datum)]
        quads=pattern.quads
        edges=np.roll(xy[quads],-1,axis=1)-xy[quads]
        if np.all(cross(edges,np.roll(edges,-1,axis=1))<0.):quads=quads[:,::-1]
        mesh=Mesh(xyz,quads);quality=mesh.validate()
        if quality['maximum_edge_ratio']>8.:raise ValueError('Angular connector exceeds edge-ratio limit')
        slopes=mesh.slope_statistics(exaggeration)
        if slopes['maximum']>.1+1e-8:raise ValueError('Angular connector exceeds physical slope limit')
        entrance=EntranceSlopeGuard().validate(mesh,front_ids,tangent,physical)
        target=domain.field.values(xy[back_ids])
        if np.max(abs(target-physical[back_ids]))>1e-8:raise ValueError('Angular body-cut profile mismatch')
        assert np.array_equal(mesh.xyz[front_ids,:2],mouth[front_ids])
        keys=[('bridge',role,end,'inside',i) for i in range(len(xyz))]
        for label,ids in (('mouth',front_ids),('body cut',back_ids)):
            for i,index in enumerate(ids):keys[index]=('bridge',role,end,label,i)
        self.interfaces=dict(keys=keys,mouth=(mesh.xyz[front_ids].copy(),[keys[i] for i in front_ids]),
                             body_cut=(mesh.xyz[back_ids].copy(),[keys[i] for i in back_ids]),
                             inner_ids=context.patterns.boundary(pattern,2))
        return mesh,dict(role=role,end=end,margin=margin,quality=quality,physical_slopes=slopes,
            entrance_longitudinal_bound=entrance,exact_original_mouth_xy=True,
            body_cut_profile_error=float(np.max(abs(target-physical[back_ids]))),
            whole_profile_intervals=across,all_direct_quads=True)
