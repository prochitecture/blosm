"""Retain source curb corners and curve accuracy on whole-profile entrance strips."""
import math
import numpy as np
from .bundle_angular_ports import AngularProfileBridge,ray_jet
from .bundle_angular_open_islands import clip_ray
from .bundle_angular_curves import polyline_distance
from .bundle_axis import _simplify
from .mesh import Mesh,cross
from .entrance_slopes import EntranceSlopeGuard


class SourcePreservingAngularBridge(AngularProfileBridge):
    def _errors(self,mesh,outer_ids,inner_ids,role):
        errors=[];frame=self.domain.frame
        for ids,ray in zip((outer_ids,inner_ids),self.domain.sides(role)):
            xy=mesh.xyz[ids,:2];station=frame.local(xy[[0,-1]])[:,0]
            reference=clip_ray(ray,*sorted(station),{})
            query=np.concatenate([(1.-t)*xy[:-1]+t*xy[1:] for t in np.linspace(0.,1.,9)])
            errors.append(max(float(polyline_distance(query,reference).max()),
                              float(polyline_distance(reference,xy).max())))
        return errors

    def build(self,role,end,margin,*,physical_output=False,extra_rows=0):
        limit=1e-7 if self.domain.ports[role,end].path.exact else .9*self.domain.context.network.settings.curve_tolerance
        failure=None
        if not getattr(self,'sampled_selected',False):
            mesh,report=super().build(role,end,margin,physical_output=physical_output,extra_rows=extra_rows)
            across=self.domain.ports[role,end].intervals
            outer=np.arange(0,len(mesh.xyz),across+1);inner=self.interfaces['inner_ids']
            errors=self._errors(mesh,outer,inner,role)
            if max(errors)<=limit:
                if not hasattr(self,'initial_parity'):self.initial_parity=(len(inner)-1)%2
                report.update(source_curb_errors=errors,source_curb_limit=limit,source_sampled=False)
                return mesh,report
            failure=errors
        self.sampled_selected=True
        desired=(self.initial_parity+extra_rows)%2 if hasattr(self,'initial_parity') else None
        mesh,report=self._sampled(role,end,margin,physical_output,desired)
        if not hasattr(self,'initial_parity'):self.initial_parity=(len(self.interfaces['inner_ids'])-1)%2
        report['previous_curb_errors']=failure
        return mesh,report

    def _sampled(self,role,end,margin,physical_output,desired):
        domain=self.domain;context=domain.context;frame=domain.frame;settings=context.network.settings
        port=domain.ports[role,end];plane=domain.head+margin if end==0 else domain.tail-margin
        curves=domain.sides(role);across=port.intervals
        front=port.coordinates(np.array((0.,1.)))[0];starts=frame.local(front)[:,0]
        sides=[]
        for index,ray in enumerate(curves):
            points=clip_ray(ray,*sorted((starts[index],plane)),{})
            if end==1:points=points[::-1]
            points[0]=front[index]
            points=_simplify(points,1e-9 if port.path.exact else .15*settings.curve_tolerance)
            dense=[]
            for a,b in zip(points,points[1:]):
                count=max(1,int(math.ceil(np.linalg.norm(b-a)/settings.spacing)))
                dense.extend(a+np.linspace(0.,1.,count+1)[:-1,None]*(b-a))
            sides.append(np.vstack((*dense,points[-1])))
        def split(index):
            points=sides[index];at=int(np.argmax(np.linalg.norm(np.diff(points,axis=0),axis=1)))
            sides[index]=np.insert(points,at+1,.5*(points[at]+points[at+1]),axis=0)
        if desired is not None and (len(sides[1])-1)%2!=desired:split(1)
        if (len(sides[0])-len(sides[1]))%2:split(0)
        if max(map(len,sides))>512:raise ValueError('Source entrance exceeds its boundary sample budget')
        pattern=domain.patterns_for(port).quad(len(sides[0])-1,across,len(sides[1])-1,across)
        if len(pattern.quads)>4096:raise ValueError('Source entrance exceeds its cell budget')
        u,v=pattern.uv.T;ids0,ids1=(context.patterns.boundary(pattern,side) for side in (0,2))
        stations=[frame.local(points)[:,0] for points in sides]
        longitudinal=[np.interp(u,pattern.uv[ids,0],values) for ids,values in zip((ids0,ids1),stations)]
        derivatives=[(values[1]-values[0])/(pattern.uv[ids[1],0]-pattern.uv[ids[0],0])
                     for ids,values in zip((ids0,ids1),stations)]
        mouth,tangent,physical=port.coordinates(v);mouth_uv=frame.local(mouth)
        correction=mouth_uv[:,0]-((1.-v)*starts[0]+v*starts[1])
        station=(1.-v)*longitudinal[0]+v*longitudinal[1]+(1.-u)*correction
        dsdu=(1.-v)*derivatives[0]+v*derivatives[1]-correction
        outer,inner=curves;y0,_=ray_jet(outer,station);y1,_=ray_jet(inner,station)
        a,da=ray_jet(outer,mouth_uv[:,0]);b,db=ray_jet(inner,mouth_uv[:,0])
        shift=mouth_uv[:,1]-((1.-v)*a+v*b)
        radial=mouth-frame.center;radial/=np.linalg.norm(radial,axis=1)[:,None]
        radius=np.linalg.norm(mouth-frame.center,axis=1);angular=np.c_[-radial[:,1],radial[:,0]]
        source_slope=-radius/frame.scale*np.sum(radial*tangent,axis=1)/np.sum(angular*tangent,axis=1)
        derivative=(source_slope-((1.-v)*da+v*db))*dsdu
        shift[(v==0.)|(v==1.)]=0.;derivative[(v==0.)|(v==1.)]=0.
        y=(1.-v)*y0+v*y1+(1.-3.*u*u+2.*u*u*u)*shift+(u-2.*u*u+u*u*u)*derivative
        xy=frame.world(np.c_[station,y]);xy[ids0]=sides[0];xy[ids1]=sides[1]
        front_ids=context.patterns.boundary(pattern,3);back_ids=context.patterns.boundary(pattern,1)
        xy[front_ids]=mouth[front_ids]
        exaggeration=1. if physical_output else settings.exaggeration
        xyz=np.c_[xy,port.profile.datum+exaggeration*(physical-port.profile.datum)]
        faces=pattern.quads;edges=np.roll(xy[faces],-1,axis=1)-xy[faces]
        if np.all(cross(edges,np.roll(edges,-1,axis=1))<0.):faces=faces[:,::-1]
        mesh=Mesh(xyz,faces);quality=mesh.validate()
        if quality['maximum_edge_ratio']>8.:raise ValueError('Source entrance exceeds edge-ratio limit')
        slopes=mesh.slope_statistics(exaggeration)
        if slopes['maximum']>.1+1e-8:raise ValueError('Source entrance exceeds physical slope limit')
        entrance=EntranceSlopeGuard().validate(mesh,front_ids,tangent,physical)
        mismatch=float(np.max(abs(domain.field.values(xy[back_ids])-physical[back_ids])))
        if mismatch>1e-8:raise ValueError('Source entrance changed the body-cut profile')
        errors=self._errors(mesh,ids0,ids1,role)
        limit=1e-7 if port.path.exact else .9*settings.curve_tolerance
        if max(errors)>limit:raise ValueError('Source entrance exceeded the curb approximation budget')
        keys=[('bridge',role,end,'inside',i) for i in range(len(xyz))]
        for label,ids in (('mouth',front_ids),('body cut',back_ids)):
            for i,index in enumerate(ids):keys[index]=('bridge',role,end,label,i)
        self.interfaces=dict(keys=keys,mouth=(mesh.xyz[front_ids].copy(),[keys[i] for i in front_ids]),
            body_cut=(mesh.xyz[back_ids].copy(),[keys[i] for i in back_ids]),inner_ids=ids1)
        return mesh,dict(role=role,end=end,margin=margin,quality=quality,physical_slopes=slopes,
            entrance_longitudinal_bound=entrance,exact_original_mouth_xy=True,body_cut_profile_error=mismatch,
            whole_profile_intervals=across,all_direct_quads=True,source_curb_errors=errors,
            source_curb_limit=limit,source_sampled=True,source_side_intervals=[len(points)-1 for points in sides])
