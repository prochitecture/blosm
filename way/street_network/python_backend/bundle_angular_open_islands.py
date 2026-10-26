"""Direct quad separator intervals bounded by computational cuts and physical tips."""
import math
import numpy as np
from .bundle_angular_curves import RayCurb
from .bundle_angular_islands import AngularIslandBuilder
from .four_side_charts import FourSidedQuadChart
from .quality_correction import QuadShapeCorrection
from .mesh import Mesh,Assembler,cross
from .validation import validate_overlaps


def clip_ray(ray,first,last,roots):
    middle=(ray.stations>first+1e-8)&(ray.stations<last-1e-8)
    points=np.vstack((ray(first),ray.points[middle],ray(last)))
    if first in roots:points[0]=roots[first]
    if last in roots:points[-1]=roots[last]
    return points


def gap_intervals(domain,first,last):
    upper,lower=domain.rays[1:3]
    nodes=np.unique(np.concatenate((np.asarray((first,last)),
        *(ray.stations[(ray.stations>first)&(ray.stations<last)] for ray in (upper,lower)))))
    values=domain.frame.local(upper(nodes))[:,1]-domain.frame.local(lower(nodes))[:,1]
    if np.any(abs(values[1:-1])<1e-9):
        raise ValueError('Prototype needs a source-knot contact policy')
    roots={}
    for i,(a,b) in enumerate(zip(nodes,nodes[1:])):
        if values[i]*values[i+1]>=0.:continue
        edges=[]
        for ray in (upper,lower):
            index=int(np.clip(np.searchsorted(ray.stations,.5*(a+b),side='right')-1,0,len(ray.points)-2))
            edges.append(ray.points[index:index+2])
        p,d=edges[0][0],edges[0][1]-edges[0][0]
        q,e=edges[1][0],edges[1][1]-edges[1][0]
        determinant=float(cross(d,e))
        if abs(determinant)<1e-12*np.linalg.norm(d)*np.linalg.norm(e):
            raise ValueError('Angular gap event has no isolated line intersection')
        point=p+float(cross(q-p,e)/determinant)*d
        station=float(domain.frame.local(point)[0])
        if not a-1e-7<station<b+1e-7:raise ValueError('Gap event is outside its source segments')
        roots[station]=point
    events=sorted((first,last,*roots));regions=[]
    for a,b in zip(events,events[1:]):
        edges=domain.field.edges(.5*(a+b))
        if edges[0]<=edges[2] or edges[1]<=edges[3]:
            raise ValueError('Prototype needs another outer-envelope ownership branch')
        positive=bool(edges[1]>edges[2])
        if regions and regions[-1][2]==positive:regions[-1]=(regions[-1][0],b,positive)
        else:regions.append((a,b,positive))
    return regions,roots


class OpenAngularIsland:
    def __init__(self,domain):self.domain=domain

    def build(self,first,last,roots,index,*,distance=.5,shear=.2):
        domain=self.domain;context=domain.context;settings=context.network.settings;frame=domain.frame
        upper,lower=(clip_ray(ray,first,last,roots) for ray in domain.rays[1:3])
        closed=(first in roots,last in roots)
        start=first+(distance if closed[0] else 0.);finish=last-(distance if closed[1] else 0.)
        if finish<=start+settings.spacing:raise ValueError('Island tips consume its positive body')
        sampler=AngularIslandBuilder(frame,context.patterns,domain.bundle.key,spacing=settings.spacing,
            curve_tolerance=settings.curve_tolerance,height=domain.field.ceiling,max_cells=5000,
            boundary_error=(0. if all(port.path.exact for port in domain.ports.values()) else .1*settings.curve_tolerance)+1e-7)
        rays=[RayCurb(frame,chain) for chain in (upper,lower)]
        error=max((sampler._cap_accuracy(ray,*interval) for ray in rays
                   for interval in ((first,start),(finish,last)) if interval[1]>interval[0]),default=0.)
        if error>sampler.cap_error_budget:raise ValueError('Open island tip exceeds the original curve budget')
        smooth=all(not port.path.exact for port in domain.ports.values())
        top,bottom=(sampler._trimmed(ray,start,finish,smooth,.1*settings.curve_tolerance) for ray in rays)
        if (len(top)-len(bottom))%2:
            at=int(np.argmax(np.linalg.norm(np.diff(top,axis=0),axis=1)))
            top=np.insert(top,at+1,.5*(top[at]+top[at+1]),axis=0)
        upper_keys=[('open island',index,'upper',i) for i in range(len(top))]
        lower_keys=[('open island',index,'lower',i) for i in range(len(bottom))]
        cuts=[];cut_keys=[]
        for end,at in ((0,0),(1,-1)):
            low,high=bottom[at],top[at]
            if closed[end]:
                uv=.5*(frame.local(low)+frame.local(high))
                uv[0]+=distance*shear*(1 if end==0 else -1)
                points=np.vstack((low,frame.world(uv),high))
            else:
                count=max(4,4*int(math.ceil(np.linalg.norm(high-low)/(4.*settings.spacing))))
                points=np.linspace(low,high,count+1)
            keys=[('open island',index,'cut',end,i) for i in range(len(points))]
            keys[0],keys[-1]=lower_keys[at],upper_keys[at]
            cuts.append(points);cut_keys.append(keys)
        rim=np.vstack((bottom,cuts[1][1:-1],top[::-1],cuts[0][-2:0:-1]))
        keys=lower_keys+cut_keys[1][1:-1]+upper_keys[::-1]+cut_keys[0][-2:0:-1]
        corners=(0,len(bottom)-1,len(bottom)+len(cuts[1])-2,len(bottom)+len(cuts[1])+len(top)-3)
        chart=FourSidedQuadChart(context.patterns,max_cells=5000)
        mapped,counts=chart.map(frame.local(rim),corners)
        xy=frame.world(mapped.xyz[:,:2]);xy[:len(rim)]=rim
        body=Mesh(np.c_[xy,np.full(len(xy),domain.field.ceiling)],mapped.quads,'grass')
        body,correction=QuadShapeCorrection(max_cells=5000,max_iterations=200).improve(body)
        if not correction['qualified']:raise ValueError('Open island needs another local body chart')
        keys.extend(('open island',index,'inside',i) for i in range(len(rim),len(body.xyz)))
        assembly=Assembler();assembly.add(body.xyz,body.quads,keys)
        exported=[(top,upper_keys),(bottom,lower_keys)]
        for end in (0,1):
            if not closed[end]:continue
            root=roots[first if end==0 else last];key=('open island',index,'tip',end)
            if end==0:
                points=np.vstack((root,bottom[0],cuts[0][1],top[0]))
                names=[key,lower_keys[0],cut_keys[0][1],upper_keys[0]]
            else:
                points=np.vstack((root,top[-1],cuts[1][1],bottom[-1]))
                names=[key,upper_keys[-1],cut_keys[1][1],lower_keys[-1]]
            assembly.add(np.c_[points,np.full(4,domain.field.ceiling)],np.array(((0,1,2,3),)),names)
            for role,(points,names) in enumerate(exported):
                exported[role]=(np.vstack((root,points)) if end==0 else np.vstack((points,root)),
                                [key]+names if end==0 else names+[key])
        mesh=assembly.mesh('grass');quality=mesh.validate()
        self.keys=[None]*len(mesh.xyz)
        for key,at in assembly.lookup.items():self.keys[at]=key
        self.cuts=tuple((points,names) for points,names in zip(cuts,cut_keys))
        if quality['maximum_edge_ratio']>8.:raise ValueError('Open angular island exceeds edge-ratio limit')
        validate_overlaps({domain.bundle.key:mesh},cell_size=2.)
        return mesh,exported,dict(quality=quality,closed_tips=closed,tip_error=error,
            rim_counts=counts,correction=correction,all_direct_quads=True,
            flat_cut_intervals=[len(cuts[end])-1 if not closed[end] else None for end in (0,1)])
