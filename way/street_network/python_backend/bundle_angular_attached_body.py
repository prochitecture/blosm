"""Physical angular pavement assembled with whole-profile mouths and shared separator rims."""
import numpy as np
from .bundle_angular_island_chain import MultiIslandSurfaceBuilder
from .mesh import Assembler
from .model import Key
from .validation import validate_overlaps
from .bundle_angular_source_strips import SourcePreservingAngularBridge as AngularProfileBridge
from .bundle_angular_open_islands import OpenAngularIsland,gap_intervals,clip_ray


class AngularAttachedBody:
    def __init__(self,domain):self.domain=domain

    def build(self,margin=2.):
        domain=self.domain;frame=domain.frame;settings=domain.context.network.settings
        first,last=domain.head+margin,domain.tail-margin
        regions,roots=gap_intervals(domain,first,last)
        if not regions[0][2] or not regions[-1][2]:
            raise ValueError('Prototype needs separated full profiles at both body cuts')
        surface=MultiIslandSurfaceBuilder(domain.context.patterns,spacing=settings.spacing,
            curve_tolerance=settings.curve_tolerance,max_cells=30000)
        grass=Assembler();rims={};islands=[];cut_interfaces={}
        for index,(a,b,positive) in enumerate(regions):
            if not positive:continue
            failures=[]
            for distance in (.5,5./7.,1.,9./7.):
                try:
                    builder=OpenAngularIsland(domain)
                    mesh,boundary,report=builder.build(a,b,roots,index,distance=distance)
                    break
                except ValueError as error:failures.append(str(error))
            else:raise ValueError('Open grass body failed: '+str(failures))
            report['prior_failures']=failures
            grass.add(mesh.xyz,mesh.quads,builder.keys);rims[index]=boundary;islands.append(report)
            if a==first:cut_interfaces[0]=builder.cuts[0]
            if b==last:cut_interfaces[1]=builder.cuts[1]
        # The same crest chain closes the overlap for both pavement halves.
        chains=[[],[]]
        for index,(a,b,positive) in enumerate(regions):
            if positive:
                for role in (0,1):chains[role].append(rims[index][role])
            else:
                if a not in roots or b not in roots:
                    raise ValueError('Prototype overlap interval needs two explicit contacts')
                ridge=surface.ridge(frame,domain.field,roots[a],roots[b],('attached crest',index))
                ridge[1][0]=rims[index-1][0][1][-1]
                ridge[1][-1]=rims[index+1][0][1][0]
                for role in (0,1):chains[role].append(ridge)
        seams=[surface.join_chains(parts) for parts in chains]
        pavement=Assembler();connectors={};bridge_reports=[];prescribed={};self.bridge_interfaces={}
        for role in (0,1):
            for end in (0,1):
                builder=AngularProfileBridge(domain)
                mesh,report=builder.build(role,end,margin,physical_output=True)
                if role==1 and hasattr(domain,'cap_parities'):
                    first_rows=len(self.bridge_interfaces[0,end]['inner_rim'])-1
                    own_rows=len(builder.interfaces['inner_ids'])-1
                    if (first_rows+own_rows+domain.cap_parities[end])%2:
                        mesh,report=builder.build(role,end,margin,physical_output=True,extra_rows=1)
                        report['coordinated_extra_row']=True
                    if (first_rows+len(builder.interfaces['inner_ids'])-1+domain.cap_parities[end])%2:
                        raise ValueError('Angular entrance samples disagree with the actual cap boundary parity')
                interface=builder.interfaces
                cut,keys=interface['body_cut'];cut=cut.copy();keys=list(keys)
                target=seams[role][0][0 if end==0 else -1]
                if np.linalg.norm(cut[-1,:2]-target)>1e-8:
                    raise ValueError('Grass and pavement have different inner cut positions')
                cut[-1,:2]=target
                interface['keys'][interface['keys'].index(keys[-1])]=seams[role][1][0 if end==0 else -1]
                keys[-1]=seams[role][1][0 if end==0 else -1]
                for point,key in zip(cut,keys):prescribed[key]=point.copy()
                for point,key in zip(mesh.xyz,interface['keys']):
                    if key in prescribed:point[:]=prescribed[key]
                pavement.add(mesh.xyz,mesh.quads,interface['keys'])
                interface['inner_rim']=mesh.xyz[interface['inner_ids']].copy()
                self.bridge_interfaces[role,end]=interface
                connectors[role,end]=(cut,keys)
                bridge_reports.append(report)
        body_reports=[]
        for role in (0,1):
            ray=domain.sides(role)[0];curb=clip_ray(ray,first,last,{})
            curb,kinds,names=surface.simplify_sample(curb,[('outer',role)]*(len(curb)-1),
                .6*settings.curve_tolerance if not domain.ports[role,0].path.exact else 1e-9,
                ('attached outer',role))
            head,hk=connectors[role,0];tail,tk=connectors[role,1]
            if max(np.linalg.norm(head[0,:2]-curb[0]),np.linalg.norm(tail[0,:2]-curb[-1]))>1e-8:
                raise ValueError('Pavement outer curb misses the prescribed profile cut')
            curb[[0,-1]]=head[0,:2],tail[0,:2];names[0],names[-1]=hk[0],tk[0]
            boundary=(np.vstack((head[::-1,:2],curb[1:],tail[1:,:2])),
                [('cut','head')]*(len(head)-1)+kinds+[('cut','tail')]*(len(tail)-1),
                hk[::-1]+names[1:]+tk[1:])
            name=f'attached {role}/body'
            boundary=surface.balanced_physical(boundary,seams[role],name)
            mesh,keys,report=surface._direct_single(frame,domain.field,boundary,seams[role],name)
            for point,key in zip(mesh.xyz,keys):
                if key in prescribed:
                    if np.linalg.norm(point-prescribed[key])>1e-8:
                        raise ValueError('Angular body changed a prescribed whole profile')
                    point[:]=prescribed[key]
            pavement.add(mesh.xyz,mesh.quads,keys);body_reports.append(report)
        pave,grassy=pavement.mesh(),grass.mesh('grass')
        # Only physical interfaces of this partial body are shared with pavement;
        # its two open grass cuts are reserved for the actual rounded end caps.
        pe=surface.boundary_edges(pave);ge=surface.boundary_edges(grassy)
        available={tuple(sorted((tuple(pave.xyz[a]),tuple(pave.xyz[b])))) for a,b in pe}
        cuts={tuple(sorted((tuple(np.r_[a,domain.field.ceiling]),tuple(np.r_[b,domain.field.ceiling]))))
            for xy,keys in cut_interfaces.values() for a,b in zip(xy,xy[1:])}
        interfaces=0
        for a,b in ge:
            edge=tuple(sorted((tuple(grassy.xyz[a]),tuple(grassy.xyz[b]))))
            if edge in cuts:continue
            if edge not in available:raise ValueError('Partial grass-body interface is not a whole pavement edge')
            interfaces+=1
        slopes=pave.slope_statistics();longitudinal=frame.longitudinal_bound(pave)
        assert slopes['maximum']<=.1+1e-8 and longitudinal<=.03+1e-8
        assert pave.xyz[:,2].min()>=domain.field.datum-1e-8 and pave.xyz[:,2].max()<=domain.field.ceiling+1e-8
        overlap=validate_overlaps({domain.bundle.key:pave,Key('separator',domain.bundle.id):grassy},cell_size=2.)
        self.pavement_keys=[None]*len(pave.xyz)
        for key,at in pavement.lookup.items():self.pavement_keys[at]=key
        self.grass_keys=[None]*len(grassy.xyz)
        for key,at in grass.lookup.items():self.grass_keys[at]=key
        self.cut_interfaces=cut_interfaces
        return pave,grassy,dict(pavement=pave.validate(),grass=grassy.validate(),islands=islands,
            connectors=bridge_reports,body=body_reports,physical_slopes=slopes,
            angular_longitudinal_bound=longitudinal,shared_grass_edges=interfaces,
            overlap_check=overlap,rounded_end_caps_included=False)
