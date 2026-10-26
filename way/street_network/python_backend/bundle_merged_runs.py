"""Original outer junctions/caps connected to the new merged run body."""
from dataclasses import replace
import numpy as np
from . import model
from .bundle_angular_ports import AngularPort
from .bundle_angular_source_strips import SourcePreservingAngularBridge
from .bundle_visible_runs import VisibleBundleRuns
from .bundle_regions import interior_crossings
from .bundle_cap_sampling import SourceRimSampling
from .bundles import BundleProvider,Member,Cap
from .separator_ports import PORTS,SeparatorPortIndex
from .routes import RouteTopology
from .mesh import Assembler,Mesh
from .validation import validate_overlaps
from .bundle_merged_body import MergedRunBody


class ExternalRunDomain:
    def __init__(self,context,bundle,approach):
        self.context,self.frame,self.field = context,approach.frame_chart,approach.field
        self.rays,self.bundle,self.ports = self.field.curves,bundle,{}
        crossing = interior_crossings(context.network,self.bundle)
        runs = VisibleBundleRuns.build(context.network,self.bundle,crossing)
        wanted = tuple(path.section.id for path in approach.paths)
        self.run, = (run for run in runs if run.source_sections == wanted)
        self.end = 1-approach.end_index
        for role,member in enumerate(self.run.members):
            visit, = member.visits
            endpoint = member.ends[self.end]
            port = context[PORTS].resolve(endpoint)
            if port is None:
                raise ValueError('The original outer run endpoint has no physical junction mouth')
            path = context[model.Key('path',visit.section_id)]
            mouth = context[port.site.layout_key].mouths[endpoint.section_id,endpoint.end]
            profile = context[model.Key('profile',visit.section_id)].oriented(0 if visit.forward else 1)
            self.ports[role,self.end] = AngularPort(role,self.end,path,visit.forward,mouth.station,mouth.transverse_intervals,profile)
        values = np.concatenate([self.frame.local(port.coordinates(np.array((0.,1.)))[0])[:,0]
                                 for port in self.ports.values()])
        if self.end == 0:
            self.head,self.tail = float(values.max()),approach.origin_station
        else:
            self.head,self.tail = approach.origin_station,float(values.min())

    def sides(self,role):
        return (self.rays[0],self.rays[1]) if role == 0 else (self.rays[3],self.rays[2])

    def patterns_for(self,port):
        return self.context.patterns


def assign_declared_samples(mesh,points,keys,prefix):
    """Match only supplied physical interfaces; these are not persistent IDs."""
    result = [(prefix,'inside',i) for i in range(len(mesh.xyz))]
    for point,key in zip(points,keys):
        distances = np.linalg.norm(mesh.xyz-point,axis=1)
        hits = np.flatnonzero(distances < 1e-8)
        if len(hits) != 1:
            raise ValueError('A complete declared cap/interface sample is missing or ambiguous')
        at = int(hits[0])
        result[at] = key
        mesh.xyz[at] = point
    return result


class MergedRunWithExternalCap:
    def build(self,context,bundle,approach,merged_cut):
        domain = ExternalRunDomain(context,bundle,approach)
        end,settings = domain.end,context.network.settings
        members = tuple(Member((),np.empty((0,3)),member.ends) for member in domain.run.members)
        sampling = SourceRimSampling(RouteTopology(context.network).attachments(context.network))
        actual = BundleProvider.cap(context,*members,end,rim_source=sampling)
        if actual is None:
            raise ValueError('The original attached end lost its rounded physical cap')
        margin = 2.*settings.spacing
        meshes,interfaces,reports = {},{},[]
        for role in (0,1):
            bridge = SourcePreservingAngularBridge(domain)
            mesh,report = bridge.build(role,end,margin,physical_output=True)
            if role == 1:
                first = len(interfaces[0]['inner_ids'])-1
                second = len(bridge.interfaces['inner_ids'])-1
                if (first+second+len(actual.rim)-1)%2:
                    mesh,report = bridge.build(role,end,margin,physical_output=True,extra_rows=1)
                    report['coordinated_extra_row'] = True
            meshes[role],interfaces[role] = mesh,bridge.interfaces
            reports.append(report)
        remote_station = domain.head+margin if end == 0 else domain.tail-margin
        body = MergedRunBody()
        paving,grass,body_report = body.build(context,bundle,approach,merged_cut,remote_station=remote_station,
            remote_cuts=tuple(interfaces[role]['body_cut'][0] for role in (0,1)))
        if grass is None or end not in body.open_grass_cuts:
            raise ValueError('An attached separated mouth needs its real grass-body cut')
        assembly = Assembler()
        assembly.add(paving.xyz,paving.quads,body.pavement_keys)
        declared_points,declared_keys,rims = [],[],[]
        for role in (0,1):
            interface,mesh = interfaces[role],meshes[role]
            keys = [('external bridge',role,i) for i in range(len(mesh.xyz))]
            cut_points,cut_keys = interface['body_cut']
            back_ids = [interface['keys'].index(key) for key in cut_keys]
            mesh.xyz[back_ids] = body.cuts[role,end]
            for i,key in zip(back_ids,body.cut_keys[role,end]):
                keys[i] = key
            assembly.add(mesh.xyz,mesh.quads,keys)
            inner = interface['inner_ids']
            rim = mesh.xyz[inner].copy()
            rims.append(rim if end == 0 else rim[::-1])
            declared_points.extend(mesh.xyz[inner])
            declared_keys.extend(keys[i] for i in inner)
        cuts,cut_names = body.open_grass_cuts[end]
        # OpenAngularIsland stores its flat cut from low to high. Cap.free is
        # member-0 to member-1; use the existing whole samples in reverse order.
        flat = Cap.free(rims[0][-1 if end == 0 else 0],rims[1][-1 if end == 0 else 0],len(cuts)-1)
        cap_members = tuple(replace(member,rim=rim) for member,rim in zip(members,rims))
        cap = BundleProvider.region(context,domain.bundle,cap_members,
            [actual,flat] if end == 0 else [flat,actual],0)
        cap_slope=cap.mesh.slope_statistics()['maximum']
        if cap_slope>.1+1e-8:
            # Retain every physical rim value and reconcile only interior
            # elevations when the plateau-distance seed is too steep.
            from .slope_projection import SlopeProjection
            from .drainage import BoundaryDrainage
            edges=np.sort(np.stack((cap.mesh.quads,np.roll(cap.mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
            unique,counts=np.unique(edges,axis=0,return_counts=True)
            fixed=np.zeros(len(cap.mesh.xyz),bool);fixed[np.unique(unique[counts==1])]=True
            boundary=cap.mesh.xyz[fixed].copy()
            floors=np.full(len(cap.mesh.xyz),domain.field.datum)
            corrected=cap.mesh;rounds=[]
            for iteration in range(8):
                corrected,projection=SlopeProjection(max_cells=body.max_cells).solve(corrected,fixed,
                    datum=domain.field.datum,ceiling=domain.field.ceiling,maximum_slope=.1,minimum_heights=floors)
                filled,drainage=BoundaryDrainage(max_cells=body.max_cells).fill(corrected,fixed)
                rounds.append(dict(projection=projection,drainage=drainage))
                if drainage['maximum_fill']<=1e-6:break
                raised=filled>corrected.xyz[:,2]+1e-10
                floors[raised]=np.maximum(floors[raised],filled[raised])
            else:raise ValueError('Rounded grass-cap reconciliation exhausted its finite iteration budget')
            if not np.array_equal(corrected.xyz[fixed],boundary):
                raise ValueError('Rounded grass-cap reconciliation changed a physical boundary')
            cap.mesh=corrected
            cap.metadata=dict(cap.metadata,interior_height_reconciliation=dict(initial_slope=cap_slope,
                rounds=rounds,exact_boundary=True,crown_ceiling=domain.field.ceiling))
        for point,key in zip(np.c_[cuts,np.full(len(cuts),domain.field.ceiling)],cut_names):
            if key not in declared_keys:
                declared_points.append(point)
                declared_keys.append(key)
        cap_keys = assign_declared_samples(cap.mesh,np.asarray(declared_points),declared_keys,'outer grass cap')
        grassy = Assembler()
        grassy.add(grass.xyz,grass.quads,body.grass_keys)
        grassy.add(cap.mesh.xyz,cap.mesh.quads,cap_keys)
        pavement,separator = assembly.mesh(),grassy.mesh('grass')
        if len(pavement.quads)+len(separator.quads) > body.max_cells:
            raise ValueError('Complete run plus actual cap exceeds the combined cell budget')
        # Every separator edge except the actual junction/crossroad-facing cap
        # must be a complete pavement edge, including the external bridge rims.
        from .bundle_angular_surfaces import AngularSurfaceBuilder
        edge_reader = AngularSurfaceBuilder(context.patterns,spacing=settings.spacing,curve_tolerance=settings.curve_tolerance)
        pave_edges = edge_reader.boundary_edges(pavement)
        grass_edges = edge_reader.boundary_edges(separator)
        available = {tuple(sorted((tuple(pavement.xyz[a]),tuple(pavement.xyz[b])))) for a,b in pave_edges}
        physical_cap = {tuple(sorted((tuple(a),tuple(b)))) for a,b in zip(actual.rim,actual.rim[1:])}
        shared,original_cap = 0,0
        for a,b in grass_edges:
            edge = tuple(sorted((tuple(separator.xyz[a]),tuple(separator.xyz[b]))))
            if edge in available:
                shared += 1
            elif edge in physical_cap:
                original_cap += 1
            else:
                raise ValueError('Complete run separator has an undeclared exterior edge')
        published_pavement = Mesh(pavement.xyz.astype(np.float32),pavement.quads)
        published_separator = Mesh(separator.xyz.astype(np.float32),separator.quads,'grass')
        if max(published_pavement.validate()['maximum_edge_ratio'],published_separator.validate()['maximum_edge_ratio']) > 8.:
            raise ValueError('Complete run publication needs better cell grading')
        slopes = dict(pavement=published_pavement.slope_statistics()['maximum'],
                      separator=published_separator.slope_statistics()['maximum'])
        if max(slopes.values()) > .1+1e-6:
            raise ValueError('Complete attached run exceeded its physical slope limits: '+str(slopes))
        overlap = validate_overlaps({model.Key('pave_probe',0):pavement,model.Key('grass_probe',0):separator},cell_size=2.)
        self.pavement_keys = [None]*len(pavement.xyz)
        for key,at in assembly.lookup.items():
            self.pavement_keys[at] = key
        self.body,self.domain = body,domain
        return pavement,separator,dict(body=body_report,external_connectors=reports,cap=cap.metadata,
            pavement=pavement.validate(),separator=separator.validate(),physical_slopes=slopes,
            shared_grass_edges=shared,actual_cap_edges=original_cap,overlap_check=overlap,
            source_sections=list(domain.run.source_sections),actual_external_rounded_cap=True,
            normal_addon_integration=True,all_direct_quads=True)
