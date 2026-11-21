"""One finite original-source body between two prepared merged mouths.

This consumer has no candidate discovery hook. Callers supply the two actual
completed junction mouths, aligned hydrated source spans and a separately
certified run frame. Junction pavement is never constructed here.
"""
from types import SimpleNamespace
import numpy as np
from .bundle_angular_curves import AngularCrown
from .bundle_merged_bridge import MergedProfileBridge
from .bundle_merged_body import MergedRunBody
from .drainage import BoundaryDrainage
from .four_side_charts import quality
from .mesh import Assembler,Mesh
from .model import End,Key
from .ownership import CompoundProduct,surface_key
from .prepared_merged_profiles import PreparedMergedMouth
from .source_mesh_containment import certify_mesh
from .validation import validate_overlaps


class PairedSourceBodyPiece:
    max_cells=16000

    @staticmethod
    def make_body():
        from .split_merged_body import BoundedMergedRunBody
        return BoundedMergedRunBody()

    @staticmethod
    def _sources(context,bundle,spans,frame):
        if len(spans)!=2 or tuple(span.role for span in spans)!=(0,1):
            raise ValueError('A paired source body needs one aligned hydrated span per original role')
        by_end={e:j for j in (*context.network.intersections,*context.network.connections) for e in j.ends}
        curves=[];identities=[];proofs=[]
        for role,span in enumerate(spans):
            if not span.visits or not span.source_intervals or len(span.curves)!=2:
                raise ValueError('A paired source body needs hydrated original visits, intervals and curb curves')
            original=tuple(bundle.members[role].path);visits=tuple(span.visits)
            starts=[i for i,v in enumerate(original) if v==visits[0]]
            if len(starts)!=1 or original[starts[0]:starts[0]+len(visits)]!=visits:
                raise ValueError('A paired body span changed its original ordered Bundle visits')
            if len(span.source_intervals)!=len(visits):
                raise ValueError('A paired body span lost an original source interval')
            if tuple((p.section.id,p.forward) for p in span.route.pieces)!=tuple((v.section_id,v.forward) for v in visits):
                raise ValueError('A paired body route does not retain its original source visits')
            for visit,row,piece in zip(visits,span.source_intervals,span.route.pieces):
                path=context[Key('path',visit.section_id)]
                if (piece.path is not path or row['section']!=visit.section_id or row['forward']!=visit.forward
                        or not 0.<=row['first']<row['last']<=path.length):
                    raise ValueError('A paired body interval changed its original finite Path or orientation')
                proofs.append(frame.certify(path,forward=visit.forward));identities.append(visit.section_id)
            for previous,following in zip(visits,visits[1:]):
                exit=End(previous.section_id,int(previous.forward));entry=End(following.section_id,int(not following.forward))
                site=by_end.get(exit)
                if site is None or entry not in site.ends:
                    raise ValueError('A paired body span has no original consecutive source incidence')
            for curve in span.curves:
                points=np.asarray(curve.points if hasattr(curve,'points') else curve,float)
                if points.ndim!=2 or points.shape[1]!=2 or len(points)<2 or not np.isfinite(points).all():
                    raise ValueError('A paired body needs two actual finite original curb chains per role')
                curves.append(points)
        if len(identities)!=len(set(identities)):
            raise ValueError('A paired body repeats an original source Section across roles')
        return tuple(curves),tuple(identities),min(proofs)

    @staticmethod
    def _mouth(context,mouth,end):
        if not isinstance(mouth,PreparedMergedMouth) or mouth.approach.end_index!=end:
            raise ValueError('A paired source piece needs ordered actual head/tail prepared mouths')
        product=context.products.get(surface_key(mouth.owner))
        if (not isinstance(product,CompoundProduct) or product.owner!=mouth.owner
                or product.group.physical_layout.approaches[mouth.approach_index] is not mouth.approach):
            raise ValueError('A paired source mouth has no authoritative completed physical junction product')
        if product.pavement.mesh is None or (product.grass is not None and product.grass.mesh is None):
            raise ValueError('A paired source piece needs hydrated actual junction pavement and grass')
        actual=np.asarray([product.physical_rims[loop][vertex] for loop,vertex in mouth.rim_references])
        if not np.array_equal(actual,mouth.xyz):
            raise ValueError('A paired source mouth changed its completed physical junction rim')
        return product

    @staticmethod
    def _span_mouths(spans,mouths):
        for end,mouth in enumerate(mouths):
            aliases=mouth.approach.source_mouths(mouth.station)
            for role,(span,alias) in enumerate(zip(spans,aliases)):
                endpoint=span.head if end==0 else span.tail
                if endpoint.owner!=mouth.owner or endpoint.source_station is None or abs(endpoint.source_station-alias.station)>1e-8:
                    raise ValueError('A paired source span lost its authoritative physical mouth owner or original station')
                if alias.original_incidence:
                    if endpoint.endpoint!=alias.endpoint:
                        raise ValueError('A paired source span changed an original physical endpoint incidence')
                else:
                    key=endpoint.registry_key
                    if (endpoint.endpoint is not None or key is None or key.owner!=mouth.owner
                            or (key.section_id,key.end)!=(alias.endpoint.section_id,alias.endpoint.end)):
                        raise ValueError('A paired source span lost its owner-scoped derived alias registry')
                row=span.source_intervals[0 if end==0 else -1]
                station=row['first'] if (end==0)==row['forward'] else row['last']
                if row['section']!=alias.endpoint.section_id or abs(station-alias.station)>1e-7:
                    raise ValueError('A paired source interval does not end at its actual prepared mouth')

    def build(self,context,bundle,spans,mouths,*,frame,namespace):
        if len(bundle.members)!=2 or len(mouths)!=2:
            raise ValueError('A finite paired source body needs two original roles and two complete mouths')
        try:hash(namespace)
        except TypeError:raise ValueError('A paired body namespace must be a stable hashable assembly identity')
        products=[self._mouth(context,mouth,end) for end,mouth in enumerate(mouths)]
        if products[0].owner==products[1].owner:
            raise ValueError('A paired source body needs two distinct already owned physical junctions')
        from .paired_source_sampling import verify_shared_source_phase
        phase=verify_shared_source_phase(tuple(mouth.approach for mouth in mouths))
        curves,sections,progress=self._sources(context,bundle,spans,frame)
        self._span_mouths(spans,mouths)
        profile=mouths[0].approach.field
        if abs(profile.datum-mouths[1].approach.field.datum)+abs(profile.rise-mouths[1].approach.field.rise)>1e-9:
            raise ValueError('A paired source body needs its original shared crown and minimum datum')
        field=AngularCrown(frame,curves,datum=profile.datum,rise=profile.rise)
        visits=tuple(tuple(span.visits) for span in spans)
        bridges=[];cuts={};reports=[]
        for end,mouth in enumerate(mouths):
            inventory=SimpleNamespace(frame=frame,field=field,source_sections=sections,
                source_visits=visits,entrance=mouth.approach)
            bridge=MergedProfileBridge()
            mesh,report=bridge.build(context,mouth.approach,mouth.station,mouth.xyz,target_inventory=inventory)
            halves=mouth.split_body_cut(bridge.body_cut,crest_column=report.get('shared_crest_column'))
            for role,points in enumerate(halves):cuts[role,end]=points
            bridges.append((bridge,mesh));reports.append(report)
        body=self.make_body()
        paths=tuple(context[Key('path',sid)] for sid in sections)
        approach=SimpleNamespace(frame_chart=frame,field=field,paths=paths,sign=1.)
        pavement,grass,body_report=body.build_between(context,bundle,approach,cuts,
            namespace=(namespace,'source body'))
        total=len(pavement.quads)+sum(len(mesh.quads) for bridge,mesh in bridges)+(len(grass.quads) if grass is not None else 0)
        limit=self.max_cells*context.network.settings.resolution
        if total>limit:raise ValueError('Paired source connectors, pavement and grass exceed their total quad budget')
        assembly=Assembler();assembly.add(pavement.xyz,pavement.quads,body.pavement_keys)
        front_keys={}
        for end,(mouth,(bridge,mesh),report) in enumerate(zip(mouths,bridges,reports)):
            keys=[(namespace,'connector',end,int(i)) for i in range(len(mesh.xyz))]
            crest=report.get('shared_crest_column',mouth.crest_column);columns=np.arange(len(bridge.back_ids))
            first,last=columns[:crest+1],columns[crest:][::-1]
            columns_by_role=(first,last) if mouth.approach.end_index==0 else (last,first)
            for role,column_ids in enumerate(columns_by_role):
                if len(column_ids)!=len(body.cut_keys[role,end]):
                    raise ValueError('A paired source body dropped a complete connector cut partition')
                for index,(column,key) in enumerate(zip(column_ids,body.cut_keys[role,end])):
                    node=int(bridge.back_ids[int(column)])
                    if not np.array_equal(mesh.xyz[node],body.cuts[role,end][index]):
                        raise ValueError('A paired source connector changed its authoritative body cut')
                    keys[node]=key
            front_keys[end]=[keys[int(i)] for i in bridge.front_ids]
            assembly.add(mesh.xyz,mesh.quads,keys)
        complete=assembly.mesh()
        edges=np.sort(np.stack((complete.quads,np.roll(complete.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
        unique,counts=np.unique(edges,axis=0,return_counts=True)
        incidence={tuple(edge):int(count) for edge,count in zip(unique,counts)}
        shared=0
        for names in body.cut_keys.values():
            nodes=[assembly.lookup[key] for key in names]
            for a,b in zip(nodes,nodes[1:]):
                if incidence.get(tuple(sorted((a,b))))!=2:
                    raise ValueError('A paired source connector and body lack one identical whole shared edge')
                shared+=1
        self.pavement_keys=[None]*len(complete.xyz)
        for key,at in assembly.lookup.items():self.pavement_keys[at]=key
        self.front_ids={end:np.asarray([assembly.lookup[key] for key in names],np.int32) for end,names in front_keys.items()}
        self.front_keys=front_keys
        for end,mouth in enumerate(mouths):
            if not np.array_equal(complete.xyz[self.front_ids[end]],mouth.xyz):
                raise ValueError('A paired source body moved its actual physical junction interface')
        self.grass_keys=body.grass_keys;self.body=body;self.field=field;self.frame=frame
        payloads={Key('paired_source_pavement',0):complete}
        if grass is not None:payloads[Key('paired_source_grass',0)]=grass
        checks=[]
        for mesh in payloads.values():
            double=quality(mesh);mesh.validate()
            stored=Mesh(mesh.xyz.astype(np.float32).astype(float),mesh.quads)
            row=quality(stored);stored.validate()
            if (max(row['maximum_edge_ratio'],double['maximum_edge_ratio'])>8.
                    or (mesh.material != 'grass'
                        and min(row['minimum_corner_sine'],double['minimum_corner_sine'])<.03)):
                raise ValueError('A paired source piece fails stored-coordinate cell quality')
            slopes=stored.slope_statistics()['maximum'];along=frame.longitudinal_bound(stored)
            if slopes>.1+1e-6 or (mesh.material != 'grass' and along>.03+1e-8):
                raise ValueError('A paired source piece exceeds physical slope bounds')
            fixed=np.zeros(len(mesh.xyz),bool);fixed[np.unique(MergedRunBody._boundary_vertices(mesh))]=True
            drainage=BoundaryDrainage(max_cells=limit).fill(stored,fixed)[1]
            if drainage['maximum_fill']>1e-6:raise ValueError('A paired source piece creates an interior drainage depression')
            checks.append(dict(material=mesh.material,quality=row,double_quality=double,
                minimum_corner_sine_limit=.03 if mesh.material != 'grass' else None,
                maximum_slope=slopes,longitudinal_bound=along,
                carriageway_longitudinal_limit=.03 if mesh.material != 'grass' else None,drainage=drainage))
        containment=certify_mesh(complete,curves,source_error=.1*context.network.settings.curve_tolerance,
            mesh_error=.6*context.network.settings.curve_tolerance,
            maximum_error=.7*context.network.settings.curve_tolerance+1e-7,max_cells=limit)
        for product in products:
            payloads[Key('prepared_junction_pavement',product.owner.id)]=product.pavement.mesh
            if product.grass is not None:payloads[Key('prepared_junction_grass',product.owner.id)]=product.grass.mesh
        overlap=validate_overlaps(payloads,cell_size=2.)
        self.cut_keys=body.cut_keys;self.open_grass_cuts=body.open_grass_cuts
        return complete,grass,dict(construction='finite paired original source piece between actual prepared junctions',
            source_sections=list(sections),source_intervals=[list(span.source_intervals) for span in spans],
            original_run_progress_lower=progress,source_containment=containment,
            shared_source_phase=phase,
            connectors=reports,body=body_report,stored_checks=checks,overlap=overlap,total_quads=total,
            exact_physical_junction_mouths=True,exact_internal_body_seams=True,
            shared_connector_body_edges=shared,junction_area_remeshed=False,all_direct_quads=True)
