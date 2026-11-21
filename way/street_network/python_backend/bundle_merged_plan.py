"""Complete merged-entrance Bundle output in the normal preparation graph."""
from dataclasses import replace
import numpy as np
from .model import Key
from .updates import BuildContext
from .ownership import PLAN, CompoundProduct, StandaloneProduct, JunctionOwnership, surface_key
from .bundle_merged_entrances import PairedAngularEntrance
from .bundle_merged_bridge import MergedProfileBridge
from .bundle_merged_runs import MergedRunWithExternalCap
from .mesh import Assembler
from .providers import Part
from .surface_query import QuadHeightField


def required_physical_owners(context, bundle, specification):
    """Resolve only original source incidences and explicitly prepared owners."""
    plan = context.products.get(PLAN)
    if plan is None:
        raise ValueError('Merged Bundle physical ownership requires the current source plan')
    sections = set(bundle.sections)
    sites = {site.key for site in context.network.intersections
             if any(end.section_id in sections for end in site.ends)}
    for prepared in (specification, *getattr(specification, 'paired_sources', ())):
        sites.update(prepared.group.members)
    return tuple(sorted({plan.owner(site) or site for site in sites}))


def physical_context(context, *, required_owners=None):
    """Reuse prepared junctions and their exact physical rims, even at E0.

    Released geometry is rebuilt only for explicitly required physical owners,
    against the current source plan. Accepted cached products remain untouched.
    """
    products = dict(context.products)
    physical = BuildContext(replace(context.network, settings=replace(context.network.settings, exaggeration=1.)),
                            products, context.patterns, context.task_digests, context.source_digests)
    plan = products.get(PLAN)
    if required_owners is None:
        required_owners = {Key('intersection', key.id) for key, product in products.items()
                           if key.kind == 'junction_surface' and isinstance(product, (CompoundProduct, StandaloneProduct))}
    declared = {site.key for site in context.network.intersections}
    owners = set()
    for site in required_owners:
        if site not in declared:
            raise ValueError('Required physical geometry has no original junction incidence')
        owners.add((plan.owner(site) if plan is not None else None) or site)
    hydrated = []
    for owner in sorted(owners):
        key = surface_key(owner)
        product = products.get(key)
        if isinstance(product, CompoundProduct):
            group = plan.by_site.get(owner) if plan is not None else None
            if (product.owner != owner or (group is not None and product.group.members != group.members)):
                raise ValueError('Required physical geometry disagrees with the current source ownership plan')
            missing = (getattr(product, 'physical_pavement', None) is None and product.pavement.mesh is None)
            missing |= (product.grass is not None and getattr(product, 'physical_grass', None) is None
                        and product.grass.mesh is None)
            # A display mesh at E0/E10 is not an authoritative physical buffer.
            missing |= product.exaggeration != 1. and (getattr(product, 'physical_pavement', None) is None
                or (product.grass is not None and getattr(product, 'physical_grass', None) is None))
        else:
            missing = (not isinstance(product, StandaloneProduct) or product.pavement.mesh is None
                       or product.exaggeration != 1.)
        if missing:
            if plan is None:
                raise ValueError('Released junction geometry needs the current source ownership plan')
            group = plan.by_site.get(owner)
            product = (CompoundProduct.build(physical, owner, group) if group is not None else
                       JunctionOwnership.surface(physical, owner))
            if product is None:
                raise ValueError('Required physical owner has no normal prepared junction surface')
            products[key] = product
            hydrated.append(owner)
        if isinstance(product, CompoundProduct):
            if product.physical_rims is None:
                raise ValueError('Merged Bundle cap requires its prepared physical rim')
            pavement=getattr(product,'physical_pavement',None)
            grass=getattr(product,'physical_grass',None)
            if product.exaggeration == 1.:
                pavement=pavement if pavement is not None else product.pavement.mesh
                if product.grass is not None:grass=grass if grass is not None else product.grass.mesh
            if pavement is None or (product.grass is not None and grass is None):
                raise ValueError('Merged Bundle needs hydrated authoritative physical junction meshes')
            products[key] = replace(product, rims=product.physical_rims, exaggeration=1.,
                pavement=replace(product.pavement,mesh=pavement),
                grass=None if product.grass is None else replace(product.grass,mesh=grass))
    physical.physical_owners = tuple(sorted(owners))
    physical.hydrated_physical_owners = tuple(hydrated)
    return physical


class MergedBundlePlan:
    mode = 'union'
    transition_owners = ()
    requires_separator_query = True

    def __init__(self, context, bundle, specification):
        if getattr(specification,'source_phase_proof',None) is not None:
            self._build_hydrated(context,bundle,specification)
            return
        physical = physical_context(context, required_owners=required_physical_owners(context,bundle,specification))
        owner = min(specification.crossing.sites)
        product = physical[surface_key(owner)]
        if not isinstance(product, CompoundProduct):
            raise ValueError('Merged Bundle needs its completed shared physical junction')
        self.bundle, self.settings = bundle, context.network.settings
        self.source_ids = tuple(visit.section_id for member in bundle.members for visit in member.path)
        self.datum = product.field.datum
        paving, grass, reports = Assembler(), Assembler(), []
        for index, approach in enumerate(specification.layout.approaches):
            if not isinstance(approach, PairedAngularEntrance):
                continue
            station = specification.layout.outward_stations[index]
            identities = {(end.section_id, end.end) for end in approach.original_ends}
            refs = sorted({(loop, node) for (sid, end, u), (loop, node, at) in product.port_points.items()
                           if (sid, end) in identities})
            mouth = np.asarray([product.rims[loop][node] for loop, node in refs])
            if len(mouth) < 3:
                raise ValueError('Merged Bundle lost its complete physical mouth')
            first = approach.point(station, -.5*approach.width)
            last = approach.point(station, .5*approach.width)
            mouth = mouth[np.argsort((mouth[:, :2]-first)@(last-first))]
            bridge = MergedProfileBridge()
            connector, connector_report = bridge.build(physical, approach, station, mouth)
            run = MergedRunWithExternalCap()
            pavement, separator, run_report = run.build(physical, bundle, approach, bridge.body_cut)
            joined = Assembler()
            keys = [('connector', i) for i in range(len(connector.xyz))]
            for row, node in enumerate(bridge.back_ids):
                keys[node] = ('merged cut', row)
            joined.add(connector.xyz, connector.quads, keys)
            joined.add(pavement.xyz, pavement.quads, run.pavement_keys)
            mesh = joined.mesh()
            paving.add(mesh.xyz, mesh.quads)
            grass.add(separator.xyz, separator.quads)
            reports.append(dict(source_sections=list(run.domain.run.source_sections),
                                connector=connector_report, body=run_report))
        if len(reports) != 2:
            raise ValueError('Merged Bundle needs its two complete source-owned runs')
        physical_paving, physical_grass = paving.mesh(), grass.mesh('grass')
        total = len(physical_paving.quads)+len(physical_grass.quads)
        if total > 16000*self.settings.resolution:
            raise ValueError('Merged Bundle exceeds its complete finite cell budget')
        query = Assembler()
        query.add(physical_paving.xyz, physical_paving.quads)
        query.add(physical_grass.xyz, physical_grass.quads)
        self.query = QuadHeightField(query.mesh(), datum=self.datum, ceiling=product.field.ceiling)
        for mesh in (physical_paving, physical_grass):
            mesh.xyz[:, 2] = self.datum+self.settings.exaggeration*(mesh.xyz[:, 2]-self.datum)
        metadata = dict(construction='source-owned merged Bundle entrances and complete rounded runs',
            source_sections=list(self.source_ids), member_junctions=[key.text() for key in specification.crossing.sites],
            consumed_sections=list(specification.crossing.sections), runs=reports,
            centerline_proof_checks=specification.centerline_proof_checks,
            complete_quads=total, all_direct_quads=True, default_provider_integration=True)
        self._pavement = Part(bundle.key, physical_paving, dict(metadata, quality=physical_paving.validate()))
        self._separator = Part(bundle.key, physical_grass, dict(metadata, quality=physical_grass.validate()))

    def _build_hydrated(self,context,bundle,specification):
        """Consume jointly prepared source leaves without remeshing junctions."""
        from .hydrated_run_assembly import HydratedRunAssembly,interval_inventory
        from .merged_run_inventory import inventory_frame,hydrate_evaluation_spans
        from .validation import validate_overlaps
        physical=physical_context(context,required_owners=required_physical_owners(context,bundle,specification))
        specifications=tuple(getattr(specification,'paired_sources',()))
        self.bundle,self.settings=bundle,context.network.settings
        self.source_ids=tuple(visit.section_id for member in bundle.members for visit in member.path)
        paving,grassy,query=Assembler(),Assembler(),Assembler()
        reports,assemblies,paired_owners=[],[],set()
        paths_used=[]
        payloads={}
        query_fields=[]
        for number,approach in enumerate(specification.layout.approaches):
            if not isinstance(approach,PairedAngularEntrance):continue
            scope=tuple(tuple((visit.section_id,visit.forward,0.,physical[Key('path',visit.section_id)].length)
                              for visit in member.visits) for member in approach.run.members)
            from .bundle_mouth_frames import candidate_frames
            frame_failures=[]
            for frame,frame_proof in candidate_frames(physical,bundle,approach,specification):
                try:
                    spans=hydrate_evaluation_spans(physical,bundle,approach,specifications,frame=frame)
                    expected=interval_inventory(spans)
                    assembly=HydratedRunAssembly()
                    pavement,grass,report=assembly.build(physical,bundle,spans,frame=frame,
                        namespace=(bundle.key,'hydrated run',number),expected_intervals=expected,scope_intervals=scope)
                except ValueError as error:
                    frame_failures.append(dict(frame=frame_proof,reason=str(error)))
                    continue
                if frame_proof['family']!='established':
                    report=dict(report,source_coordinate_selection=dict(frame_proof,prior_failures=frame_failures))
                break
            else:
                raise ValueError('No qualified source-certified frame for the complete hydrated run: '+str(frame_failures))
            paving.add(pavement.xyz,pavement.quads,assembly.pavement_keys)
            if grass is not None:grassy.add(grass.xyz,grass.quads,assembly.grass_keys)
            query.add(assembly.body_query_mesh.xyz,assembly.body_query_mesh.quads)
            payloads[Key('hydrated_bundle_pavement',number)]=pavement
            if grass is not None:payloads[Key('hydrated_bundle_grass',number)]=grass
            paired_owners.update(assembly.paired_owners)
            paths_used.extend(visit.section_id for member in approach.run.members for visit in member.visits)
            reports.append(report);assemblies.append(assembly)
        expected_runs=getattr(specification,'expected_run_count',2)
        if (expected_runs not in (1,2) or len(reports)!=expected_runs
                or len(paths_used)!=len(set(paths_used)) or set(paths_used)!=set(self.source_ids)):
            raise ValueError('Hydrated Bundle runs do not retain every original member Section exactly once')
        physical_pavement=paving.mesh()
        physical_grass=grassy.mesh('grass') if grassy.quads else None
        total=len(physical_pavement.quads)+(len(physical_grass.quads) if physical_grass is not None else 0)
        if total>50000*self.settings.resolution:
            raise ValueError('Complete hydrated Bundle exceeds its finite cell budget')
        for owner in sorted(paired_owners):
            product=physical[surface_key(owner)]
            query_fields.append(product.field)
            for material,part in (('pavement',product.pavement),('grass',product.grass)):
                if part is None:continue
                if part.mesh is None:raise ValueError('Hydrated Bundle query lost an accepted physical junction mesh')
                query.add(part.mesh.xyz,part.mesh.quads)
                payloads[Key('owned_bundle_junction_'+material,owner.id)]=part.mesh
        overlap=validate_overlaps(payloads,cell_size=2.)
        root_owner=(specification.root_owner if hasattr(specification,'root_owner')
                    else min(specification.crossing.sites))
        first=physical[surface_key(root_owner)].field
        self.datum=first.datum
        if any(abs(field.datum-self.datum)>1e-9 for field in query_fields):
            raise ValueError('Hydrated Bundle junction queries need their original common minimum height')
        ceiling=max(field.ceiling for field in query_fields)
        self.query=QuadHeightField(query.mesh(),datum=first.datum,ceiling=ceiling)
        self.transition_owners=tuple(sorted({owner for assembly in assemblies for owner in assembly.transition_owners}))
        for mesh in (physical_pavement,physical_grass):
            if mesh is not None:mesh.xyz[:,2]=self.datum+self.settings.exaggeration*(mesh.xyz[:,2]-self.datum)
        metadata=dict(construction='jointly prepared original-source Bundle runs and rounded physical junctions',
            source_sections=list(self.source_ids),member_junctions=[key.text() for key in sorted(paired_owners)],
            consumed_sections=list(specification.consumed_sections if hasattr(specification,'consumed_sections')
                                   else specification.crossing.sections),runs=reports,
            complete_quads=total,all_direct_quads=True,default_provider_integration=True,
            shared_source_phase=specification.source_phase_proof,complete_selected_source_coverage=True,
            query_only_junction_owners=[key.text() for key in sorted(paired_owners)],
            junction_geometry_republished=False,overlap=overlap)
        self._pavement=Part(bundle.key,physical_pavement,dict(metadata,quality=physical_pavement.validate()))
        self._separator=(Part(bundle.key,physical_grass,dict(metadata,quality=physical_grass.validate()))
                         if physical_grass is not None else None)

    def pavement(self, context):
        return self._pavement

    def separator(self, context):
        return self._separator

    def height(self, xy):
        xy = np.asarray(xy)
        shape = xy.shape[:-1]
        value = self.query.evaluate(xy.reshape(-1, 2), default=self.datum)
        return (self.datum+self.settings.exaggeration*(value-self.datum)).reshape(shape)

    def release_geometry(self):
        for part in (self._pavement, self._separator):
            if part is not None and part.mesh is not None:
                part.metadata['bounds'] = part.mesh.bounds
                part.mesh = None
