"""Extensible rendering orchestration over a single prepared street snapshot."""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import json
import bpy

from ...app.asset_library import BlenderAssetLoader
from . import assets
from .adapter import PreparedNetworkAdapter, StreetDesign
from .authoring import EditableNetworkBuilder, SourceMeshBuilder
from ..gn.editable_input import model, runtime
from .transaction import ImportTransaction
from .coverage import CoverageReport
from .boundary import UnresolvedBoundaryPolicy
from .evidence import GeometryEvidence
import time


@dataclass
class RenderContext:
    app: object
    snapshot: object
    collection: object
    design: StreetDesign
    loader: BlenderAssetLoader
    outputs: list = field(default_factory=list)
    results: list = field(default_factory=list)
    timings: dict = field(default_factory=dict)


class StreetFeature(ABC):
    """Future sidewalks/markings/track providers share topology and asset services."""
    required_libraries = ()
    ground_passages = False
    @abstractmethod
    def render(self, context):
        pass


class CarriagewayFeature(StreetFeature):
    required_libraries = (assets.LIBRARY,)
    def __init__(self, *, junction_ids=None, street_ids=(), boundary_clearance=0., open_bundle_clearance=0.):
        # Optional qualification/selection boundary, never sample-specific logic.
        self.junction_ids = None if junction_ids is None else frozenset(junction_ids)
        self.street_ids = frozenset(street_ids)
        if not 0. <= boundary_clearance < float('inf'):
            raise ValueError('Selection boundary clearance must be finite and nonnegative')
        self.boundary_clearance = boundary_clearance
        if not 0. <= open_bundle_clearance < float('inf'):
            raise ValueError('Open Bundle boundary clearance must be finite and nonnegative')
        # Explicit selection/inspection policy only. Ordinary free endpoints
        # retain their source positions unless a caller requests this cut.
        self.open_bundle_clearance = open_bundle_clearance

    def render(self, context):
        load_started = time.perf_counter()
        loaded = context.loader.load_many((assets.SOLVER, assets.PAVEMENT, assets.GRASS))
        group = loaded[assets.SOLVER]
        assets.validate_solver(group)
        context.timings['asset_loading_seconds'] = time.perf_counter()-load_started
        builder = EditableNetworkBuilder(context.design, dict(solver=group,
            pavement=loaded[assets.PAVEMENT], grass=loaded[assets.GRASS]))
        snapshot = context.snapshot
        selected = self.junction_ids
        if selected is not None:
            missing = selected - snapshot.junctions.keys()
            if missing:
                raise ValueError('Selected intersections are outside the supported domain: '+str(sorted(missing)))
            wanted = self.street_ids | {s.source_id for jid in selected for s in snapshot.junctions[jid].streets}
            unknown = wanted - snapshot.streets.keys()
            if unknown:
                raise ValueError('Unknown selected Streets: '+str(sorted(unknown)))
        else:
            wanted = set(snapshot.streets)
        if self.open_bundle_clearance:
            for sid in wanted:
                street = snapshot.streets[sid]
                if street['bundle_id'] is None:
                    continue
                sections = [s for s in street['items'] if s['kind']=='Section']
                for end, jid in enumerate(street['endpoints']):
                    if sections and jid is None:
                        section = sections[-1 if end else 0]
                        section['clip_end' if end else 'clip_start'] = self.open_bundle_clearance
                        snapshot.issue('Street', sid, 'EXPLICIT_OPEN_BUNDLE_CLIP',
                            f'Endpoint {end} explicitly clipped by {self.open_bundle_clearance:g} m for this selection')
        if selected is not None:
            for sid in wanted:
                street = snapshot.streets[sid]
                sections = [s for s in street['items'] if s['kind']=='Section']
                for end, jid in enumerate(street['endpoints']):
                    if sections and jid is not None and jid not in selected:
                        section = sections[-1 if end else 0]
                        section['clip_end' if end else 'clip_start'] = self.boundary_clearance
        # Keep complete eligible Bundle ownership; unsupported junctions receive
        # an explicit excluded area instead of removing every connected Street.
        candidates = {sid for sid, s in snapshot.streets.items() if s['status']=='SUPPORTED' and sid in wanted}
        changed = True
        while changed:
            rejected = set()
            for sid in candidates:
                street = snapshot.streets[sid]
                if street['bundle_id'] is not None:
                    bundle = snapshot.bundles[street['bundle_id']]
                    if bundle.get('status') != 'SUPPORTED' or not set(bundle['head']) <= candidates:
                        rejected.add(sid)
            changed = bool(rejected)
            candidates -= rejected
        adjacency = {sid:set() for sid in candidates}
        for bundle in snapshot.bundles.values():
            members = set(bundle['head'])
            if members <= candidates:
                for sid in members:
                    adjacency[sid].update(members-{sid})
        active_junctions = {jid:j for jid,j in snapshot.junctions.items()
            if (selected is None or jid in selected) and {s.source_id for s in j.streets} <= candidates}
        if selected is None:
            UnresolvedBoundaryPolicy().apply(snapshot, candidates, active_junctions)
        for junction in active_junctions.values():
            members = {s.source_id for s in junction.streets}
            if members <= candidates:
                for sid in members:
                    adjacency[sid].update(members-{sid})
        remaining = set(candidates)
        components = []
        while remaining:
            component, todo = set(), [min(remaining)]
            while todo:
                sid = todo.pop()
                if sid not in component:
                    component.add(sid); todo.extend(adjacency[sid]-component)
            remaining -= component
            components.append(component)
        # Compiling later bindings dirties Blender's dependency graph. Build the
        # largest network last instead of reevaluating it for every tiny island.
        for component in sorted(components, key=lambda c:(len(c),min(c))):
            junction_ids = [jid for jid,j in active_junctions.items()
                            if {s.source_id for s in j.streets} <= component]
            try:
                started = time.perf_counter()
                with ImportTransaction():
                    obj = builder.create('Streets / network '+str(min(component)), snapshot,
                                         component, junction_ids, context.collection)
                context.outputs.append(obj)
                context.results.append(dict(streets=sorted(component), intersections=junction_ids,
                    status=obj.st_network.status, reason=obj.st_network.diagnostic,
                    evidence=GeometryEvidence().read(obj),
                    transitions={r.transition_id:r.street_id for r in obj.st_network.transitions},
                    native_timings=runtime._metrics.get(obj.name, {}).copy(),
                    seconds=time.perf_counter()-started))
            except Exception as exc:
                context.results.append(dict(streets=sorted(component), intersections=junction_ids,
                                            status='FAILED', reason=str(exc)))
        for sid in set(snapshot.streets)-candidates:
            street = snapshot.streets[sid]
            excluded = selected is not None and sid not in wanted
            deferred = street['reason'] and street['reason'].startswith('Reserved component/category:')
            street['surface_status'] = 'DELIBERATELY_EXCLUDED' if excluded or deferred else 'UNSUPPORTED'
            street['surface_reason'] = ('Outside explicit qualification selection' if excluded else
                street['reason'] or 'Connected network contains an unsupported required mouth or Bundle')
        for jid, row in snapshot.intersections.items():
            row['surface_status'] = 'DELIBERATELY_EXCLUDED' if selected is not None and jid not in selected else 'UNSUPPORTED'
        for result in context.results:
            evidence = result.get('evidence', {})
            present = set(evidence.get('section_ids', ()))
            present_junctions = set(evidence.get('intersection_ids', ()))
            present_transitions = set(evidence.get('transition_ids', ()))
            diagnostic = evidence.get('diagnostics', ())
            for sid in result['streets']:
                street = snapshot.streets[sid]
                required = {s['id'] for s in street['items'] if s['kind']=='Section'}
                errors = [d for d in diagnostic if d['section_id'] in required
                    or street['bundle_id'] is not None and d['bundle_id']==street['bundle_id']
                    or d['section_id'] < 0 and d['intersection_id'] < 0 and d['bundle_id'] < 0]
                attached = {j for j in street['endpoints'] if j in active_junctions}
                attached.update(i['id'] for i in street['items'] if i.get('surface_treatment')=='NATIVE_JUNCTION' and i['id'] in active_junctions)
                required_transitions = {tid for tid,owner in result.get('transitions',{}).items() if owner==sid}
                expected_transitions = sum(bool(i.get('transition')) for i in street['items'])
                complete = (required <= present and attached <= present_junctions
                    and len(required_transitions)==expected_transitions and required_transitions <= present_transitions
                    and not errors and not street.get('unresolved_boundaries'))
                street['surface_status'] = ('FAILED' if result['status']=='FAILED' else 'GENERATED' if complete else 'PARTIAL')
                street['surface_reason'] = ('Surface stops before unresolved junctions; see unresolved_boundaries'
                    if street.get('unresolved_boundaries') else None if complete else result.get('reason'))
                street['native_diagnostics'] = errors
            for jid in result['intersections']:
                snapshot.intersections[jid]['surface_status'] = 'GENERATED' if jid in present_junctions else result['status']


class GeometryNodesStreetRenderer:
    def __init__(self, app, *, design=None, features=None, asset_loader=None):
        self.app = app
        self.design = design or StreetDesign.from_app(app)
        if features is None and getattr(app,'streetConstruction','current')=='original_spans':
            from .original_span_feature import OriginalSpanCarriagewayFeature
            features=(OriginalSpanCarriagewayFeature(),)
        self.features = tuple(features) if features is not None else (CarriagewayFeature(),)
        self.loader = asset_loader or BlenderAssetLoader.for_app(app,
            assets.catalog(getattr(app, 'streetAssetFile', None)))
        self.adapter = PreparedNetworkAdapter(self.design,
            ground_passages=any(getattr(feature, 'ground_passages', False) for feature in self.features))
        self.context = None

    def prepare(self):
        if bpy.app.version < (5, 2, 0):
            raise ValueError('The experimental street Geometry Nodes backend requires Blender 5.2 or newer')
        model.register()
        # Validate library availability before creating output collections.
        for feature in self.features:
            for library in feature.required_libraries:
                self.loader.resolve(library.key)

    def render(self, manager, data):
        with ImportTransaction():
            return self._render(manager, data)

    def _render(self, manager, data):
        adaptation_started = time.perf_counter()
        snapshot = self.adapter.read(manager)
        adaptation_seconds = time.perf_counter()-adaptation_started
        collection = bpy.data.collections.new('Blosm / Geometry Nodes streets')
        bpy.context.scene.collection.children.link(collection)
        self.context = context = RenderContext(self.app, snapshot, collection, self.design, self.loader)
        context.timings['adaptation_seconds'] = adaptation_seconds
        collection['blosm_street_input_schema'] = '0.1'
        collection['blosm_projection_origin'] = (float(getattr(data,'lat',0.)), float(getattr(data,'lon',0.)))
        collection['blosm_input_file'] = str(getattr(self.app,'osmFilepath',''))
        collection['blosm_street_asset_revision'] = ', '.join(dict.fromkeys(
            library.key+':'+library.revision for feature in self.features for library in feature.required_libraries))
        collection['blosm_connectivity_review_events'] = len(snapshot.preparation_audit.get('events', ()))
        for feature in self.features:
            feature.render(context)
        generated = {sid for row in context.results if row['status'] != 'FAILED' for sid in row['streets']}
        inspection = bpy.data.collections.new('Streets / unsupported source context')
        collection.children.link(inspection)
        meshes = SourceMeshBuilder()
        for sid, street in snapshot.streets.items():
            if sid in generated:
                continue
            sections = [i for i in street['items'] if i['kind']=='Section']
            if not sections:
                continue
            obj = meshes.create('Unresolved / Street '+str(sid), sections, inspection,
                                interpolation=self.design.interpolation)
            obj['blosm_street_status'] = street['reason'] or 'Connected network contains an unsupported mouth'
            obj.color = (1., .15, .06, 1.)
        report = bpy.data.texts.new('Blosm street import report')
        report.write(json.dumps(dict(network=snapshot.report(), surfaces=context.results,
            coverage=CoverageReport().build(snapshot), timings=context.timings),indent=2))
        collection['blosm_street_report'] = report.name
        collection['blosm_street_coverage'] = 'PARTIAL' if any(s.get('surface_status') in ('UNSUPPORTED','PARTIAL','FAILED') for s in snapshot.streets.values()) or any(r['status']!='READY' for r in context.results) else 'COMPLETE'
        if getattr(self.app, 'op', None):
            ready = sum(r['status']=='READY' for r in context.results)
            self.app.op.report({'INFO' if collection['blosm_street_coverage']=='COMPLETE' else 'WARNING'},
                f"Street import: {ready} ready networks. Coverage: {collection['blosm_street_coverage']}. See Text Editor: {report.name}")
        return context

    def finalize(self):
        pass

    def cleanup(self):
        # Keep authored sources, registry and native lookup bindings after import.
        pass
