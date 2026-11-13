"""Normal Blosm import adapter for the direct numerical mesh backend."""
import json
import time
from types import SimpleNamespace
from uuid import uuid4
import bpy
from ....app.asset_library import BlenderAssetLoader
from .. import assets
from ..backends import StreetRenderBackend
from ..adapter import PreparedNetworkAdapter
from ..design import StreetDesign
from ..physical_sources import PhysicalStreetInventory
from ..physical_ports import GroundSurfacePolicy
from ...way_properties import wayCategoryProps
from ..transaction import ImportTransaction
from .adapter import PreparedPythonAdapter, components
from .model import Settings
from .scene import NetworkScene, StreetInputChanged
from .source_support import ClippedApproachSupport
from . import ui
from ....progressive.core import YieldToUI
from ....progressive.ownership import alive


class PythonStreetRenderer(StreetRenderBackend):
    def __init__(self, app, *, design=None, features=None, asset_loader=None):
        self.app = app
        self.design = design or StreetDesign.from_app(app)
        self.loader = asset_loader or BlenderAssetLoader.for_app(app,
            assets.material_catalog(getattr(app, 'streetMaterialFile', None)))
        self.adapter = PreparedNetworkAdapter(self.design, ground_passages=True)
        self.converter = PreparedPythonAdapter(self.design.categories, GroundSurfacePolicy.reason, wayCategoryProps)
        self.context = None
        self.features = None if features is None else tuple(features)

    def prepare(self):
        if bpy.app.version < (5, 2, 0):
            raise ValueError('The Python street backend targets Blender 5.2 or newer')
        ui.register()
        self.loader.resolve(assets.MATERIAL_LIBRARY.key)

    def render(self, manager, data):
        with ImportTransaction():
            return self._drain(self._render_steps(manager, data))

    @staticmethod
    def _linked(parent, target):
        return parent == target or any(PythonStreetRenderer._linked(child, target) for child in parent.children)

    @staticmethod
    def _drain(producer):
        while True:
            try:
                next(producer)
            except StopIteration as done:
                return done.value

    @staticmethod
    def _drive(producer, session):
        # Numerical steps allocate no Blender data; the scene publisher
        # explicitly forwards each allocation to the session ownership ledger.
        def advance():
            try:
                return False, next(producer)
            except StopIteration as done:
                return True, done.value
        try:
            while True:
                done, value = advance()
                if done:
                    return value
                yield value
        finally:
            producer.close()

    def render_steps(self, manager, data, session):
        try:
            return (yield from self._render_steps(manager, data, session))
        finally:
            # Discard process-local engines for any pending component that the
            # session will remove; accepted controllers keep their warm cache.
            from .scene import _ENGINES
            for kind, item in tuple(session.ownership.pending.values()):
                if kind == 'objects' and alive(item) and item.get('st_py_network_uuid'):
                    _ENGINES.pop(item['st_py_network_uuid'], None)
            # Accepted reports survive cancellation alongside accepted networks.
            # Never allocate new report data from teardown/history handlers.
            if self.context is not None and self.context.report.get('coverage') == 'IN_PROGRESS':
                text = self.context.text
                collection = self.context.collection
                if alive(text) and alive(collection):
                    self.context.report['coverage'] = 'PARTIAL'
                    self.context.report['interrupted'] = True
                    text.clear()
                    text.write(json.dumps(self.context.report, indent=2))
                    collection['blosm_street_coverage'] = 'PARTIAL'
                    session.revision += 1

    def load_materials(self):
        """Explicit import override, with the existing library defaults otherwise."""
        pavement = getattr(self.app, 'streetCarriagewayMaterial', None)
        refs = (assets.GRASS,) if pavement is not None else (assets.PAVEMENT, assets.GRASS)
        loaded = self.loader.load_many(refs)
        return {'pavement': pavement if pavement is not None else loaded[assets.PAVEMENT],
                'grass': loaded[assets.GRASS]}

    def _render_steps(self, manager, data, session=None):
        def call(callback):
            return session.call(callback)[0] if session is not None else callback()

        def progress(label):
            if session is not None:
                session.phase = label

        start = time.perf_counter()
        progress('Preparing street network snapshot')
        yield None
        snapshot = self.adapter.read(manager)
        yield None
        physical = PhysicalStreetInventory.view(snapshot)
        settings = Settings(self.design.cell_size, self.design.resolution, self.design.exaggeration)
        adaptation = self.converter.read(physical, str(uuid4()), settings,
            rise=self.design.profile_rise, interpolation=self.design.interpolation)
        yield None
        materials = call(self.load_materials)
        def create_collection():
            collection = bpy.data.collections.new('Blosm / Python streets')
            bpy.context.scene.collection.children.link(collection)
            return collection
        collection = call(create_collection)
        def validate():
            if session is not None:
                session.validate()
                if not alive(collection) or not self._linked(session.scene.collection, collection):
                    raise StreetInputChanged('The street import collection was deleted or unlinked')
        results = []
        support = ClippedApproachSupport(physical, snapshot.clipped_source_context)
        for number, network in enumerate(components(adaptation.network)):
            progress(f'Preparing street component {number + 1}: boundary support')
            yield YieldToUI() if session is not None else None
            preparation = yield from support.prepare_steps(network, on_progress=progress)
            yield None
            network = preparation.network
            controller = call(lambda: NetworkScene.create(network, collection, name=f'Python streets / Network {number}',
                                             materials=materials, generate=False, features=self.features))
            from .terrain.blender import configure
            configure(controller, self.app)
            controller['st_py_source_support'] = json.dumps(preparation.report, separators=(',', ':'), allow_nan=False)
            ids = {s.id for s in network.sections}
            result = dict(network_uuid=network.uuid, sections=sorted(ids),
                          streets=sorted({s.street_id for s in network.sections}), controller=controller.name,
                          source_support=preparation.report)
            try:
                if ids.intersection(adaptation.blocked_sections):
                    raise ValueError('Required source context or Bundle membership is unsupported; see preparation issues')
                if preparation.error:
                    raise ValueError(preparation.error)
                if session is None:
                    update = NetworkScene.update(controller)
                else:
                    update = yield from self._drive(NetworkScene.update_steps(controller,
                        validate=validate, on_progress=progress, parent_ownership=session.ownership), session)
                result.update(status='READY', quads=sum(len(o.data.polygons)
                    for o in controller['st_py_outputs'].objects if o.type == 'MESH'), timings=update.timings)
                result['fitted_corner_radii'] = [fit for obj in controller['st_py_outputs'].objects
                    if obj.get('st_py_metadata') for fit in
                    json.loads(obj['st_py_metadata']).get('fitted_corner_radii', ())]
                if session is not None:
                    for warning in update.timings.get('publication_cleanup_errors', ()):
                        session.report({'WARNING'}, 'Street cleanup: ' + warning)
            except StreetInputChanged:
                raise
            except (ValueError, RuntimeError) as exc:
                result.update(status='FAILED', reason=str(exc))
                if alive(controller):
                    controller['st_py_error'] = str(exc)
                if session is not None:
                    session.report({'ERROR'}, 'Street component %s: %s' % (number + 1, exc))
            results.append(result)
            if session is not None:
                validate()
                call(lambda: self._write_report(snapshot, adaptation, collection, results, start, complete=False))
                session._write_origin()
                session.ownership.accept()
                session.revision += 1
                session.work_done += 1
                if result['status'] == 'READY':
                    session.committed_units += 1
                    if session.first_result is None:
                        session.first_result = time.monotonic() - session.started
                yield YieldToUI()
        context = call(lambda: self._write_report(snapshot, adaptation, collection, results, start, complete=True))
        report = context.report
        text = context.text
        if getattr(self.app, 'op', None):
            ready = sum(r['status'] == 'READY' for r in results)
            self.app.op.report({'INFO' if report['coverage'] == 'COMPLETE' else 'WARNING'},
                f'Python streets: {ready}/{len(results)} networks ready. Coverage: {report["coverage"]}. See {text.name}')
        return context

    def _write_report(self, snapshot, adaptation, collection, results, start, *, complete):
        report = dict(backend='python', version=1, source=snapshot.report(),
            preparation_issues=adaptation.issues, surfaces=results,
            coverage='IN_PROGRESS' if not complete else 'COMPLETE' if results and all(r['status'] == 'READY' for r in results)
                     and not adaptation.issues else 'PARTIAL', seconds=time.perf_counter() - start)
        text = self.context.text if self.context is not None and alive(self.context.text) else bpy.data.texts.new('Blosm Python street import report')
        text.clear()
        text.write(json.dumps(report, indent=2))
        collection['blosm_street_report'] = text.name
        collection['blosm_street_coverage'] = report['coverage']
        collection['blosm_street_backend'] = 'python'
        collection['blosm_input_file'] = str(getattr(self.app, 'osmFilepath', ''))
        self.context = SimpleNamespace(app=self.app, snapshot=snapshot, collection=collection,
                                       results=results, report=report, text=text)
        return self.context

    def finalize(self):
        pass

    def cleanup(self):
        pass
