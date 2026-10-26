"""Normal Blosm import adapter for the direct numerical mesh backend."""
import json
import time
from types import SimpleNamespace
from uuid import uuid4
import bpy
from ....app.asset_library import BlenderAssetLoader
from .. import assets
from ..adapter import PreparedNetworkAdapter
from ..design import StreetDesign
from ..physical_sources import PhysicalStreetInventory
from ..physical_ports import GroundSurfacePolicy
from ...way_properties import wayCategoryProps
from ..transaction import ImportTransaction
from .adapter import PreparedPythonAdapter, components
from .model import Settings
from .scene import NetworkScene
from .source_support import ClippedApproachSupport
from . import ui


class PythonStreetRenderer:
    def __init__(self, app, *, design=None, features=None, asset_loader=None):
        self.app = app
        self.design = design or StreetDesign.from_app(app)
        self.loader = asset_loader or BlenderAssetLoader.for_app(app,
            assets.catalog(getattr(app, 'streetAssetFile', None)))
        self.adapter = PreparedNetworkAdapter(self.design, ground_passages=True)
        self.converter = PreparedPythonAdapter(self.design.categories, GroundSurfacePolicy.reason, wayCategoryProps)
        self.context = None
        self.features = None if features is None else tuple(features)

    def prepare(self):
        if bpy.app.version < (5, 2, 0):
            raise ValueError('The Python street backend targets Blender 5.2 or newer')
        ui.register()
        self.loader.resolve(assets.LIBRARY.key)

    def render(self, manager, data):
        with ImportTransaction():
            return self._render(manager, data)

    def _render(self, manager, data):
        start = time.perf_counter()
        snapshot = self.adapter.read(manager)
        physical = PhysicalStreetInventory.view(snapshot)
        settings = Settings(self.design.cell_size, self.design.resolution, self.design.exaggeration)
        adaptation = self.converter.read(physical, str(uuid4()), settings,
            rise=self.design.profile_rise, interpolation=self.design.interpolation)
        loaded = self.loader.load_many((assets.PAVEMENT, assets.GRASS))
        materials = {'pavement': loaded[assets.PAVEMENT], 'grass': loaded[assets.GRASS]}
        collection = bpy.data.collections.new('Blosm / Python streets')
        bpy.context.scene.collection.children.link(collection)
        results = []
        support = ClippedApproachSupport(physical, snapshot.clipped_source_context)
        for number, network in enumerate(components(adaptation.network)):
            preparation = support.prepare(network)
            network = preparation.network
            controller = NetworkScene.create(network, collection, name=f'Python streets / Network {number}',
                                             materials=materials, generate=False, features=self.features)
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
                update = NetworkScene.update(controller)
                result.update(status='READY', quads=sum(len(o.data.polygons)
                    for o in controller['st_py_outputs'].objects if o.type == 'MESH'), timings=update.timings)
                result['fitted_corner_radii'] = [fit for obj in controller['st_py_outputs'].objects
                    if obj.get('st_py_metadata') for fit in
                    json.loads(obj['st_py_metadata']).get('fitted_corner_radii', ())]
            except (ValueError, RuntimeError) as exc:
                result.update(status='FAILED', reason=str(exc))
                controller['st_py_error'] = str(exc)
            results.append(result)
        report = dict(backend='python', version=1, source=snapshot.report(),
            preparation_issues=adaptation.issues, surfaces=results,
            coverage='COMPLETE' if results and all(r['status'] == 'READY' for r in results)
                     and not adaptation.issues else 'PARTIAL', seconds=time.perf_counter() - start)
        text = bpy.data.texts.new('Blosm Python street import report')
        text.write(json.dumps(report, indent=2))
        collection['blosm_street_report'] = text.name
        collection['blosm_street_coverage'] = report['coverage']
        collection['blosm_street_backend'] = 'python'
        collection['blosm_input_file'] = str(getattr(self.app, 'osmFilepath', ''))
        self.context = SimpleNamespace(app=self.app, snapshot=snapshot, collection=collection,
                                       results=results, report=report)
        if getattr(self.app, 'op', None):
            ready = sum(r['status'] == 'READY' for r in results)
            self.app.op.report({'INFO' if report['coverage'] == 'COMPLETE' else 'WARNING'},
                f'Python streets: {ready}/{len(results)} networks ready. Coverage: {report["coverage"]}. See {text.name}')
        return self.context

    def finalize(self):
        pass

    def cleanup(self):
        pass
