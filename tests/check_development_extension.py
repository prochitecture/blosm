"""Check the buildings_gn checkout through a real Blender extension repository.

Link the checkout as REPOSITORY/blosm and pass --repository and --artifacts.
Launch Blender 5.2 with an isolated profile, --factory-startup and --disable-autoexec.
"""
import argparse
import importlib
import json
from pathlib import Path
import sys
from types import ModuleType
from unittest.mock import patch

import addon_utils
import bpy

sys.dont_write_bytecode = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--artifacts', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.artifacts.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    assert (args.repository / 'blosm').resolve() == root
    repo = bpy.context.preferences.extensions.repos.new(
        name='Blosm development test', module='blosm_development',
        custom_directory=str(args.repository.resolve()))
    package = 'bl_ext.' + repo.module + '.blosm'
    checks = []
    assert package in {module.__name__ for module in addon_utils.modules(refresh=True)}
    checks.append('extension discovery in the add-on list')

    sentinels = {name: ModuleType(name) for name in
                 ('app', 'util', 'building', 'pml', 'antlr4', 'parse', 'defs', 'grammar', 'way', 'padl', 'lib', 'debug', 'pmlParser')}
    sys.modules.update(sentinels)
    paths_before = list(sys.path)
    addon = addon_utils.enable(package, default_set=True)
    assert addon is not None, 'Extension registration failed'
    assert addon_utils.check(package)[1]
    assert addon.BlosmPreferences.bl_idname == package
    assert Path(addon.__file__).resolve().parent == root
    assert sys.path == paths_before
    assert all(sys.modules[name] is value for name, value in sentinels.items())
    checks.append('registration, preferences and namespace isolation')

    def module(name):
        return importlib.import_module(package + '.' + name)

    for name in ('setup.base', 'setup.premium_default', 'setup.realistic_dev',
                 'style', 'threed_tiles.blender', 'threed_tiles.progressive',
                 'building.padl_integration.backends', 'way.renderer_streets',
                 'way.waymap.waymap', 'action.generate_streets', 'lib.pygeos.index_chain'):
        module(name)
    checks.append('building, street, realistic and 3D Tiles module imports')
    nx = module('lib.blosm_networkx')
    graph = module('way.waymap.waymap').WayMap()
    graph.add_edge(1, 2)
    graph.add_edge(2, 3)
    assert nx.shortest_path(graph, 1, 3) == [1, 2, 3]
    assert len(nx.graph_atlas_g()) == 1253
    checks.append('namespaced WayMap graph operations and bundled graph-atlas resource')

    pml_file = args.artifacts / 'extension_style.pml'
    pml_file.write_text('footprint { height: 10; }', encoding='utf-8')
    store = module('style').StyleStore(str(pml_file), str(args.artifacts))
    assert 'extension_style' in store.styles
    checks.append('ANTLR 4.13.1 PML parsing and style compilation')

    prefs = bpy.context.preferences.addons[package].preferences
    prefs.dataDir = ''
    addon.blenderApp.app.setDataDir(bpy.context, str(root), package)
    assert not Path(addon.blenderApp.app.dataDir).resolve().is_relative_to(root)
    setup_file = args.artifacts / 'custom_setup.py'
    setup_file.write_text('from ..app.blender import app\ndef setup(*args):\n    return app\n', encoding='utf-8')
    assert module('util.setup_script').load_setup_script(setup_file)() is addon.blenderApp.app
    checks.append('external data storage and namespaced custom setup scripts')

    props = bpy.types.Scene.bl_rna.properties['blosm'].fixed_type.properties
    assert 'buildingBackend' in props and 'streetBackend' in props
    assert hasattr(bpy.types.Object, 'blosm_street_pavement')
    assert hasattr(bpy.types.Object, 'blosm_street_terrain')
    assert bpy.context.scene.blosm.streetBackend == 'python'
    checks.append('Geometry Nodes buildings and Python street controls')

    osm_file = args.artifacts / 'local.osm'
    osm_file.write_text('''<osm version="0.6">
      <node id="1" lat="60.1700" lon="24.9400"/>
      <node id="2" lat="60.1700" lon="24.9402"/>
      <node id="3" lat="60.1701" lon="24.9402"/>
      <node id="4" lat="60.1701" lon="24.9400"/>
      <way id="1"><nd ref="1"/><nd ref="2"/><nd ref="3"/><nd ref="4"/><nd ref="1"/>
        <tag k="building" v="yes"/><tag k="height" v="10"/><tag k="roof:shape" v="flat"/>
      </way></osm>''', encoding='utf-8')
    cases = ((mode, backend, single) for mode, backend in
             (('2D', 'geometry_nodes'), ('3Dsimple', 'geometry_nodes'), ('3Dsimple', 'python'))
             for single in (True, False))
    for mode, backend, single in cases:
        for obj in list(bpy.data.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        settings = bpy.context.scene.blosm
        settings.dataType, settings.osmSource, settings.mode = 'osm', 'file', mode
        settings.buildingBackend = backend
        settings.singleObject = single
        settings.osmFilepath = str(osm_file)
        settings.loadMissingMembers = False
        settings.forests = settings.vegetation = settings.highways = settings.railways = settings.water = False
        settings.buildings = True
        module('gui').addDefaultLevels()
        with patch('urllib.request.urlopen', side_effect=AssertionError('Local OSM attempted network')):
            assert bpy.ops.blosm.import_data() == {'FINISHED'}
        bpy.context.view_layer.update()
        geometry = []
        for obj in bpy.context.scene.objects:
            if obj.type != 'MESH' or obj.hide_get():
                continue
            evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
            mesh = evaluated.to_mesh()
            try:
                if mesh.polygons:
                    geometry.extend((evaluated.matrix_world @ vertex.co).z for vertex in mesh.vertices)
            finally:
                evaluated.to_mesh_clear()
        assert geometry, (mode, backend, 'No rendered building geometry')
        if mode == '3Dsimple':
            assert abs(max(geometry) - min(geometry) - 10.0) < .01, (backend, min(geometry), max(geometry))
        checks.append(f'offline local OSM import: {mode} / {backend} / single={single}')

    assert sys.path == paths_before
    assert all(sys.modules[name] is value for name, value in sentinels.items())
    addon_utils.disable(package, default_set=True)
    assert not hasattr(bpy.types.Scene, 'blosm')
    assert not hasattr(bpy.types.Object, 'blosm_street_pavement')
    assert addon_utils.enable(package, default_set=True) is not None
    addon_utils.disable(package, default_set=True)
    checks.append('disable, re-enable and cleanup')
    result = dict(status='passed', blender=bpy.app.version_string, package=package, checks=checks)
    (args.artifacts / 'extension_validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
