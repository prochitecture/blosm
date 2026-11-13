"""Native undo, rollback, material-only Update and bounded cache lifecycle.

Run only in isolated factory Blender. Existing independent native regression
checks are reused against this checkout, never the sibling main Blosm folder.
"""
import argparse
from dataclasses import replace
import importlib
import importlib.util
import json
from pathlib import Path
import sys

import bpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
                                            submodule_search_locations=[str(ROOT)])
package = importlib.util.module_from_spec(spec)
sys.modules['blosm'] = package
spec.loader.exec_module(package)
from blosm.way.street_network.python_backend import model, ui
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime
from _python_street_fixtures import fixture_networks
from check_street_terrain import make_terrain
from check_street_uv import mesh_snapshot, compare


def invoke(name, *arguments):
    previous = sys.argv
    try:
        sys.argv = [str(ROOT / 'tests' / (name + '.py')), '--', *map(str, arguments)]
        importlib.import_module(name).main()
    finally:
        sys.argv = previous


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    invoke('check_python_street_lifecycle', '--output', args.output / 'cache_lifetime')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    clear_runtime()
    ui.register()
    network = next(n for n in fixture_networks(model) if n.uuid == 'bent_two_intersections')
    network = replace(network, settings=replace(network.settings, resolution=1, exaggeration=1.))
    controller = NetworkScene.create(network)
    controller.blosm_street_terrain.target = make_terrain('Lifecycle original terrain', (0., 0., 0.))
    controller.blosm_street_terrain.mode = 'PROFILE'
    NetworkScene.update(controller)
    before = mesh_snapshot(controller)
    identities = {o['st_py_key']: (o.as_pointer(), o.data.as_pointer())
                  for o in controller['st_py_outputs'].objects if o.type == 'MESH'}
    # A material-only edit must retain actual mesh datablocks and all maps.
    material = bpy.data.materials.new('Performance / material-only check')
    controller.blosm_street_pavement = material
    NetworkScene.update(controller)
    compare(before, mesh_snapshot(controller))
    assert identities == {o['st_py_key']: (o.as_pointer(), o.data.as_pointer())
                          for o in controller['st_py_outputs'].objects if o.type == 'MESH'}
    scene = args.output / 'lifecycle_fixture.blend'
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    bpy.ops.wm.save_as_mainfile(filepath=str(scene.resolve()))
    invoke('check_python_street_undo', '--scene', scene, '--output', args.output / 'undo')
    invoke('check_python_street_publication', '--scene', scene, '--output', args.output / 'rollback')
    checks = dict(passed=True, checkout=str(ROOT), native_cache_lifetime=True,
        material_only_retains_mesh_and_uv=True, undo_redo=True, publication_rollback=True,
        terrain_enabled=True)
    (args.output / 'checks.json').write_text(json.dumps(checks, indent=2), encoding='utf8')
    print('PERFORMANCE_LIFECYCLE', json.dumps(checks), flush=True)


if __name__ == '__main__':
    from blosm.app.import_performance import ImportPerformance
    with ImportPerformance():
        main()
