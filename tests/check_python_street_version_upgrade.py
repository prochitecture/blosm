"""Qualify a geometry-preserving version upgrade before same-version reopening.

This is an explicit test entry point, not a permissive migration tool. It rejects
any changed source record, mesh position, connectivity or semantic output key.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import bpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from blosm.way.street_network.python_backend import VERSION
from blosm.way.street_network.python_backend.model import fingerprint
from blosm.way.street_network.python_backend.scene import NetworkScene, SourceObjects, clear_runtime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    controllers = [o for o in bpy.data.objects if o.get('st_py_source_data')]
    assert controllers
    records = []
    clear_runtime()
    for controller in controllers:
        previous = json.loads(controller['st_py_manifest'])['implementation']
        source = fingerprint(SourceObjects.read(controller))
        before = mesh_buffers(controller)
        keys = {o['st_py_key'] for o in controller['st_py_outputs'].objects}
        start = time.perf_counter()
        NetworkScene.update(controller)
        seconds = time.perf_counter() - start
        assert source == fingerprint(SourceObjects.read(controller)), 'Upgrade changed authored source data'
        assert equal(before, mesh_buffers(controller)), 'Upgrade changed published mesh buffers'
        assert keys == {o['st_py_key'] for o in controller['st_py_outputs'].objects}, 'Upgrade changed semantic output identity'
        assert json.loads(controller['st_py_manifest'])['implementation'] == VERSION
        records.append(dict(network=controller['st_py_network_uuid'], previous=previous,
                            implementation=VERSION, seconds=seconds,
                            exact_source=True, exact_mesh_buffers=True, exact_semantic_keys=True))
        print(json.dumps(dict(phase='version_upgrade', **records[-1])), flush=True)
    scene = args.output / 'python_street_import_upgraded.blend'
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    bpy.ops.wm.save_as_mainfile(filepath=str(scene.resolve()))
    report = dict(passed=True, source_scene=str(args.scene.resolve()), scene=str(scene.resolve()),
                  implementation=VERSION, controllers=len(controllers), migrations=records,
                  quads=sum(len(o.data.polygons) for c in controllers
                            for o in c['st_py_outputs'].objects if o.type == 'MESH'))
    (args.output / 'upgrade_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(dict(passed=True, scene=report['scene'], quads=report['quads'])), flush=True)


if __name__ == '__main__':
    main()
