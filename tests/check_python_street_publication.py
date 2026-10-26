"""Inject a late publication error and verify last-good scene/engine ownership."""
import argparse
import json
from pathlib import Path
import sys
import bpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from blosm.way.street_network.python_backend import scene
from blosm.way.street_network.python_backend.providers import Part


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    controller = next(c for c in bpy.data.objects if c.get('st_py_network_uuid') == 'bent_two_intersections')
    scene.NetworkScene.update(controller)
    before = mesh_buffers(controller)
    manifest = controller['st_py_manifest']
    objects, meshes = set(bpy.data.objects.keys()), set(bpy.data.meshes.keys())
    source = next(o for o in controller['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    original = source.data.vertices[1].co.copy()
    source.data.vertices[1].co.y += .2
    source.data.update()
    publish = scene.SceneWriter.publish

    def fail_after_object_creation(writer, controller, network, engine, update):
        part = next(update.products[k] for k in update.rebuilt if isinstance(update.products[k], Part))
        # An extension supplied invalid serializable metadata. The object and
        # mesh have already been allocated when normal publication detects it.
        part.metadata['injected_invalid_metadata'] = object()
        return publish(writer, controller, network, engine, update)

    scene.SceneWriter.publish = fail_after_object_creation
    try:
        try:
            scene.NetworkScene.update(controller)
        except TypeError:
            pass
        else:
            raise AssertionError('Injected late publication error was not raised')
    finally:
        scene.SceneWriter.publish = publish
    assert set(bpy.data.objects.keys()) == objects, 'Failed publication leaked a Blender object'
    assert set(bpy.data.meshes.keys()) == meshes, 'Failed publication leaked a mesh datablock'
    assert controller['st_py_manifest'] == manifest
    assert controller['st_py_error']
    assert equal(before, mesh_buffers(controller))
    engine = scene._ENGINES[controller['st_py_network_uuid']]
    assert engine.manifest() == json.loads(manifest), 'Failed publication advanced the engine'
    scene.NetworkScene.update(controller)
    assert not equal(before, mesh_buffers(controller))
    source.data.vertices[1].co = original
    source.data.update()
    scene.NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    report = dict(passed=True, injected_late_publication_error=True, no_orphan_objects=True,
        no_orphan_meshes=True, prior_output_preserved=True, prior_manifest_preserved=True,
        engine_revision_preserved=True, retry_and_exact_restoration=True)
    (args.output / 'publication_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
