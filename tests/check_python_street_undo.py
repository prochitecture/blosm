"""Native global Undo/Redo checkpoints around source edits and manual Update.

The undo pushes provide deterministic background-Blender history boundaries;
the same registered Update operator has Blender's UNDO flag in the UI.
All Blender object references are reacquired after a history restoration.
"""
import argparse
import json
from pathlib import Path
import sys
import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from blosm.way.street_network.python_backend.scene import NetworkScene
from blosm.way.street_network.python_backend import ui


def controller():
    return next(c for c in bpy.data.objects if c.get('st_py_network_uuid') == 'bent_two_intersections')


def source():
    return next(o for o in controller()['st_py_sources'].objects
                if 1 in {v.value for v in o.data.attributes['st_section_id'].data})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    ui.register()
    bpy.context.preferences.edit.use_global_undo = True
    NetworkScene.update(controller())
    before = mesh_buffers(controller())
    before_manifest = controller()['st_py_manifest']
    original = tuple(source().data.vertices[1].co)
    bpy.ops.object.select_all(action='DESELECT')
    source().select_set(True)
    bpy.context.view_layer.objects.active = source()
    assert bpy.ops.ed.undo_push(message='Python streets initial') == {'FINISHED'}
    source().data.vertices[1].co.y += .3
    source().data.update()
    edited = tuple(source().data.vertices[1].co)
    assert bpy.ops.ed.undo_push(message='Python streets edited source') == {'FINISHED'}
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    after = mesh_buffers(controller())
    assert not equal(before, after)
    after_manifest = controller()['st_py_manifest']
    assert bpy.ops.ed.undo_push(message='Python streets generated output') == {'FINISHED'}
    assert bpy.ops.ed.undo() == {'FINISHED'}
    assert tuple(source().data.vertices[1].co) == edited
    assert controller()['st_py_manifest'] == before_manifest
    assert equal(before, mesh_buffers(controller()))
    assert bpy.ops.ed.redo() == {'FINISHED'}
    assert controller()['st_py_manifest'] == after_manifest
    assert equal(after, mesh_buffers(controller()))
    # Move back through both publication and source edit, then reuse an engine
    # whose last Python runtime state came from the later history branch.
    assert bpy.ops.ed.undo() == {'FINISHED'}
    assert bpy.ops.ed.undo() == {'FINISHED'}
    assert tuple(source().data.vertices[1].co) == original
    assert equal(before, mesh_buffers(controller()))
    NetworkScene.update(controller())
    assert equal(before, mesh_buffers(controller()))
    source().data.vertices[1].co = edited
    source().data.update()
    NetworkScene.update(controller())
    assert equal(after, mesh_buffers(controller()))
    report = dict(passed=True, source_undo=True, output_undo=True, output_redo=True,
        update_after_history_restore=True, deterministic_regeneration=True,
        history='native global undo checkpoints in background Blender')
    (args.output / 'undo_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
