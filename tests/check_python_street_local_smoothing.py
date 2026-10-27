"""Native manual-Update and saved-mesh checks for localized centerline bends."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import bpy
import bmesh
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from blosm.way.street_network.python_backend.scene import NetworkScene, SourceObjects
from blosm.way.street_network.python_backend import model, ui


def hashes(controller):
    result = {}
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH':
            continue
        xyz = np.empty(len(obj.data.vertices) * 3, dtype=np.float32)
        indices = np.empty(len(obj.data.loops), dtype=np.int32)
        obj.data.vertices.foreach_get('co', xyz)
        obj.data.loops.foreach_get('vertex_index', indices)
        result[obj['st_py_key']] = hashlib.sha256(xyz.tobytes() + indices.tobytes()).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path)
    parser.add_argument('--section', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    if args.scene:
        bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
        controller = next(o for o in bpy.data.objects if o.get('st_py_source_data')
            and args.section in model.Network.from_dict(json.loads(o['st_py_source_data'])).section_index)
    else:
        for obj in list(bpy.data.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        source = model.Network('local-courtyard', (
            model.Section(1, 1, ((0., 0., 0.), (0., -80., 0.), (120., -80., 0.), (120., 0., 0.)), width=5.),
            model.Section(99, 99, ((0., 80., 0.), (120., 80., 0.)))),
            settings=model.Settings(exaggeration=1.), furniture=(model.Placement(81, 1, 45.),))
        controller = NetworkScene.create(source)
    if not ui.BLOSM_PT_PythonStreets.is_registered:
        ui.register()
    before = hashes(controller)
    pointers = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    cold = NetworkScene.update(controller)
    assert hashes(controller) == before
    assert cold.timings.get('published_objects', 0) == 0
    source = SourceObjects.read(controller)
    section = source.section_index[args.section]
    obj = next(o for o in controller['st_py_sources'].objects
        if any(d.value == args.section for d in o.data.attributes['st_section_id'].data))
    degree = {}
    for edge in obj.data.edges:
        if obj.data.attributes['st_section_id'].data[edge.index].value == args.section:
            for vertex in edge.vertices:
                degree[vertex] = degree.get(vertex, 0) + 1
    controls = sorted(i for i, count in degree.items() if count == 2)
    vertex = controls[len(controls) // 2] if controls else min(degree)
    position = obj.data.vertices[vertex].co.copy()
    # A Street can contain several two-control Sections. Move both endpoint
    # copies at their explicit seam, preserving the existing editing contract.
    vertices = [vertex] if controls else [v.index for v in obj.data.vertices
        if (v.co - position).length < 1e-6]
    original = {i: obj.data.vertices[i].co.copy() for i in vertices}
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bm.verts.ensure_lookup_table()
    for index in vertices:
        bm.verts[index].co.x += .1
    bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
    bpy.ops.object.mode_set(mode='OBJECT')
    assert hashes(controller) == before, 'Mesh changed before the Update button'
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    changed = hashes(controller)
    key = f'road:{args.section}'
    assert changed[key] != before[key]
    assert 'centerline_smoothing' in json.loads(next(o for o in controller['st_py_outputs'].objects
        if o.get('st_py_key') == key)['st_py_metadata'])
    rebuilt = set(json.loads(controller['st_py_last_update'])['rebuilt'])
    # Related output can legitimately be republished with byte-identical
    # rounded coordinates, for example after a tangent/dependency change.
    # Only outputs outside the actual rebuild set must retain their objects.
    remote = {key: pointer for key, pointer in pointers.items()
              if key in changed and key not in rebuilt and before.get(key) == changed[key]}
    assert remote
    current = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    assert all(current[key] == pointer for key, pointer in remote.items())
    widths = obj.data.attributes['st_width']
    old_widths = [d.value for d in widths.data]
    for edge in obj.data.edges:
        if obj.data.attributes['st_section_id'].data[edge.index].value == args.section:
            widths.data[edge.index].value += .1
    obj.data.update()
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    assert hashes(controller)[key] != changed[key]
    widths.data.foreach_set('value', old_widths)
    for index, position in original.items():
        obj.data.vertices[index].co = position
    obj.data.update()
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    assert hashes(controller) == before, 'Restoration differs from the saved current geometry'
    assert all(m.type != 'NODES' for o in bpy.data.objects for m in o.modifiers)
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    scene = args.output / 'local_smoothing_editable.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(scene.resolve()))
    report = dict(passed=True, section=args.section, saved_cold_identity=True,
        real_edit_mode=True, manual_update=True, width_update=True, exact_restoration=True,
        unchanged_output_identity=True, retained_mesh_objects=len(remote),
        no_geometry_nodes=True, scene=str(scene.resolve()))
    (args.output / 'native_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
