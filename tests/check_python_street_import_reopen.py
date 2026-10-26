"""Cold reconstruction and local-versus-full Update on a real imported scene."""
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
from _python_street_edit import SectionProbeEdit, SectionSeamProbeEdit
from blosm.way.street_network.python_backend import VERSION
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime
from blosm.way.street_network.python_backend import ui


def identities(controller):
    return {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--edit-kind', choices=('terminal', 'seam'), default='terminal')
    parser.add_argument('--connection', type=int)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    ui.register()
    controllers = [o for o in bpy.data.objects if o.get('st_py_source_data')]
    assert controllers
    for controller in controllers:
        previous = json.loads(controller['st_py_manifest'])['implementation']
        if previous != VERSION:
            raise ValueError(f'Scene uses {previous}, current backend is {VERSION}; qualify the explicit '
                             'version upgrade before same-version reopening/identity checks')
    saved = {c.name: (mesh_buffers(c), identities(c)) for c in controllers}
    timings = {}
    clear_runtime()
    for controller in controllers:
        start = time.perf_counter()
        NetworkScene.update(controller)
        timings[controller.name] = dict(cold_seconds=time.perf_counter() - start)
        assert equal(saved[controller.name][0], mesh_buffers(controller))
        assert saved[controller.name][1] == identities(controller)
    controller = max(controllers, key=lambda c: len(json.loads(c['st_py_source_data'])['sections']))
    if args.connection is not None and args.edit_kind != 'seam':
        raise ValueError('A connection selection requires the seam edit kind')
    edit = (SectionSeamProbeEdit(controller, args.connection) if args.edit_kind == 'seam'
            else SectionProbeEdit(controller))
    network, section, free = edit.network, edit.section, edit.free
    edit.apply()
    before = saved[controller.name][0]
    initial_ids = identities(controller)
    assert equal(before, mesh_buffers(controller)), 'Output changed before explicit Update'
    start = time.perf_counter()
    NetworkScene.update(controller)
    local_seconds = time.perf_counter() - start
    edited = mesh_buffers(controller)
    assert not equal(before, edited)
    after_ids = identities(controller)
    retained = sum(after_ids.get(key) == value for key, value in initial_ids.items())
    assert retained > len(initial_ids) // 2, 'A local edit replaced most network outputs'
    start = time.perf_counter()
    NetworkScene.update(controller, force=True)
    force_seconds = time.perf_counter() - start
    assert equal(edited, mesh_buffers(controller)), 'Incremental and full real-network builds differ'
    edit.restore()
    start = time.perf_counter()
    NetworkScene.update(controller)
    restore_seconds = time.perf_counter() - start
    assert equal(before, mesh_buffers(controller)), 'Restoration changed the imported buffers'
    for other in controllers:
        if other != controller:
            assert equal(saved[other.name][0], mesh_buffers(other))
            assert saved[other.name][1] == identities(other)
    report = dict(passed=True, scene=str(args.scene.resolve()), controllers=len(controllers),
        sections=len(network.sections), edited_section=section.id, edited_street=section.street_id,
        free_end=free, terminal_kind=edit.terminal_kind,
        cold_updates=timings, local_seconds=local_seconds,
        forced_full_seconds=force_seconds, restore_seconds=restore_seconds,
        retained_outputs=retained, initial_outputs=len(initial_ids),
        explicit_update=True, incremental_equals_full=True, exact_restoration=True,
        unchanged_other_components=True, gn_modifiers=sum(m.type == 'NODES' for o in bpy.data.objects for m in o.modifiers))
    report['edit_kind'] = args.edit_kind
    if args.edit_kind == 'seam':
        report.update(edited_connection=edit.connection.id,
                      edited_sections=[s.id for s in edit.sections],
                      seam_vertex_copies=len(edit.positions),
                      coordinated_width_sections=len(edit.sections))
    (args.output / 'import_reopen_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
