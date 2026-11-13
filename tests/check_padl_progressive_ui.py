"""Native building event-loop, accepted-data editing and undo qualification.

Launch in a dedicated Blender GUI process with the release runner arguments.
The test drives the real import operator and its registered controller timer.
"""
from pathlib import Path
import importlib
import json
import sys
import time
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
from release_support import arguments, enable, configure_osm, catalog, PACKAGE

args = arguments()
enable(args)
driver = importlib.import_module(PACKAGE + '.progressive.blender')
runtime = importlib.import_module(PACKAGE + '.app.osm_import')
building_assets = importlib.import_module(PACKAGE + '.building.padl_integration.assets')
read_footprints = importlib.import_module(PACKAGE + '.padl.persistence').read_footprints
checks, popups, heartbeats, progress_runs = [], [], [], []
cancel_durations = []
driver.show_result = lambda level, message: popups.append((level, message))
started = time.monotonic()
current_progress = None
asset_checks = []
original_tree = building_assets.BuildingProgramLibrary.tree


def verified_tree(library, role):
    asset_checks.append(role)
    return original_tree(library, role)


building_assets.BuildingProgramLibrary.tree = verified_tree


class PadlTestUserEdit(bpy.types.Operator):
    bl_idname = 'blosm.test_padl_user_edit'
    bl_label = 'Native building test user edit'
    bl_options = {'UNDO'}

    def execute(self, context):
        mesh = bpy.data.meshes.new('User mesh between building groups')
        mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0, 1, 0)], [], [(0, 1, 2)])
        obj = bpy.data.objects.new('User edit between building groups', mesh)
        context.scene.collection.objects.link(obj)
        return {'FINISHED'}


bpy.utils.register_class(PadlTestUserEdit)


def owners(scene_only=True):
    source = bpy.context.scene.objects if scene_only else bpy.data.objects
    return [obj for obj in source if obj.get('padl_source_keys')]


def assert_no_temporary_readers():
    assert not [item.name for item in bpy.data.node_groups
                if item.name.startswith('PADL temporary status reader')]
    assert not [item.name for item in bpy.data.objects
                if item.name.startswith('PADL temporary status owner')]


def report_state(owner_count, status_count=None, *, scene_only=True):
    """Discover restored RNA after Undo/Redo; never reuse old Python result IDs."""
    found = owners(scene_only)
    assert len(found) == owner_count, ('editable input count', owner_count, len(found))
    assert_no_temporary_readers()
    if not found:
        return None
    report_text = found[0].get('blosm_building_report')
    assert isinstance(report_text, bpy.types.Text), 'Building lost its persisted report'
    report = json.loads(report_text.as_string())
    assert report['input_count'] == 3, report
    if status_count is not None:
        assert len(report['statuses']) == status_count, report
    assert all(not item.get('error') for item in report['statuses']), report
    reported_ids = {item['footprint_id'] for item in report['statuses']}
    for obj in found:
        assert obj.get('blosm_building_report') == report_text
        assert '_blosm_pending_import' not in obj
        footprints = read_footprints(obj)
        assert len(footprints) == 1
        assert footprints[0].footprint_id in reported_ids
        if bpy.context.scene.blosm.buildingBackend == 'python':
            assert obj.hide_render
            if scene_only:
                assert obj.hide_get()
        else:
            assert obj.modifiers.get('PADL Buildings') is not None
    return report


def setup(backend):
    entry = next(item for item in catalog() if item['id'] == 'shared_edges')
    settings = configure_osm(args, entry, '3Dsimple', True)
    settings.progressiveImport = True
    settings.buildingBackend = backend
    settings.buildingSeed = 113
    return settings


def start(backend, case):
    global current_progress
    assert bpy.ops.blosm.import_data('EXEC_DEFAULT', True) == {'FINISHED'}
    assert driver.active is not None
    current_progress = dict(backend=backend, case=case, commits=[], asset_check_start=len(asset_checks))
    progress_runs.append(current_progress)


def observe_progress():
    assert_no_temporary_readers()
    if driver.active is None or current_progress is None:
        return
    count = driver.active.adapter.committed_units
    commits = current_progress['commits']
    if count and (not commits or count != commits[-1]):
        assert count == (commits[-1] + 1 if commits else 1), ('Missed UI commit', commits, count)
        commits.append(count)
        # Deleted accepted buildings are intentionally still counted in the
        # cumulative import report; all surviving owners must retain that report.
        report_state(len(owners()), count)


def wait_count(count):
    while driver.active and driver.active.adapter.committed_units < count:
        yield
    assert driver.active is not None, ('Import ended before expected group', count)
    assert driver.active.adapter.committed_units == count
    observe_progress()


def wait_finished(outcome='success'):
    while driver.active:
        yield
    result = runtime.last_result
    assert result['outcome'] == outcome, result
    assert_no_temporary_readers()
    if outcome == 'success':
        assert not result['errors'], result
        assert result['committed'] == 3, result
        assert current_progress['commits'] == [1, 2, 3], current_progress
        # Real UI redraws and undo checkpoints run between these commits. A
        # headless loop cannot expose deferred graph updates from publication.
        checks = len(asset_checks) - current_progress['asset_check_start']
        assert checks == (1 if current_progress['backend'] == 'geometry_nodes' else 0), (
            'Generator revalidated without an asset edit', current_progress, checks)
    return result


def history(action):
    assert getattr(bpy.ops.ed, action)() == {'FINISHED'}
    yield
    assert_no_temporary_readers()


def cancel():
    before = time.monotonic()
    assert bpy.ops.blosm.cancel_progressive() == {'FINISHED'}
    cancel_durations.append(time.monotonic() - before)
    assert driver.active is None


def scenario():
    user_name = 'User edit between building groups'
    for backend in ('geometry_nodes', 'python'):
        settings = setup(backend)
        before_popups = len(popups)
        start(backend, 'uninterrupted_undo_redo')
        yield from wait_count(1)
        report_state(1, 1)
        yield from wait_finished()
        finished_report = report_state(3, 3)
        assert settings.singleObject, 'Saved blocking preference changed'
        assert len(popups) == before_popups + 1
        assert popups[-1][1] == runtime.last_result['message']
        yield from history('undo')
        report_state(0)
        assert not [obj for obj in bpy.context.scene.objects if obj.type == 'MESH']
        yield from history('redo')
        assert report_state(3, 3) == finished_report
        assert sum(obj.type == 'MESH' for obj in bpy.context.scene.objects) == (
            6 if backend == 'python' else 3)
        checks.append(backend + ': event-loop commits [1, 2, 3], one Undo removes import, Redo restores inputs/output/report')

        setup(backend)
        start(backend, 'user_edit_undo_boundaries')
        yield from wait_count(1)
        first_report = report_state(1, 1)
        assert bpy.ops.blosm.test_padl_user_edit('EXEC_DEFAULT', True) == {'FINISHED'}
        yield from wait_finished()
        report_state(3, 3)
        yield from history('undo')
        assert user_name in bpy.data.objects
        assert report_state(1, 1) == first_report
        yield from history('undo')
        assert user_name not in bpy.data.objects
        assert report_state(1, 1) == first_report
        yield from history('undo')
        report_state(0)
        yield from history('redo')
        assert report_state(1, 1) == first_report
        yield from history('redo')
        assert user_name in bpy.data.objects
        assert report_state(1, 1) == first_report
        yield from history('redo')
        report_state(3, 3)
        checks.append(backend + ': intervening user edit has its own Undo/Redo step and reports follow restored groups')

        setup(backend)
        start(backend, 'cancel_after_first')
        yield from wait_count(1)
        before = report_state(1, 1)
        cancel()
        assert runtime.last_result['outcome'] == 'cancelled', runtime.last_result
        assert runtime.last_result['committed'] == 1
        assert report_state(1, 1) == before
        yield from history('undo')
        report_state(0)
        yield from history('redo')
        assert report_state(1, 1) == before
        checks.append(backend + ': cancellation retains one complete building and its report through Undo/Redo')

        setup(backend)
        start(backend, 'rename_and_delete_accepted')
        yield from wait_count(1)
        first = runtime.app.padlBuildingResults[0]
        first.owner.name = 'User renamed PADL input'
        first.output.name = 'User renamed PADL output'
        preserved_names = {first.owner.name, first.output.name}
        yield from wait_count(2)
        second = runtime.app.padlBuildingResults[1]
        removed = (second.owner,) if second.owner == second.output else (second.output, second.owner)
        removed_names = [obj.name for obj in removed]
        for obj in removed:
            bpy.data.objects.remove(obj, do_unlink=True)
        yield from wait_finished()
        report_state(2, 3)
        assert all(name in bpy.data.objects for name in preserved_names)
        assert all(name not in bpy.data.objects for name in removed_names)
        checks.append(backend + ': renamed accepted input/output survive; deleting another accepted group does not interrupt remaining groups')

        setup(backend)
        start(backend, 'destination_deleted')
        yield from wait_count(1)
        session = driver.active.adapter
        bpy.data.collections.remove(session.destination)
        yield from wait_finished('partial')
        assert runtime.last_result['committed'] == 1
        assert any('collection' in message for message in runtime.last_result['errors'])
        report_state(1, 1, scene_only=False)
        assert not session.ownership.pending
        checks.append(backend + ': deleting the destination stops with a partial-result error and no temporary status reader')


steps = scenario()


def tick():
    heartbeats.append(time.monotonic())
    try:
        if time.monotonic() - started > 300:
            raise TimeoutError('Native building GUI qualification exceeded five minutes')
        observe_progress()
        next(steps)
        return .01
    except StopIteration:
        status, error = 'passed', None
    except Exception:
        status, error = 'failed', traceback.format_exc()
    if driver.active:
        driver.active.close(silent=True)
    report = dict(status=status, error=error, checks=checks, progress_runs=progress_runs,
                  heartbeats=len(heartbeats),
                  max_heartbeat_gap_seconds=max((b-a for a, b in zip(heartbeats, heartbeats[1:])), default=0),
                  max_cancel_seconds=max(cancel_durations, default=0),
                  elapsed_seconds=time.monotonic() - started)
    (args.artifacts / 'padl_ui_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PADL_UI_CHECKS', json.dumps(report), flush=True)
    bpy.ops.wm.quit_blender()
    return None


bpy.app.timers.register(tick, first_interval=.5)
