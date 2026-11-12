"""Real event-loop/undo qualification in a dedicated factory-profile GUI process."""
from pathlib import Path
import importlib
import json
import sys
import time
import traceback
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
import addon_utils
from release_support import arguments, enable, configure_osm, catalog, PACKAGE

args = arguments()
enable(args)
driver = importlib.import_module(PACKAGE+'.progressive.blender')
osm_runtime = importlib.import_module(PACKAGE+'.app.osm_import')
checks, popups, heartbeats = [], [], []
cancel_durations = []
driver.show_result = lambda level, message: popups.append((level, message))
started = time.monotonic()


class TestUserEdit(bpy.types.Operator):
    bl_idname = 'blosm.test_osm_user_edit'
    bl_label = 'User scene edit'
    bl_options = {'UNDO'}
    name: bpy.props.StringProperty(default='User cube')

    def execute(self, context):
        mesh = bpy.data.meshes.new(self.name)
        mesh.from_pydata([(0, 0, 0), (1, 0, 0), (0, 1, 0)], [], [(0, 1, 2)])
        obj = bpy.data.objects.new(self.name, mesh)
        context.scene.collection.objects.link(obj)
        return {'FINISHED'}


bpy.utils.register_class(TestUserEdit)


def mesh_count():
    return sum(obj.type == 'MESH' for obj in bpy.context.scene.objects)


def setup(fixture='karl_marx_allee'):
    entry = next(e for e in catalog() if e['id'] == fixture)
    settings = configure_osm(args, entry, '3Dsimple', True)
    settings.progressiveImport = True
    return settings


def start():
    assert bpy.ops.blosm.import_data('EXEC_DEFAULT', True) == {'FINISHED'}
    assert driver.active is not None


def wait_finished():
    while driver.active is not None:
        yield
    assert osm_runtime.last_result['outcome'] == 'success', osm_runtime.last_result


def wait_partial():
    while driver.active is not None and not driver.active.adapter.committed_units:
        yield
    assert driver.active is not None
    assert driver.active.adapter.committed_units > 0


def history(action):
    getattr(bpy.ops.ed, action)()
    yield


def cancel():
    processes = []
    if driver.active:
        for renderer in osm_runtime.app.renderers or []:
            queue = getattr(getattr(renderer, 'textureExporter', None), 'bakeQueue', None)
            if queue and queue.process:
                processes.append(queue.process)
    before = time.monotonic()
    bpy.ops.blosm.cancel_progressive()
    cancel_durations.append(time.monotonic()-before)
    assert all(process.poll() is not None for process in processes), 'Render child survived cancellation'


def scenario():
    prefs = bpy.context.preferences.addons[PACKAGE].preferences
    saved_directory = prefs.dataDir
    try:
        for failure in ('data_directory', 'osm_file', 'setup_script', 'corrupt_osm'):
            settings = setup('rectangle')
            if failure == 'data_directory':
                prefs.dataDir = ''
            elif failure == 'osm_file':
                settings.osmFilepath = str(args.artifacts/'missing_original.osm')
            elif failure == 'setup_script':
                settings.setupScript = str(args.artifacts/'missing_original_setup.py')
            else:
                corrupt = args.artifacts/'corrupt_original.osm'
                corrupt.write_text('<osm><node')
                settings.osmFilepath = str(corrupt)
            before = len(popups)
            start()
            while driver.active:
                yield
            result = osm_runtime.last_result
            assert result['outcome'] == 'failed' and result['errors']
            assert len(popups) == before+1
            assert popups[-1] == ('ERROR', result['message'])
            assert result['errors'][0] in popups[-1][1]
            if failure == 'corrupt_osm':
                assert not result['preliminary_failure']
                assert 'OSM import failed: 0 completed,' in popups[-1][1]
            else:
                assert result['preliminary_failure']
                assert popups[-1][1] == result['errors'][0]
            prefs.dataDir = saved_directory
        checks.append('preflight errors show only the original message; parsing failures retain the import summary')
    finally:
        prefs.dataDir = saved_directory

    setup()
    start()
    yield from wait_partial()
    first = mesh_count()
    assert first > 0
    yield from wait_finished()
    total = mesh_count()
    assert total > first
    assert bpy.context.scene.blosm.singleObject
    assert popups and 'success' in popups[-1][1]
    checks.append('actual operator, timer visibility, separate objects and aggregate result')
    yield from history('undo')
    assert mesh_count() == 0, ('one Undo', mesh_count())
    yield from history('redo')
    assert mesh_count() == total
    checks.append('uninterrupted OSM import is one undo step, redo restores geometry')

    setup()
    start()
    yield from wait_partial()
    retained = mesh_count()
    bpy.ops.blosm.test_osm_user_edit('EXEC_DEFAULT', True)
    yield
    yield from wait_finished()
    yield from history('undo')
    assert 'User cube' in bpy.data.objects
    assert mesh_count() >= retained+1
    yield from history('undo')
    assert 'User cube' not in bpy.data.objects
    assert mesh_count() >= retained
    yield from history('undo')
    assert mesh_count() == 0
    checks.append('intervening user edits keep their own undo step')

    setup()
    start()
    yield from wait_partial()
    session = driver.active.adapter
    obj = next(obj for obj in bpy.context.scene.objects if obj.type == 'MESH')
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    count = session.committed_units
    until = time.monotonic()+.3
    while time.monotonic() < until:
        yield
    assert driver.active.paused and session.committed_units == count
    bpy.ops.object.mode_set(mode='OBJECT')
    yield
    session.destination.name = 'User renamed import'
    obj.name = 'User renamed building'
    cancel()
    assert driver.active is None
    assert osm_runtime.last_result['outcome'] == 'cancelled'
    assert mesh_count() > 0
    checks.append('Edit Mode pause, renamed data and explicit cancellation retain committed objects')

    setup()
    start()
    yield from wait_partial()
    session = driver.active.adapter
    bpy.data.collections.remove(session.destination)
    while driver.active:
        yield
    assert osm_runtime.last_result['outcome'] == 'partial'
    assert any('collection' in e for e in osm_runtime.last_result['errors'])
    checks.append('destination deletion stops with one partial-result error')

    setup()
    start()
    other = bpy.data.scenes.new('Other scene')
    original = bpy.context.window.scene
    bpy.context.window.scene = other
    until = time.monotonic()+.2
    while time.monotonic() < until:
        yield
    assert driver.active.paused
    bpy.context.window.scene = original
    bpy.data.scenes.remove(other)
    cancel()
    checks.append('scene switching pauses and cancellation before first work is safe')

    if args.assets:
        settings = setup('rectangle')
        settings.mode = '3Drealistic'
        settings.assetsDir = str(args.assets)
        settings.importForExport = True
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        prefs.enableExperimentalFeatures = True
        folder = args.artifacts/'cancel_bake'
        folder.mkdir(exist_ok=True)
        prefs.dataDir = str(folder)
        start()
        while driver.active and driver.active.adapter.phase != 'Rendering export textures':
            yield
        assert driver.active
        # A user undo checkpoint can capture staged RNA. Restored pending IDs
        # must be discarded after cancellation, including subsequent Redo.
        bpy.ops.blosm.test_osm_user_edit('EXEC_DEFAULT', True, name='User during bake')
        yield
        cancel()
        assert driver.active is None
        for action in ('undo', 'redo', 'undo', 'redo'):
            yield from history(action)
            assert not [obj for obj in bpy.data.objects if obj.get('_blosm_pending_import')]
        assert 'User during bake' in bpy.data.objects
        checks.append('uncached bake yields; cancel/repeated history removes pending data and preserves user edit')

        settings = setup('karl_marx_allee')
        settings.mode, settings.assetsDir = '3Drealistic', str(args.assets)
        settings.importForExport = True
        folder = args.artifacts/'committed_bake'
        folder.mkdir(exist_ok=True)
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        prefs.enableExperimentalFeatures = True
        prefs.dataDir = str(folder)
        start()
        yield from wait_partial()
        bpy.ops.blosm.test_osm_user_edit('EXEC_DEFAULT', True, name='User after baked building')
        cancel()
        for action in ('undo', 'redo'):
            yield from history(action)
            assert not [scene for scene in bpy.data.scenes if scene.get('_blosm_transient_import')]
        assert 'User after baked building' in bpy.data.objects
        checks.append('history after an accepted baked building restores no private template scenes')

    setup()
    start()
    addon_utils.disable(PACKAGE, default_set=True)
    assert driver.active is None
    checks.append('disable before first step cleans the unstarted session')


steps = scenario()


def tick():
    heartbeats.append(time.monotonic())
    try:
        if time.monotonic()-started > 240:
            raise TimeoutError('OSM GUI qualification exceeded four minutes')
        next(steps)
        return .03
    except StopIteration:
        status, error = 'passed', None
    except Exception:
        status, error = 'failed', traceback.format_exc()
    if driver.active:
        driver.active.close(silent=True)
    report = dict(status=status, error=error, checks=checks, heartbeats=len(heartbeats),
                  max_heartbeat_gap_seconds=max((b-a for a, b in zip(heartbeats, heartbeats[1:])), default=0),
                  max_cancel_seconds=max(cancel_durations, default=0),
                  elapsed_seconds=time.monotonic()-started)
    (args.artifacts/'osm_ui_checks.json').write_text(json.dumps(report, indent=2))
    print('OSM_UI_CHECKS', json.dumps(report), flush=True)
    bpy.ops.wm.quit_blender()
    return None


bpy.app.timers.register(tick, first_interval=.5)
