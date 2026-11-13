"""Native street timers, early previews, pause, cancellation and Undo/Redo.

Run in a newly launched factory-profile Blender GUI (never --background) with
--factory-startup --disable-autoexec --python this_file -- --artifacts <folder>.
The first import uses the unchanged production driver and records heartbeat
gaps. A second import holds one genuine early visible preview checkpoint to exercise
pending-data history restoration deterministically; its timing is separate.
--osm runs only the unchanged production timer on the selected real extract;
--expected-hash optionally checks its final exact scene digest against the
background parity report. Visibility timing means visible scene objects sampled
by the UI heartbeat, not pixel-level presentation timing.
The isolated test window closes itself. No preferences or source .blend saves.
"""
import argparse
from functools import wraps
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

import bpy

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_street_progressive as support
import check_street_previews as preview_support

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--artifacts', required=True, type=Path)
parser.add_argument('--osm', type=Path,
                    help='Qualify normal production timers on one real local extract only')
parser.add_argument('--expected-hash',
                    help='Expected exact_geometry_sha256 from check_street_progressive.py')
parser.add_argument('--performance-request', action='store_true',
                    help='Enable the normal scoped Import performance preference for this run')
parser.add_argument('--no-navigation', action='store_true',
                    help='Keep the observer heartbeat but disable its simulated viewport navigation')
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
if args.expected_hash and not args.osm:
    parser.error('--expected-hash requires --osm')
if args.osm and not args.osm.is_file():
    parser.error('The selected local OSM extract does not exist')
assert not bpy.app.background, 'This check requires its own interactive Blender process'
args.artifacts.mkdir(parents=True, exist_ok=True)
(args.artifacts / 'data').mkdir(exist_ok=True)
support.ARGS = args
support.FIXTURE = args.osm.resolve() if args.osm else args.artifacts / 'three_disconnected_networks.osm'
if not args.osm:
    support.synthetic_osm(support.FIXTURE)
support.enable()
from blosm.app import osm_import as runtime
from blosm.progressive import blender as driver, ownership, core
from blosm.way.street_network.python_backend import scene as scene_runtime
support.runtime, support.driver = runtime, driver
support.ownership, support.scene_runtime = ownership, scene_runtime
checks = support.StreetProgressiveChecks()
checks.prepare()
bpy.context.preferences.edit.use_global_undo = True

results, heartbeats, cases = [], [], []
driver.show_result = lambda level, message: results.append((level, message))
started = time.monotonic()
phase = 'normal completion'
original_update = scene_runtime.NetworkScene.update_steps
stage_held = False
other_scene = None
first_state = complete_state = None
normal_result = None
normal_started = None
normal_completed = None
normal_geometry_digest = None
last_progress = started - 5.
original_job_step = runtime.OsmJob.step
original_job_init = runtime.OsmJob.__init__
original_core_step = core.Job.step
original_controller_tick = driver.Controller.tick
controller_tick_depth = 0
source_hashes = {str(Path(path).resolve()): hashlib.sha256(Path(path).read_bytes()).hexdigest()
                 for path in (__file__, support.__file__, preview_support.__file__, runtime.__file__,
                    driver.__file__, core.__file__, scene_runtime.__file__,
                    *(Path(scene_runtime.__file__).with_name(name) for name in
                      ('previews.py', 'backend.py', 'updates.py', 'renderer.py')))}


def empty_scope():
    return dict(count=0, cumulative_seconds=0., maximum_seconds=0., longest_callbacks=[])


SCOPE_NAMES = ('core_job_step', 'osm_job_step', 'controller_tick', 'osm_job_step_inside_controller')
wall_scopes = {group: {name: empty_scope() for name in SCOPE_NAMES}
               for group in ('all', 'normal_import')}
scope_gaps = {group: {name: dict(last_end=None, samples=[]) for name in SCOPE_NAMES}
              for group in wall_scopes}
yield_counts = {group: dict(total=0, commands={}) for group in wall_scopes}


def record_scope(name, before, after, session, harness_phase, before_phase, before_committed):
    elapsed = after - before
    for group in ('all', 'normal_import') if harness_phase == 'normal completion' else ('all',):
        scope = wall_scopes[group][name]
        gaps = scope_gaps[group][name]
        if gaps['last_end'] is not None:
            gaps['samples'].append(max(0., before - gaps['last_end']))
        gaps['last_end'] = after
        scope['count'] += 1
        scope['cumulative_seconds'] += elapsed
        scope['maximum_seconds'] = max(scope['maximum_seconds'], elapsed)
        longest = scope['longest_callbacks']
        if len(longest) < 10 or elapsed > longest[-1]['seconds']:
            longest.append(dict(seconds=elapsed, harness_phase=harness_phase,
                session_phase_before=before_phase, session_phase_after=getattr(session, 'phase', None),
                committed_before=before_committed, committed_after=getattr(session, 'committed_units', None)))
            longest.sort(key=lambda item: item['seconds'], reverse=True)
            del longest[10:]


class CountedSteps:
    """Count only commands returned to Job; retain send/throw/close semantics."""
    def __init__(self, source):
        self.source = source

    @staticmethod
    def counted(command):
        for group in ('all', 'normal_import') if phase == 'normal completion' else ('all',):
            counts = yield_counts[group]
            counts['total'] += 1
            name = type(command).__name__
            counts['commands'][name] = counts['commands'].get(name, 0) + 1
        return command

    def send(self, value):
        return self.counted(self.source.send(value))

    def throw(self, *arguments):
        return self.counted(self.source.throw(*arguments))

    def close(self):
        return self.source.close()

    def __iter__(self):
        return self

    def __next__(self):
        return self.send(None)


@wraps(original_job_init)
def measured_job_init(self, *positional, **keywords):
    original_job_init(self, *positional, **keywords)
    self.steps = CountedSteps(self.steps)


@wraps(original_core_step)
def measured_core_step(self, *positional, **keywords):
    before, harness_phase = time.monotonic(), phase
    session = getattr(self, 'session', None)
    before_phase = getattr(session, 'phase', None)
    before_committed = getattr(session, 'committed_units', None)
    try:
        return original_core_step(self, *positional, **keywords)
    finally:
        record_scope('core_job_step', before, time.monotonic(), session, harness_phase,
                     before_phase, before_committed)


@wraps(original_job_step)
def measured_job_step(self, *positional, **keywords):
    before, harness_phase = time.monotonic(), phase
    session = self.session
    before_phase, before_committed = session.phase, session.committed_units
    nested = controller_tick_depth > 0
    try:
        return original_job_step(self, *positional, **keywords)
    finally:
        after = time.monotonic()
        record_scope('osm_job_step', before, after, session, harness_phase, before_phase, before_committed)
        if nested:
            record_scope('osm_job_step_inside_controller', before, after, session, harness_phase,
                         before_phase, before_committed)


@wraps(original_controller_tick)
def measured_controller_tick(self, *positional, **keywords):
    global controller_tick_depth
    before, harness_phase = time.monotonic(), phase
    session = self.adapter
    before_phase = getattr(session, 'phase', None)
    before_committed = getattr(session, 'committed_units', None)
    controller_tick_depth += 1
    try:
        return original_controller_tick(self, *positional, **keywords)
    finally:
        controller_tick_depth -= 1
        record_scope('controller_tick', before, time.monotonic(), session, harness_phase,
                     before_phase, before_committed)


def wall_scope_report(*, compact=False):
    now = time.monotonic()
    elapsed = now - started
    normal_elapsed = ((normal_completed or now) - normal_started if normal_started is not None else 0.)
    derived = {}
    for group, seconds in (('all', elapsed), ('normal_import', normal_elapsed)):
        scopes = wall_scopes[group]
        driver_seconds = scopes['controller_tick']['cumulative_seconds']
        nested_seconds = scopes['osm_job_step_inside_controller']['cumulative_seconds']
        derived[group] = dict(observed_wall_seconds=seconds,
            controller_excluding_nested_job_seconds=max(0., driver_seconds - nested_seconds),
            outside_controller_seconds=max(0., seconds - driver_seconds),
            osm_wrapper_excluding_core_seconds=max(0., scopes['osm_job_step']['cumulative_seconds']
                                                       - scopes['core_job_step']['cumulative_seconds']),
            osm_job_outside_controller_seconds=max(0., scopes['osm_job_step']['cumulative_seconds'] - nested_seconds))
    scopes_report = {}
    for group, scopes in wall_scopes.items():
        scopes_report[group] = {}
        for name, scope in scopes.items():
            values = sorted(scope_gaps[group][name]['samples'])
            length = len(values)
            gaps = dict(count=length, sum_seconds=sum(values),
                median_seconds=((values[(length - 1) // 2] + values[length // 2]) / 2 if length else None),
                p95_seconds=values[max(0, (95 * length + 99) // 100 - 1)] if length else None)
            scopes_report[group][name] = dict(
                (key, value) for key, value in scope.items() if not compact or key != 'longest_callbacks')
            scopes_report[group][name]['start_after_previous_end_gaps'] = gaps
    return dict(measurement='Elapsed wall scopes, not CPU time; outside includes UI idle/redraw and harness work',
                scopes=scopes_report, derived=derived, generator_yields=yield_counts)


runtime.OsmJob.step = measured_job_step
runtime.OsmJob.__init__ = measured_job_init
core.Job.step = measured_core_step
driver.Controller.tick = measured_controller_tick


def pending():
    return [(kind, item) for kind in ownership.KINDS for item in getattr(bpy.data, kind)
            if item.get(ownership.PENDING)]


def coherent():
    assert not pending(), 'Undo/Redo restored pending street data'
    assert driver.active is None, 'Undo/Redo restarted a completed import'
    controllers = support.controllers()
    assert all(controller.get('st_py_manifest') for controller in controllers), \
        'Undo/Redo exposed a half-published source controller'
    for controller in controllers:
        outputs = tuple(controller['st_py_outputs'].objects)
        assert outputs, 'Accepted manifest has no generated outputs'
        assert all(obj.get('st_py_controller') == controller and obj.get('st_py_key')
                   for obj in outputs)
    return support.scene_state()


def held_update(controller, **kwargs):
    """Hold the first visible preview of component two, before validation."""
    global stage_held
    producer = original_update(controller, **kwargs)
    try:
        while True:
            try:
                command = next(producer)
            except StopIteration as done:
                return done.value
            active = driver.active
            visible = preview_support.previews(controller, visible=True)
            if active is not None and active.adapter.committed_units == 1 and visible:
                stage_held = True
                while stage_held:
                    yield core.YieldToUI()
            yield command
    finally:
        producer.close()


def start_import():
    global normal_started
    bpy.context.preferences.addons['blosm'].preferences.preferPerformanceDuringImport = args.performance_request
    if phase == 'normal completion':
        normal_started = time.monotonic()
    assert bpy.ops.blosm.import_data('EXEC_DEFAULT', True) == {'FINISHED'}
    assert driver.active is not None


def undo_redo(expected, *, cycles=2):
    for _ in range(cycles):
        assert bpy.ops.ed.undo() == {'FINISHED'}
        yield
        state = coherent()
        assert not state or state == expected, 'Undo restored altered accepted street geometry'
        assert bpy.ops.ed.redo() == {'FINISHED'}
        yield
        assert coherent() == expected, 'Redo changed accepted street geometry'


def scenario():
    global phase, other_scene, first_state, complete_state, normal_result, stage_held, normal_geometry_digest
    global normal_completed
    start_import()
    while driver.active is not None:
        yield
    normal_completed = time.monotonic()
    normal_result = dict(runtime.last_result)
    assert normal_result['outcome'] == 'success', normal_result
    complete_state = coherent()
    component_count = len(complete_state)
    assert component_count > 0, 'The extract produced no accepted street network'
    assert normal_result['committed'] == component_count, normal_result
    normal_geometry_digest = support.digest(complete_state)
    if args.expected_hash:
        assert normal_geometry_digest == args.expected_hash, \
            'GUI geometry differs from the background parity result'
    complete_beats = [row for row in heartbeats if row['phase'] == 'normal completion']
    assert len(complete_beats) >= 5, 'No meaningful timer heartbeat observation'
    assert any(row['visible_components'] > 0 for row in complete_beats), \
        'No visible published street geometry was observed'
    preview_beats = [row for row in complete_beats if row['visible_previews']]
    assert preview_beats, 'No actual visible early street preview was observed'
    first_accepted = next(row for row in complete_beats if row['visible_components'])
    assert preview_beats[0]['import_seconds'] < first_accepted['import_seconds'], \
        'First preview appeared only after component commit'
    assert preview_beats[0]['committed'] == 0
    assert normal_result['first_preview_seconds'] < normal_result['first_result_seconds']
    assert normal_result['preview_batches'] > 1
    if component_count > 1:
        assert any(0 < row['visible_components'] < component_count for row in complete_beats), \
            'No partially completed import was presented between components'
    cases.append('unmodified production timer: early visible previews before component commit')
    if args.osm:
        cases.append(dict(real_extract=str(args.osm.resolve()), accepted_components=component_count,
                          exact_geometry_sha256=normal_geometry_digest,
                          expected_hash_checked=bool(args.expected_hash)))
        return
    assert component_count == 3, 'Synthetic fixture component count changed'
    assert any(row['committed'] == 1 for row in complete_beats), 'First component was never presented'
    assert any(row['committed'] == 2 for row in complete_beats), 'Second component was never presented'
    phase = 'complete undo redo'
    yield from undo_redo(complete_state)
    cases.append('completed import: repeated native Undo/Redo preserves exact geometry')

    phase = 'early preview cancellation setup'
    checks.prepare()
    scene_runtime.NetworkScene.update_steps = staticmethod(held_update)
    start_import()
    while not stage_held:
        assert driver.active is not None, 'Import ended before early preview cancellation was exercised'
        yield
    phase = 'paused with early street previews'
    session = driver.active.adapter
    assert session.committed_units == 1
    first_state = support.scene_state()
    assert len(first_state) == 1
    visible_previews = preview_support.previews(visible=True)
    assert visible_previews, 'No real visible preview mesh was held'
    assert session.preview_count == len(visible_previews)
    assert all(obj.get(ownership.PENDING) and obj.data.get(ownership.PENDING)
               for obj in visible_previews)
    preview_controller = visible_previews[0]['st_py_controller']
    assert not preview_controller.get('st_py_manifest')
    assert not preview_controller['st_py_outputs'].objects
    assert pending(), 'Pending component has no explicit ownership'
    sentinel, note = checks.sentinel()
    sentinel_name, text_name = sentinel.name, note.name
    assert bpy.ops.ed.undo_push(message='User edit during pending street generation') == {'FINISHED'}
    target_scene = bpy.context.window.scene
    other_scene = bpy.data.scenes.new('User navigation scene')
    bpy.context.window.scene = other_scene
    before, revision = support.inventory(), session.revision
    deadline = time.monotonic() + .3
    while time.monotonic() < deadline:
        yield
    assert driver.active.paused
    assert support.inventory() == before, 'Paused import mutated Blender data'
    assert session.revision == revision
    bpy.context.window.scene = target_scene
    # Resume through the actual event loop before pressing the Cancel operator.
    deadline = time.monotonic() + .1
    while time.monotonic() < deadline:
        yield
    assert not driver.active.paused
    assert bpy.context.view_layer.objects.active == sentinel
    assert sentinel.select_get()
    cases.append('scene-switch pause/resume preserves pending allocation and user selection')
    phase = 'cancel early previews'
    before_cancel = time.monotonic()
    assert bpy.ops.blosm.cancel_progressive() == {'FINISHED'}
    latency = time.monotonic() - before_cancel
    assert driver.active is None
    assert runtime.last_result['outcome'] == 'cancelled', runtime.last_result
    assert runtime.last_result['committed'] == 1
    assert coherent() == first_state
    assert bpy.data.objects.get(sentinel_name) is not None
    assert bpy.data.texts.get(text_name) is not None
    assert not session.ownership.pending
    assert session.finish(silent=True) is None
    cases.append('early preview cancel preserves prior network and unrelated user object/text')
    cases.append(dict(cancel_latency_seconds=latency))
    # The user edit's undo snapshot captured pending source/staged IDs. Repeated
    # history restoration must remove those IDs before they become visible.
    phase = 'cancelled undo redo'
    yield from undo_redo(first_state, cycles=3)
    cases.append('cancelled pending snapshot: repeated native Undo/Redo has no orphan previews')


steps = scenario()


def complete(error=None):
    scene_runtime.NetworkScene.update_steps = staticmethod(original_update)
    try:
        if driver.active is not None:
            driver.active.close(silent=True)
    finally:
        runtime.OsmJob.step = original_job_step
        runtime.OsmJob.__init__ = original_job_init
        core.Job.step = original_core_step
        driver.Controller.tick = original_controller_tick
    normal_beats = [row for row in heartbeats if row['phase'] == 'normal completion']
    gaps = [second['seconds'] - first['seconds'] for first, second in zip(normal_beats, normal_beats[1:])]
    first_visible = next((row['import_seconds'] for row in normal_beats
                          if row['visible_components'] and row['import_seconds'] is not None), None)
    first_preview = next((row['import_seconds'] for row in normal_beats
                          if row['visible_previews'] and row['import_seconds'] is not None), None)
    report = dict(passed=error is None, blender=bpy.app.version_string, checkout=str(support.ROOT),
                  checks=cases, error=error, normal_result=normal_result,
                  source_path=str(support.FIXTURE.resolve()),
                  source_sha256=hashlib.sha256(support.FIXTURE.read_bytes()).hexdigest(),
                  performance_request=args.performance_request,
                  simulated_navigation=not args.no_navigation,
                  source_hashes_at_start=source_hashes, wall_timing=wall_scope_report(),
                  real_extract_only=bool(args.osm), expected_geometry_sha256=args.expected_hash,
                  normal_exact_geometry_sha256=normal_geometry_digest,
                  normal_first_visible_scene_seconds=first_visible,
                  normal_first_visible_preview_seconds=first_preview,
                  normal_preview_to_commit_seconds=(first_visible - first_preview
                    if first_visible is not None and first_preview is not None else None),
                  normal_peak_visible_preview_meshes=max((row['visible_previews'] for row in normal_beats), default=0),
                  visibility_measurement='Visible preview meshes and accepted components sampled from the UI heartbeat; not pixel timing',
                  normal_worst_heartbeat_gap_seconds=max(gaps, default=None),
                  heartbeats=heartbeats, result_notifications=results)
    (args.artifacts / 'ui_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('STREET_PROGRESSIVE_UI_CHECKS', json.dumps(dict(passed=report['passed'], error=error)), flush=True)
    bpy.ops.wm.quit_blender()


def heartbeat():
    global last_progress
    try:
        assert time.monotonic() - started < 600, 'Street GUI check timed out'
        if not args.no_navigation:
            for area in bpy.context.screen.areas:
                if area.type == 'VIEW_3D':
                    area.spaces.active.region_3d.view_distance += .001
                    area.tag_redraw()
        active = driver.active
        visible = 0
        if phase == 'normal completion':
            layer = bpy.context.view_layer
            visible = sum(any(obj.type == 'MESH' and len(obj.data.polygons)
                              and obj.name in layer.objects and obj.visible_get(view_layer=layer)
                              for obj in controller['st_py_outputs'].objects)
                          for controller in support.controllers(committed=True))
        heartbeats.append(dict(phase=phase, seconds=time.monotonic() - started,
            import_seconds=(time.monotonic() - normal_started if normal_started is not None else None),
            visible_components=visible,
            visible_previews=len(preview_support.previews(visible=True)),
            session_phase=active.adapter.phase if active else None,
            committed=active.adapter.committed_units if active else
            (runtime.last_result['committed'] if runtime.last_result else 0),
            paused=bool(active and active.paused)))
        now = time.monotonic()
        if now - last_progress >= 5.:
            progress = dict(phase=phase, session_phase=active.adapter.phase if active else None,
                            committed=heartbeats[-1]['committed'], paused=heartbeats[-1]['paused'],
                            elapsed_seconds=now - started,
                            import_elapsed_seconds=now - normal_started if normal_started is not None else None,
                            performance_request=args.performance_request,
                            simulated_navigation=not args.no_navigation,
                            wall_timing=wall_scope_report(compact=True))
            (args.artifacts / 'ui_progress.json').write_text(json.dumps(progress, indent=2), encoding='utf-8')
            print('STREET_PROGRESSIVE_UI_PROGRESS', json.dumps(progress), flush=True)
            last_progress = now
        next(steps)
    except StopIteration:
        complete()
        return None
    except Exception:
        complete(traceback.format_exc())
        return None
    return .01


bpy.app.timers.register(heartbeat, first_interval=.1)
