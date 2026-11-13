"""Offline native visibility, import parity, persistence and scale checks.

Launch with factory startup, disabled autoexec and the chat's isolated profile.
Optional --ui exercises real progressive timers and Undo/Redo as well.
"""
import argparse
import cProfile
import importlib
import json
from pathlib import Path
import sys
import time
import traceback

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT.parent), str(ROOT/'tests')]
import bpy
from mapillary_test_support import MemoryTransport, BOUNDS

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--count', type=int, default=10000)
parser.add_argument('--ui', action='store_true')
args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type = 'NONE'
addon = importlib.import_module(ROOT.name)
addon.register()
filters = importlib.import_module(ROOT.name+'.mapillary.filters')
adapter = importlib.import_module(ROOT.name+'.mapillary.progressive')
runtime = importlib.import_module(ROOT.name+'.progressive.blender')
Projection = importlib.import_module(ROOT.name+'.util.transverse_mercator').TransverseMercator
runtime.show_result = lambda *_: None
adapter.notify_result = lambda *_: None
checks, metrics = [], {}


def check(label, condition):
    assert condition, label
    checks.append(label)


def settings():
    return bpy.context.view_layer.blosm_mapillary_filters


def task(count=96, reverse=True):
    return adapter.MapillaryImport(bpy.context, BOUNDS, Projection(lat=60.17, lon=24.94),
        (60.17,24.94), False, 'fixture-token', transport=MemoryTransport(count, reverse=reverse))


def matches(obj):
    return filters.compiled(settings().applied).matches(
        obj.get('mapillary_captured_at'), obj.get('mapillary_camera_model'))


def apply():
    return bpy.ops.blosm.apply_mapillary_filters('EXEC_DEFAULT', True)


def reset():
    return bpy.ops.blosm.reset_mapillary_filters('EXEC_DEFAULT', True)


def parity(t):
    return all(obj.hide_get() == (not matches(obj)) for obj in t.collection.objects)


def report(status, error=None):
    (args.artifacts/'filters.json').write_text(json.dumps(dict(status=status, error=error,
        checks=checks, metrics=metrics, blender=bpy.app.version_string), indent=2), encoding='utf-8')
    print('FILTER_CHECKS', status, len(checks), metrics, flush=True)


def basic():
    scene, layer = bpy.context.scene, bpy.context.view_layer
    original = (layer.objects.active, layer.active_layer_collection, scene.camera,
                tuple(scene.cursor.location), tuple(bpy.context.selected_objects))
    t = task()
    t.run_blocking()
    objects = list(t.collection.objects)
    manual = next(o for o in objects if o['mapillary_camera_model']=='perspective')
    manual.hide_set(True)
    s = settings()
    s.ordinary = s.fisheye = s.other = False
    s.date_from = s.date_to = '2023-11-14'
    apply()
    check('date/type filter hides existing nonmatching cameras', parity(t))
    check('ordinary scene camera and render flags are unchanged', not scene.camera.hide_get()
          and all(not o.hide_render and not o.hide_viewport for o in objects))
    check('unrelated active/selection/collection/cursor/render camera are preserved', original ==
          (layer.objects.active, layer.active_layer_collection, scene.camera,
           tuple(scene.cursor.location), tuple(bpy.context.selected_objects)))
    owned = next(item.object for item in s.hidden)
    owned.name = 'Renamed filtered camera'
    removed = next(item.object for item in s.hidden if item.object != owned)
    bpy.data.objects.remove(removed, do_unlink=True)
    reset()
    check('reset tolerates rename/deletion and preserves manual hiding',
          manual.hide_get() and not owned.hide_get() and len(s.hidden)==0)
    manual.hide_set(False)
    # Exercise the optimized native path with unrelated hidden objects,
    # disabled cameras, non-selectable cameras and a hidden collection.
    light = bpy.data.objects['Light']
    light.hide_set(True)
    disabled = next(o for o in t.collection.objects if o != manual and o['mapillary_camera_model']=='perspective')
    locked = next(o for o in t.collection.objects if o not in (manual,disabled) and o['mapillary_camera_model']=='perspective')
    disabled.hide_viewport = True
    locked.hide_select = True
    selected_camera = next(o for o in t.collection.objects if o['mapillary_camera_model']=='spherical')
    selected_camera.select_set(True)
    selection_flag = selected_camera.select_get()
    s.ordinary = s.fisheye = s.other = False
    apply()
    check('batch filtering respects disabled and non-selectable cameras', disabled.hide_get() and locked.hide_get()
          and disabled.hide_viewport and locked.hide_select)
    check('batch filtering preserves selection of matching cameras', selected_camera.select_get()==selection_flag)
    reset()
    check('batch reset keeps unrelated manual hiding and global flags', light.hide_get()
          and disabled.hide_viewport and locked.hide_select and not disabled.hide_get() and not locked.hide_get())
    selected_camera.select_set(False)
    disabled.hide_viewport = locked.hide_select = False
    light.hide_set(False)
    s.ordinary = s.fisheye = s.other = False
    apply()
    previous = s.applied
    before = {o.name:o.hide_get() for o in t.collection.objects}
    s.date_from = '2024-02-30'
    try:
        result = apply()
    except RuntimeError:  # Blender raises for an ERROR report in Python calls.
        result = {'CANCELLED'}
    check('invalid date is rejected atomically', result=={'CANCELLED'} and s.applied==previous
          and before=={o.name:o.hide_get() for o in t.collection.objects})
    s.date_from = ''
    # Pending GUI edits must not change the predicate for new cameras.
    s.ordinary, s.panoramic = True, False
    pending = task(reverse=False)
    pending.run_blocking()
    check('unapplied GUI edits do not leak into incoming cameras', parity(pending))
    apply()
    progressive = task()
    while not progressive.job.done:
        progressive.job.step(max_steps=7)
        if progressive.collection:
            assert parity(progressive)
    check('progressive commits obey the same predicate immediately', progressive.count==96 and parity(progressive))
    block = task(reverse=False)
    block.run_blocking()
    snapshot = lambda t: sorted((o['mapillary_image_id'], o.hide_get()) for o in t.collection.objects)
    check('blocking/progressive visibility parity with out-of-order metadata', snapshot(progressive)==snapshot(block))
    other = scene.view_layers.new('Unfiltered view layer')
    check('filters do not affect other view layers', all(not o.hide_get(view_layer=other) for o in block.collection.objects)
          and not other.blosm_mapillary_filters.applied)
    other.name = 'Renamed view layer'
    # Save/reload exercises RNA settings and owned object references together.
    s.expanded = True
    scene.blosm.dataType = 'mapillary'
    path = args.artifacts/'filtered_cameras.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(path))
    bpy.ops.wm.open_mainfile(filepath=str(path))
    check('applied snapshot and ownership survive save/reload', settings().applied
          and len(settings().hidden)>0 and all(item.object.hide_get() for item in settings().hidden))
    reset()
    check('reset after reload restores owned visibility', not settings().applied and not len(settings().hidden)
          and all(not o.hide_get() for o in bpy.context.view_layer.objects if filters.is_imported_camera(o)))
    if args.ui:
        return  # Let Blender rebuild its window context after opening the file.
    # Scale cameras are created by the production path, not test-only objects.
    settings().ordinary = settings().fisheye = settings().other = False
    apply()
    scale = task(args.count)
    started = time.perf_counter()
    max_turn = 0.
    while not scale.job.done:
        turn = time.perf_counter()
        scale.job.step()
        max_turn = max(max_turn, time.perf_counter()-turn)
    assert scale.job.error is None, scale.job.error
    metrics['scale_import_seconds'] = time.perf_counter()-started
    metrics['max_filtered_import_turn_seconds'] = max_turn
    check('large incoming camera set is filtered at creation', parity(scale))
    settings().ordinary, settings().panoramic = True, False
    profiler = cProfile.Profile()
    started = time.perf_counter()
    profiler.enable()
    apply()
    profiler.disable()
    metrics['apply_seconds'] = time.perf_counter()-started
    profiler.dump_stats(str(args.artifacts/'apply.prof'))
    metrics['camera_count'] = sum(filters.is_imported_camera(o) for o in bpy.context.view_layer.objects)
    check('large camera set has exact visibility without deleting objects', scale.count==args.count and parity(scale))
    started = time.perf_counter()
    reset()
    metrics['reset_seconds'] = time.perf_counter()-started
    check('large camera set restores visibility', all(not o.hide_get() for o in scale.collection.objects))


state = {'phase':'start', 'started':time.monotonic(), 'ticks':0}


def ui_tick():
    try:
        assert time.monotonic()-state['started'] < 180, 'Filter UI check timeout'
        state['ticks'] += 1
        phase = state['phase']
        if phase == 'start':
            bpy.context.window.event_simulate(type='ESC', value='PRESS')
            bpy.context.window.event_simulate(type='ESC', value='RELEASE')
            bpy.context.preferences.edit.use_global_undo = True
            settings().ordinary = settings().fisheye = settings().other = False
            bpy.ops.ed.undo_push(message='Before camera filters')
            apply()
            check('native Apply owns hidden flags', len(settings().hidden)>0)
            bpy.ops.ed.undo()
            check('native Undo restores visibility and applied settings', not settings().applied and not len(settings().hidden))
            bpy.ops.ed.redo()
            check('native Redo restores filter and visibility', bool(settings().applied)
                  and all(item.object.hide_get() for item in settings().hidden))
            state['task'] = task(640)
            runtime.start(state['task'], bpy.context)
            state['phase'] = 'import'
        elif phase == 'import':
            t = state['task']
            if t.count and not t.job.done:
                check('native timer creates filtered cameras progressively', parity(t))
                settings().ordinary, settings().panoramic = True, False
                apply()
                check('filter change updates already imported cameras during import', parity(t))
                state['phase'] = 'changed'
        elif phase == 'changed':
            t = state['task']
            assert parity(t)
            if runtime.active is None:
                check('later timer commits use updated filter', t.count==640 and t.job.error is None)
                state['phase'], state['ticks'] = 'screenshot', 0
                for area in bpy.context.screen.areas:
                    if area.type == 'VIEW_3D':
                        area.spaces.active.show_region_ui = True
        elif phase == 'screenshot' and state['ticks'] == 4:
            area = next(a for a in bpy.context.screen.areas if a.type=='VIEW_3D')
            region = next(r for r in area.regions if r.type=='UI')
            region.active_panel_category = 'Blosm'
            area.tag_redraw()
        elif phase == 'screenshot' and state['ticks'] > 9:
            area = next(a for a in bpy.context.screen.areas if a.type=='VIEW_3D')
            check('filter controls are displayed in the actual Blosm sidebar', settings().expanded
                  and next(r for r in area.regions if r.type=='UI').active_panel_category=='Blosm')
            bpy.ops.screen.screenshot(filepath=str(args.artifacts/'filters.png'))
            # Import checkpoints after Apply are a separate undo group; Undo
            # must stop further work and preserve the intervening filter edit.
            bpy.ops.ed.undo()
            check('undo after progressive completion keeps intervening filter edit',
                  filters.compiled(settings().applied).categories==frozenset({'ordinary'}))
            t = task(640)
            state['task'] = t
            runtime.start(t, bpy.context)
            state['phase'] = 'cancel'
        elif phase == 'cancel' and state['task'].count:
            t = state['task']
            runtime.active.close()
            check('cancellation retains completed filtered cameras', 0<t.count<640 and parity(t))
            reset()
            check('reset also reveals completed cancelled-import cameras', all(not o.hide_get() for o in t.collection.objects))
            report('passed')
            addon.unregister()
            bpy.ops.wm.quit_blender()
            return None
        return .03
    except Exception:
        report('failed', traceback.format_exc())
        if runtime.active is not None:
            runtime.active.close(silent=True)
        addon.unregister()
        bpy.ops.wm.quit_blender()
        return None


try:
    basic()
    if args.ui:
        assert not bpy.app.background
        state['started'] = time.monotonic()
        bpy.app.timers.register(ui_tick, first_interval=.3)
    else:
        report('passed')
        addon.unregister()
except Exception:
    report('failed', traceback.format_exc())
    raise
