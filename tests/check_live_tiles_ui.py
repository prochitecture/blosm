"""Actual UI/timer/undo/cancel checks with live Helsinki and real GLB failures."""
import argparse
import importlib
import json
from pathlib import Path
import sys
import time
import traceback

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gltf_test_support import bpy, REPO, reset_scene
from live_tile_datasets import DATASETS
import addon_utils

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)
package = REPO.name
addon_utils.enable(package, default_set=True)
runtime = importlib.import_module(package + '.progressive.blender')
transport = importlib.import_module(package + '.downloads.transport')
popups, processes, heartbeats, checks = [], [], [], []
show_result = runtime.show_result
def record_popup(level, message):
    popups.append((level, message))
    show_result(level, message)
runtime.show_result = record_popup
phase = 'start'
started = time.monotonic()
saved = None
cancelled_count = 0


def meshes():
    objects = [obj for obj in bpy.context.scene.objects if obj.type == 'MESH']
    return len(objects), sum(len(obj.data.vertices) for obj in objects)


def begin(name, join=False):
    config = DATASETS[name]
    settings = bpy.context.scene.blosm
    settings.dataType = '3d-tiles'
    settings.threedTilesSource = 'custom'
    settings.threedTilesUrl = config['url']
    settings.minLon, settings.minLat, settings.maxLon, settings.maxLat = config['bounds']
    settings.geometricError = config['geometric_error']
    settings.progressiveImport = True
    settings.concurrentDownloads = True
    settings.downloadWorkers = 8
    settings.join3dTilesObjects = join
    assert bpy.ops.blosm.import_data('EXEC_DEFAULT', True) == {'FINISHED'}
    assert runtime.active is not None
    processes.append(runtime.active.adapter.job.transport.process)


def finish(status, error=None):
    if runtime.active:
        runtime.active.close(silent=True)
    for process in processes:
        process.wait(timeout=5)
    assert not transport._live and not transport.reap()
    report = {'blender': bpy.app.version_string, 'status': status, 'live_network': True,
              'workers': 8, 'checks': checks, 'popups': popups, 'heartbeats': heartbeats,
              'helper_exit_codes': [process.returncode for process in processes], 'error': error}
    (args.artifacts / 'live_ui_checks.json').write_text(json.dumps(report, indent=2))
    print('LIVE_UI_CHECKS', status, error or '', flush=True)
    addon_utils.disable(package, default_set=True)
    bpy.ops.wm.quit_blender()


def tick():
    global phase, saved, cancelled_count
    try:
        assert time.monotonic() - started < 180, 'Live UI check timed out'
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                area.spaces.active.region_3d.view_distance += .001
                area.tag_redraw()
        count, vertices = meshes()
        heartbeats.append((phase, count, time.monotonic() - started))
        if phase == 'start':
            reset_scene()
            bpy.ops.ed.undo_push(message='Live import baseline')
            begin('helsinki_dense')
            phase = 'complete'
        elif phase == 'complete' and runtime.active is None:
            assert count > 1 and vertices > 0, popups
            observed = {n for name, n, _ in heartbeats if name == 'complete'}
            assert any(0 < n < count for n in observed), 'No intermediate tiles were observed'
            assert len(popups) == 1 and popups[-1][0] == 'INFO'
            saved = meshes()
            bpy.ops.ed.undo()
            assert meshes() == (0, 0)
            bpy.ops.ed.redo()
            assert meshes() == saved and runtime.active is None
            checks.append('live Helsinki tiles appear incrementally while viewport timer runs; one Undo/Redo covers import')
            reset_scene()
            begin('helsinki_dense', join=True)
            phase = 'join'
        elif phase == 'join' and runtime.active is None:
            assert count == 1 and vertices == saved[1], popups
            bpy.ops.ed.undo()
            assert meshes() == (0, 0)
            bpy.ops.ed.redo()
            assert meshes() == (1, saved[1]) and runtime.active is None
            checks.append('live joined import shares one Undo/Redo step')
            reset_scene()
            begin('helsinki_dense')
            phase = 'cancel'
        elif phase == 'cancel' and count > 0:
            assert runtime.active is not None, 'Import finished before cancellation was exercised'
            cancelled_count = count
            assert bpy.ops.blosm.cancel_progressive() == {'FINISHED'}
            assert runtime.active is None and meshes()[0] == cancelled_count
            checks.append('Cancel with live requests retains committed tiles and terminates helper')
            reset_scene()
            begin('svirnas')
            phase = 'unsupported'
        elif phase == 'unsupported' and runtime.active is None:
            assert meshes() == (0, 0)
            assert len(popups) == 4, popups
            assert popups[-1][0] == 'ERROR' and '18 errors' in popups[-1][1], popups[-1]
            checks.append('18 real unsupported GLB tiles produce one aggregate UI error and no partial objects')
            finish('passed')
            return None
    except Exception:
        finish('failed', traceback.format_exc())
        return None
    return .025


bpy.app.timers.register(tick, first_interval=.3)
