"""Run in a separate factory-startup Blender GUI process, not --background.

Real timers, download subprocess, undo and scene editing are exercised while
the viewport changes. No preferences or .blend files are saved.
"""
import argparse
import base64
import importlib
import json
from math import radians
from pathlib import Path
import subprocess
import sys
import time
import traceback

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gltf_test_support import bpy, REPO, fixture, make_glb, make_b3dm, reset_scene
import addon_utils

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)
sys.path.insert(0, str(REPO.parent))
package = REPO.name
addon_utils.enable(package, default_set=True)
manager_module = importlib.import_module(package + '.threed_tiles.manager')
runtime = importlib.import_module(package + '.progressive.blender')
transport_module = importlib.import_module(package + '.progressive.transport')
from io_scene_gltf2.io.imp.gltf2_io_gltf import glTFImporter

center = manager_module.BaseManager.fromGeographic(60.17, 24.94, 0)
doc, binary, glb = fixture()
b3dm = make_b3dm(glb, list(center))
doc['nodes'][0]['translation'] = [center.x, center.z, -center.y]
glb = make_glb(doc, binary)
region = [radians(24.93), radians(60.16), radians(24.95), radians(60.18), -100, 100]
tileset = {'root': {'geometricError': 0, 'boundingVolume': {'region': region}, 'children': [
    {'geometricError': 0, 'boundingVolume': {'region': region}, 'content': {'uri': 'one.glb'}},
    {'geometricError': 0, 'boundingVolume': {'region': region}, 'content': {'uri': 'two.b3dm'}}]}}
payloads = {'root.json': json.dumps(tileset).encode(), 'one.glb': glb, 'two.b3dm': b3dm}
server = subprocess.Popen([transport_module.python_executable(), '-I', '-B',
    str(Path(__file__).with_name('progressive_http_fixture.py'))],
    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
server.stdin.write((json.dumps({'payloads': {key: base64.b64encode(value).decode()
                                           for key, value in payloads.items()},
                               'delays': {'root.json': .3, 'one.glb': .4, 'two.b3dm': .8}}) + '\n').encode())
server.stdin.flush()
port = int(server.stdout.readline())
server.stdin.close()
server.stdout.close()
scene = bpy.context.scene
settings = scene.blosm
settings.dataType = '3d-tiles'
settings.threedTilesSource = 'custom'
settings.threedTilesUrl = 'http://127.0.0.1:%s/root.json' % port
settings.minLon, settings.maxLon = 24.93, 24.95
settings.minLat, settings.maxLat = 60.16, 60.18
settings.geometricError = 1
settings.progressiveImport = True
settings.join3dTilesObjects = False
phase = 'start'
checks = []
heartbeats = []
started = time.monotonic()
stage_started = started
processes = []
hold_until = 0
paused_once = False
extra_scene = None


def start_import():
    global stage_started
    assert bpy.ops.blosm.import_data() == {'FINISHED'}
    assert runtime.active is not None
    processes.append(runtime.active.adapter.job.transport.process)
    stage_started = time.monotonic()


def complete(status, error=None):
    if runtime.active is not None:
        runtime.active.close(silent=True)
    server.kill()
    server.wait(timeout=5)
    for process in processes:
        process.wait(timeout=5)
    report = {'blender': bpy.app.version_string, 'status': status, 'checks': checks,
              'heartbeat_count': len(heartbeats), 'error': error}
    (args.artifacts / 'ui_checks.json').write_text(json.dumps(report, indent=2))
    print('PROGRESSIVE_UI_CHECKS', json.dumps(report), flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    global phase, hold_until, paused_once, extra_scene, scene
    try:
        assert time.monotonic() - started < 50, 'UI test timed out'
        assert not hasattr(glTFImporter, '_patch_convert_functions')
        # Every heartbeat runs through Blender's real event loop and moves the
        # viewport while the helper waits on deliberately delayed HTTP replies.
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                area.spaces.active.region_3d.view_distance += .001
                area.tag_redraw()
        count = sum(obj.type == 'MESH' for obj in scene.objects)
        heartbeats.append((phase, count, time.monotonic() - started))
        if phase == 'start':
            reset_scene()
            start_import()
            phase = 'complete'
        elif phase == 'complete' and runtime.active is None:
            assert count == 2
            assert sum(p == 'complete' and n == 0 for p, n, t in heartbeats) > 3
            assert sum(p == 'complete' and n == 1 for p, n, t in heartbeats) > 3
            checks.append('real HTTP helper + viewport timer heartbeats + incremental tiles')
            reset_scene()
            start_import()
            phase = 'cancel'
        elif phase == 'cancel' and count == 1:
            assert bpy.ops.blosm.cancel_progressive() == {'FINISHED'}
            assert runtime.active is None
            assert sum(obj.type == 'MESH' for obj in scene.objects) == 1
            checks.append('Cancel button keeps first tile and stops helper')
            reset_scene()
            start_import()
            phase = 'undo'
        elif phase == 'undo' and count == 1:
            bpy.ops.ed.undo()
            assert runtime.active is None
            scene = bpy.context.scene  # Undo replaces Blender data-block wrappers.
            assert not [obj for obj in scene.objects if obj.type == 'MESH']
            checks.append('Undo cancels active operation before restoring scene')
            start_import()
            phase = 'edit'
        elif phase == 'edit' and count >= 1 and not paused_once:
            obj = next(obj for obj in scene.objects if obj.type == 'MESH')
            bpy.context.view_layer.objects.active = obj
            obj.select_set(True)
            bpy.ops.object.mode_set(mode='EDIT')
            hold_until = time.monotonic() + .5
            paused_once = True
        elif phase == 'edit' and paused_once and time.monotonic() >= hold_until:
            assert runtime.active is not None and runtime.active.paused
            bpy.ops.object.mode_set(mode='OBJECT')
            checks.append('Edit Mode pauses scene commits without blocking the UI')
            phase = 'resume'
        elif phase == 'resume' and runtime.active is None:
            checks.append('Object Mode resumes and completes import')
            reset_scene()
            start_import()
            phase = 'scene'
            extra_scene = bpy.data.scenes.new('Other scene')
            bpy.context.window.scene = extra_scene
            hold_until = time.monotonic() + .3
        elif phase == 'scene' and time.monotonic() >= hold_until:
            assert runtime.active is not None and runtime.active.paused
            bpy.context.window.scene = scene
            runtime.stop_for_scene_change()  # same registered load_pre callback
            assert runtime.active is None
            bpy.data.scenes.remove(extra_scene)
            checks.append('scene switch pauses; load/undo handler closes pending work')
            start_import()
            addon_utils.disable(package, default_set=True)
            assert runtime.active is None
            checks.append('disabling addon cancels helper and unregisters timer')
            complete('passed')
            return None
        return .03
    except Exception:
        complete('failed', traceback.format_exc())
        return None


bpy.app.timers.register(tick, first_interval=.5)
