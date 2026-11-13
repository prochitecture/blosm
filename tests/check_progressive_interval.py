"""Native Controller interval/lifecycle checks in an isolated Blender process.

Launch with --background --factory-startup --disable-autoexec and the chat's
five isolated resource overrides. The required artifact directory also becomes
Blender's temporary directory; the user's normal profile is never used.
"""
import argparse
import importlib
import json
from pathlib import Path
import sys
import traceback
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))

import bpy

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)

addon = importlib.import_module(ROOT.name)
addon.register()
runtime = importlib.import_module(ROOT.name+'.progressive.blender')
window = bpy.context.window
assert window is not None, 'These checks require the native factory-startup window context'
scene = window.scene
view_layer = window.view_layer
checks = []


class FakeJob:
    def __init__(self):
        self.done = False
        self.error = None
        self.steps = self.polls = 0
        self.action = None

    def step(self):
        self.steps += 1
        if self.action is not None:
            self.action()

    def poll(self):
        self.polls += 1


class FakeAdapter:
    def __init__(self):
        self.scene, self.view_layer = scene, view_layer
        self.revision = 0
        self.job = FakeJob()
        self.cancelled = 0
        self.finished = []
        self.status = 'Synthetic processing'

    def cancel(self):
        self.cancelled += 1
        self.job.done = True

    def finish(self,silent=False):
        self.finished.append(silent)


def start(adapter):
    assert runtime.active is None
    controller = runtime.start(adapter,bpy.context)
    assert runtime.active is controller
    assert bpy.app.timers.is_registered(runtime._tick)
    return controller


def stop():
    if runtime.active is not None:
        runtime.active.close(silent=True)
    assert not bpy.app.timers.is_registered(runtime._tick)


def check_hint(value,expected):
    adapter = FakeAdapter()
    adapter.processing_interval = value
    controller = start(adapter)
    try:
        assert runtime._tick() == expected, (value,controller.next_interval)
        assert controller.next_interval == expected
        assert adapter.job.steps == 1 and adapter.job.polls == 0
    finally:
        stop()


report = {}
try:
    adapter = FakeAdapter()  # Legacy adapters do not expose the optional hint.
    controller = start(adapter)
    assert runtime._tick() == .03
    assert adapter.job.steps == 1
    stop()
    checks.append('missing adapter hint preserves the 30 ms default')

    for value, expected in ((.005,.005),(.02,.02),(.03,.03),(.001,.005),
                            (0,.005),(-1,.005),(1,.03),(10**1000,.03),(-10**1000,.005)):
        check_hint(value,expected)
    checks.append('finite numeric hints clamp to 5–30 ms, including huge integers')

    for value in (True,False,None,'0.005',object(),float('nan'),float('inf'),float('-inf')):
        check_hint(value,.03)
    checks.append('invalid, boolean and nonfinite hints fall back to 30 ms')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    assert runtime._tick() == .005
    def revise():
        adapter.revision += 1
    adapter.job.action = revise
    with patch.object(controller,'checkpoint',wraps=controller.checkpoint) as checkpoint:
        assert runtime._tick() == .03
        assert checkpoint.call_count == 1
    adapter.job.action = None
    assert runtime._tick() == .005
    stop()
    checks.append('a scene revision uses the default interval and normal checkpoint')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    other_scene = bpy.data.scenes.new('Interval pause scene')
    other_layer = scene.view_layers.new('Interval pause layer')
    try:
        assert runtime._tick() == .005
        before = adapter.job.steps
        window.scene = other_scene
        assert runtime._tick() == .03 and controller.paused
        assert adapter.job.steps == before and adapter.job.polls == 1
        window.scene, window.view_layer = scene, view_layer
        assert runtime._tick() == .005 and not controller.paused
        before = adapter.job.steps
        window.view_layer = other_layer
        assert runtime._tick() == .03 and controller.paused
        assert adapter.job.steps == before and adapter.job.polls == 2
        window.view_layer = view_layer
        assert runtime._tick() == .005 and not controller.paused
        mesh_object = next(obj for obj in scene.objects if obj.type=='MESH')
        bpy.ops.object.select_all(action='DESELECT')
        mesh_object.select_set(True)
        view_layer.objects.active = mesh_object
        bpy.ops.object.mode_set(mode='EDIT')
        before = adapter.job.steps
        assert runtime._tick() == .03 and controller.paused
        assert adapter.job.steps == before and adapter.job.polls == 3
        bpy.ops.object.mode_set(mode='OBJECT')
        assert runtime._tick() == .005 and not controller.paused
    finally:
        if bpy.context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        window.scene, window.view_layer = scene, view_layer
        stop()
        scene.view_layers.remove(other_layer)
        bpy.data.scenes.remove(other_scene)
    checks.append('scene, view-layer and Edit Mode pauses reset the interval and only poll')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    adapter.job.action = lambda: setattr(adapter.job,'done',True)
    assert runtime._tick() is None
    assert runtime.active is None and not bpy.app.timers.is_registered(runtime._tick)
    assert adapter.finished == [False] and adapter.cancelled == 0
    checks.append('completion closes the native timer without cancellation')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    assert runtime._tick() == .005
    assert bpy.ops.blosm.cancel_progressive() == {'FINISHED'}
    assert runtime.active is None and not bpy.app.timers.is_registered(runtime._tick)
    assert adapter.cancelled == 1 and adapter.finished == [False]
    assert runtime._tick() is None
    checks.append('native cancel operator closes the timer and retains normal finalization')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    def fail():
        raise RuntimeError('Synthetic interval failure')
    adapter.job.action = fail
    assert runtime._tick() is None
    assert runtime.active is None and not bpy.app.timers.is_registered(runtime._tick)
    assert isinstance(adapter.job.error,RuntimeError)
    assert adapter.cancelled == 1 and adapter.finished == [False]
    checks.append('step failure records the error and removes the timer')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    controller.window = object()  # An owning window no longer present in Blender.
    assert runtime._tick() is None
    assert runtime.active is None and not bpy.app.timers.is_registered(runtime._tick)
    assert adapter.cancelled == 1 and adapter.finished == [True]
    checks.append('missing owning window silently cancels and removes the timer')

    adapter = FakeAdapter()
    adapter.processing_interval = .005
    controller = start(adapter)
    adapter.job.action = lambda: controller.close(silent=True)
    assert runtime._tick() is None
    assert runtime.active is None and not bpy.app.timers.is_registered(runtime._tick)
    assert adapter.cancelled == 1 and adapter.finished == [True]
    checks.append('closure during a step cannot revive the accelerated timer')

    report = {'status':'passed','checks':checks,'blender':bpy.app.version_string}
except Exception:
    report = {'status':'failed','checks':checks,'blender':bpy.app.version_string,
              'error':traceback.format_exc()}
    raise
finally:
    stop()
    (args.artifacts/'interval.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PROGRESSIVE_INTERVAL',json.dumps(report),flush=True)
    addon.unregister()
