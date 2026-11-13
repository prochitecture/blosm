"""Real isolated Blender UI timers, HTTP helper, grouped undo and scale."""
import argparse
import ctypes
import importlib
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback
from unittest.mock import patch

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT.parent),str(ROOT/'tests')]
import bpy
addon=importlib.import_module(ROOT.name)
addon.register()
module=importlib.import_module(ROOT.name+'.mapillary.progressive')
runtime=importlib.import_module(ROOT.name+'.progressive.blender')
transport_module=importlib.import_module(ROOT.name+'.downloads.transport')
Projection=importlib.import_module(ROOT.name+'.util.transverse_mercator').TransverseMercator
from mapillary_test_support import BOUNDS

parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
parser.add_argument('--count',type=int,default=5000)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
gate=args.artifacts/'continue_metadata'
gate.touch()
checks=[]
messages=[]
processes=[]
servers=[]
tasks=[]
progress_examples={}
heartbeat=[]
import_gaps=[]
step_times=[]
checkpoint_times=[]
history_times=[]
previous_active=False
adapter=None
server=None
started=time.monotonic()
runtime.show_result=lambda *args: messages.append(args)
original_checkpoint=runtime.Controller.checkpoint
def timed_checkpoint(self,*a,**kw):
    start=time.monotonic()
    try: return original_checkpoint(self,*a,**kw)
    finally: checkpoint_times.append(time.monotonic()-start)
runtime.Controller.checkpoint=timed_checkpoint


def memory_bytes():
    if sys.platform!='win32': return None
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong),
                  ('peak',ctypes.c_size_t),('working',ctypes.c_size_t),
                  ('poolpeak',ctypes.c_size_t),('pool',ctypes.c_size_t),
                  ('nonpoolpeak',ctypes.c_size_t),('nonpool',ctypes.c_size_t),
                  ('page',ctypes.c_size_t),('pagepeak',ctypes.c_size_t)]
    counters=Counters()
    counters.cb=ctypes.sizeof(counters)
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.GetCurrentProcess.restype=ctypes.c_void_p
    psapi=ctypes.WinDLL('psapi',use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
    if psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(counters),counters.cb):
        return counters.working


def launch_server(count, **options):
    process=subprocess.Popen([transport_module.python_executable(),'-I','-B',
        str(ROOT/'tests/mapillary_http_fixture.py')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    config={'count':count,'gate':str(gate),'stats':str(args.artifacts/('http-%s.json'%len(servers))),**options}
    process.stdin.write((json.dumps(config)+'\n').encode())
    process.stdin.close()
    port=int(process.stdout.readline())
    process.stdout.close()
    servers.append(process)
    return port


def settings():
    scene=bpy.context.scene
    scene.blosm.minLon,scene.blosm.minLat,scene.blosm.maxLon,scene.blosm.maxLat=BOUNDS
    scene.blosm.dataType='mapillary'
    scene.blosm.concurrentDownloads=True
    scene.blosm.downloadWorkers=4
    scene['height_offset']=20.
    return scene


def clear():
    for obj in list(bpy.data.objects):
        if 'mapillary_image_id' in obj:
            data=obj.data
            bpy.data.objects.remove(obj,do_unlink=True)
            if data.users==0: bpy.data.cameras.remove(data)
    for collection in list(bpy.data.collections):
        if collection.get('source')=='Mapillary': bpy.data.collections.remove(collection)
    for obj in list(bpy.context.scene.objects):
        if obj.name.startswith('User edit'):
            bpy.data.objects.remove(obj,do_unlink=True)


def count():
    return sum('mapillary_image_id' in obj for obj in bpy.context.scene.objects)


def begin(port, *, progressive=True):
    global adapter
    settings()
    bpy.context.scene.blosm.progressiveImport=progressive
    factory=module.MapillaryImport
    def create(*a,**kw):
        global adapter
        kw.update(tile_url='http://127.0.0.1:%s/tiles/{z}/{x}/{y}'%port,
                  graph_url='http://127.0.0.1:%s/images'%port)
        # Exercise the actual import operator without persisting a credential in
        # Blender preferences or the normal user profile.
        a=list(a)
        a[5]='fixture-client-token'
        adapter=factory(*a,**kw)
        assert 'Estimating…' in adapter.status and '0 cameras created' in adapter.status
        progress_examples.setdefault('initial',adapter.status)
        original_step=adapter.job.step
        def timed_step(*a,**kw):
            start=time.monotonic()
            try: return original_step(*a,**kw)
            finally: step_times.append(time.monotonic()-start)
        adapter.job.step=timed_step
        tasks.append(adapter)
        processes.append(adapter.transport.process)
        return adapter
    with patch.object(module,'MapillaryImport',create):
        result=bpy.ops.blosm.import_data('EXEC_DEFAULT',True)
    assert result=={'FINISHED'},result
    if progressive: assert runtime.active is not None


def wait_done():
    while runtime.active is not None:
        yield


def history(action):
    start=time.monotonic()
    getattr(bpy.ops.ed,action)()
    history_times.append((action,time.monotonic()-start))
    yield


def scenario():
    port=launch_server(96)
    settings()
    clear()
    begin(port,progressive=False)
    yield
    assert count()==96
    assert '100%' in adapter.status and '96 cameras created' in adapter.status
    progress_examples['complete']=adapter.status
    yield from history('undo')
    assert count()==0
    yield from history('redo')
    assert count()==96
    checks.append('blocking operator: one Undo/Redo step')
    clear()
    begin(port)
    yield from wait_done()
    assert count()==96
    yield from history('undo')
    assert count()==0
    yield from history('redo')
    assert count()==96 and runtime.active is None
    checks.append('progressive completion: grouped Undo/Redo; no download restart')

    clear()
    terrain_mesh=bpy.data.meshes.new('UI terrain surface')
    terrain_mesh.from_pydata([(-100,-100,0),(100,-100,0),(100,100,0),(-100,100,0)],[],[(0,1,2,3)])
    terrain=bpy.data.objects.new('UI terrain surface',terrain_mesh)
    bpy.context.scene.collection.objects.link(terrain)
    terrain.location.z=7.
    terrain.hide_viewport=True
    terrain.hide_render=True
    terrain.hide_set(True)
    settings().blosm.terrainObject=terrain.name
    bpy.context.view_layer.update()
    bpy.ops.ed.undo_push(message='Terrain fixture setup')
    for progressive in (False,True):
        begin(port,progressive=progressive)
        if progressive: yield from wait_done()
        else: yield
        assert count()==96 and all(abs(o.location.z-9.)<1e-5
            for o in bpy.context.scene.objects if 'mapillary_image_id' in o)
        yield from history('undo')
        assert count()==0 and bpy.context.scene.objects.get('UI terrain surface') is not None
        yield from history('redo')
        assert count()==96 and all(abs(o.location.z-9.)<1e-5
            for o in bpy.context.scene.objects if 'mapillary_image_id' in o)
        terrain=bpy.context.scene.objects['UI terrain surface']
        assert terrain.location.z==7. and all(v.co.z==0 for v in terrain.data.vertices)
        assert terrain.hide_viewport and terrain.hide_render and terrain.hide_get()
        assert not any(c.name.startswith('.Mapillary terrain evaluation') for c in bpy.data.collections)
        clear()
    terrain_mesh=terrain.data
    bpy.data.objects.remove(terrain,do_unlink=True)
    bpy.data.meshes.remove(terrain_mesh)
    settings().blosm.terrainObject=''
    checks.append('terrain-relative blocking/progressive placement and Undo/Redo preserve terrain')

    clear()
    gate.unlink(missing_ok=True)
    begin(port)
    while not count(): yield
    before_edit=count()
    bpy.ops.object.empty_add('EXEC_DEFAULT',True,location=(10,0,0))
    bpy.context.object.name='User edit kept'
    gate.touch()
    yield from wait_done()
    assert count()==96
    yield from history('undo')
    assert count()==before_edit and any(o.type=='EMPTY' for o in bpy.context.scene.objects), (count(),before_edit)
    yield from history('undo')
    assert count()==before_edit and not any(o.type=='EMPTY' for o in bpy.context.scene.objects), (count(),before_edit)
    yield from history('undo')
    assert count()==0
    checks.append('intervening user edit splits import undo groups')

    clear()
    gate.unlink(missing_ok=True)
    begin(port)
    while not count(): yield
    kept=count()
    before=len(messages)
    bpy.ops.blosm.cancel_progressive()
    yield
    assert runtime.active is None and count()==kept and len(messages)==before+1
    yield from history('undo')
    assert count()==0
    yield from history('redo')
    assert count()==kept
    checks.append('cancel retains cameras; partial import undo/redo and one result')

    clear()
    gate.touch()
    begin(port)
    while not count(): yield
    bpy.data.collections.remove(adapter.collection)
    yield from wait_done()
    assert adapter.job.error is not None and messages[-1][0]=='ERROR'
    checks.append('deleted destination stops real timer import')

    clear()
    gate.unlink(missing_ok=True)
    begin(port)
    while not count(): yield
    original=bpy.context.scene
    other=bpy.data.scenes.new('Isolated pause test')
    bpy.context.window.scene=other
    for _ in range(8): yield
    assert runtime.active.paused
    bpy.context.window.scene=original
    bpy.data.scenes.remove(other)
    cube=next(o for o in original.objects if o.type=='MESH')
    bpy.context.view_layer.objects.active=cube
    cube.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    gate.touch()
    paused_count=count()
    for _ in range(8): yield
    assert runtime.active.paused and count()==paused_count
    bpy.ops.object.mode_set(mode='OBJECT')
    yield from wait_done()
    checks.append('scene switch and Edit Mode pause/resume without scene commits')

    clear()
    bad_port=launch_server(96,fail_metadata=True)
    begin(bad_port)
    yield from wait_done()
    assert count()==96 and adapter.errors['metadata_request']>0
    assert all(o['mapillary_pose_approximate'] for o in bpy.context.scene.objects if 'mapillary_image_id' in o)
    checks.append('real HTTP failure: explicit approximate fallback and aggregate errors')

    clear()
    scale_port=launch_server(args.count)
    before_rss=memory_bytes()
    scale_started=time.monotonic()
    before_beats=len(heartbeat)
    initial_gaps,initial_steps,initial_checkpoints=len(import_gaps),len(step_times),len(checkpoint_times)
    begin(scale_port)
    scale_task=adapter
    yield from wait_done()
    scale_elapsed=time.monotonic()-scale_started
    assert count()==args.count and scale_task.metrics['peak_metadata_records']<=32
    scale_rss=memory_bytes()
    yield from history('undo')
    assert count()==0
    after_undo_rss=memory_bytes()
    yield from history('redo')
    assert count()==args.count
    checks.append('scale: real HTTP, live timers, Global Undo enabled, grouped Undo/Redo')
    scale={'count':args.count,'elapsed':scale_elapsed,'metrics':scale_task.metrics,
           'rss_before':before_rss,'rss_with_cameras_and_undo':scale_rss,'rss_after_undo':after_undo_rss,
           'heartbeat_count':len(heartbeat)-before_beats,'global_undo':bpy.context.preferences.edit.use_global_undo,
           'max_import_heartbeat_gap':max(import_gaps[initial_gaps:],default=0),
           'max_job_step':max(step_times[initial_steps:],default=0),
           'max_undo_checkpoint':max(checkpoint_times[initial_checkpoints:],default=0),
           'checkpoint_seconds':sum(checkpoint_times[initial_checkpoints:]),
           'history_operations':history_times[-2:]}

    clear()
    gate.unlink(missing_ok=True)
    begin(port)
    while not count(): yield
    addon.unregister()
    assert runtime.active is None
    checks.append('addon disable stops helper and timers')
    return scale


steps=scenario()
last=time.monotonic()


def finish(status, error=None, scale=None):
    if runtime.active is not None: runtime.active.close(silent=True)
    for process in servers:
        process.kill()
        process.wait(timeout=5)
    for process in processes:
        process.wait(timeout=5)
    report={'status':status,'checks':checks,'error':error,'scale':scale,
            'progress_examples':progress_examples,
            'blender':bpy.app.version_string,'elapsed':time.monotonic()-started,
            'max_heartbeat_gap':max(heartbeat,default=0),'messages':messages}
    (args.artifacts/'ui.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    global last,previous_active
    now=time.monotonic()
    heartbeat.append(now-last)
    if previous_active and runtime.active is not None:
        import_gaps.append(now-last)
    previous_active=runtime.active is not None
    last=now
    try:
        assert now-started<550,'UI check deadline'
        if adapter is not None and not adapter.job.done:
            assert '100%' not in adapter.status and 'cameras created' in adapter.status
            if adapter.count and adapter.progress.percent() is not None:
                progress_examples.setdefault('during_import',adapter.status)
        for area in bpy.context.screen.areas:
            if area.type=='VIEW_3D':
                area.spaces.active.region_3d.view_distance+=.001
                area.tag_redraw()
        next(steps)
        return .03
    except StopIteration as result:
        finish('passed',scale=result.value)
    except Exception:
        finish('failed',traceback.format_exc())


bpy.app.timers.register(tick,first_interval=.5)
