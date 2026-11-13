"""Credentialed live Mapillary test; secrets are read from a local file only.

Run via run_mapillary_checks.py --live [--ui]. No token is written into reports,
preferences, command lines, fixtures or scenes. HTTP uses production transport.
"""
import argparse
from collections import Counter
import importlib
import json
import math
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
addon=importlib.import_module(ROOT.name)
addon.register()
module=importlib.import_module(ROOT.name+'.mapillary.progressive')
runtime=importlib.import_module(ROOT.name+'.progressive.blender')
from importlib import import_module
Projection=import_module(ROOT.name+'.util.transverse_mercator').TransverseMercator
sanitize=import_module(ROOT.name+'.downloads.diagnostics').redact_text

parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
parser.add_argument('--token-file',type=Path,required=True)
parser.add_argument('--bounds',type=float,nargs=4,default=(24.9398,60.1698,24.9402,60.1702))
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
token=args.token_file.read_text(encoding='utf-8-sig').strip()
if not token or any(c.isspace() for c in token):
    raise RuntimeError('The local token file must contain only the client access token')
redaction_headers={'Authorization':'OAuth '+token}
scene=bpy.context.scene
settings=scene.blosm
settings.minLon,settings.minLat,settings.maxLon,settings.maxLat=args.bounds
settings.dataType='mapillary'
settings.downloadWorkers=4
settings.concurrentDownloads=True
origin=((args.bounds[1]+args.bounds[3])/2,(args.bounds[0]+args.bounds[2])/2)
projection=Projection(lat=origin[0],lon=origin[1])
messages=[]
runtime.show_result=lambda *a:messages.append(a)
report={'provider':'Mapillary live API','bounds':args.bounds,'blender':bpy.app.version_string,
        'token_source':'local file; not retained','modes':{},'checks':[]}
state={'task':None,'started':None,'blocking':None,'turns':[],'heartbeats':[],'last':None}


def create():
    return module.MapillaryImport(bpy.context,args.bounds,projection,origin,True,token,max_images=2000)


def snapshot(task):
    bpy.context.view_layer.update()
    objects=list(task.collection.objects) if task.collection else []
    result={}
    for obj in objects:
        identity=obj['mapillary_image_id']
        assert identity.isdecimal()
        assert obj.type=='CAMERA' and all(math.isfinite(v) for v in (*obj.location,*obj.rotation_quaternion))
        result[identity]={'location':list(obj.location),'quaternion':list(obj.rotation_quaternion),
                          'model':obj.data.type,'lens':obj.data.lens}
    return result


def summarize(task,elapsed):
    objects=list(task.collection.objects) if task.collection else []
    fields=Counter()
    heading_errors=[]
    heights=[]
    negative_provider_altitudes=0
    for obj in objects:
        # Check every persisted property against the exact secret, without ever
        # emitting it or writing request URLs/raw HTTP bodies into a report.
        for block in (obj,obj.data,task.collection,scene):
            assert token not in json.dumps(dict(block.items()),default=str),'Credential found in scene metadata'
        fields.update([obj['mapillary_orientation_source'],obj['mapillary_altitude_source'],obj['mapillary_camera_model']])
        metadata=json.loads(obj['mapillary_metadata'])
        heights.append(obj.location.z)
        negative_provider_altitudes+=int(obj.get('mapillary_source_altitude',0.)<0)
        expected=metadata.get('computed_compass_angle')
        if isinstance(expected,(int,float)) and obj['mapillary_orientation_source']=='computed_rotation':
            forward=obj.rotation_quaternion.to_matrix()@__import__('mathutils').Vector((0,0,-1))
            lon,lat=obj['mapillary_source_longitude'],obj['mapillary_source_latitude']
            grid=module.pose.projected_basis(projection,lon,lat)
            east=grid[0][0]*forward.x+grid[1][0]*forward.y
            north=grid[0][1]*forward.x+grid[1][1]*forward.y
            if math.hypot(east,north)>.05:
                angle=math.degrees(math.atan2(east,north))
                heading_errors.append(abs((angle-expected+180)%360-180))
    heading_errors.sort()
    return {'cameras':task.count,'approximate':task.approximate,'errors':dict(task.errors),
            'progress':{'counted_tiles':task.progress.counted_tiles,'total_tiles':task.progress.tile_count,
                        'candidates':task.progress.candidates,'processed':task.progress.processed,
                        'estimated_total':task.progress.estimated_total,'status':task.status},
            'fatal_error':sanitize(task.job.error,headers=redaction_headers) if task.job.error else None,
            'elapsed':elapsed,'metrics':task.metrics,'sources_and_models':dict(fields),
            'with_metadata':sum(bool(o['mapillary_metadata_available']) for o in objects),
            'vertical_placement':{'minimum_z':min(heights,default=None),'maximum_z':max(heights,default=None),
                                  'negative_provider_altitudes':negative_provider_altitudes,
                                  'estimated_height_above_ground':2.},
            'heading_comparison':{'count':len(heading_errors),'median_degrees':heading_errors[len(heading_errors)//2] if heading_errors else None,
                                  'maximum_degrees':max(heading_errors,default=None)}}


def finish(status,error=None):
    task=state['task']
    if runtime.active is not None:
        runtime.active.close(silent=True)
    elif task and not task.job.done:
        task.cancel()
    report.update(status=status,error=sanitize(error,headers=redaction_headers) if error else None)
    text=json.dumps(report,indent=2)
    assert token not in text,'Credential found in report'
    (args.artifacts/'live.json').write_text(text,encoding='utf-8')
    print(text,flush=True)
    if not bpy.app.background:
        bpy.ops.wm.quit_blender()
    elif status!='passed':
        raise RuntimeError('Live Mapillary check failed; see sanitized live.json')


def blocking():
    task=create()
    state['task']=task
    selected=tuple(bpy.context.selected_objects)
    active=bpy.context.view_layer.objects.active
    camera=scene.camera
    started=time.monotonic()
    task.run_blocking()
    report['modes']['blocking']=summarize(task,time.monotonic()-started)
    assert task.job.error is None,'Live blocking import stopped'
    assert task.count>0,'No cameras found in selected live area'
    assert task.progress.processed==task.progress.candidates,'Candidate processing count mismatch'
    assert '100%' in task.status and 'cameras created' in task.status,'Incorrect completion progress'
    assert report['modes']['blocking']['with_metadata']>0,'No live camera metadata was retrieved'
    assert all(abs(o.location.z-2.)<1e-5 for o in task.collection.objects),'Flat scene cameras must be above Z=0'
    assert tuple(bpy.context.selected_objects)==selected and bpy.context.view_layer.objects.active==active
    assert scene.camera==camera
    report['checks'].append('Real vector coverage, exact IDs and Graph metadata imported as cameras')
    state['blocking']=snapshot(task)
    for obj in list(task.collection.objects):
        data=obj.data
        bpy.data.objects.remove(obj,do_unlink=True)
        bpy.data.cameras.remove(data)
    bpy.data.collections.remove(task.collection)


def progressive_complete():
    task=state['task']
    report['modes']['progressive']=summarize(task,time.monotonic()-state['started'])
    assert task.job.error is None,'Live progressive import stopped'
    after=snapshot(task)
    before=state['blocking']
    assert set(before)==set(after),'Live dataset IDs changed between blocking and progressive checks'
    for identity in before:
        for key in ('location','quaternion'):
            assert max(abs(a-b) for a,b in zip(before[identity][key],after[identity][key]))<1e-4, 'Live camera pose mismatch'
    report['checks'].append('Blocking/progressive camera ID and pose parity')
    report['checks'].append('No credential present in persisted camera/scene properties')
    report['max_job_step']=max(state['turns'],default=0)
    report['max_import_heartbeat_gap']=max(state['heartbeats'],default=0)
    report['heartbeat_count']=len(state['heartbeats'])
    if not bpy.app.background:
        report['checks'].append('Real timer controller and viewport heartbeats during live import')
    # Save only scene data and redacted metadata, never the in-memory token.
    bpy.ops.wm.save_as_mainfile(filepath=str(args.artifacts/'live_cameras.blend'),check_existing=False)
    finish('passed')


def start_progressive():
    task=create()
    state.update(task=task,started=time.monotonic(),last=time.monotonic())
    original_step=task.job.step
    def step(*a,**kw):
        started=time.monotonic()
        try: return original_step(*a,**kw)
        finally: state['turns'].append(time.monotonic()-started)
    task.job.step=step
    runtime.start(task,bpy.context)


def tick():
    try:
        if state['started'] is None:
            blocking()
            start_progressive()
        else:
            now=time.monotonic()
            state['heartbeats'].append(now-state['last'])
            state['last']=now
            assert now-state['started']<240,'Live check exceeded test deadline'
            for area in bpy.context.screen.areas:
                if area.type=='VIEW_3D':
                    area.spaces.active.region_3d.view_distance+=.001
                    area.tag_redraw()
            if runtime.active is None:
                progressive_complete()
                return None
        return .03
    except Exception as error:
        finish('failed',error)
        return None


if bpy.app.background:
    try:
        blocking()
        # Stepped path is useful for a quick credential probe; --ui is required
        # to validate the real controller, viewport responsiveness and undo.
        task=create()
        state.update(task=task,started=time.monotonic())
        while not task.job.done:
            started=time.monotonic()
            task.job.step()
            state['turns'].append(time.monotonic()-started)
            if task.job.waiting: time.sleep(.001)
        task.finish()
        progressive_complete()
    except Exception as error:
        finish('failed',error)
else:
    bpy.app.timers.register(tick,first_interval=.5)
