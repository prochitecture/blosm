"""Factory-profile background Blender integration and scale checks."""
import argparse
import importlib
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT.parent), str(ROOT / 'tests')]
import bpy
from mapillary_test_support import MemoryTransport, BOUNDS, BASE_ID
addon = importlib.import_module(ROOT.name)
addon.register()
adapter = importlib.import_module(ROOT.name + '.mapillary.progressive')
from mapillary.pose import number
Projection = importlib.import_module(ROOT.name + '.util.transverse_mercator').TransverseMercator

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--count', type=int, default=5000)
args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)
scene = bpy.context.scene
settings = scene.blosm
settings.minLon,settings.minLat,settings.maxLon,settings.maxLat = BOUNDS
scene['height_offset'] = 20.
origin=(60.17,24.94)
checks=[]


def check(label, condition):
    assert condition, label
    checks.append(label)


def task(count=96, **kwargs):
    transport=MemoryTransport(count, **kwargs)
    return adapter.MapillaryImport(bpy.context,BOUNDS,Projection(lat=origin[0],lon=origin[1]),
                                   origin,True,'fixture-client-token',transport=transport)


def drive(t):
    turns=[]
    while not t.job.done:
        started=time.monotonic()
        t.job.step()
        turns.append(time.monotonic()-started)
    assert t.job.error is None, t.job.error
    return turns


def cameras(t):
    return sorted((o['mapillary_image_id'],tuple(round(x,6) for x in o.location),
                   tuple(round(x,6) for x in o.rotation_quaternion),o.data.type,round(o.data.lens,6))
                  for o in t.collection.objects)


selected=tuple(bpy.context.selected_objects)
active=bpy.context.view_layer.objects.active
layer=bpy.context.view_layer.active_layer_collection
cursor=tuple(scene.cursor.location)
render_camera=scene.camera
messages=[]
adapter.notify_result=lambda *args: messages.append(args)
t=task()
turns=drive(t)
baseline=cameras(t)
check('96 cameras, duplicate removed, exact IDs',t.count==96 and str(BASE_ID)==baseline[0][0])
check('selection preserved',tuple(bpy.context.selected_objects)==selected)
check('active object / collection / cursor / render camera preserved',
      bpy.context.view_layer.objects.active==active and bpy.context.view_layer.active_layer_collection==layer
      and tuple(scene.cursor.location)==cursor and scene.camera==render_camera)
check('metadata provenance and geographic origin',scene['lat']==origin[0] and all('mapillary_orientation_source' in o for o in t.collection.objects))
bpy.context.view_layer.update()
geographic_alignment=True
for obj in t.collection.objects:
    expected=t.projection.fromGeographic(obj['mapillary_source_latitude'],obj['mapillary_source_longitude'])
    geographic_alignment &= abs(obj.matrix_world.translation.x-expected[0])<1e-4 and abs(obj.matrix_world.translation.y-expected[1])<1e-4
check('camera world positions agree with standard Blosm geographic projection',geographic_alignment)
check('panoramas are actual camera objects',any(o.data.type=='PANO' for o in t.collection.objects) and all(o.type=='CAMERA' for o in t.collection.objects))
check('no credentials persisted',all('fixture-client-token' not in str(dict(o.items())) for o in t.collection.objects))
t.finish()
check('one aggregate message',len(messages)==1)
t2=task(reverse=False)
t2.run_blocking()
check('blocking/progressive and out-of-order parity',cameras(t2)==baseline)

t3=task()
while t3.count<4: t3.job.step()
count=t3.count
obj=next(iter(t3.collection.objects))
obj.name='user renamed camera'
bpy.data.objects.remove(obj,do_unlink=True)
while not t3.job.done and t3.count<count+16: t3.job.step()
check('deleting or renaming cameras does not stop loading',t3.job.error is None and t3.count>count)
kept=t3.count
t3.cancel()
check('cancellation retains committed cameras',t3.count==kept and len(t3.collection.objects)==kept-1 and t3.transport.closed)

for operation in ('unlink','delete','rename'):
    t4=task()
    while not t4.count: t4.job.step()
    if operation=='unlink': scene.collection.children.unlink(t4.collection)
    elif operation=='delete': bpy.data.collections.remove(t4.collection)
    else: t4.collection.name='renamed destination'
    if operation=='rename':
        drive(t4)
        check('renamed destination remains valid',t4.count==96)
    else:
        t4.job.step()
        check(operation+' destination stops with error',t4.job.done and t4.job.error is not None)

t5=task(missing=True)
drive(t5)
check('missing metadata falls back explicitly',t5.count==96 and t5.errors['missing_metadata']>0 and
      any(not o['mapillary_metadata_available'] and o['mapillary_pose_approximate'] for o in t5.collection.objects))

wide_bounds=(24.92,60.16,24.96,60.18)
cross=adapter.MapillaryImport(bpy.context,wide_bounds,Projection(lat=origin[0],lon=origin[1]),
                              origin,False,'fixture-client-token',transport=MemoryTransport(8))
drive(cross)
check('deduplicate buffered image IDs across multiple tiles',cross.count==8 and cross.transport.requests>2)

limited=adapter.MapillaryImport(bpy.context,BOUNDS,Projection(lat=origin[0],lon=origin[1]),
                                origin,False,'fixture-client-token',transport=MemoryTransport(96),max_images=10)
while not limited.job.done: limited.job.step()
check('explicit bounded deduplication failure',limited.job.error is not None and 'image limit' in str(limited.job.error))

large=task(args.count)
started=time.monotonic()
scale_turns=drive(large)
scale=dict(count=large.count,elapsed=time.monotonic()-started,max_turn=max(scale_turns),
           metrics=large.metrics,source='synthetic; MemoryTransport; background stepped Job; no GUI undo')
check('bounded decoded metadata',large.metrics['peak_metadata_records']<=32)
check('scale count and processing',large.count==args.count)
report={'checks':checks,'scale':scale,'small_max_turn':max(turns),'blender':bpy.app.version_string}
(args.artifacts/'background.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
addon.unregister()
