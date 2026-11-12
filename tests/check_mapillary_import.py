"""Factory-profile background Blender integration and scale checks."""
import argparse
import importlib
import json
import math
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
check('without terrain all cameras are 2 m above Z=0 despite provider altitude and scene offset',
      all(abs(o.location.z-2.)<1e-6 and o['mapillary_vertical_position_source']=='scene_plane'
          and o['mapillary_camera_height_source']=='assumed' and o['mapillary_pose_approximate'] for o in t.collection.objects))
check('negative provider altitude is preserved without lowering the camera',
      any(o.get('mapillary_source_altitude')==-3. and o.location.z==2. for o in t.collection.objects))
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

# Ground projection uses evaluated mesh geometry in world coordinates, while
# preserving the terrain's mesh and transforms. A translated, tilted, scaled
# surface catches incorrect local-space rays and double height-offset handling.
mesh=bpy.data.meshes.new('Terrain test surface')
mesh.from_pydata([(-100,-100,0),(100,-100,0),(100,100,0),(-100,100,0)],[],[(0,1,2,3)])
terrain=bpy.data.objects.new('Mapillary test terrain',mesh)
scene.collection.objects.link(terrain)
terrain.location=(5.,-3.,8.)
terrain.rotation_euler.y=.2
terrain.scale=(2.,3.,1.5)
terrain['height_offset']=999.
settings.terrainObject=terrain.name
bpy.context.view_layer.update()
vertices=[tuple(v.co) for v in mesh.vertices]
matrix=terrain.matrix_world.copy()
grounded=task(8)
drive(grounded)
check('cameras follow translated tilted scaled terrain in world coordinates',
      all(abs(o.location.z-(8.-math.tan(.2)*(o.location.x-5.)+2.))<1e-4
          and o['mapillary_vertical_position_source']=='terrain_surface' for o in grounded.collection.objects))
check('terrain geometry and transform preserved',vertices==[tuple(v.co) for v in mesh.vertices]
      and terrain.matrix_world==matrix)
terrain.rotation_euler.y=0.
terrain.location.z=10.
terrain.scale.z=2.
modifier=terrain.modifiers.new('Evaluated top surface','SOLIDIFY')
modifier.thickness=3.
modifier.offset=1.
bpy.context.view_layer.update()
evaluated=task(8)
drive(evaluated)
check('terrain modifiers evaluated and upper surface sampled',
      all(abs(o.location.z-18.)<1e-4 for o in evaluated.collection.objects))
terrain.modifiers.remove(modifier)
terrain.location.z=-10.
bpy.context.view_layer.update()
below_plane=task(8)
drive(below_plane)
check('negative terrain elevation remains terrain-relative',
      all(abs(o.location.z+8.)<1e-4 for o in below_plane.collection.objects))
terrain.location.x=1000.
bpy.context.view_layer.update()
outside=task(8)
drive(outside)
check('outside terrain uses marked plane fallback',outside.terrain_misses==8 and
      all(o.location.z==2. and o['mapillary_vertical_position_source']=='scene_plane_outside_terrain'
          for o in outside.collection.objects))
terrain.location.x=5.
bpy.context.view_layer.update()
renamed=task(96)
while not renamed.count: renamed.job.step()
terrain.name='Renamed terrain while importing'
drive(renamed)
check('renaming terrain preserves sampling identity',renamed.count==96 and renamed.job.error is None)
settings.terrainObject=terrain.name
for operation in ('unlink','delete'):
    interrupted=task(96)
    while not interrupted.count: interrupted.job.step()
    kept=interrupted.count
    if operation=='unlink': scene.collection.objects.unlink(terrain)
    else: bpy.data.objects.remove(terrain,do_unlink=True)
    interrupted.job.step()
    check(operation+' terrain stops import and retains cameras',interrupted.job.done
          and interrupted.job.error is not None and interrupted.count==kept and interrupted.transport.closed)
    if operation=='unlink': scene.collection.objects.link(terrain)
settings.terrainObject=''

large=task(args.count)
started=time.monotonic()
scale_turns=drive(large)
scale=dict(count=large.count,elapsed=time.monotonic()-started,max_turn=max(scale_turns),
           metrics=large.metrics,source='synthetic; MemoryTransport; background stepped Job; no GUI undo')
check('bounded decoded metadata',large.metrics['peak_metadata_records']<=32)
check('scale count and processing',large.count==args.count)

# Exercise thousands of native cached rays on a moderately detailed mesh. This
# is separate from the no-terrain baseline and excludes GUI undo overhead.
grid=100
verts=[(-100.+2*x,-100.+2*y,35.) for y in range(grid+1) for x in range(grid+1)]
faces=[(y*(grid+1)+x,y*(grid+1)+x+1,(y+1)*(grid+1)+x+1,(y+1)*(grid+1)+x)
       for y in range(grid) for x in range(grid)]
dense_mesh=bpy.data.meshes.new('Detailed terrain fixture')
dense_mesh.from_pydata(verts,[],faces)
dense_terrain=bpy.data.objects.new('Detailed terrain fixture',dense_mesh)
scene.collection.objects.link(dense_terrain)
settings.terrainObject=dense_terrain.name
bpy.context.view_layer.update()
ground_large=task(args.count)
started=time.monotonic()
ground_turns=drive(ground_large)
ground_scale=dict(count=ground_large.count,elapsed=time.monotonic()-started,max_turn=max(ground_turns),
                  metrics=ground_large.metrics,terrain_quads=len(faces),source='synthetic; no GUI undo')
check('scale terrain rays place every camera above the mesh',ground_large.count==args.count
      and ground_large.terrain_misses==0 and all(abs(o.location.z-37.)<1e-5 for o in ground_large.collection.objects))
check('ground sampling does not create extra terrain geometry',len(dense_mesh.vertices)==len(verts)
      and len(dense_mesh.polygons)==len(faces) and len(dense_terrain.modifiers)==0)
report={'checks':checks,'scale':scale,'ground_scale':ground_scale,
        'small_max_turn':max(turns),'blender':bpy.app.version_string}
(args.artifacts/'background.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
addon.unregister()
