"""Read-only user-scene reproduction; synthetic cameras, no HTTP or credentials."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT.parent),str(ROOT/'tests')]
import bpy
from mathutils import Vector
from mapillary_test_support import MemoryTransport,BASE_ID,mvt
addon=importlib.import_module(ROOT.name)
addon.register()
adapter=importlib.import_module(ROOT.name+'.mapillary.progressive')
Projection=importlib.import_module(ROOT.name+'.util.transverse_mercator').TransverseMercator
parser=argparse.ArgumentParser()
parser.add_argument('--scene',type=Path,required=True)
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
digest=hashlib.sha256(args.scene.read_bytes()).digest()
bpy.ops.wm.open_mainfile(filepath=str(args.scene),load_ui=False)
scene=bpy.context.scene
terrain=scene.objects['Terrain']
scene.blosm.terrainObject=terrain.name
original_vertices=[tuple(v.co) for v in terrain.data.vertices]
original_matrix=terrain.matrix_world.copy()
flags=(terrain.hide_viewport,terrain.hide_render,terrain.hide_get())
selected=tuple(bpy.context.selected_objects)
active=bpy.context.view_layer.objects.active
layer=bpy.context.view_layer.active_layer_collection
cursor=tuple(scene.cursor.location)
render_camera=scene.camera
try:
    terrain.evaluated_get(bpy.context.evaluated_depsgraph_get()).ray_cast(Vector((0,0,10000)),Vector((0,0,-1)))
except RuntimeError as error:
    native_error=str(error).strip()
else:
    native_error=None
# This supplied fixture has no terrain modifiers; triangle centroids provide
# independently known positions directly on its original world-space surface.
assert not terrain.modifiers,'This scene regression expects the supplied unmodified terrain'
terrain.data.calc_loop_triangles()
triangles=list(terrain.data.loop_triangles)
points=[]
for i in range(0,len(triangles),max(1,len(triangles)//16)):
    tri=triangles[i]
    points.append(original_matrix @ (sum((terrain.data.vertices[v].co for v in tri.vertices),Vector())/3))
    if len(points)==16:break
# Older OSM setup files may not store scene lat/lon. The regression only needs
# a reversible XY/geographic mapping; use the known fixture-area centre then.
projection=Projection(lat=float(scene.get('lat',52.51816)),lon=float(scene.get('lon',13.434145)))
coords=[projection.toGeographic(p.x,p.y) for p in points]
bounds=(min(lon for lat,lon in coords)-.00001,min(lat for lat,lon in coords)-.00001,
        max(lon for lat,lon in coords)+.00001,max(lat for lat,lon in coords)+.00001)
def task():
    transport=MemoryTransport(16)
    for i,(lat,lon) in enumerate(coords):
        metadata=transport.entries[str(BASE_ID+i)]
        metadata['geometry']={'coordinates':[lon,lat]}
        metadata['computed_geometry']={'coordinates':[lon,lat]}
    transport.tile=next(adapter.tile_coordinates(bounds))
    transport.content=mvt([(BASE_ID+i,lon,lat) for i,(lat,lon) in enumerate(coords)],transport.tile)
    return adapter.MapillaryImport(bpy.context,bounds,projection,(projection.lat,projection.lon),False,
                                   'offline-fixture-token',transport=transport)
adapter.notify_result=lambda *a:None
tasks=[]
for blocking in (True,False):
    current=task()
    if blocking:current.run_blocking()
    else:
        while not current.job.done:current.job.step()
    assert current.job.error is None,str(current.job.error)
    assert current.count==len(points)
    for obj in current.collection.objects:
        point=points[int(obj['mapillary_image_id'])-BASE_ID]
        assert abs(obj.location.z-point.z-2.)<.002,(tuple(obj.location),tuple(point))
        assert obj['mapillary_vertical_position_source']=='terrain_surface'
    tasks.append(current)
assert original_vertices==[tuple(v.co) for v in terrain.data.vertices]
assert terrain.matrix_world==original_matrix
assert flags==(terrain.hide_viewport,terrain.hide_render,terrain.hide_get())
assert selected==tuple(bpy.context.selected_objects) and active==bpy.context.view_layer.objects.active
assert layer==bpy.context.view_layer.active_layer_collection and cursor==tuple(scene.cursor.location)
assert render_camera==scene.camera
assert not any(c.name.startswith('.Mapillary terrain evaluation') for c in bpy.data.collections)
assert digest==hashlib.sha256(args.scene.read_bytes()).digest(),'Source scene file changed'
report={'status':'passed','native_error_reproduced':native_error,'terrain_vertices':len(original_vertices),
        'source_hidden_flags':flags,'cameras_per_mode':[t.count for t in tasks],
        'mesh_transform_visibility_context_preserved':True,'temporary_ids_removed':True,
        'source_file_unchanged':True,'blender':bpy.app.version_string}
(args.artifacts/'scene_regression.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
