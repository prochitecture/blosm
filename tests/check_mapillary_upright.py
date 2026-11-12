"""Inspect saved real cameras against Mapillary's displayed up-vector convention."""
import argparse
from collections import Counter
import importlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
from mathutils import Vector
pose=importlib.import_module(ROOT.name+'.mapillary.pose')
parser=argparse.ArgumentParser()
source=parser.add_mutually_exclusive_group(required=True)
source.add_argument('--scene',type=Path)
source.add_argument('--pose-report',type=Path,
                    help='Replay bounded samples from check_mapillary_area_pose.py without HTTP')
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
if args.scene:
    bpy.ops.wm.open_mainfile(filepath=str(args.scene))
else:
    # Reuse the actual production commit path for the bounded real metadata
    # samples. A dummy transport guarantees this replay makes no HTTP requests
    # and needs no development credential. The original large scene is untouched.
    addon=importlib.import_module(ROOT.name)
    addon.register()
    adapter=importlib.import_module(ROOT.name+'.mapillary.progressive')
    Projection=importlib.import_module(ROOT.name+'.util.transverse_mercator').TransverseMercator
    document=json.loads(args.pose_report.read_text(encoding='utf-8'))
    bounds=document['bounds']
    origin=((bounds[1]+bounds[3])/2,(bounds[0]+bounds[2])/2)
    projection=Projection(lat=origin[0],lon=origin[1])
    transport=SimpleNamespace(max_requests=1,close=lambda:None)
    task=adapter.MapillaryImport(bpy.context,bounds,projection,origin,True,
                                'offline-fixture-token',transport=transport)
    for group in document['samples'].values():
        for sample in group:
            metadata,record=sample['metadata'],sample['record']
            camera_pose=pose.interpret(metadata,record,projection)
            task.commit(record,metadata,camera_pose)
    task.cancel()
    bpy.context.view_layer.update()
    bpy.ops.wm.save_as_mainfile(filepath=str(args.artifacts/'sample_cameras.blend'),check_existing=False)
rows=[]
counts=Counter()
samples={}
max_up_z_error=0.
max_forward_z_error=0.
skipped=0
for obj in bpy.context.scene.objects:
    if 'mapillary_image_id' not in obj: continue
    metadata=json.loads(obj['mapillary_metadata'])
    try:
        raw=pose.rotation_vector(metadata.get('computed_rotation'))
    except ValueError:
        skipped+=1
        continue
    exif=metadata.get('exif_orientation',1)
    matrix=obj.matrix_world.to_3x3().normalized()
    up=matrix@Vector((0,1,0))
    forward=matrix@Vector((0,0,-1))
    # Transform.upVector(): column-major matrix elements correspond to rows
    # of the raw world-to-camera rotation.
    raw_axis={1:(0,-1,0),3:(0,1,0),6:(-1,0,0),8:(1,0,0)}.get(exif,(0,-1,0))
    viewer_up=Vector([sum(raw[j][i]*raw_axis[j] for j in range(3)) for i in range(3)])
    # Geographic grid convergence rotates around world Z, so these vertical
    # components are directly comparable without using the adapter's grid math.
    # Compare actual saved Blender transforms, including their quaternion
    # conversion, with Mapillary's independently specified viewer directions.
    max_up_z_error=max(max_up_z_error,abs(up.z-viewer_up.z))
    max_forward_z_error=max(max_forward_z_error,abs(forward.z-raw[2][2]))
    counts[(str(exif),'up' if up.z>=0 else 'down')]+=1
    row={'id':obj['mapillary_image_id'],'exif':exif,
                 'model':metadata.get('camera_type'),'blender_up':list(up),
                 'viewer_up':list(viewer_up),'forward':list(forward),
                 'rotation':metadata['computed_rotation'],
                 'scene_z':obj.matrix_world.translation.z,
                 'source_altitude':obj.get('mapillary_source_altitude')}
    rows.append(row)
    samples.setdefault(str(exif),row)
report={'count':len(rows),'skipped_without_computed_rotation':skipped,
        'orientation_counts':{'/'.join(k):v for k,v in counts.items()},
        'minimum_up_z':min((row['blender_up'][2] for row in rows),default=None),
        'maximum_up_z_error_vs_mapillary':max_up_z_error,
        'maximum_forward_z_error_vs_mapillary':max_forward_z_error,
        'minimum_scene_z':min((row['scene_z'] for row in rows),default=None),
        'maximum_scene_z':max((row['scene_z'] for row in rows),default=None),
        'sample_by_exif':samples,
        'downward_ids':[row['id'] for row in rows if row['blender_up'][2]<0],
        'downward_samples':[row for row in rows if row['blender_up'][2]<0]}
(args.artifacts/'upright.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
assert rows,'No cameras with computed rotation found'
assert max_up_z_error<1e-5,'Saved camera up direction differs from Mapillary'
assert max_forward_z_error<1e-5,'Saved camera viewing direction differs from Mapillary'
if args.pose_report:
    assert all(abs(row['scene_z']-2.)<1e-5 for row in rows),'Offline sample cameras must be above the scene plane'
