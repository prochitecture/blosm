"""Opt-in real-photo probe. Read the local token once; persist no request URLs."""
import argparse
import importlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy

parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
parser.add_argument('--token-file',type=Path,required=True)
parser.add_argument('--scene',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
addon=importlib.import_module(ROOT.name)
addon.register()
photo=importlib.import_module(ROOT.name+'.mapillary.photo')
token=args.token_file.read_text(encoding='utf-8-sig').strip()
photo.runtime.notify_result=lambda *_args:None
report={'status':'failed','photos':[],'blender':bpy.app.version_string}
try:
    bpy.ops.wm.open_mainfile(filepath=str(args.scene))
    bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
    for model in ('perspective','fisheye','spherical'):
        camera=next(obj for obj in bpy.context.scene.objects if obj.type=='CAMERA'
                    and obj.get('mapillary_camera_model')==model)
        bpy.context.view_layer.objects.active=camera
        camera.select_set(True)
        task=photo.PhotoImport(bpy.context,camera,token)
        assert task.run_blocking(),'Live photo download failed'
        image=photo.cached_image(camera['mapillary_image_id'])
        assert image is not None and image.packed_file is not None
        assert bpy.context.scene.camera==camera
        assert any(bg.image==image and bg.alpha==1. and bg.frame_method=='FIT'
                   for bg in camera.data.background_images)
        assert photo.cached_image(camera['mapillary_image_id'])==image
        report['photos'].append({'image_id':camera['mapillary_image_id'],'camera_model':model,
                                  'dimensions':list(image.size),'packed':True,
                                  'creator_available':bool(image.get('mapillary_creator'))})
    for block in (*bpy.data.images,*bpy.data.objects,*bpy.data.scenes):
        assert token not in json.dumps(dict(block.items()),default=str)
    bpy.ops.wm.save_as_mainfile(filepath=str(args.artifacts/'live_photos.blend'),check_existing=False)
    report['status']='passed'
except Exception as error:
    report['error']=photo.redact_text(error,headers={'Authorization':'OAuth '+token})
finally:
    text=json.dumps(report,indent=2)
    assert token not in text
    (args.artifacts/'live.json').write_text(text,encoding='utf-8')
    print(text,flush=True)
    addon.unregister()
if report['status']!='passed':
    raise RuntimeError('Live photo check failed; see sanitized report')
