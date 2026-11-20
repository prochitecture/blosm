"""Native transient source-camera visibility, history, filters and save checks.

Only run in an isolated UI Blender with factory startup, disabled autoexec,
event simulation and the five chat-specific user-resource paths. No network.
"""
import argparse
import importlib
import json
from pathlib import Path
import shutil
import sys
import time
import traceback
from unittest.mock import patch

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
from mathutils import Quaternion
parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
addon=importlib.import_module(ROOT.name);addon.register()
viewport=importlib.import_module(ROOT.name+'.viewport')
pano=importlib.import_module(ROOT.name+'.mapillary.panorama')
photo=importlib.import_module(ROOT.name+'.mapillary.photo')
visibility=importlib.import_module(ROOT.name+'.viewport.visibility')
host=importlib.import_module(ROOT.name+'.viewport.session')
checks=[]
state={'phase':0,'ticks':0,'start':time.monotonic(),'registered':True}
SOURCE='9007199254740993'
OTHER='9007199254740994'
HIDDEN='9007199254740995'


def check(label,condition):
    assert condition,label
    checks.append(label)


def camera(identity):
    return next(o for o in bpy.context.scene.objects if o.get('mapillary_image_id')==identity)


def create(identity,photo_cache=True):
    obj=bpy.data.objects.new('Visibility camera '+identity,bpy.data.cameras.new('Visibility camera'))
    bpy.context.scene.collection.objects.link(obj)
    obj.location=(0,0,2);obj.rotation_mode='QUATERNION'
    obj.rotation_quaternion=Quaternion((1,0,0),1.5707963267948966)
    for key,value in dict(image_id=identity,camera_model='spherical',metadata='{}',
                          orientation_source='computed_rotation',exif_orientation=1).items():
        obj['mapillary_'+key]=value
    if photo_cache and photo.cached_image(identity) is None:
        image=bpy.data.images.new('Visibility cached '+identity,width=128,height=64)
        image.pixels.foreach_set([.05,.1,.2,1]*(128*64));image.pack();image.use_fake_user=True
        image['mapillary_photo_source']='Mapillary';image['mapillary_image_id']=identity
    bpy.context.view_layer.update()
    return obj


def activate(obj):
    for selected in tuple(bpy.context.selected_objects):selected.select_set(False)
    bpy.context.view_layer.objects.active=obj
    obj.select_set(True)


def open_camera(identity=SOURCE):
    state['session']=viewport.open(bpy.context,pano.Panorama(camera(identity)))
    return state['session']


def empty_records():
    return not visibility._leases and all(not layer.blosm_temporary_hidden
        for scene in bpy.data.scenes for layer in scene.view_layers)


def next_phase(phase):state.update(phase=phase,ticks=0)


def finish(error=None):
    if viewport.current():viewport.current().close('USER')
    if state['registered']:addon.unregister()
    result=dict(status='failed' if error else 'passed',error=error,checks=checks,
                blender=bpy.app.version_string,elapsed=time.monotonic()-state['start'])
    (args.artifacts/'visibility.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        assert time.monotonic()-state['start']<120,'Visibility test timeout'
        state['ticks']+=1
        if state['ticks']<5:return .04
        state['ticks']=0
        window=bpy.context.window
        area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
        region=next(r for r in area.regions if r.type=='WINDOW')
        with bpy.context.temp_override(window=window,area=area,region=region):
            p=state['phase']
            layer=bpy.context.view_layer
            if p==0:
                window.event_simulate(type='ESC',value='PRESS');window.event_simulate(type='ESC',value='RELEASE')
                for identity in (SOURCE,OTHER,HIDDEN):create(identity)
                camera(HIDDEN).hide_set(True)
                bpy.context.scene.view_layers.new('Other view layer')
                activate(camera(SOURCE));state['pose']=camera(SOURCE).matrix_world.copy()
                state['render_camera']=bpy.context.scene.camera
                next_phase(1)
            elif p==1:
                if 'session' not in state:
                    open_camera();return .04
                obj=camera(SOURCE)
                check('source camera alone is hidden in the owner view layer',obj.hide_get() and not camera(OTHER).hide_get())
                check('temporary hiding does not affect another view layer',not obj.hide_get(view_layer=bpy.context.scene.view_layers['Other view layer']))
                check('global viewport/render visibility and camera pose stay unchanged',not obj.hide_viewport
                      and not obj.hide_render and obj.matrix_world==state['pose'] and bpy.context.scene.camera==state['render_camera'])
                obj.name='Renamed hidden source'
                state['session'].close('USER')
                check('exit restores source visibility and its untouched active selection',not obj.hide_get() and obj.select_get()
                      and layer.objects.active==obj and empty_records())
                open_camera(HIDDEN);state['session'].close('USER')
                check('a previously hidden source remains hidden after exit',camera(HIDDEN).hide_get())
                activate(obj);first=open_camera();activate(camera(OTHER));second=open_camera(OTHER)
                check('viewer replacement restores the old source and hides the new source',first.closed
                      and not obj.hide_get() and camera(OTHER).hide_get())
                second.close('USER')
                check('restoring the old source does not steal the new selection',layer.objects.active==camera(OTHER)
                      and camera(OTHER).select_get() and not obj.select_get())
                mesh=bpy.data.objects['Cube'];activate(mesh);open_camera();bpy.ops.object.mode_set(mode='EDIT')
                next_phase(2)
            elif p==2:
                check('Edit Mode retains panorama and temporarily hidden source',viewport.current() is state['session']
                      and camera(SOURCE).hide_get() and not empty_records())
                bpy.ops.object.mode_set(mode='OBJECT')
                open_camera();state['session'].error='Fixture viewport failure';next_phase(3)
            elif p==3:
                check('error exit restores source visibility',viewport.current() is None and not camera(SOURCE).hide_get() and empty_records())
                open_camera();bpy.data.objects.remove(camera(SOURCE),do_unlink=True);next_phase(4)
            elif p==4:
                check('deleting the hidden source leaves no visibility lease or marker',viewport.current() is None and empty_records())
                create(SOURCE);activate(camera(SOURCE))
                try:
                    with patch.object(pano,'shader_create',side_effect=RuntimeError('Fixture shader failure')):
                        open_camera()
                except RuntimeError:pass
                check('failed startup restores visibility and selection',viewport.current() is None and not camera(SOURCE).hide_get()
                      and camera(SOURCE).select_get() and empty_records())
                for i in range(48):create(str(8000000000000000+i),False)
                open_camera()
                settings=layer.blosm_mapillary_filters;settings.panoramic=False
                bpy.ops.blosm.apply_mapillary_filters()
                check('filter ownership records the temporarily hidden source',any(item.object==camera(SOURCE) for item in settings.hidden)
                      and state['session'].surface.visibility.hidden)
                state['session'].close('USER')
                check('exit respects a filter applied during panorama mode',camera(SOURCE).hide_get())
                bpy.ops.blosm.reset_mapillary_filters()
                check('filter reset reveals the source while preserving manual hiding',not camera(SOURCE).hide_get() and camera(HIDDEN).hide_get())
                activate(camera(SOURCE));open_camera()
                settings.panoramic=False;bpy.ops.blosm.apply_mapillary_filters();bpy.ops.blosm.reset_mapillary_filters()
                check('filter reset during viewing keeps the source temporarily hidden',camera(SOURCE).hide_get()
                      and not state['session'].surface.visibility.hidden)
                state['session'].close('USER')
                check('exit after filter reset reveals the source',not camera(SOURCE).hide_get() and empty_records())
                activate(camera(SOURCE));open_camera()
                state['path']=str(args.artifacts/'saved_during_panorama.blend')
                bpy.ops.wm.save_as_mainfile(filepath=state['path'])
                check('saving resumes temporary hiding in the still-open viewer',not state['session'].closed and camera(SOURCE).hide_get())
                visibility.save_pre()
                check('save preparation suspends temporary hiding and its markers',not camera(SOURCE).hide_get()
                      and not layer.blosm_temporary_hidden)
                visibility.save_post()
                check('save completion/failure callback resumes the hide',camera(SOURCE).hide_get()
                      and len(layer.blosm_temporary_hidden)==1)
                inspection=args.artifacts/'saved_for_inspection.blend'
                shutil.copyfile(state['path'],inspection)
                with bpy.data.libraries.load(str(inspection),link=False) as (source,dest):dest.scenes=[bpy.context.scene.name]
                saved=dest.scenes[0];saved_layer=saved.view_layers[layer.name]
                saved_obj=next(o for o in saved.objects if o.get('mapillary_image_id')==SOURCE)
                check('saved file contains visible source and no temporary markers',not saved_obj.hide_get(view_layer=saved_layer)
                      and not saved_layer.blosm_temporary_hidden)
                bpy.data.scenes.remove(saved)
                next_phase(5)
            elif p==5:
                # load_pre must clean the running session before IDs disappear;
                # this test timer is persistent only inside this isolated process.
                bpy.ops.wm.open_mainfile(filepath=state['path']);next_phase(6)
            elif p==6:
                check('opening a file exits the viewer with a visible source',viewport.current() is None
                      and not camera(SOURCE).hide_get() and empty_records())
                activate(camera(SOURCE));bpy.ops.ed.undo_push(message='Before panorama snapshot')
                open_camera();bpy.ops.ed.undo_push(message='Snapshot with hidden panorama source')
                state['session'].close('USER');bpy.ops.ed.undo_push(message='After panorama exit')
                bpy.ops.ed.undo()
                check('undo into an old closed-viewer snapshot cannot hide the camera permanently',viewport.current() is None
                      and not camera(SOURCE).hide_get() and empty_records())
                bpy.ops.ed.redo()
                check('redo also leaves source visibility restored',not camera(SOURCE).hide_get() and empty_records())
                activate(camera(SOURCE));bpy.ops.ed.undo_push(message='Before live viewer rename')
                open_camera();camera(SOURCE).name='Hidden source renamed in history'
                bpy.ops.ed.undo_push(message='Rename during panorama');bpy.ops.ed.undo()
                check('undo while viewing restores source and closes session',viewport.current() is None
                      and not camera(SOURCE).hide_get() and empty_records())
                bpy.ops.ed.redo()
                check('history recovery follows object pointers through camera renames',camera(SOURCE).name=='Hidden source renamed in history'
                      and not camera(SOURCE).hide_get() and empty_records())
                open_camera();addon.unregister();state['registered']=False
                check('addon disable restores the source visibility',not camera(SOURCE).hide_get() and not visibility._leases)
                addon.register();state['registered']=True
                check('re-enabling addon leaves no temporary markers',empty_records())
                finish();return None
        return .04
    except Exception:
        finish(traceback.format_exc());return None


bpy.app.timers.register(tick,first_interval=.5,persistent=True)
