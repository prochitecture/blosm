"""Native click/drag selection and mode-switch checks with a cached panorama.

Run only in an isolated Blender UI process, with factory startup, disabled
autoexec, --enable-event-simulate, and this chat's five user-resource paths.
This script uses no token, network, installed addon or existing Blender session.
"""
import argparse
import importlib
import json
import math
from pathlib import Path
import sys
import time
import traceback

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
from mathutils import Quaternion
from bpy_extras.view3d_utils import location_3d_to_region_2d

parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
addon=importlib.import_module(ROOT.name)
addon.register()
viewport=importlib.import_module(ROOT.name+'.viewport')
pano=importlib.import_module(ROOT.name+'.mapillary.panorama')
photo=importlib.import_module(ROOT.name+'.mapillary.photo')
ui=importlib.import_module(ROOT.name+'.viewport.ui')
runtime=importlib.import_module(ROOT.name+'.progressive.blender')
runtime.show_result=lambda *_:None
window=bpy.context.window
area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
region=next(r for r in area.regions if r.type=='WINDOW')
checks=[]
state={'step':0,'ticks':0,'started':time.monotonic()}


def check(label,condition):
    assert condition,label
    checks.append(label)


mesh=bpy.data.objects['Cube']
mesh.name='Selectable mesh'
mesh.location=(-2.2,11,2)
other=bpy.data.objects.new('Second selectable mesh',mesh.data.copy())
bpy.context.scene.collection.objects.link(other)
other.location=(2.2,11,2)
camera=bpy.data.objects.new('Other Mapillary camera',bpy.data.cameras.new('Other camera'))
bpy.context.scene.collection.objects.link(camera)
camera.location=(0,8,4.2)
camera.rotation_euler=(math.pi/2,0,0)
camera.data.display_size=1.2
for key,value in dict(image_id='9007199254740994',camera_model='perspective',
                      metadata='{}',orientation_source='computed_rotation').items():
    camera['mapillary_'+key]=value
obj=bpy.data.objects.new('Panorama eye',bpy.data.cameras.new('Panorama eye'))
bpy.context.scene.collection.objects.link(obj)
obj.location=(0,0,2)
obj.rotation_mode='QUATERNION'
obj.rotation_quaternion=Quaternion((1,0,0),math.pi/2)
for key,value in dict(image_id='9007199254740993',camera_model='spherical',
                      metadata='{}',orientation_source='computed_rotation',exif_orientation=1).items():
    obj['mapillary_'+key]=value
image=bpy.data.images.new('Selection panorama fixture',width=512,height=256)
image.pixels.foreach_set([.06,.12,.2,1]*(512*256))
image.pack();image.use_fake_user=True
image['mapillary_photo_source']='Mapillary'
image['mapillary_image_id']=obj['mapillary_image_id']
with bpy.context.temp_override(window=window,area=area,region=region):
    bpy.ops.object.armature_add(location=(5,12,2))
rig=bpy.context.object
rig.name='Pose mode fixture'
for selected in tuple(bpy.context.selected_objects):selected.select_set(False)
bpy.context.view_layer.objects.active=mesh
mesh.select_set(True)
bpy.context.view_layer.update()
initial_camera=bpy.context.scene.camera
pose=obj.matrix_world.copy()
mesh_pose=mesh.matrix_world.copy()


def selected():return set(bpy.context.selected_objects)


def native_event(kind,value='PRESS',point=None,shift=False):
    if point is None:
        point=(region.width//2,region.height//2)
    x,y=point
    window.event_simulate(type=kind,value=value,x=round(region.x+x),y=round(region.y+y),shift=shift)


def point(target):
    bpy.context.view_layer.update()
    result=location_3d_to_region_2d(region,area.spaces.active.region_3d,target.matrix_world.translation)
    assert result is not None and 80<result.x<region.width-80 and 130<result.y<region.height-80, (target.name,result)
    return tuple(result)


def click(target=None,shift=False,jitter=False):
    where=point(target) if target else (100,region.height*.58)
    if shift:native_event('LEFT_SHIFT',point=where)
    native_event('MOUSEMOVE','NOTHING',where,shift=shift)
    native_event('LEFTMOUSE',point=where,shift=shift)
    end=(where[0]+1,where[1]+1) if jitter else where
    if jitter:native_event('MOUSEMOVE','NOTHING',end,shift=shift)
    native_event('LEFTMOUSE','RELEASE',end,shift=shift)
    if shift:native_event('LEFT_SHIFT','RELEASE',end)


def open_view():
    state['session']=viewport.open(bpy.context,pano.Panorama(obj))
    state['surface']=state['session'].surface
    state['surface'].set_opacity(.35)
    state['session'].invalidate()


def stable(label):
    check(label,viewport.current() is state['session'] and not state['session'].closed
          and state['surface'].obj==obj and obj.matrix_world==pose
          and bpy.context.scene.camera==initial_camera)


def finish(error=None):
    if viewport.current():viewport.current().close('USER')
    addon.unregister()
    report=dict(status='failed' if error else 'passed',error=error,checks=checks,
                blender=bpy.app.version_string,elapsed=time.monotonic()-state['started'])
    (args.artifacts/'selection.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        assert time.monotonic()-state['started']<90,'Panorama selection check timed out'
        state['ticks']+=1
        if state['ticks']<4:return .08
        state['ticks']=0
        step=state['step']
        with bpy.context.temp_override(window=window,area=area,region=region):
            if step==0:
                native_event('ESC');native_event('ESC','RELEASE')
            elif step==1:
                check('shared photo opacity defaults to 30 percent',bpy.context.scene.blosm.mapillaryPhotoOpacity==30.)
                open_view()
                state['look']=state['surface'].look.copy()
            elif step==2:
                check('cached panorama ready without network: '+state['surface'].status,state['surface'].ready and runtime.active is None)
                click(mesh,jitter=True)
            elif step==3:
                check('native click selects a mesh despite a small pointer jitter',selected()=={mesh}
                      and bpy.context.object==mesh)
                check('click jitter does not rotate the panorama',state['surface'].look.rotation_difference(state['look']).angle<1e-5)
                stable('selection preserves panorama camera and render camera')
                click(other,shift=True)
            elif step==4:
                check('native Shift-click extends selection and changes active object',selected()=={mesh,other}
                      and bpy.context.object==other)
                click(other,shift=True)
            elif step==5:
                check('native Shift-click toggles the active selected object off',selected()=={mesh})
                click(shift=True)
            elif step==6:
                check('Shift-click empty space keeps the selection',selected()=={mesh})
                click(camera)
            elif step==7:
                if state.pop('toolbar_pending',False):
                    check('blank overlay toolbar does not pick through or clear selection',selected()=={camera})
                    click()
                else:
                    check('native click selects another Mapillary camera',selected()=={camera} and bpy.context.object==camera)
                    stable('another camera becomes active without replacing the panorama')
                    bpy.ops.screen.screenshot(filepath=str(args.artifacts/'camera_selected.png'))
                    toolbar=state['surface'].toolbar
                    scale=state['session'].frame.scale
                    where=((toolbar.x+toolbar.width-10)*scale,(toolbar.y+toolbar.height*.5)*scale)
                    native_event('LEFTMOUSE',point=where)
                    native_event('LEFTMOUSE','RELEASE',where)
                    state.update(step=6,toolbar_pending=True)
            elif step==8:
                check('plain click empty space deselects objects',not selected())
                state['before_drag']=state['surface'].look.copy()
                where=point(mesh)
                native_event('LEFTMOUSE',point=where)
                native_event('MOUSEMOVE','NOTHING',(where[0]+60,where[1]))
                native_event('LEFTMOUSE','RELEASE',(where[0]+60,where[1]))
            elif step==9:
                check('native drag rotates without selecting or moving geometry',not selected()
                      and state['surface'].look.rotation_difference(state['before_drag']).angle>.01
                      and state['surface'].navigation.yaw<0
                      and mesh.matrix_world==mesh_pose)
                check('drag keeps viewpoint at the panorama eye',
                      (area.spaces.active.region_3d.view_matrix.inverted().translation-obj.location).length<1e-4)
                slider=state['surface'].slider.rect
                where=(slider.x+slider.width*.5,slider.y+slider.height*.5)
                native_event('LEFTMOUSE',point=(where[0]*state['session'].frame.scale,where[1]*state['session'].frame.scale))
                native_event('LEFTMOUSE','RELEASE',(where[0]*state['session'].frame.scale,where[1]*state['session'].frame.scale))
            elif step==10:
                check('overlay slider handles clicks without selecting objects',not selected() and abs(state['surface'].opacity-.5)<.02)
                check('panorama slider updates the shared percentage',abs(bpy.context.scene.blosm.mapillaryPhotoOpacity-50)<2)
                bpy.context.scene.blosm.mapillaryPhotoOpacity=27.
                check('ordinary photo percentage updates a live panorama',abs(state['surface'].opacity-.27)<1e-6)
                click(mesh)
            elif step==11:
                check('mesh can be selected after looking around',bpy.context.object==mesh and selected()=={mesh})
                state['final_look']=area.spaces.active.region_3d.view_rotation.copy()
                native_event('TAB');native_event('TAB','RELEASE')
            elif step==12:
                check('native Tab enters Edit Mode while keeping panorama',viewport.current() is state['session']
                      and not state['session'].closed and bpy.context.mode=='EDIT_MESH')
                check('mode change retains final viewport orientation',area.spaces.active.region_3d.view_rotation.rotation_difference(state['final_look']).angle<1e-5)
                check('Edit Mode clicks and transform keys pass through to Blender',all(
                    state['surface'].handle(ui.InputEvent(key,'PRESS',100,region.height*.58))==ui.PASS_THROUGH
                    for key in ('LEFTMOUSE','G','R','S')) and not state['session'].closed)
                bpy.ops.transform.translate(value=(.1,0,0))
                state['edit_look']=state['surface'].look.copy()
                where=(150,region.height*.58)
                native_event('MIDDLEMOUSE',point=where)
                native_event('MOUSEMOVE','NOTHING',(where[0]+60,where[1]))
                native_event('MIDDLEMOUSE','RELEASE',(where[0]+60,where[1]))
            elif step==13:
                check('middle-drag looks around in Edit Mode without closing panorama',viewport.current() is state['session']
                      and state['surface'].look.rotation_difference(state['edit_look']).angle>.01)
                bpy.ops.object.mode_set(mode='OBJECT')
                bpy.ops.object.mode_set(mode='SCULPT')
            elif step==14:
                check('operator/UI-style Sculpt Mode change keeps panorama',bpy.context.mode=='SCULPT'
                      and not state['session'].closed and viewport.current() is state['session'])
                check('Sculpt brush strength shortcut stays native',
                      state['surface'].handle(ui.InputEvent('F','PRESS',100,region.height*.58,shift=True))==ui.PASS_THROUGH
                      and not state['session'].closed)
                open_view()
                check('a cached panorama can attach in another mode',viewport.current() is state['session'])
                bpy.ops.object.mode_set(mode='OBJECT')
                open_view()
                bpy.ops.object.mode_set(mode='WEIGHT_PAINT')
            elif step==15:
                check('Weight Paint Mode keeps the panorama',bpy.context.mode=='PAINT_WEIGHT'
                      and not state['session'].closed and viewport.current() is state['session'])
                bpy.ops.object.mode_set(mode='OBJECT')
                open_view()
                bpy.ops.object.mode_set(mode='VERTEX_PAINT')
            elif step==16:
                check('Vertex Paint Mode keeps the panorama',bpy.context.mode=='PAINT_VERTEX'
                      and not state['session'].closed and viewport.current() is state['session'])
                bpy.ops.object.mode_set(mode='OBJECT')
                open_view()
                bpy.ops.object.mode_set(mode='TEXTURE_PAINT')
            elif step==17:
                check('Texture Paint Mode keeps the panorama',bpy.context.mode=='PAINT_TEXTURE'
                      and not state['session'].closed and viewport.current() is state['session'])
                bpy.ops.object.mode_set(mode='OBJECT')
                for item in tuple(selected()):item.select_set(False)
                rig.select_set(True)
                bpy.context.view_layer.objects.active=rig
                open_view()
                bpy.ops.object.mode_set(mode='POSE')
            elif step==18:
                check('Pose Mode keeps the panorama',bpy.context.mode=='POSE'
                      and not state['session'].closed and viewport.current() is state['session'])
                bpy.ops.object.mode_set(mode='OBJECT')
                check('all selection/mode checks preserve the source pose and scene camera',
                      obj.matrix_world==pose and bpy.context.scene.camera==initial_camera)
                state['surface'].set_opacity(.42)
                state['session'].close('USER')
                photo.show_photo(bpy.context,camera,photo.cached_image(obj['mapillary_image_id']))
                check('ordinary photo inherits opacity changed on a panorama',
                      all(abs(bg.alpha-.42)<1e-6 for bg in camera.data.background_images if bg.image))
                bpy.context.scene.blosm.mapillaryPhotoOpacity=17.
                open_session=viewport.open(bpy.context,pano.Panorama(obj))
                check('reopened panorama inherits opacity changed on an ordinary photo',abs(open_session.surface.opacity-.17)<1e-6)
                finish()
                return None
        state['step']+=1
        return .08
    except Exception:
        finish(traceback.format_exc())
        return None


bpy.app.timers.register(tick,first_interval=.5)
