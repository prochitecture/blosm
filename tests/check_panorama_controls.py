"""Isolated native upright-view, zoom-animation and navigation lifecycle checks.

No network or credentials. Use factory startup, disable-autoexec, event simulation
and all five isolated resource paths. Every output belongs in --artifacts.
"""
import argparse
from contextlib import nullcontext
import importlib
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time
import traceback
from types import SimpleNamespace

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
from mathutils import Quaternion,Vector
from blosm.viewport.ui import InputEvent
parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
addon=importlib.import_module(ROOT.name);addon.register()
viewport=importlib.import_module(ROOT.name+'.viewport')
pano=importlib.import_module(ROOT.name+'.mapillary.panorama')
host=importlib.import_module(ROOT.name+'.viewport.session')
window=bpy.context.window
area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
region=next(r for r in area.regions if r.type=='WINDOW')
checks=[];metrics={}
state={'phase':0,'ticks':0,'started':time.monotonic(),'frames':[],'polls':0}


def check(label,condition):
    assert condition,label
    checks.append(label)


obj=bpy.data.objects.new('Rolled capture pose',bpy.data.cameras.new('Panorama'))
bpy.context.scene.collection.objects.link(obj)
obj.location=(1,2,3)
obj.rotation_mode='QUATERNION'
obj.rotation_quaternion=(Quaternion((0,0,1),.7) @ Quaternion((1,0,0),math.pi/2+.2)
                         @ Quaternion((0,0,1),2.6))
for key,value in dict(image_id='9007199254740993',camera_model='spherical',
        metadata='{}',orientation_source='computed_rotation',exif_orientation=1).items():
    obj['mapillary_'+key]=value
image=bpy.data.images.new('Cached control fixture',width=512,height=256)
image.pixels.foreach_set([.08,.2,.3,1]*(512*256));image.pack();image.use_fake_user=True
image['mapillary_photo_source']='Mapillary';image['mapillary_image_id']=obj['mapillary_image_id']
bpy.context.view_layer.update()
source_pose=obj.matrix_world.copy()
scene_camera=bpy.context.scene.camera


def next_phase(number):state.update(phase=number,ticks=0)


def open_view():
    state['session']=viewport.open(bpy.context,pano.Panorama(obj))
    surface=state['surface']=state['session'].surface
    original_animate,original_tick=surface.animate,surface.tick
    def animate(seconds):
        state['frames'].append(seconds)
        return original_animate(seconds)
    def poll(budget):
        state['polls']+=1
        return original_tick(budget)
    surface.animate,surface.tick=animate,poll
    return surface


def upright():
    rv=area.spaces.active.region_3d
    up=rv.view_rotation @ Vector((0,1,0))
    right=rv.view_rotation @ Vector((1,0,0))
    return up.z>0 and abs(right.z)<1e-5


def wheel(steps):
    state['surface'].handle(InputEvent('WHEELUPMOUSE' if steps>0 else 'WHEELDOWNMOUSE',
                                     'PRESS',400,300))


def drag(dx,dy):
    surface=state['surface']
    surface.handle(InputEvent('LEFTMOUSE','PRESS',300,300))
    surface.handle(InputEvent('MOUSEMOVE','NOTHING',300+dx,300+dy))
    surface.handle(InputEvent('LEFTMOUSE','RELEASE',300+dx,300+dy))


def finish(error=None):
    if viewport.current():viewport.current().close('USER')
    addon.unregister()
    result=dict(status='failed' if error else 'passed',error=error,checks=checks,metrics=metrics,
                blender=bpy.app.version_string,elapsed=time.monotonic()-state['started'])
    (args.artifacts/'controls.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        assert time.monotonic()-state['started']<120,'Control test timed out'
        state['ticks']+=1
        with bpy.context.temp_override(window=window,area=area,region=region):
            p=state['phase']
            rv=area.spaces.active.region_3d
            if p==0:
                window.event_simulate(type='ESC',value='PRESS');window.event_simulate(type='ESC',value='RELEASE')
                next_phase(1)
            elif p==1 and state['ticks']>5:
                open_view();next_phase(2)
            elif p==2 and state['ticks']>5:
                surface=state['surface']
                source_forward=obj.matrix_world.to_quaternion() @ Vector((0,0,-1))
                view_forward=rv.view_rotation @ Vector((0,0,-1))
                check('rolled source pose opens with upright view axes',upright()
                      and (obj.matrix_world.to_quaternion() @ Vector((0,1,0))).z<0)
                check('opening keeps source viewing direction without copying capture roll',
                      (source_forward-view_forward).length<1e-5)
                state['opening_fov']=surface.fov
                wheel(1)
                check('wheel changes target without an immediate projection jump',surface.animating and surface.fov==state['opening_fov'])
                state['frames']=[];state['polls']=0;state['fovs']=[]
                next_phase(3)
            elif p==3:
                surface=state['surface'];state['fovs'].append(surface.fov)
                assert upright(),'Zoom introduced roll'
                if not surface.animating:
                    values=state['fovs']
                    check('wheel zoom eases monotonically and reaches its target',len(values)>3
                          and all(b<=a for a,b in zip(values,values[1:])) and surface.fov==surface.navigation.target_fov)
                    check('animation frames do not run background preparation every frame',len(state['frames'])>2*state['polls'])
                    check('native projection agrees with zoomed vertical FOV',abs(rv.window_matrix[1][1]-1/math.tan(surface.fov/2))<.002)
                    metrics['zoom_frames']=len(state['frames']);metrics['work_polls_during_zoom']=state['polls']
                    metrics['median_animation_interval_ms']=statistics.median(state['frames'])*1000
                    for _ in range(80):
                        window.event_simulate(type='WHEELUPMOUSE',value='PRESS',x=region.x+300,y=region.y+300)
                    next_phase(4)
            elif p==4 and state['ticks']>5 and not state['surface'].animating:
                surface=state['surface']
                check('native wheel burst stops at minimum FOV',abs(surface.fov-math.radians(30))<1e-6)
                for _ in range(80):
                    window.event_simulate(type='WHEELDOWNMOUSE',value='PRESS',x=region.x+300,y=region.y+300)
                next_phase(5)
            elif p==5 and state['ticks']>5 and not state['surface'].animating:
                surface=state['surface']
                hfov=2*math.atan(1/rv.window_matrix[0][0])
                check('zoom-out respects both horizontal and vertical projection bounds',hfov<=math.radians(120)+1e-5
                      and surface.fov<=math.radians(90) and surface.fov==surface.navigation.maximum_fov)
                drag(0,-100000)
                check('huge upward drag stops below the north pole',upright() and abs(surface.navigation.pitch-math.radians(89))<1e-8)
                drag(0,10)
                check('reversing at the pole responds immediately',surface.navigation.pitch<math.radians(89))
                drag(0,100000)
                check('huge downward drag stops above the south pole',upright() and abs(surface.navigation.pitch+math.radians(89))<1e-8)
                rng=random.Random(83)
                for _ in range(100):
                    drag(rng.uniform(-3000,3000),rng.uniform(-2000,2000))
                    assert upright(),'Mixed dragging inverted the native view'
                check('mixed repeated native viewport rotations cannot accumulate roll',upright())
                check('rotating and zooming keep the eye fixed',(rv.view_matrix.inverted().translation-obj.location).length<1e-4)
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'bounded_view.png'))
                wheel(1)
                current=surface.fov
                state['session'].event(bpy.context,SimpleNamespace(type='WINDOW_DEACTIVATE'))
                check('focus loss stops zoom without moving the current view',not surface.animating and surface.fov==current
                      and not state['session'].closed)
                wheel(1);state['old_width']=region.width
                bpy.ops.screen.area_split(direction='VERTICAL',factor=.5)
                next_phase(6)
            elif p==6 and state['ticks']>8:
                surface=state['surface']
                check('resize cancels animation while preserving upright navigation',region.width!=state['old_width']
                      and not surface.animating and upright() and surface.owns_view())
                check('resized native projection matches the effective FOV',abs(rv.window_matrix[1][1]-1/math.tan(surface.fov/2))<.002)
                wheel(1);current=surface.fov
                surface.handle(InputEvent('LEFTMOUSE','PRESS',300,300))
                check('selection press stops zoom at the displayed projection',not surface.animating and surface.fov==current)
                surface.cancel_gesture()
                wheel(1);next_phase(7)
            elif p==7 and state['ticks']>3:
                state['exit_view']=state['surface'].view_state()
                state['session'].close('USER');next_phase(8)
            elif p==8 and state['ticks']>8:
                check('exit during zoom keeps the displayed orientation and lens',state['surface'].view_state()==state['exit_view'])
                check('animation stops and the shared timer becomes idle after exit',viewport.current() is None
                      and not bpy.app.timers.is_registered(host._tick))
                surface=open_view()
                rv.view_rotation=Quaternion((0,1,0),math.pi)
                state['external_view']=surface.view_state()
                check('external native navigation is detected before sampling the image',not surface.owns_view())
                class Painter:
                    def rect(self,*_):pass
                    def text(self,*_,**__):pass
                    def clip(self,*_):return nullcontext()
                    def pixel_clip(self):raise AssertionError('An unowned projection sampled the panorama')
                surface.draw(state['session'].frame,Painter())
                check('photo drawing skips an externally inverted projection',True)
                next_phase(9)
            elif p==9 and state['ticks']>12:
                check('external navigation closes the viewer without restoring an old view',viewport.current() is None
                      and state['surface'].view_state()==state['external_view'])
                mesh=bpy.data.objects['Cube'];mesh.select_set(True);bpy.context.view_layer.objects.active=mesh
                open_view();wheel(1)
                bpy.ops.object.mode_set(mode='EDIT');next_phase(10)
            elif p==10 and state['ticks']>8:
                check('mode change keeps panorama and stops the zoom gesture',viewport.current() is state['session']
                      and not state['session'].closed and not state['surface'].animating)
                bpy.ops.object.mode_set(mode='OBJECT')
                check('all navigation preserves imported pose and rendering camera',obj.matrix_world==source_pose
                      and bpy.context.scene.camera==scene_camera)
                finish();return None
        return .02
    except Exception:
        finish(traceback.format_exc());return None


bpy.app.timers.register(tick,first_interval=.5)
