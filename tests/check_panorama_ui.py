"""Native isolated UI smoke/regression harness; no credentials or user profile."""
import argparse
import importlib
import json
from pathlib import Path
import sys
import time
import traceback

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
from mathutils import Quaternion, Vector
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
runtime=importlib.import_module(ROOT.name+'.progressive.blender')
runtime.show_result=lambda *_:None
window=bpy.context.window
area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
region=next(r for r in area.regions if r.type=='WINDOW')
checks=[]
state={'phase':0,'started':time.monotonic(),'frames':0,'last':None,'gaps':[]}


def check(label,condition):
    assert condition,label
    checks.append(label)


def camera(identity):
    obj=bpy.data.objects.new('Pano fixture '+identity,bpy.data.cameras.new('Pano fixture'))
    bpy.context.scene.collection.objects.link(obj)
    obj.location=(0,0,2)
    obj.rotation_mode='QUATERNION'
    obj.rotation_quaternion=Quaternion((1,0,0),1.5707963267948966)
    obj['mapillary_image_id']=identity
    obj['mapillary_camera_model']='spherical'
    obj['mapillary_metadata']='{}'
    obj['mapillary_orientation_source']='computed_rotation'
    obj['mapillary_exif_orientation']=1
    bpy.context.view_layer.update()
    return obj


obj=camera('9007199254740993')
image=bpy.data.images.new('Directional panorama fixture',width=1024,height=512)
import numpy as np
pixels=np.zeros((512,1024,4),dtype=np.float32)
pixels[:,:,3]=1
pixels[:,:,:3]=(.03,.1,.2)
pixels[:,256:512,:3]=(.8,.05,.04)
pixels[:,512:768,:3]=(.03,.7,.07)
pixels[::32,:,:3]=1
pixels[:,::64,:3]=1
image.pixels.foreach_set(pixels.ravel())
image.pack();image.use_fake_user=True
image['mapillary_photo_source']='Mapillary'
image['mapillary_image_id']=obj['mapillary_image_id']
image['mapillary_creator']='Synthetic direction fixture'
initial_camera=bpy.context.scene.camera
initial_render=(bpy.context.scene.render.resolution_x,bpy.context.scene.render.resolution_y)
initial_active=bpy.context.view_layer.objects.active
pose=obj.matrix_world.copy()


def check_gpu_projection(surface):
    """Check actual GPU colors against direction/EXIF landmarks, not a screenshot."""
    import gpu
    from mathutils import Matrix,Vector
    from blosm.viewport.assets import ImageRef
    from blosm.mapillary.panorama_math import image_uv
    gradient=bpy.data.images.new('GPU test gradient',width=256,height=128,float_buffer=True)
    pixels=np.zeros((128,256,4),dtype=np.float32)
    pixels[:,:,0]=(np.arange(256)+.5)[None,:]/256
    pixels[:,:,1]=(np.arange(128)+.5)[:,None]/128
    pixels[:,:,3]=1
    gradient.pixels.foreach_set(pixels.ravel())
    ref=ImageRef(gradient)
    target=gpu.types.GPUOffScreen(64,64,format='RGBA16F')
    old=(gpu.state.blend_get(),gpu.state.depth_test_get(),gpu.state.depth_mask_get(),gpu.state.viewport_get())
    errors=[]
    def display(value):return 12.92*value if value<=.0031308 else 1.055*value**(1/2.4)-.055
    try:
        with target.bind():
            gpu.state.viewport_set(0,0,64,64)
            gpu.state.blend_set('NONE');gpu.state.depth_test_set('NONE');gpu.state.depth_mask_set(False)
            for orientation in (1,3,6,8):
                for tilt in (False,True):
                    matrix=area.spaces.active.region_3d.window_matrix.inverted()
                    if tilt:matrix=Matrix.Rotation(.2,4,'Z') @ Matrix.Rotation(.3,4,'X') @ matrix
                    shader=surface.shader
                    shader.bind();shader.uniform_float('ray_matrix',matrix)
                    shader.uniform_float('opacity',1.)
                    shader.uniform_int('orientation',orientation)
                    shader.uniform_int('premultiplied',0)
                    shader.uniform_float('clip_rect',(0,0,64,64))
                    shader.uniform_sampler('photo_texture',ref.texture)
                    surface.batch.draw(shader)
                    for x,y in ((17,23),(46,38)):
                        ray=(matrix @ Vector(((x+.5)/32-1,(y+.5)/32-1,-1,1))).xyz
                        expected=image_uv(ray,orientation)
                        pixel=gpu.state.active_framebuffer_get().read_color(x,y,1,1,4,0,'FLOAT')
                        pixel.dimensions=4
                        errors.extend(abs(pixel[i]-display(expected[i])) for i in (0,1))
        check('GPU sampling matches raw/basic EXIF and pitched/rolled rays',max(errors)<.012)
    finally:
        gpu.state.blend_set(old[0]);gpu.state.depth_test_set(old[1]);gpu.state.depth_mask_set(old[2]);gpu.state.viewport_set(*old[3])
        target.free();ref.release();bpy.data.images.remove(gradient)


def finish(error=None):
    current=viewport.current()
    if current: current.close('USER')
    addon.unregister()
    report={'status':'failed' if error else 'passed','error':error,'checks':checks,
            'elapsed':time.monotonic()-state['started'],'heartbeat_max':max(state['gaps'],default=0),
            'blender':bpy.app.version_string}
    (args.artifacts/'ui.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        now=time.monotonic()
        if state['last'] is not None:state['gaps'].append(now-state['last'])
        state['last']=now
        assert now-state['started']<60,'UI test timed out'
        state['frames']+=1
        with bpy.context.temp_override(window=window,area=area,region=region):
            phase=state['phase']
            if phase==0:
                window.event_simulate(type='ESC',value='PRESS')
                window.event_simulate(type='ESC',value='RELEASE')
                state['phase']=1
            elif phase==1:
                surface=pano.Panorama(obj)
                session=viewport.open(bpy.context,surface)
                state.update(surface=surface,session=session,phase=2,frames=0)
            elif phase==2 and state['frames']>=8:
                session,surface=state['session'],state['surface']
                check('cached image opens without import/network',surface.ready and runtime.active is None)
                check('GPU draw succeeds',session.error is None and not session.closed)
                check('camera and scene state untouched',obj.matrix_world==pose and bpy.context.scene.camera==initial_camera
                      and initial_render==(bpy.context.scene.render.resolution_x,bpy.context.scene.render.resolution_y)
                      and bpy.context.view_layer.objects.active==initial_active)
                rv=area.spaces.active.region_3d
                eye=rv.view_matrix.inverted().translation
                check('viewport eye matches capture position',(eye-obj.location).length<1e-4)
                check('vertical field of view matches request',abs(rv.window_matrix[1][1]-1/np.tan(surface.fov/2))<.01)
                check_gpu_projection(surface)
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'panorama.png'))
                from blosm.viewport.ui import InputEvent
                surface.handle(InputEvent('LEFTMOUSE','PRESS',400,300))
                surface.handle(InputEvent('MOUSEMOVE','NOTHING',470,320))
                surface.handle(InputEvent('LEFTMOUSE','RELEASE',470,320))
                surface.handle(InputEvent('LEFTMOUSE','PRESS',400,300))
                surface.handle(InputEvent('LEFTMOUSE','RELEASE',surface.close.rect.x+10,surface.close.rect.y+10))
                check('canvas release over a control cannot leave a stuck drag',surface.drag is None and not session.closed)
                surface.handle(InputEvent('WHEELUPMOUSE','PRESS',400,300))
                surface.set_opacity(.35)
                session.invalidate()
                state.update(phase=3,frames=0)
            elif phase==3 and state['frames']>=5:
                surface=state['surface']
                check('look and zoom keep eye fixed',(area.spaces.active.region_3d.view_matrix.inverted().translation-obj.location).length<1e-4)
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'panorama_transparent.png'))
                state['old_width']=region.width
                state['old_fov']=surface.fov
                bpy.ops.screen.area_split(direction='VERTICAL',factor=.45)
                state.update(phase=31,frames=0)
            elif phase==31 and state['frames']>=5:
                surface=state['surface']
                check('real editor resize updates layout without ending the viewer',not state['session'].closed and region.width!=state['old_width']
                      and abs(state['session'].frame.region.width*state['session'].frame.scale-region.width)<1)
                check('resize preserves vertical FOV and capture position',surface.fov==state['old_fov']
                      and abs(area.spaces.active.region_3d.window_matrix[1][1]-1/np.tan(surface.fov/2))<.01
                      and (area.spaces.active.region_3d.view_matrix.inverted().translation-obj.location).length<1e-4)
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'resized.png'))
                final=area.spaces.active.region_3d.view_rotation.copy()
                state['session'].close('USER')
                check('exit preserves final orientation',area.spaces.active.region_3d.view_rotation.rotation_difference(final).angle<1e-5)
                check('exit releases handler without deleting photo',viewport.current() is None and state['session'].handler is None and image.name in bpy.data.images)
                session=viewport.open(bpy.context,pano.Panorama(obj))
                obj.name='Renamed panorama'
                session.surface.validate()
                check('rename keeps identity',not session.closed)
                bpy.data.objects.remove(obj,do_unlink=True)
                state.update(session=session,phase=4,frames=0)
            elif phase==4 and state['frames']>=5:
                check('deletion cleans viewer',state['session'].closed and viewport.current() is None)
                finish()
                return None
        return .1
    except Exception:
        finish(traceback.format_exc())
        return None

bpy.app.timers.register(tick,first_interval=1.)
