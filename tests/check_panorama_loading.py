"""Isolated native UI loading/lifecycle checks using delayed HTTP fixtures."""
import argparse
import importlib
import json
from pathlib import Path
import struct
import sys
import time
import traceback
from urllib.parse import parse_qs,urlsplit
import zlib

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
runtime=importlib.import_module(ROOT.name+'.progressive.blender')
example=importlib.import_module(ROOT.name+'.viewport.example')
runtime.show_result=lambda *_:None
photo.get_token=lambda *_:'fixture-token'
window=bpy.context.window
area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
region=next(r for r in area.regions if r.type=='WINDOW')
checks=[];transports=[]
state={'phase':0,'frames':0,'start':time.monotonic()}


def check(label,condition):
    assert condition,label
    checks.append(label)


def png():
    def chunk(kind,data):return struct.pack('>I',len(data))+kind+data+struct.pack('>I',zlib.crc32(kind+data)&0xffffffff)
    data=(b'\0'+bytes((30,140,80,255))*256)*128
    return b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',256,128,8,6,0,0,0))+chunk(b'IDAT',zlib.compress(data))+chunk(b'IEND',b'')


class Transport:
    max_requests=1
    def __init__(self,**kwargs):
        self.pending=None;self.requests=[];self.closed=False;self.release=False;self.is_photo=False;self.failure=False
        transports.append(self)
    def submit(self,request):
        self.requests.append(request)
        query=parse_qs(urlsplit(request.url).query)
        if 'image_ids' in query:
            data=json.dumps({'data':[{'id':query['image_ids'][0],'thumb_2048_url':'https://scontent.example.fbcdn.net/photo.jpg'}]}).encode()
            assert request.headers.get('Authorization')
        else:
            assert not request.headers
            self.is_photo=True;data=png()
        self.pending=(len(self.requests),data,None)
        return len(self.requests)
    def poll(self):
        if self.is_photo and not self.release:return None
        result,self.pending=self.pending,None
        if result and self.failure and self.is_photo:return result[0],None,RuntimeError('Fixture HTTP failure')
        return result
    def close(self):self.closed=True;self.pending=None


photo.DownloadTransport=Transport


def camera(identity):
    obj=bpy.data.objects.new('Load pano '+identity,bpy.data.cameras.new('Pano'))
    bpy.context.scene.collection.objects.link(obj)
    obj.location=(1,2,3);obj.rotation_mode='QUATERNION';obj.rotation_quaternion=Quaternion((1,0,0),1.57079632679)
    for key,value in {'image_id':identity,'camera_model':'spherical','metadata':'{}','orientation_source':'computed_rotation','exif_orientation':1}.items():
        obj['mapillary_'+key]=value
    bpy.context.view_layer.update()
    return obj


def start(identity):
    obj=camera(identity)
    surface=pano.Panorama(obj)
    session=viewport.open(bpy.context,surface)
    state.update(surface=surface,session=session,obj=obj)
    return session


def phase(value):state.update(phase=value,frames=0)


def finish(error=None):
    if viewport.current():viewport.current().close('USER')
    if runtime.active:runtime.active.close(silent=True)
    addon.unregister()
    result={'status':'failed' if error else 'passed','checks':checks,'error':error,'blender':bpy.app.version_string,'transports':len(transports)}
    (args.artifacts/'loading.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        assert time.monotonic()-state['start']<90,'deadline'
        state['frames']+=1
        with bpy.context.temp_override(window=window,area=area,region=region):
            p=state['phase']
            if p==0:
                window.event_simulate(type='ESC',value='PRESS');window.event_simulate(type='ESC',value='RELEASE')
                phase(1)
            elif p==1:
                state['old_camera']=bpy.context.scene.camera
                state['old_render']=(bpy.context.scene.render.resolution_x,bpy.context.scene.render.resolution_y)
                start('9100000000000001');phase(2)
            elif p==2 and transports[-1].is_photo:
                check('uncached viewer uses finite shared photo task',runtime.active is state['surface'].controller and not state['surface'].ready)
                transports[-1].release=True;phase(3)
            elif p==3 and state['surface'].ready:
                check('publication completes before display',runtime.active is None and transports[-1].closed)
                check('acquisition does not activate render camera or resize render',bpy.context.scene.camera==state['old_camera'] and state['old_render']==(bpy.context.scene.render.resolution_x,bpy.context.scene.render.resolution_y))
                check('published image is packed',photo.cached_image('9100000000000001').packed_file is not None)
                bpy.ops.wm.save_as_mainfile(filepath=str(args.artifacts/'panorama_open.blend'))
                check('saving keeps live surface usable',not state['session'].closed)
                bpy.ops.ed.undo()
                check('undo closes viewer before restoring IDs',viewport.current() is None)
                bpy.ops.ed.redo()
                check('redo restores image without restarting viewer',viewport.current() is None and photo.cached_image('9100000000000001') is not None)
                phase(4)
            elif p==4:
                start('9100000000000002');phase(5)
            elif p==5 and transports[-1].is_photo:
                t=transports[-1]
                state['session'].close('USER')
                check('close cancels only the owned task',runtime.active is None and t.closed and photo.cached_image('9100000000000002') is None)
                start('9100000000000003');phase(6)
            elif p==6 and transports[-1].is_photo:
                transports[-1].failure=True;transports[-1].release=True;phase(7)
            elif p==7 and state['surface'].failed:
                check('network failure retains usable exit/retry UI',not state['session'].closed and runtime.active is None)
                state['surface'].retry_load();phase(8)
            elif p==8 and transports[-1].is_photo and not transports[-1].closed:
                transports[-1].release=True;phase(9)
            elif p==9 and state['surface'].ready:
                check('retry can recover',photo.cached_image('9100000000000003') is not None)
                state['session'].close('USER')
                # A different photo/import owns the busy slot. Closing a queued
                # panorama must not cancel that task.
                obj=camera('9100000000000004')
                task=photo.PhotoImport(bpy.context,obj,'fixture-token',activate=False,report_result=False)
                owner=runtime.start(task,bpy.context)
                count=len(transports)
                start('9100000000000005')
                check('queued viewer starts no second helper',len(transports)==count and 'Queued' in state['surface'].status)
                state['session'].close('USER')
                check('queued close preserves unrelated import',runtime.active is owner)
                owner.close(silent=True)
                phase(10)
            elif p==10:
                demo=example.show(bpy.context,photo.cached_image('9100000000000001'))
                state['demo']=demo
                state['old_sidebar']=area.spaces.active.show_region_ui
                area.spaces.active.show_region_ui=True
                phase(11)
            elif p==11 and state['frames']>5:
                check('generic surface draws with a sidebar',not state['demo'].closed and state['demo'].error is None)
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'framework_example.png'))
                example_state=state['demo'].surface
                from blosm.viewport.ui import Frame,Rect
                example_state.layout(Frame(Rect(0,0,300,200),Rect(0,0,300,200),1))
                check('side card collapses on narrow layout',example_state.card is None and not example_state.slider.visible)
                state['demo'].close('USER')
                area.spaces.active.show_region_ui=state['old_sidebar']
                obj=camera('9100000000000001')
                state['session']=viewport.open(bpy.context,pano.Panorama(obj))
                state['scene']=window.scene
                window.scene=bpy.data.scenes.new('Viewer switch test')
                phase(12)
            elif p==12 and state['frames']>4:
                check('scene change cleans the session',viewport.current() is None and state['session'].closed)
                window.scene=state['scene']
                obj=camera('9100000000000001')
                bpy.context.view_layer.objects.active=obj
                result=bpy.ops.blosm.view_mapillary_panorama()
                check('public operator opens the selected panorama',result=={'FINISHED'} and viewport.current() is not None)
                state['final_rotation']=area.spaces.active.region_3d.view_rotation.copy()
                window.event_simulate(type='ESC',value='PRESS',x=region.x+300,y=region.y+300)
                window.event_simulate(type='ESC',value='RELEASE',x=region.x+300,y=region.y+300)
                phase(13)
            elif p==13 and state['frames']>4:
                check('real Escape event exits without rotation reset',viewport.current() is None and area.spaces.active.region_3d.view_rotation.rotation_difference(state['final_rotation']).angle<1e-5)
                finish();return None
        return .1
    except Exception:
        finish(traceback.format_exc());return None

bpy.app.timers.register(tick,first_interval=1.)
