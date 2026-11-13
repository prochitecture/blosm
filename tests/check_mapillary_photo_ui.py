"""Native UI timers, undo and packed-photo reuse with offline HTTP fixtures.

Run only in a separate factory-startup, disable-autoexec Blender UI process with
the chat's five isolated resource directories. No credentials/network are used.
"""
import argparse
import importlib
import json
from pathlib import Path
import struct
import sys
import time
import traceback
from urllib.parse import parse_qs, urlsplit
import zlib

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type = 'NONE'
assert not bpy.app.background and bpy.context.window is not None
addon = importlib.import_module(ROOT.name)
addon.register()
photo = importlib.import_module(ROOT.name+'.mapillary.photo')
runtime = importlib.import_module(ROOT.name+'.progressive.blender')
runtime.show_result = lambda *_args: None
photo.get_token = lambda *_args:'fixture-photo-token'


def png(width=320,height=180):
    def chunk(kind,payload):
        return struct.pack('>I',len(payload))+kind+payload+struct.pack('>I',zlib.crc32(kind+payload)&0xffffffff)
    rows = b''.join(b'\0'+bytes((40+y//2,110,210,255))*width for y in range(height))
    return b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',width,height,8,6,0,0,0))+chunk(b'IDAT',zlib.compress(rows))+chunk(b'IEND',b'')


payload = png()
transports = []
checks = []
heartbeats = []
state = {'phase':'start','ticks':0,'started':time.monotonic(),'previous':None,'registered':True}
scene_name = bpy.context.scene.name
window = bpy.context.window


def check(label,condition):
    assert condition,label
    checks.append(label)


class Transport:
    max_requests = 1
    def __init__(self,**_options):
        self.requests = []
        self.pending = None
        self.photo_pending = False
        self.release = len(transports)==0
        self.closed = False
        transports.append(self)

    def submit(self,request):
        assert self.pending is None
        self.requests.append(request)
        query = parse_qs(urlsplit(request.url).query)
        if 'image_ids' in query:
            identity = query['image_ids'][0]
            content = json.dumps({'data':[{'id':identity,'creator':{'username':'Offline fixture'},
                'thumb_2048_url':'https://scontent.example.fbcdn.net/photo.jpg'}]}).encode()
        else:
            assert not request.headers,'OAuth must not be sent to the photo CDN'
            content = payload
            self.photo_pending = True
        self.pending = (len(self.requests),content,None)
        return len(self.requests)

    def poll(self):
        if self.photo_pending and not self.release:
            return None
        result,self.pending = self.pending,None
        return result

    def close(self):
        self.closed = True
        self.pending = None


photo.DownloadTransport = Transport


def scene():
    return bpy.data.scenes[scene_name]


def camera(identity):
    name = 'Mapillary UI '+identity
    obj = bpy.data.objects.get(name)
    if obj is None:
        obj = bpy.data.objects.new(name,bpy.data.cameras.new(name))
        scene().collection.objects.link(obj)
        from mathutils import Vector
        obj.location = (7.,-7.,5.)
        obj.rotation_euler = (Vector((0.,0.,0.))-obj.location).to_track_quat('-Z','Y').to_euler()
        obj['mapillary_image_id'] = identity
        obj['mapillary_metadata'] = '{}'
        obj['mapillary_orientation_source'] = 'computed_rotation'
        obj['mapillary_image_url'] = 'https://www.mapillary.com/app/?pKey='+identity
    for selected in tuple(bpy.context.selected_objects):
        selected.select_set(False)
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    return obj


def start(identity):
    camera(identity)
    assert bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT',True)=={'FINISHED'}
    assert runtime.active is not None
    state['task'] = runtime.active.adapter
    state['ticks'] = 0


def displayed(identity):
    obj = bpy.data.objects['Mapillary UI '+identity]
    image = photo.cached_image(identity)
    area = next(area for area in window.screen.areas if area.type=='VIEW_3D')
    background = next((bg for bg in obj.data.background_images if bg.image==image),None)
    return (image is not None and image.packed_file is not None and scene().camera==obj
            and background is not None and background.frame_method=='FIT'
            and background.alpha==1. and background.display_depth=='FRONT'
            and obj.data.show_background_images and area.spaces.active.region_3d.view_perspective=='CAMERA'
            and tuple(image.size)==(320,180) and scene().render.resolution_x==320
            and scene().render.resolution_y==180)


def finish(status,error=None):
    if runtime.active is not None:
        runtime.active.close(silent=True)
    if state['registered']:
        addon.unregister()
        state['registered'] = False
    report = {'status':status,'error':error,'checks':checks,'heartbeats':len(heartbeats),
              'max_heartbeat_gap':max(heartbeats,default=0),'elapsed':time.monotonic()-state['started'],
              'transport_count':len(transports),'blender':bpy.app.version_string}
    (args.artifacts/'ui.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        now = time.monotonic()
        if state['previous'] is not None: heartbeats.append(now-state['previous'])
        state['previous'] = now
        assert now-state['started']<90,'UI photo fixture exceeded its deadline'
        phase = state['phase']
        state['ticks'] += 1
        if phase=='start':
            # Requires --enable-event-simulate, and affects only this isolated
            # test window. Dismiss the startup splash before visual capture.
            window.event_simulate(type='ESC',value='PRESS')
            window.event_simulate(type='ESC',value='RELEASE')
            scene().camera = None
            bpy.context.preferences.edit.use_global_undo = True
            start('9007199254740993')
            state['phase'] = 'download'
        elif phase=='download' and runtime.active is None:
            check('native timer download activates camera and fills its frame with a packed photo',displayed('9007199254740993'))
            check('metadata and CDN requests share the bounded helper',len(transports)==1 and len(transports[0].requests)==2 and transports[0].closed)
            state['phase'],state['ticks'] = 'screenshot',0
        elif phase=='screenshot' and state['ticks']>=8:
            bpy.ops.screen.screenshot(filepath=str(args.artifacts/'photo_view.png'))
            bpy.ops.ed.undo()
            check('one native Undo restores the pre-photo camera state',scene().camera is None)
            bpy.ops.ed.redo()
            check('native Redo restores the camera and packed image without downloads',displayed('9007199254740993') and len(transports)==1)
            camera('9007199254740993')
            assert bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT',True)=={'FINISHED'}
            check('cached operator reuses Blender image data without a transport',displayed('9007199254740993') and len(transports)==1 and runtime.active is None)
            bg = next(bg for bg in scene().camera.data.background_images if bg.image)
            bg.alpha = .35
            assert bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT',True)=={'FINISHED'}
            check('cached camera viewing retains adjustable photo opacity',abs(bg.alpha-.35)<1e-6)
            state['phase'],state['ticks'] = 'opacity_screenshot',0
        elif phase=='opacity_screenshot' and state['ticks']>=8:
            bpy.ops.screen.screenshot(filepath=str(args.artifacts/'photo_transparent.png'))
            start('9007199254740994')
            state['phase'] = 'before_pause'
        elif phase=='before_pause' and transports[-1].photo_pending:
            other = bpy.data.scenes.new('Photo UI paused scene')
            window.scene = other
            state['phase'],state['ticks'] = 'paused',0
        elif phase=='paused' and state['ticks']>=5:
            check('scene switch pauses the actual photo controller before image creation',runtime.active is not None and runtime.active.paused and photo.cached_image('9007199254740994') is None)
            transports[-1].release = True
            state['phase'],state['ticks'] = 'paused_reply',0
        elif phase=='paused_reply' and state['ticks']>=5:
            check('a received photo stays unapplied while paused',photo.cached_image('9007199254740994') is None and state['task'].revision==0)
            window.scene = scene()
            state['phase'] = 'resumed'
        elif phase=='resumed' and runtime.active is None:
            check('returning to the import scene resumes photo display',displayed('9007199254740994'))
            state['images_before_cancel'] = len(bpy.data.images)
            start('9007199254740995')
            state['phase'] = 'before_cancel'
        elif phase=='before_cancel' and transports[-1].photo_pending:
            bpy.ops.blosm.cancel_progressive()
            check('cancel during HTTP wait leaves no image and closes the helper',runtime.active is None and transports[-1].closed
                  and photo.cached_image('9007199254740995') is None and len(bpy.data.images)==state['images_before_cancel'])
            start('9007199254740997')
            state['phase'] = 'before_queue'
        elif phase=='before_queue' and transports[-1].photo_pending:
            state['queue_a_transport'] = transports[-1]
            obj = camera('9007199254740998')
            assert bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT',True)=={'FINISHED'}
            check('a second explicit camera request queues behind the first photo',photo.pending is not None
                  and photo.pending[-1]==obj and runtime.active.adapter.image_id=='9007199254740997')
            state['queue_a_transport'].release = True
            state['phase'] = 'queued_download'
        elif phase=='queued_download' and transports[-1] is not state['queue_a_transport'] and transports[-1].photo_pending:
            check('the first photo is cached without stealing the queued camera selection',
                  photo.cached_image('9007199254740997') is not None
                  and bpy.context.view_layer.objects.active.get('mapillary_image_id')=='9007199254740998'
                  and runtime.active.adapter.image_id=='9007199254740998')
            transports[-1].release = True
            state['phase'] = 'queued_complete'
        elif phase=='queued_complete' and runtime.active is None:
            check('the latest queued photo loads and becomes the camera view',displayed('9007199254740998')
                  and photo.pending is None)
            start('9007199254740996')
            state['phase'] = 'before_disable'
        elif phase=='before_disable' and transports[-1].photo_pending:
            addon.unregister()
            state['registered'] = False
            check('addon disable removes the active timer and helper',runtime.active is None and transports[-1].closed
                  and not bpy.app.timers.is_registered(runtime._tick))
            finish('passed')
            return None
        return .03
    except Exception:
        finish('failed',traceback.format_exc())
        return None


bpy.app.timers.register(tick,first_interval=.3)
