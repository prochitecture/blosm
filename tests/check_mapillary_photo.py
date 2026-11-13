"""Mapillary photo download/cache/view checks in an isolated Blender process.

Use --factory-startup --disable-autoexec and all five chat-specific resource
overrides. This test supplies PNG bytes and Graph metadata in memory; it never
reads credentials or contacts Mapillary. All generated files stay in --artifacts.
"""
import argparse
import importlib
import json
from pathlib import Path
import struct
import sys
import traceback
from types import SimpleNamespace
from unittest.mock import patch
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

addon = importlib.import_module(ROOT.name)
addon.register()
photo = importlib.import_module(ROOT.name+'.mapillary.photo')
runtime = importlib.import_module(ROOT.name+'.progressive.blender')
runtime.show_result = lambda *_args: None
TOKEN = 'fixture-mapillary-photo-secret'
SIGNATURE = 'fixture-private-photo-signature'
checks = []
portrait_jpeg = None


def check(label,condition):
    assert condition,label
    checks.append(label)


def png(width=12,height=8):
    def chunk(kind,data):
        return (struct.pack('>I',len(data))+kind+data+
                struct.pack('>I',zlib.crc32(kind+data)&0xffffffff))
    rows = b''.join(b'\x00'+bytes((40+y*10,120,200,255))*width for y in range(height))
    return (b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',width,height,8,6,0,0,0))+
            chunk(b'IDAT',zlib.compress(rows))+chunk(b'IEND',b''))


def exif_jpeg():
    """Use Blender's own encoder; no Pillow or bundled-Python installation."""
    settings = bpy.context.scene.render.image_settings
    old_format,old_mode = settings.file_format,settings.color_mode
    source = bpy.data.images.new('Photo EXIF fixture',width=64,height=32)
    pixels = []
    for y in range(32):
        for x in range(64):
            # Top-left red, top-right green, bottom-left blue, bottom-right yellow.
            color = ((1.,0.,0.,1.) if x<32 else (0.,1.,0.,1.)) if y>=16 else (
                    (0.,0.,1.,1.) if x<32 else (1.,1.,0.,1.))
            pixels.extend(color)
    source.pixels.foreach_set(pixels)
    source.update()
    path = args.artifacts/'raw_exif_photo.jpg'
    try:
        settings.file_format = 'JPEG'
        settings.color_mode = 'RGB'
        source.save_render(str(path),scene=bpy.context.scene)
        content = path.read_bytes()
    finally:
        settings.file_format = old_format
        settings.color_mode = old_mode
        bpy.data.images.remove(source)
    assert content[:2] == b'\xff\xd8'
    tiff = (b'II'+struct.pack('<HI',42,8)+struct.pack('<H',1)+
            struct.pack('<HHI',0x0112,3,1)+struct.pack('<H',6)+b'\0\0'+b'\0'*4)
    app1 = b'Exif\0\0'+tiff
    return content[:2]+b'\xff\xe1'+struct.pack('>H',len(app1)+2)+app1+content[2:]


def portrait_colors_correct(image):
    """EXIF6 is clockwise: blue/red across the top, yellow/green at bottom."""
    width,height = image.size
    def sample(x,y):
        start = (y*width+x)*4
        return tuple(image.pixels[start:start+3])
    tl,tr = sample(4,height-5),sample(width-5,height-5)
    bl,br = sample(4,4),sample(width-5,4)
    return (tl[2]>max(tl[:2])+.2 and tr[0]>max(tr[1:])+.2
            and min(bl[:2])>bl[2]+.2 and br[1]>max(br[0],br[2])+.2)


def camera(identity):
    data = bpy.data.cameras.new('Photo camera '+identity)
    obj = bpy.data.objects.new('Photo camera '+identity,data)
    bpy.context.scene.collection.objects.link(obj)
    obj['mapillary_image_id'] = identity
    obj['mapillary_metadata'] = '{}'
    obj['mapillary_camera_model'] = 'perspective'
    obj['mapillary_orientation_source'] = 'computed_rotation'
    obj['mapillary_exif_orientation'] = 1
    obj['mapillary_image_url'] = 'https://www.mapillary.com/app/?pKey='+identity
    return obj


def activate(obj):
    for current in tuple(bpy.context.selected_objects):
        current.select_set(False)
    if obj is not None:
        obj.select_set(True)
    bpy.context.view_layer.objects.active = obj


def panel_operators():
    """Run the actual panel draw without depending on screen pixel positions."""
    class Layout:
        def __init__(self):
            self.operators = []
        def operator(self,name,*_args,**_kwargs):
            self.operators.append(name)
            return SimpleNamespace()
        def __getattr__(self,_name):
            return lambda *_args,**_kwargs:self
    layout = Layout()
    gui = importlib.import_module(ROOT.name+'.gui')
    gui.BLOSM_PT_Extent.draw(SimpleNamespace(layout=layout),bpy.context)
    return layout.operators


class MemoryTransport:
    max_requests = 1

    def __init__(self,identity,mode='valid'):
        self.identity = identity
        self.mode = mode
        self.requests = []
        self.pending = None
        self.closed = False

    def submit(self,request):
        self.requests.append(request)
        identity = len(self.requests)
        url = urlsplit(request.url)
        if url.hostname == 'graph.mapillary.com':
            query = parse_qs(url.query)
            assert query['image_ids'] == [self.identity]
            assert 'thumb_2048_url' in query['fields'][0]
            assert request.headers.get('Authorization') == 'OAuth '+TOKEN
            entity = {'id':self.identity}
            if self.mode != 'missing_url':
                entity['thumb_2048_url'] = 'https://scontent.mapillary.com/photo.png?sig='+SIGNATURE
            if self.mode == 'wrong_id':
                entity['id'] = '123'
            content = json.dumps({'data':[entity]}).encode()
            error = RuntimeError('Fixture network failure') if self.mode == 'network' else None
        else:
            assert url.scheme == 'https'
            assert not request.headers.get('Authorization'),'Graph credential sent to photo host'
            content = (b'not a photograph' if self.mode == 'bad_image'
                       else portrait_jpeg if self.mode == 'portrait' else png())
            error = None
        self.pending = (identity,content,error)
        return identity

    def poll(self):
        result,self.pending = self.pending,None
        return result

    def close(self):
        self.closed = True
        self.pending = None


def task(identity,mode='valid'):
    obj = camera(identity)
    activate(obj)
    transport = MemoryTransport(identity,mode)
    return obj,transport,photo.PhotoImport(bpy.context,obj,TOKEN,transport=transport)


def drive(job):
    for _ in range(1000):
        if job.done:
            return
        job.step()
    raise AssertionError('Photo job did not finish')


report = {}
try:
    activate(None)
    check('no active object hides the photo action',photo.imported_camera(bpy.context) is None)
    check('the actual panel omits the action without an imported camera',
          'blosm.view_mapillary_photo' not in panel_operators())
    cube = bpy.data.objects.get('Cube')
    activate(cube)
    check('ordinary active object hides the photo action',photo.imported_camera(bpy.context) is None)
    ordinary = camera('9007199254741000')
    del ordinary['mapillary_image_id']
    activate(ordinary)
    check('ordinary camera hides the photo action',photo.imported_camera(bpy.context) is None)

    identity = '9007199254740993'
    obj,transport,operation = task(identity)
    check('imported camera exposes the photo action',photo.imported_camera(bpy.context) is obj)
    check('the actual panel draws the imported camera photo action',
          'blosm.view_mapillary_photo' in panel_operators())
    window = bpy.context.window
    original_scene = window.scene
    other_scene = bpy.data.scenes.new('Photo wrong-scene fixture')
    try:
        window.scene = other_scene
        check('switching to another scene hides the imported camera action',
              photo.imported_camera(bpy.context) is None
              and 'blosm.view_mapillary_photo' not in panel_operators())
    finally:
        window.scene = original_scene
        bpy.data.scenes.remove(other_scene)
    # Existing user backgrounds survive, and a camera sharing the original
    # datablock must not acquire the selected camera's Mapillary photograph.
    user_image = bpy.data.images.new('User camera background',width=12,height=8)
    obj.data.background_images.new().image = user_image
    twin = bpy.data.objects.new('Other camera sharing data',obj.data)
    bpy.context.scene.collection.objects.link(twin)
    cursor = tuple(bpy.context.scene.cursor.location)
    drive(operation.job)
    check('Graph and photo download complete through the shared job',
          operation.job.error is None and len(transport.requests)==2 and transport.closed)
    image = photo.cached_image(identity)
    check('downloaded image is tagged packed and decoded',image is not None and image.packed_file is not None
          and tuple(image.size)==(12,8))
    check('camera becomes the scene camera with opaque front photo background',
          bpy.context.scene.camera is obj and obj.data.show_background_images
          and any(bg.image is image and bg.alpha==1 and bg.display_depth=='FRONT'
                  and bg.frame_method=='FIT' for bg in obj.data.background_images))
    check('shared camera data is copied without changing the other camera photo',
          obj.data is not twin.data and all(bg.image is not image for bg in twin.data.background_images))
    check('an existing user background is retained',
          any(bg.image is user_image for bg in obj.data.background_images))
    render = bpy.context.scene.render
    check('camera frame matches the photo aspect without stretching',
          abs((render.resolution_x*render.pixel_aspect_x)/(render.resolution_y*render.pixel_aspect_y)-1.5)<1e-6)
    check('photo activation preserves the cursor',tuple(bpy.context.scene.cursor.location)==cursor)
    check('HTTP credentials and signed photo URL are absent from scene metadata',
          all(TOKEN not in json.dumps(dict(block.items()),default=str) and
              SIGNATURE not in json.dumps(dict(block.items()),default=str)
              for block in (obj,obj.data,image,bpy.context.scene)))

    image.name = 'User renamed the photo'
    obj.name = 'User renamed the camera'
    before_images = len(bpy.data.images)
    before_backgrounds = len(obj.data.background_images)
    obj.data.use_fake_user = True
    camera_data = obj.data
    background = next(bg for bg in obj.data.background_images if bg.image is image)
    background.alpha = .35
    photo.show_photo(bpy.context,obj,image)
    check('renamed image is reused without duplicating background entries',
          photo.cached_image(identity) is image and len(bpy.data.images)==before_images
          and len(obj.data.background_images)==before_backgrounds)
    check('a camera fake user does not create duplicate camera data on reopening',obj.data is camera_data)
    activate(obj)
    result = bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
    check('cached operator works without another download',result=={'FINISHED'}
          and len(transport.requests)==2 and photo.cached_image(identity) is image)
    check('photo opacity is retained by direct and cached activation',abs(background.alpha-.35)<1e-6)

    saved = args.artifacts/'packed_photo.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(saved),check_existing=False)
    bpy.ops.wm.open_mainfile(filepath=str(saved))
    image = photo.cached_image(identity)
    check('packed photo survives saving and reopening without its download file',
          image is not None and image.packed_file is not None and tuple(image.size)==(12,8))
    activate(next(o for o in bpy.context.scene.objects if o.get('mapillary_image_id')==identity))
    check('saved cached image remains usable by the operator',
          bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')=={'FINISHED'})
    check('photo opacity survives saving reopening and cached activation',
          any(bg.image is image and abs(bg.alpha-.35)<1e-6
              for bg in bpy.context.scene.camera.data.background_images))

    portrait_jpeg = exif_jpeg()
    portrait_id = str(int(identity)+30)
    portrait_obj,portrait_transport,portrait = task(portrait_id,'portrait')
    portrait_obj['mapillary_exif_orientation'] = 6
    drive(portrait.job)
    portrait_image = photo.cached_image(portrait_id)
    check('embedded JPEG orientation is normalized exactly once',portrait.job.error is None
          and portrait_image is not None and tuple(portrait_image.size)==(32,64)
          and portrait_colors_correct(portrait_image))
    bpy.ops.wm.save_as_mainfile(filepath=str(args.artifacts/'packed_portrait.blend'),check_existing=False)
    bpy.ops.wm.open_mainfile(filepath=str(args.artifacts/'packed_portrait.blend'))
    portrait_image = photo.cached_image(portrait_id)
    check('normalized portrait pixels survive packing and reopening',portrait_image is not None
          and tuple(portrait_image.size)==(32,64) and portrait_colors_correct(portrait_image))
    upright_id = str(int(identity)+31)
    upright_obj,upright_transport,upright = task(upright_id)
    upright_obj['mapillary_exif_orientation'] = 6
    drive(upright.job)
    check('camera Graph EXIF is not applied again to upright downloaded pixels',
          upright.job.error is None and tuple(photo.cached_image(upright_id).size)==(12,8))

    for offset,mode in enumerate(('missing_url','wrong_id','network','bad_image'),10):
        bad_id = str(int(identity)+offset)
        before_images = set(bpy.data.images.keys())
        bad_obj,bad_transport,bad = task(bad_id,mode)
        drive(bad.job)
        check(mode+' fails without retaining a partial image',bad.job.error is not None
              and photo.cached_image(bad_id) is None and set(bpy.data.images.keys())==before_images
              and bad_transport.closed)

    stale_id = str(int(identity)+60)
    stale_obj,stale_transport,stale = task(stale_id)
    while not stale_transport.requests:
        stale.job.step(max_steps=1)
    new_active = camera(str(int(identity)+61))
    activate(new_active)
    bpy.context.scene.camera = new_active
    drive(stale.job)
    check('a completed old photo is cached without stealing the newly selected camera view',
          stale.job.error is None and photo.cached_image(stale_id) is not None
          and bpy.context.view_layer.objects.active is new_active
          and bpy.context.scene.camera is new_active
          and not stale_obj.data.show_background_images)

    queued_one = camera(str(int(identity)+50))
    queued_two = camera(str(int(identity)+51))
    with patch.object(photo.runtime,'is_busy',return_value=True) as busy, \
            patch.object(photo,'get_token',side_effect=AssertionError('Queued photo read a credential')) as token_read, \
            patch.object(photo,'DownloadTransport',side_effect=AssertionError('Queued photo started another helper')) as helper, \
            patch.object(photo,'start_download',return_value=True) as start_download:
        activate(queued_one)
        first = bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
        activate(queued_two)
        second = bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
        check('busy import queues only the latest selected photo without another helper',
              first==second=={'FINISHED'} and photo.pending[-1] is queued_two
              and bpy.app.timers.is_registered(photo._pending_tick)
              and not start_download.called and not helper.called and not token_read.called)
        check('queued request waits while the shared import is busy',
              photo._pending_tick()==.1 and photo.pending[-1] is queued_two and not start_download.called)
        activate(queued_one)
        check('changing the selected camera clears a stale queued request',
              photo._pending_tick() is None and photo.pending is None
              and not bpy.app.timers.is_registered(photo._pending_tick))

        bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
        check('cancel queued photo removes its timer and request',
              bpy.ops.blosm.cancel_mapillary_photo('EXEC_DEFAULT')=={'FINISHED'}
              and photo.pending is None and not bpy.app.timers.is_registered(photo._pending_tick))

        bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
        busy.return_value = False
        check('queued photo starts the original target after the shared import finishes',
              photo._pending_tick() is None and photo.pending is None
              and start_download.call_count==1 and start_download.call_args.args[1] is queued_one
              and not bpy.app.timers.is_registered(photo._pending_tick))

        busy.return_value = True
        bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
        activate(next(o for o in bpy.context.scene.objects if o.get('mapillary_image_id')==identity))
        cached_result = bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT')
        check('cached photo works during another import and clears the queued download',
              cached_result=={'FINISHED'} and photo.pending is None
              and not bpy.app.timers.is_registered(photo._pending_tick)
              and start_download.call_count==1 and not helper.called and not token_read.called)

    cancelled_id = str(int(identity)+20)
    _,cancelled_transport,cancelled = task(cancelled_id)
    while not cancelled_transport.requests:
        cancelled.job.step(max_steps=1)
    cancelled.cancel()
    check('cancellation closes the transport without caching a photo',
          cancelled.job.cancelled and cancelled_transport.closed and photo.cached_image(cancelled_id) is None)

    deleted_id = str(int(identity)+21)
    deleted_obj,deleted_transport,deleted = task(deleted_id)
    bpy.data.objects.remove(deleted_obj,do_unlink=True)
    drive(deleted.job)
    check('deleting the target camera stops the request safely',deleted.job.error is not None
          and deleted_transport.closed and photo.cached_image(deleted_id) is None)
    report = {'status':'passed','checks':checks,'blender':bpy.app.version_string}
except Exception:
    report = {'status':'failed','checks':checks,'error':traceback.format_exc(),
              'blender':bpy.app.version_string}
    raise
finally:
    (args.artifacts/'photo.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('MAPILLARY_PHOTO',json.dumps(report),flush=True)
    if runtime.active is not None:
        runtime.active.close(silent=True)
    addon.unregister()
