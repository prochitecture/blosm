"""On-demand, packed Mapillary photos through the shared import controller."""
import hashlib
from pathlib import Path
import tempfile
import time
from urllib.parse import urlencode

import bpy
from bpy.app.handlers import persistent

from ..downloads import Request
from ..downloads.diagnostics import redact_text
from ..downloads.transport import DownloadTransport
from ..progressive import blender as runtime
from ..progressive.core import Job
from .credentials import get_token
from .tiles import identifier
from .photo_data import photo_url, image_info

MAX_BYTES = 32 * 1024 * 1024
pending = None


def imported_camera(context):
    obj = context.view_layer.objects.active
    try:
        if (obj is not None and obj.type == 'CAMERA' and obj.name in context.scene.objects
                and 'mapillary_metadata' in obj and 'mapillary_orientation_source' in obj):
            identifier(obj.get('mapillary_image_id'))
            return obj
    except (ValueError, ReferenceError):
        pass
    return None


def cached_image(image_id):
    image_id = identifier(image_id)
    for image in bpy.data.images:
        if (image.get('mapillary_photo_source') == 'Mapillary'
                and image.get('mapillary_image_id') == image_id
                and (image.has_data or image.packed_file is not None)
                and min(image.size) > 0):
            return image
    return None


def view_area(context):
    if context.area is not None and context.area.type == 'VIEW_3D':
        return context.area
    if context.window is not None:
        return next((area for area in context.window.screen.areas if area.type == 'VIEW_3D'),None)
    return None


def show_photo(context, obj, image, area=None):
    """Set the requested render camera and fill its frame without stretching."""
    area = area or view_area(context)
    if obj.name not in context.scene.objects or obj.type != 'CAMERA':
        raise RuntimeError('The Mapillary camera was deleted or unlinked')
    if area is None or area.type != 'VIEW_3D':
        raise RuntimeError('Open a 3D View to view the Mapillary photo')
    width, height = image.size
    if not width or not height:
        raise RuntimeError('The Mapillary photo has no image data')
    # Backgrounds are camera-data settings. Copy shared data so a second camera
    # using the same datablock does not acquire this camera's photograph.
    if obj.data.users-int(obj.data.use_fake_user) > 1:
        obj.data = obj.data.copy()
    camera = obj.data
    background = next((bg for bg in camera.background_images
                       if bg.image and bg.image.get('mapillary_photo_source') == 'Mapillary'),None)
    if background is None:
        background = camera.background_images.new()
        background.alpha = 1.
    background.image = image
    background.show_background_image = True
    background.display_depth = 'FRONT'
    background.frame_method = 'FIT'
    background.offset = (0.,0.)
    background.rotation = 0.
    background.scale = 1.
    background.use_flip_x = background.use_flip_y = False
    camera.show_background_images = True
    scene = context.scene
    scene.camera = obj
    scale = max(1,(4+min(width,height)-1)//min(width,height))
    scene.render.resolution_x, scene.render.resolution_y = width*scale, height*scale
    scene.render.pixel_aspect_x = scene.render.pixel_aspect_y = 1.
    context.view_layer.objects.active = obj
    obj.select_set(True)
    space = area.spaces.active
    space.use_local_camera = False
    space.overlay.show_overlays = True
    for region in space.region_quadviews or (space.region_3d,):
        region.view_perspective = 'CAMERA'
        region.view_camera_offset = (0.,0.)
        region.view_camera_zoom = 0.
    area.tag_redraw()


def load_photo(content, image_id, creator):
    width, height, orientation = image_info(content)
    image = None
    try:
        # Blender loads files, not byte buffers. Keep the transient file in its
        # configured temporary directory, pack it, then remove the directory.
        with tempfile.TemporaryDirectory(prefix='blosm-mapillary-',dir=bpy.app.tempdir) as folder:
            path = Path(folder)/'photo'
            path.write_bytes(content)
            image = bpy.data.images.load(str(path),check_existing=False)
            if tuple(image.size) != (width,height):
                raise ValueError('Mapillary photo dimensions do not match its header')
            if orientation != 1:
                # Blender JPEG decoding ignores EXIF. Normalize the downloaded
                # JPEG's own orientation once, never the original API EXIF tag.
                import numpy as np
                pixels = np.empty(width*height*4,dtype=np.float32)
                image.pixels.foreach_get(pixels)
                top = pixels.reshape(height,width,4)[::-1]
                transforms = {2:lambda a:a[:,::-1],3:lambda a:a[::-1,::-1],
                              4:lambda a:a[::-1],5:lambda a:a.transpose(1,0,2),
                              6:lambda a:np.rot90(a,-1),
                              7:lambda a:a.transpose(1,0,2)[::-1,::-1],
                              8:lambda a:np.rot90(a,1)}
                upright = transforms[orientation](top)
                corrected = bpy.data.images.new('Mapillary photo '+image_id,
                                                width=upright.shape[1],height=upright.shape[0])
                try:
                    corrected.colorspace_settings.name = image.colorspace_settings.name
                    corrected.alpha_mode = image.alpha_mode
                    corrected.pixels.foreach_set(np.ascontiguousarray(upright[::-1]).ravel())
                    corrected.pack()
                except Exception:
                    bpy.data.images.remove(corrected)
                    raise
                bpy.data.images.remove(image)
                image = corrected
            else:
                image.pack()
        image.name = 'Mapillary photo '+image_id
        image.filepath_raw = ''
        image.use_fake_user = True
        image['mapillary_photo_source'] = 'Mapillary'
        image['mapillary_image_id'] = image_id
        image['mapillary_creator'] = creator
        image['mapillary_license'] = 'CC-BY-SA'
        return image
    except Exception:
        if image is not None and image.name in bpy.data.images:
            bpy.data.images.remove(image)
        raise


class PhotoImport:
    def __init__(self,context,obj,token,*,transport=None,activate=True,report_result=True):
        if not token or any(c.isspace() for c in token):
            raise ValueError('Set a valid Mapillary client access token in Blosm preferences')
        self.scene, self.view_layer = context.scene, context.view_layer
        self.obj, self.area, self.window = obj, view_area(context), context.window
        self.image_id = identifier(obj.get('mapillary_image_id'))
        self.headers = {'Authorization':'OAuth '+token}
        self.auth_scope = hashlib.sha256(token.encode()).hexdigest()
        self.revision = 0
        self.phase = 'Downloading photo metadata…'
        self.image = None
        self.result_error = None
        self.activate,self.report_result = activate,report_result
        self.transport = transport if transport is not None else DownloadTransport(
            max_requests=1,max_bytes=MAX_BYTES,max_buffer_bytes=MAX_BYTES)
        self.job = Job(self.steps(),self.transport,validate=self.validate)

    def validate(self):
        try:
            if (self.obj.type == 'CAMERA' and self.scene.objects.get(self.obj.name) == self.obj
                    and identifier(self.obj.get('mapillary_image_id')) == self.image_id):
                return
        except (ReferenceError,ValueError):
            pass
        raise RuntimeError('The Mapillary camera was deleted, unlinked or changed')

    def steps(self):
        query = urlencode({'image_ids':self.image_id,'fields':'id,thumb_2048_url,creator'})
        content = yield Request('https://graph.mapillary.com/images?'+query,
                                headers=self.headers,auth_scope=self.auth_scope,decode_gzip=True)
        url, creator = photo_url(content,self.image_id)
        content = None
        self.phase = 'Downloading photo…'
        # Signed CDN URLs carry their own authorization. Never forward OAuth.
        content = yield Request(url,decode_gzip=True)
        self.validate()
        image = load_photo(content,self.image_id,redact_text(creator,headers=self.headers))
        try:
            if self.window is None or self.window not in list(bpy.context.window_manager.windows):
                raise RuntimeError('The photo view window was closed')
            with bpy.context.temp_override(window=self.window,scene=self.scene,view_layer=self.view_layer):
                # A download must not steal the view back after the user has
                # selected a different camera (possibly queuing its photo).
                if self.activate and self.view_layer.objects.active == self.obj:
                    show_photo(bpy.context,self.obj,image,self.area)
        except Exception:
            bpy.data.images.remove(image)
            raise
        self.image = image
        self.revision += 1

    @property
    def status(self):
        return 'Mapillary: '+self.phase

    def cancel(self):
        self.job.cancel()

    def finish(self,silent=False):
        self.result_error = redact_text(self.job.error,headers=self.headers) if self.job.error else None
        if not silent and self.report_result:
            if self.job.error:
                runtime.notify_result('ERROR','Mapillary photo failed: '+redact_text(self.job.error,headers=self.headers))
            elif self.job.cancelled:
                runtime.notify_result('INFO','Mapillary photo download cancelled')
            else:
                runtime.notify_result('INFO','Mapillary photo loaded')
        self.headers.clear()

    def run_blocking(self):
        while not self.job.done:
            self.job.step()
            if self.job.waiting:
                time.sleep(.001)
        self.finish()
        return self.job.error is None


def start_download(context,obj):
    from ..app.blender import app
    prefs = context.preferences.addons.get(app.addonName)
    token = get_token(prefs.preferences.mapillaryAccessToken if prefs else '')
    task = None
    try:
        task = PhotoImport(context,obj,token)
        if bpy.app.background:
            return task.run_blocking()
        runtime.start(task,context)
        return True
    except Exception as error:
        if task is not None:
            task.cancel()
        raise RuntimeError(redact_text(error,headers={'Authorization':'OAuth '+token} if token else {})) from None


@persistent
def clear_pending(*_args):
    global pending
    pending = None
    if bpy.app.timers.is_registered(_pending_tick):
        bpy.app.timers.unregister(_pending_tick)


def _pending_tick():
    item = pending
    if item is None:
        return None
    window,scene,layer,obj = item
    try:
        if (window not in list(bpy.context.window_manager.windows)
                or scene.objects.get(obj.name) != obj or layer.objects.active != obj):
            clear_pending()
            runtime.redraw()
            return None
        if runtime.is_busy() or window.scene != scene or window.view_layer != layer or bpy.context.mode != 'OBJECT':
            return .1
        clear_pending()
        with bpy.context.temp_override(window=window,scene=scene,view_layer=layer):
            # Cached activation is a scene change too: preserve the operator's
            # ordinary undo step instead of mutating the camera from this timer.
            bpy.ops.blosm.view_mapillary_photo('EXEC_DEFAULT',True)
    except Exception as error:
        clear_pending()
        runtime.notify_result('ERROR','Mapillary photo failed: '+str(error))
    runtime.redraw()
    return None


class BLOSM_OT_ViewMapillaryPhoto(bpy.types.Operator):
    bl_idname = 'blosm.view_mapillary_photo'
    bl_label = 'View Mapillary photo'
    bl_description = 'View this Mapillary camera with its photo; download once and reuse the packed image'
    bl_options = {'REGISTER','UNDO'}

    @classmethod
    def poll(cls,context):
        return imported_camera(context) is not None and context.mode == 'OBJECT'

    def execute(self,context):
        global pending
        from .. import viewport
        session=viewport.current()
        if session is not None: session.request_close('REPLACED')
        obj = imported_camera(context)
        try:
            image = cached_image(obj['mapillary_image_id'])
            if image is not None:
                clear_pending()
                show_photo(context,obj,image)
                return {'FINISHED'}
            if runtime.is_busy():
                clear_pending()
                pending = (context.window,context.scene,context.view_layer,obj)
                bpy.app.timers.register(_pending_tick,first_interval=.1)
                self.report({'INFO'},'Photo queued until the current Blosm import finishes or is cancelled')
                runtime.redraw()
                return {'FINISHED'}
            return {'FINISHED'} if start_download(context,obj) else {'CANCELLED'}
        except Exception as error:
            self.report({'ERROR'},str(error))
            return {'CANCELLED'}


class BLOSM_OT_CancelMapillaryPhoto(bpy.types.Operator):
    bl_idname = 'blosm.cancel_mapillary_photo'
    bl_label = 'Cancel queued photo'

    def execute(self,context):
        clear_pending()
        runtime.redraw()
        return {'FINISHED'}


def register():
    bpy.utils.register_class(BLOSM_OT_ViewMapillaryPhoto)
    bpy.utils.register_class(BLOSM_OT_CancelMapillaryPhoto)
    for handlers in (bpy.app.handlers.load_pre,bpy.app.handlers.undo_pre,bpy.app.handlers.redo_pre):
        handlers.append(clear_pending)


def unregister():
    clear_pending()
    for handlers in (bpy.app.handlers.load_pre,bpy.app.handlers.undo_pre,bpy.app.handlers.redo_pre):
        if clear_pending in handlers:
            handlers.remove(clear_pending)
    bpy.utils.unregister_class(BLOSM_OT_CancelMapillaryPhoto)
    bpy.utils.unregister_class(BLOSM_OT_ViewMapillaryPhoto)
