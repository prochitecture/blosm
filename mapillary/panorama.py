"""Immersive Mapillary surface. HTTP and photo publication remain finite jobs."""
import math
import bpy
import gpu
from mathutils import Quaternion, Vector
from gpu_extras.batch import batch_for_shader
from .. import viewport
from ..viewport.session import Surface
from ..viewport.ui import Button, Slider, Container, Router, Rect, HANDLED, PASS_THROUGH, shell_layout
from ..viewport.assets import ImageRef
from ..viewport.draw import DISPLAY_RGB_GLSL
from ..progressive import blender as runtime
from ..downloads.diagnostics import redact_text
from . import photo
from .panorama_math import valid_image_size
from .photo_data import MAX_PIXELS, MAX_SIDE


def is_panorama(obj):
    try:
        return (obj is not None and obj.type=='CAMERA'
                and isinstance(obj.get('mapillary_camera_model'),str)
                and obj['mapillary_camera_model'] in {'spherical','equirectangular'})
    except ReferenceError:
        return False


def shader_create():
    info=gpu.types.GPUShaderCreateInfo()
    info.vertex_in(0,'VEC2','pos')
    info.push_constant('MAT4','ray_matrix')
    info.push_constant('FLOAT','opacity')
    info.push_constant('INT','orientation')
    info.push_constant('INT','premultiplied')
    info.push_constant('VEC4','clip_rect')
    info.sampler(0,'FLOAT_2D','photo_texture')
    interface=gpu.types.GPUStageInterfaceInfo('blosm_panorama_interface')
    interface.smooth('VEC2','screen_uv')
    info.vertex_out(interface)
    info.fragment_out(0,'VEC4','fragColor')
    info.vertex_source('void main(){screen_uv=pos;gl_Position=vec4(pos,0.0,1.0);}')
    info.fragment_source(DISPLAY_RGB_GLSL+'''
    void main(){
        if(any(lessThan(gl_FragCoord.xy,clip_rect.xy)) || any(greaterThanEqual(gl_FragCoord.xy,clip_rect.zw))) discard;
        vec3 local=(ray_matrix*vec4(screen_uv,-1.0,1.0)).xyz;
        vec3 d=normalize(vec3(local.x,-local.y,-local.z));
        if(orientation==3) d.xy=-d.xy;
        else if(orientation==6) d.xy=vec2(d.y,-d.x);
        else if(orientation==8) d.xy=vec2(-d.y,d.x);
        float u=fract(0.5+atan(d.x,d.z)/6.28318530718);
        float v=0.5+atan(d.y,length(d.xz))/3.14159265359;
        vec2 uv=vec2(u,v);
        if(orientation==3) uv=1.0-uv;
        else if(orientation==6) uv=vec2(1.0-v,u);
        else if(orientation==8) uv=vec2(v,1.0-u);
        uv.y=1.0-uv.y;
        // Pole clamp prevents filtering the opposite pole. Wrap only longitude.
        vec2 half_texel=0.5/vec2(textureSize(photo_texture,0));
        if(orientation==6 || orientation==8) uv.x=clamp(uv.x,half_texel.x,1.0-half_texel.x);
        else uv.y=clamp(uv.y,half_texel.y,1.0-half_texel.y);
        // Explicit longitude wrap works with Blender-owned sampler settings,
        // including at the seam where a clamped sampler alone would leave a line.
        bool portrait=orientation==6 || orientation==8;
        float n=float(portrait?textureSize(photo_texture,0).y:textureSize(photo_texture,0).x);
        float pixel=(portrait?uv.y:uv.x)*n-0.5;
        float a=(mod(floor(pixel),n)+0.5)/n;
        float b=(mod(floor(pixel)+1.0,n)+0.5)/n;
        vec2 uv_a=portrait?vec2(uv.x,a):vec2(a,uv.y);
        vec2 uv_b=portrait?vec2(uv.x,b):vec2(b,uv.y);
        vec4 c=mix(texture(photo_texture,uv_a),texture(photo_texture,uv_b),fract(pixel));
        if(premultiplied==1 && c.a>0.00001) c.rgb/=c.a;
        fragColor=vec4(display_rgb(c.rgb),c.a*opacity);
    }''')
    return gpu.shader.create_from_info(info)


class Panorama(Surface):
    requires_object_mode=True
    close_label='Exit panorama'

    def __init__(self,obj):
        self.obj=obj
        self.identity=str(obj['mapillary_image_id'])
        self.image_ref=self.shader=self.batch=self.task=self.controller=None
        self.drag=None
        self.status='Loading panorama…'
        self.failed=False
        self.ready=False
        self.opacity=1.
        self.fov=math.radians(60.)
        self.snapshot=None
        self.last_shape=None
        self.image_ticket=None

    @property
    def capturing(self): return self.drag is not None or self.router.capture is not None

    def attach(self,session):
        super().attach(session)
        space=session.space
        if space.lock_camera or space.lock_object or space.region_3d.lock_rotation or space.local_view:
            raise RuntimeError('Unlock the viewport and leave Local View before opening a panorama')
        self.orientation=int(self.obj.get('mapillary_exif_orientation',1))
        if self.orientation not in (1,3,6,8):
            raise RuntimeError('Mirrored EXIF panoramas are not supported')
        if self.orientation!=1 and self.obj.get('mapillary_orientation_source')!='computed_rotation':
            raise RuntimeError('This oriented panorama has no complete rotation metadata')
        self.pose=self.obj.matrix_world.copy()
        if self.pose.to_3x3().determinant()<=0:
            raise RuntimeError('A mirrored or singular camera transform cannot define a panorama')
        self.rotation=self.pose.to_quaternion().normalized()
        self.position=self.pose.translation.copy()
        self.look=self.rotation.copy()
        rv=space.region_3d
        self.snapshot=(rv.view_perspective,rv.view_rotation.copy(),rv.view_location.copy(),rv.view_distance,space.lens)
        self.shader=shader_create()
        self.batch=batch_for_shader(self.shader,'TRI_FAN',{'pos':((-1,-1),(1,-1),(1,1),(-1,1))})
        photo.clear_pending()  # supersede an older queued flat-photo action
        self.close=Button('close','Exit panorama',session.request_close)
        self.slider=Slider('opacity','Photo opacity',lambda:self.opacity,self.set_opacity)
        self.retry=Button('retry','Retry',self.retry_load)
        self.credit=Button('credit','Mapillary · CC-BY-SA',self.open_credit)
        self.root=Container('panorama_controls',[self.close,self.slider,self.retry,self.credit])
        self.router=Router(self.root)
        self.image_ticket=session.ticket(self.identity)
        for bg in self.obj.data.background_images:
            if bg.image and bg.image.get('mapillary_image_id')==self.identity:
                self.opacity=bg.alpha
                break
        # The view stays in place until a usable cached/downloaded image exists.
        self.tick(.004)

    def validate(self):
        try:
            if self.obj.name not in self.session.scene.objects or not is_panorama(self.obj):
                raise RuntimeError('The panorama camera was removed')
            if str(self.obj.get('mapillary_image_id'))!=self.identity or self.obj.matrix_world!=self.pose:
                raise RuntimeError('The panorama camera changed; reopen its viewer')
        except ReferenceError:
            raise RuntimeError('The panorama camera was removed') from None

    def retry_load(self):
        self.failed=False
        self.task=self.controller=None
        self.status='Loading panorama…'

    def open_credit(self):
        bpy.ops.wm.url_open(url='https://www.mapillary.com/app/?pKey='+self.identity)

    def set_opacity(self,value): self.opacity=value

    def cancel_gesture(self):
        self.drag=None
        if hasattr(self,'router'): self.router.reset_capture()

    def tick(self,budget_seconds):
        self.validate()
        if self.image_ref is not None:
            if not self.image_ref.image.has_data: raise RuntimeError('The panorama image was removed')
            return .2
        if self.failed: return .2
        previous=self.status
        if self.task is not None:
            # Controller completion and its undo checkpoint precede publication.
            if self.controller is runtime.active: return .1
            if self.task.job.error or self.task.job.cancelled:
                self.failed=True
                self.status=('Panorama download cancelled' if self.task.job.cancelled else
                             'Panorama download failed: '+(self.task.result_error or 'Retry.'))
            elif self.task.image is not None and self.session.accepts(self.image_ticket):
                self.use_image(self.task.image)
            else:
                self.failed=True
                self.status='The panorama could not be loaded. Retry.'
        else:
            image=photo.cached_image(self.identity)
            if image is not None:
                self.use_image(image)
            elif runtime.is_busy(): self.status='Queued until the current Blosm import finishes…'
            else:
                from ..app.blender import app
                prefs=bpy.context.preferences.addons.get(app.addonName)
                token=photo.get_token(prefs.preferences.mapillaryAccessToken if prefs else '')
                try:
                    self.task=photo.PhotoImport(bpy.context,self.obj,token,activate=False,report_result=False)
                    self.controller=runtime.start(self.task,bpy.context)
                    self.status='Downloading panorama…'
                except Exception as exc:
                    if self.task is not None: self.task.cancel()
                    self.failed=True
                    self.status=redact_text(exc,headers={'Authorization':'OAuth '+token} if token else {})
        if previous!=self.status:
            self.layout(self.session.frame) if self.session.frame else None
            self.session.invalidate()
        return .1 if not self.ready and not self.failed else .2

    def use_image(self,image):
        width,height=image.size
        if max(width,height)>MAX_SIDE or width*height>MAX_PIXELS:
            self.failed=True
            self.status='Panorama exceeds the supported image size'
            return
        if not valid_image_size(*image.size,self.orientation):
            self.failed=True
            self.status='Unsupported panorama dimensions/orientation'
            return
        self.image_ref=ImageRef(image)
        self.status='Drag to look around · Wheel to zoom · Esc to exit'
        self.ready=True
        self.apply_view()

    def apply_view(self):
        space=self.session.space
        rv=space.region_3d
        rv.view_perspective='PERSP'
        rv.view_rotation=self.look
        rv.view_distance=10.
        rv.view_location=self.position+self.look @ Vector((0,0,-10.))
        # Blender's viewport lens uses a 32mm sensor over the larger dimension.
        # Read actual projection below to calibrate its vertical FOV on resize.
        rv.update()
        actual=rv.window_matrix[1][1]
        wanted=1/math.tan(self.fov/2)
        if actual>0: space.lens=max(1.,min(250.,space.lens*wanted/actual))
        rv.update()
        self.session.invalidate()

    def layout(self,frame):
        if frame is None:return
        self.root.rect=frame.available
        self.toolbar,_=shell_layout(frame.available)
        r=self.toolbar
        self.close.rect=Rect(r.x+6,r.y+6,min(126.,max(20.,r.width-12)),30.)
        self.close.label='Exit panorama' if r.width>=145 else 'Exit'
        self.compact_slider=self.ready and 155<=r.width<330 and r.height>=76
        self.slider.visible=self.ready and (r.width>=330 and r.height>=38 or self.compact_slider)
        self.slider.rect=(Rect(r.x+6,r.y+42,r.width-12,30) if self.compact_slider else
                          Rect(r.x+142,r.y+6,min(220.,max(0.,r.width-155)),30.))
        self.retry.visible=self.failed and r.width>=155
        self.retry.rect=Rect(r.x+144,r.y+6,80,30)
        if self.failed and r.width<240:
            self.close.label='Exit'
            self.close.rect=Rect(r.x+6,r.y+6,60,30)
            self.retry.rect=Rect(r.x+76,r.y+6,min(80,r.width-82),30)
        self.credit.visible=frame.available.width>=100 and frame.available.height>145
        self.credit.rect=Rect(frame.available.x+12,frame.available.y+frame.available.height-36,min(205,frame.available.width-24),28)
        shape=(frame.region.width,frame.region.height)
        if self.ready and shape!=self.last_shape:
            self.last_shape=shape
            self.apply_view()

    def handle(self,event):
        # A canvas drag owns its release even when crossing a toolbar button.
        if self.drag is not None and event.type=='LEFTMOUSE' and event.value=='RELEASE':
            self.drag=None
            return HANDLED
        if self.drag is not None and event.type=='MOUSEMOVE':
            x,y=self.drag
            self.drag=(event.x,event.y)
            self.look=(self.look @ Quaternion((0,1,0),-(event.x-x)*.004)
                       @ Quaternion((1,0,0),(event.y-y)*.004)).normalized()
            self.apply_view()
            return HANDLED
        if self.router.handle(event)==HANDLED: return HANDLED
        if event.type=='LEFTMOUSE':
            self.drag=(event.x,event.y) if event.value=='PRESS' and self.ready else None
            return HANDLED
        if event.type in {'WHEELUPMOUSE','WHEELDOWNMOUSE'}:
            self.fov=max(math.radians(20),min(math.radians(110),self.fov*(.9 if event.type=='WHEELUPMOUSE' else 1/.9)))
            if self.ready:self.apply_view()
            return HANDLED
        # Navigation that would move the eye is owned by this surface. Ordinary
        # modelling keys exit first and continue through Blender's keymap.
        if event.type in {'MIDDLEMOUSE','TRACKPADPAN','TRACKPADZOOM','NDOF_MOTION'}: return HANDLED
        if event.value=='PRESS' and event.type in {'G','R','S','TAB','NUMPAD_0','NUMPAD_1','NUMPAD_3','NUMPAD_7'}:
            self.session.request_close()
        return PASS_THROUGH

    def draw(self,frame,painter):
        if self.image_ref is not None:
            rv=self.session.space.region_3d
            # Rotation-only capture inverse prevents scale from distorting rays.
            world_to_capture=self.rotation.inverted().to_matrix().to_4x4()
            view_to_world=rv.view_matrix.inverted().to_3x3().to_4x4()
            matrix=world_to_capture @ view_to_world @ rv.window_matrix.inverted()
            self.shader.bind()
            self.shader.uniform_float('ray_matrix',matrix)
            self.shader.uniform_float('opacity',self.opacity)
            self.shader.uniform_int('orientation',self.orientation)
            self.shader.uniform_int('premultiplied',int(self.image_ref.premultiplied))
            self.shader.uniform_float('clip_rect',painter.pixel_clip())
            self.shader.uniform_sampler('photo_texture',self.image_ref.texture)
            self.batch.draw(self.shader)
        painter.rect(self.toolbar,(.025,.035,.05,.94))
        if self.toolbar.height>56 and not self.compact_slider:
            painter.text(self.status,self.toolbar.x+8,self.toolbar.y+51,max_width=self.toolbar.width-16)
        self.root.draw(painter,self.router)
        if self.image_ref is not None and frame.available.height>180:
            creator=self.image_ref.image.get('mapillary_creator','')
            if creator:
                painter.text('Photo by '+creator,frame.available.x+12,frame.available.y+frame.available.height-57,
                             max_width=frame.available.width-24)

    def detach(self,reason):
        self.cancel_gesture()
        if self.controller is not None and runtime.active is self.controller:
            self.controller.close(silent=True)
        if self.image_ref is not None: self.image_ref.release()
        self.image_ref=self.shader=self.batch=None
        if hasattr(self,'root'): self.root.dispose()
        if reason=='START_FAILED' and self.snapshot and self.session.valid():
            rv=self.session.space.region_3d
            rv.view_perspective,rv.view_rotation,rv.view_location,rv.view_distance,self.session.space.lens=self.snapshot


class BLOSM_OT_ViewMapillaryPanorama(bpy.types.Operator):
    bl_idname='blosm.view_mapillary_panorama'
    bl_label='View Mapillary panorama'
    bl_description='Look around a 360-degree photo; Esc exits while keeping your view'
    @classmethod
    def poll(cls,context):
        return context.mode=='OBJECT' and is_panorama(photo.imported_camera(context))
    def execute(self,context):
        try:
            viewport.open(context,Panorama(photo.imported_camera(context)))
            return {'FINISHED'}
        except Exception as exc:
            self.report({'ERROR'},str(exc))
            return {'CANCELLED'}


def register(): bpy.utils.register_class(BLOSM_OT_ViewMapillaryPanorama)
def unregister(): bpy.utils.unregister_class(BLOSM_OT_ViewMapillaryPanorama)
