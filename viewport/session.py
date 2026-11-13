"""One owned viewport session, independent of scene-import transactions."""
from contextlib import contextmanager
import itertools
import bpy
from bpy.app.handlers import persistent
from .ui import Frame, Rect, InputEvent, HANDLED, PASS_THROUGH
from .draw import Painter
from .assets import AssetTicket

active = None
_launch = None
_ids = itertools.count(1)


class Surface:
    """Override only the hooks needed by the consumer. No provider dispatch."""
    requires_object_mode = False
    close_label = 'Close view'
    def attach(self,session): self.session = session
    def layout(self,frame): pass
    def handle(self,event): return PASS_THROUGH
    def tick(self,budget_seconds): return .2
    def draw(self,frame,painter): pass
    def detach(self,reason): pass
    def cancel_gesture(self): pass


class Session:
    def __init__(self,context,surface):
        if bpy.app.background or context.window is None or context.area is None or context.area.type!='VIEW_3D':
            raise RuntimeError('Open this viewer from a 3D Viewport')
        self.window,self.screen,self.area = context.window,context.window.screen,context.area
        self.space = self.area.spaces.active
        if self.space.region_quadviews:
            raise RuntimeError('Leave Quad View before opening this viewer')
        self.region = next(r for r in self.area.regions if r.type=='WINDOW')
        self.scene,self.layer = context.scene,context.view_layer
        self.surface,self.id = surface,next(_ids)
        self.generation,self.layout_revision = 0,0
        self.frame,self._layout_key,self.painter,self.handler = None,None,None,None
        self.closed,self.attached,self.error = False,False,None

    @contextmanager
    def context(self):
        with bpy.context.temp_override(window=self.window,area=self.area,region=self.region):
            yield bpy.context

    def valid(self):
        try:
            return (self.window in list(bpy.context.window_manager.windows)
                and self.window.screen==self.screen and self.area in list(self.screen.areas)
                and self.area.type=='VIEW_3D' and self.area.spaces.active==self.space
                and not self.space.region_quadviews and self.region in list(self.area.regions)
                and self.window.scene==self.scene and self.window.view_layer==self.layer
                and (not self.surface.requires_object_mode or self.layer.objects.active is None
                     or self.layer.objects.active.mode=='OBJECT'))
        except (ReferenceError,RuntimeError):
            return False

    def ticket(self,key): return AssetTicket(self.id,self.generation,str(key))

    def accepts(self,ticket):
        return not self.closed and self.valid() and (ticket.session_id,ticket.generation)==(self.id,self.generation)

    def invalidate(self):
        if not self.closed:
            try: self.area.tag_redraw()
            except ReferenceError: pass

    def update_layout(self):
        scale = max(.5,bpy.context.preferences.system.ui_scale)
        r = self.region
        left,bottom,right,top = 0.,0.,float(r.width),float(r.height)
        # In overlap mode Blender's N-panel/tools can cover the WINDOW region.
        for panel in self.area.regions:
            if panel.type in {'HEADER','TOOL_HEADER'} and panel.height>1:
                y=panel.y-r.y
                if y<=1: bottom=max(bottom,y+panel.height)
                else: top=min(top,y)
                continue
            if panel.type not in {'UI','TOOLS'} or panel.width<=1: continue
            x = panel.x-r.x
            if x <= 1: left=max(left,x+panel.width)
            else: right=min(right,x)
        key = (r.width,r.height,left,bottom,right,top,scale)
        if key==self._layout_key: return
        self._layout_key = key
        self.layout_revision += 1
        self.frame = Frame(Rect(0,0,r.width/scale,r.height/scale),
                           Rect(left/scale,bottom/scale,max(0,right-left)/scale,max(0,top-bottom)/scale),
                           scale,self.layout_revision)
        self.surface.cancel_gesture()
        self.surface.layout(self.frame)
        self.invalidate()

    def start(self):
        self.painter = Painter()
        self.attached = True  # partial attach also receives rollback
        self.surface.attach(self)
        self.update_layout()
        self.handler = bpy.types.SpaceView3D.draw_handler_add(self.draw,(),'WINDOW','POST_PIXEL')
        self.invalidate()

    def draw(self):
        if self.closed: return
        c=bpy.context
        if c.window!=self.window or c.area!=self.area or c.region!=self.region: return
        if not self.valid(): return
        try:
            # Size-dependent preparation is normally performed on events/ticks.
            # A draw after resize waits for that update; it never loads assets.
            if self.frame is not None:
                with self.painter.state(self.frame): self.surface.draw(self.frame,self.painter)
        except Exception:
            # Avoid rendering/logging the same failure every frame. Main timer
            # performs teardown outside the draw callback.
            self.error = 'The viewport graphics could not be drawn'

    def event(self,context,event):
        if self.closed: return {'FINISHED'}
        if not self.valid():
            self.close('OWNER_CHANGED')
            return {'FINISHED'}
        if event.type=='WINDOW_DEACTIVATE':
            self.surface.cancel_gesture()
            return {'PASS_THROUGH'}
        captured=getattr(self.surface,'capturing',False)
        if context.area!=self.area and not captured: return {'PASS_THROUGH'}
        self.update_layout()
        if event.type=='ESC' and event.value=='PRESS':
            self.close('USER')
            return {'FINISHED'}
        r,s=self.region,self.frame.scale
        e=InputEvent(event.type,event.value,(event.mouse_x-r.x)/s,(event.mouse_y-r.y)/s,
                     event.shift,event.ctrl,event.alt)
        # Sidebar/header events stay with Blender, except a captured drag.
        inside = self.frame.available.contains(e.x,e.y)
        if not inside and not captured: return {'PASS_THROUGH'}
        handled = self.surface.handle(e)==HANDLED
        if handled: self.invalidate()
        return {'RUNNING_MODAL'} if handled else {'PASS_THROUGH'}

    def request_close(self,reason='USER'): self.close(reason)

    def close(self,reason='USER'):
        global active
        if self.closed: return
        self.closed = True
        self.generation += 1
        if active is self: active=None
        try:
            if self.handler is not None:
                bpy.types.SpaceView3D.draw_handler_remove(self.handler,'WINDOW')
                self.handler=None
        finally:
            try:
                if self.attached:
                    self.attached=False
                    self.surface.detach(reason)
            finally:
                if self.painter: self.painter.release()
                self.painter=None
                try: self.area.tag_redraw()
                except ReferenceError: pass


def open(context,surface):
    """Open in the invoking area. Returns the session; raises on failed startup."""
    global active,_launch
    candidate=Session(context,surface)
    if active is not None: active.close('REPLACED')
    active=candidate
    try:
        with candidate.context():
            candidate.start()
            _launch=candidate
            bpy.ops.blosm.viewport_session('INVOKE_DEFAULT')
        if not bpy.app.timers.is_registered(_tick): bpy.app.timers.register(_tick,first_interval=.05)
        return candidate
    except Exception:
        candidate.close('START_FAILED')
        raise
    finally:
        _launch=None


def _tick():
    current=active
    if current is None: return None
    try:
        if not current.valid(): current.close('OWNER_CHANGED')
        elif current.error:
            current.close('ERROR')
            report_error(current.error)
        else:
            with current.context():
                current.update_layout()
                delay=current.surface.tick(.004)
            if not current.closed: return max(.05,min(.25,float(delay or .2)))
    except Exception as exc:
        # Surface contract: provider exceptions must be sanitized before raising.
        current.error=str(exc) or 'The viewport session stopped'
        current.close('ERROR')
        report_error(current.error)
    return .05 if active is not None else None


class BLOSM_OT_ViewportSession(bpy.types.Operator):
    bl_idname='blosm.viewport_session'
    bl_label='Blosm viewport'
    bl_options={'INTERNAL'}
    def invoke(self,context,event):
        if _launch is None: return {'CANCELLED'}
        self.session=_launch
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    def modal(self,context,event):
        try: return self.session.event(context,event)
        except Exception:
            self.session.close('ERROR')
            self.report({'ERROR'},'The viewport interaction stopped')
            return {'CANCELLED'}


class BLOSM_OT_CloseViewport(bpy.types.Operator):
    bl_idname='blosm.close_viewport'
    bl_label='Close Blosm view'
    def execute(self,context):
        if active is not None: active.close('USER')
        return {'FINISHED'}


class BLOSM_OT_ViewportReport(bpy.types.Operator):
    bl_idname='blosm.viewport_report'
    bl_label='Blosm viewport'
    bl_options={'INTERNAL'}
    message:bpy.props.StringProperty()
    def execute(self,context):
        self.report({'WARNING'},self.message)
        return {'FINISHED'}


def report_error(message):
    try:bpy.ops.blosm.viewport_report(message=message)
    except RuntimeError:pass


@persistent
def stop(*_):
    if active is not None: active.close('LIFECYCLE')
    if bpy.app.timers.is_registered(_tick): bpy.app.timers.unregister(_tick)


def register():
    for cls in (BLOSM_OT_ViewportSession,BLOSM_OT_CloseViewport,BLOSM_OT_ViewportReport): bpy.utils.register_class(cls)
    for handlers in (bpy.app.handlers.load_pre,bpy.app.handlers.undo_pre,bpy.app.handlers.redo_pre):
        handlers.append(stop)


def unregister():
    stop()
    for handlers in (bpy.app.handlers.load_pre,bpy.app.handlers.undo_pre,bpy.app.handlers.redo_pre):
        if stop in handlers: handlers.remove(stop)
    for cls in (BLOSM_OT_ViewportReport,BLOSM_OT_CloseViewport,BLOSM_OT_ViewportSession): bpy.utils.unregister_class(cls)
