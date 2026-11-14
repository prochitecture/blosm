"""Temporary per-view-layer hiding with save/undo-safe ownership.

The RNA records travel with undo snapshots, so undoing into an old viewer state
can restore visibility even after that viewer has closed. Normal saves suspend
the hiding and omit its records; the live session resumes after saving.
"""
import uuid
import bpy
from bpy.app.handlers import persistent
from bpy.props import BoolProperty, CollectionProperty, PointerProperty

_leases={}


def key(obj,layer):return (obj.as_pointer(),layer.as_pointer())


def hidden_get(obj,layer):
    """Underlying visibility, excluding an active temporary hide."""
    lease=_leases.get(key(obj,layer))
    return lease.hidden if lease is not None else obj.hide_get(view_layer=layer)


def underlying_changes(layer,changes):
    """Update base visibility while keeping temporarily hidden objects hidden.

    Used by native/batched visibility clients (currently camera filters). They
    must use hidden_get for ownership tests and this before changing flags.
    """
    result=dict(changes)
    for obj,hidden in changes.items():
        lease=_leases.get(key(obj,layer))
        if lease is not None:
            lease.hidden=bool(hidden)
            if not lease.suspended:lease.record().hidden=lease.hidden
            result[obj]=lease.hidden if lease.suspended else True
    return result


def restore(obj,layer,hidden,selected):
    if obj is None:return
    try:
        obj.hide_set(hidden,view_layer=layer)
        if not hidden and selected and layer.objects.active==obj:
            obj.select_set(True,view_layer=layer)
    except (ReferenceError,RuntimeError):
        # Deleted/unlinked objects and removed view layers have no flag to own.
        pass


class TemporaryHide:
    def __init__(self,obj,layer):
        self.obj,self.layer=obj,layer
        self.key=key(obj,layer)
        if self.key in _leases:raise RuntimeError('This object already has a temporary visibility owner')
        self.token=uuid.uuid4().hex
        self.hidden=obj.hide_get(view_layer=layer)
        self.selected=obj.select_get(view_layer=layer)
        self.closed=self.suspended=False
        _leases[self.key]=self
        try:self.resume()
        except Exception:
            self.release()
            raise

    def record(self):
        records=self.layer.blosm_temporary_hidden
        item=records.get(self.token)
        if item is None:
            item=records.add();item.name=self.token
        item.object=self.obj
        item.hidden,item.selected=self.hidden,self.selected
        return item

    def remove_record(self):
        records=self.layer.blosm_temporary_hidden
        index=records.find(self.token)
        if index>=0:records.remove(index)

    def enforce(self):
        if not self.closed and not self.suspended and not self.obj.hide_get(view_layer=self.layer):
            self.obj.hide_set(True,view_layer=self.layer)

    def resume(self):
        if self.closed:return
        self.record()
        self.obj.hide_set(True,view_layer=self.layer)
        self.suspended=False

    def suspend(self):
        if self.closed:return
        restore(self.obj,self.layer,self.hidden,self.selected)
        self.remove_record()
        self.suspended=True

    def release(self):
        if self.closed:return
        self.closed=True
        _leases.pop(self.key,None)
        restore(self.obj,self.layer,self.hidden,self.selected)
        try:self.remove_record()
        except (ReferenceError,RuntimeError):pass


class BLOSM_PG_TemporaryHidden(bpy.types.PropertyGroup):
    object:PointerProperty(type=bpy.types.Object)
    hidden:BoolProperty()
    selected:BoolProperty()


@persistent
def recover(*_):
    """Release stale records restored by load/Undo/Redo, not live sessions."""
    for scene in bpy.data.scenes:
        for layer in scene.view_layers:
            records=layer.blosm_temporary_hidden
            for index in range(len(records)-1,-1,-1):
                item=records[index]
                if item.object is not None:
                    lease=_leases.get(key(item.object,layer))
                    if lease is not None and lease.token==item.name:continue
                    restore(item.object,layer,item.hidden,item.selected)
                records.remove(index)


def _resume_timer():
    for lease in tuple(_leases.values()):
        if lease.suspended:
            try:lease.resume()
            except (ReferenceError,RuntimeError):lease.release()
    return None


@persistent
def save_pre(*_):
    recover()
    for lease in tuple(_leases.values()):
        try:lease.suspend()
        except (ReferenceError,RuntimeError):lease.release()
    # Also recover after a failed save on Blender versions without save_post_fail.
    if _leases and not bpy.app.timers.is_registered(_resume_timer):
        bpy.app.timers.register(_resume_timer,first_interval=0.)


@persistent
def save_post(*_):
    _resume_timer()
    if bpy.app.timers.is_registered(_resume_timer):bpy.app.timers.unregister(_resume_timer)


def handlers():
    return ((bpy.app.handlers.load_post,recover),(bpy.app.handlers.undo_post,recover),
            (bpy.app.handlers.redo_post,recover),(bpy.app.handlers.save_pre,save_pre),
            (bpy.app.handlers.save_post,save_post)) + (
            ((bpy.app.handlers.save_post_fail,save_post),) if hasattr(bpy.app.handlers,'save_post_fail') else ())


def register():
    bpy.utils.register_class(BLOSM_PG_TemporaryHidden)
    bpy.types.ViewLayer.blosm_temporary_hidden=CollectionProperty(type=BLOSM_PG_TemporaryHidden,options={'HIDDEN'})
    for group,callback in handlers():group.append(callback)
    recover()


def unregister():
    for lease in tuple(_leases.values()):lease.release()
    recover()
    for group,callback in handlers():
        if callback in group:group.remove(callback)
    if bpy.app.timers.is_registered(_resume_timer):bpy.app.timers.unregister(_resume_timer)
    del bpy.types.ViewLayer.blosm_temporary_hidden
    bpy.utils.unregister_class(BLOSM_PG_TemporaryHidden)
