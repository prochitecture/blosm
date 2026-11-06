"""Blender lifecycle shared by progressive data-type adapters.

Adapters expose job, scene, view_layer, revision, status, cancel(), and finish().
Consecutive import revisions share an undo step; user edits split the group.
Undo/file loading/add-on disable cancels outstanding work before Blender
invalidates scene references.
"""
from contextlib import contextmanager
import textwrap
import bpy
from bpy.app.handlers import persistent
from .transport import reap

active = None


@contextmanager
def preserve_context():
    context = bpy.context
    layer = context.view_layer
    selected = list(context.selected_objects)
    current = layer.objects.active
    collection = layer.active_layer_collection
    cursor = context.scene.cursor.location.copy()
    try:
        yield
    finally:
        for obj in list(context.selected_objects):
            obj.select_set(False)
        for obj in selected:
            try:
                if obj.name in layer.objects:
                    obj.select_set(True)
            except ReferenceError:
                pass
        try:
            layer.objects.active = current
        except (ReferenceError, RuntimeError):
            pass
        try:
            layer.active_layer_collection = collection
        except (ReferenceError, RuntimeError):
            layer.active_layer_collection = layer.layer_collection
        context.scene.cursor.location = cursor


def notify(level, message):
    try:
        bpy.ops.blosm.progressive_report(level=level, message=message)
    except RuntimeError:
        # ERROR reports make bpy.ops raise after displaying the report.
        pass


def show_result(level, message):
    """Timer-driven operator reports alone do not reliably show a UI banner."""
    def draw(menu, context):
        for paragraph in message.splitlines():
            for line in textwrap.wrap(paragraph, width=80) or ['']:
                menu.layout.label(text=line)
    bpy.context.window_manager.popup_menu(
        draw, title='Blosm import', icon='ERROR' if level == 'ERROR' else 'INFO')


def notify_result(level, message):
    notify(level, message)
    if not bpy.app.background and bpy.context.window is not None:
        show_result(level, message)


def redraw():
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == 'VIEW_3D':
                area.tag_redraw()


def start(adapter, context):
    global active
    if active is not None:
        raise RuntimeError('A Blosm import is already running')
    active = Controller(adapter, context.window)
    try:
        active.checkpoint(initial=True)
        bpy.app.timers.register(_tick, first_interval=0.05)
    except Exception:
        active.close(silent=True)
        raise
    redraw()
    return active


class Controller:
    def __init__(self, adapter, window):
        self.adapter = adapter
        self.window = window
        self.paused = False

    @property
    def status(self):
        return ('Paused: return to the import scene in Object Mode' if self.paused
                else self.adapter.status)

    def tick(self):
        adapter = self.adapter
        if self.window not in list(bpy.context.window_manager.windows):
            self.close(silent=True)
            return
        if adapter.scene.name not in bpy.data.scenes:
            self.close(silent=True)
            return
        self.paused = (self.window.scene != adapter.scene or
                       self.window.view_layer != adapter.view_layer or
                       adapter.scene != bpy.context.scene or bpy.context.mode != 'OBJECT')
        if self.paused:
            adapter.job.poll()
            redraw()
            return
        revision = adapter.revision
        with bpy.context.temp_override(window=self.window, scene=adapter.scene,
                                       view_layer=adapter.view_layer):
            adapter.job.step()
            if adapter.revision != revision:
                self.checkpoint()
        if adapter.job.done:
            self.close()
        redraw()

    def checkpoint(self, initial=False):
        if not bpy.app.background and bpy.context.preferences.edit.use_global_undo:
            if initial:
                # Keep a baseline and separate consecutive imports, including
                # when start() is called without the outer import operator.
                bpy.ops.ed.undo_push(message='Before Blosm progressive import')
            else:
                # Explicitly enable undo for this Python operator call.
                bpy.ops.blosm.progressive_checkpoint('EXEC_DEFAULT', True)

    def close(self, silent=False):
        global active
        if active is not self:
            return
        active = None
        try:
            if not self.adapter.job.done:
                self.adapter.cancel()
            if self.window in list(bpy.context.window_manager.windows):
                with bpy.context.temp_override(window=self.window):
                    self.adapter.finish(silent=silent)
            else:
                self.adapter.finish(silent=True)
        finally:
            if bpy.app.timers.is_registered(_tick):
                bpy.app.timers.unregister(_tick)
            if reap() and not bpy.app.timers.is_registered(_reap):
                bpy.app.timers.register(_reap, first_interval=0.1)
            redraw()


def _tick():
    if active is None:
        return None
    try:
        active.tick()
    except Exception as exc:
        controller = active
        if controller is not None:
            controller.adapter.job.error = exc
            controller.close()
    return 0.03 if active is not None else None


def _reap():
    return 0.1 if reap() else None


@persistent
def stop_for_scene_change(*_):
    if active is not None:
        active.close(silent=True)


class BLOSM_OT_CancelProgressive(bpy.types.Operator):
    bl_idname = 'blosm.cancel_progressive'
    bl_label = 'Cancel import'
    bl_description = 'Stop loading and keep the tiles already imported'

    def execute(self, context):
        if active is not None:
            active.close()
        return {'FINISHED'}


class BLOSM_OT_ProgressiveCheckpoint(bpy.types.Operator):
    bl_idname = 'blosm.progressive_checkpoint'
    bl_label = 'Blosm progressive import'
    bl_options = {'INTERNAL', 'UNDO_GROUPED'}

    @classmethod
    def poll(cls, context):
        return active is not None

    def execute(self, context):
        # Blender replaces only a consecutive checkpoint with this label.
        # A user edit has its own label and therefore keeps its undo boundary.
        return {'FINISHED'}


class BLOSM_OT_ProgressiveReport(bpy.types.Operator):
    bl_idname = 'blosm.progressive_report'
    bl_label = 'Blosm import'
    bl_options = {'INTERNAL'}
    level: bpy.props.StringProperty(default='INFO')
    message: bpy.props.StringProperty()

    def execute(self, context):
        self.report({self.level}, self.message)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(BLOSM_OT_CancelProgressive)
    bpy.utils.register_class(BLOSM_OT_ProgressiveCheckpoint)
    bpy.utils.register_class(BLOSM_OT_ProgressiveReport)
    for handlers in (bpy.app.handlers.load_pre, bpy.app.handlers.undo_pre, bpy.app.handlers.redo_pre):
        handlers.append(stop_for_scene_change)


def unregister():
    stop_for_scene_change()
    for handlers in (bpy.app.handlers.load_pre, bpy.app.handlers.undo_pre, bpy.app.handlers.redo_pre):
        if stop_for_scene_change in handlers:
            handlers.remove(stop_for_scene_change)
    bpy.utils.unregister_class(BLOSM_OT_ProgressiveReport)
    bpy.utils.unregister_class(BLOSM_OT_ProgressiveCheckpoint)
    bpy.utils.unregister_class(BLOSM_OT_CancelProgressive)
