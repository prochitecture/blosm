"""Manual Python Update and authored Section properties; no edit observers."""
import bpy
import bmesh
from .scene import controller_for, NetworkScene, register_runtime, unregister_runtime


class BLOSM_OT_UpdatePythonStreets(bpy.types.Operator):
    bl_idname = 'blosm.update_python_streets'
    bl_label = 'Update'
    bl_description = 'Regenerate changed streets, dependent surfaces and attached furniture'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return controller_for(context.active_object) is not None

    def execute(self, context):
        controller = controller_for(context.active_object)
        try:
            update = NetworkScene.update(controller)
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        count = update.timings.get('published_objects', 0) + update.timings.get('removed_objects', 0)
        self.report({'INFO'}, f'Updated {count} street objects' if count else 'Already up to date')
        return {'FINISHED'}


def selected_section(obj):
    if obj is None or obj.type != 'MESH' or not obj.get('st_py_controller'):
        raise ValueError('Select a Section in its source mesh')
    temporary = obj.mode != 'EDIT'
    bm = bmesh.new() if temporary else bmesh.from_edit_mesh(obj.data)
    if temporary:
        bm.from_mesh(obj.data)
    layer = bm.edges.layers.int.get('st_section_id')
    ids = {e[layer] for e in bm.edges if e.select or any(v.select for v in e.verts)} if layer else set()
    if len(ids) != 1:
        if temporary:
            bm.free()
        raise ValueError('Select controls from exactly one Section')
    return bm, temporary, ids.pop()


class BLOSM_OT_PythonStreetSection(bpy.types.Operator):
    bl_idname = 'blosm.python_street_section'
    bl_label = 'Section properties'
    bl_description = 'Set authored Section attributes; press Update to regenerate'
    bl_options = {'REGISTER', 'UNDO'}

    width: bpy.props.FloatProperty(name='Width', min=.01, default=6.)
    offset: bpy.props.FloatProperty(name='Offset', default=0.)
    rise: bpy.props.FloatProperty(name='Crown rise', min=0., default=.08)
    interpolation: bpy.props.EnumProperty(name='Interpolation', items=(
        ('0', 'Exact polyline', 'Preserve polyline segments and mitered joins'),
        ('3', 'Smooth interpolation', 'Cubic smoothing with bounded local corner transitions')), default='3')

    def invoke(self, context, event):
        try:
            bm, temporary, identity = selected_section(context.active_object)
            try:
                layer = bm.edges.layers.int['st_section_id']
                edge = next(e for e in bm.edges if e[layer] == identity)
                for name in ('width', 'offset', 'rise'):
                    setattr(self, name, edge[bm.edges.layers.float['st_' + name]])
                self.interpolation = str(edge[bm.edges.layers.int['st_interpolation']])
            finally:
                if temporary:
                    bm.free()
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        obj = context.active_object
        try:
            bm, temporary, identity = selected_section(obj)
            try:
                identity_layer = bm.edges.layers.int['st_section_id']
                for edge in bm.edges:
                    if edge[identity_layer] != identity:
                        continue
                    for name in ('width', 'offset', 'rise'):
                        edge[bm.edges.layers.float['st_' + name]] = getattr(self, name)
                    edge[bm.edges.layers.int['st_interpolation']] = int(self.interpolation)
                if temporary:
                    bm.to_mesh(obj.data)
                    obj.data.update()
                else:
                    bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
            finally:
                if temporary:
                    bm.free()
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}


class BLOSM_PT_PythonStreets(bpy.types.Panel):
    bl_label = 'Python streets'
    bl_idname = 'BLOSM_PT_python_streets'
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Streets'

    @classmethod
    def poll(cls, context):
        return controller_for(context.active_object) is not None

    def draw(self, context):
        controller = controller_for(context.active_object)
        self.layout.operator('blosm.update_python_streets', icon='FILE_REFRESH')
        self.layout.operator('blosm.python_street_section')
        for key, label in (('st_py_cell_size', 'Cell size (m)'), ('st_py_resolution', 'Resolution'),
                           ('st_py_exaggeration', 'Exaggeration')):
            self.layout.prop(controller, '["' + key + '"]', text=label)
        self.layout.label(text='Edit source controls, then press Update')
        if controller.get('st_py_error'):
            self.layout.label(text='Output is outdated: see error', icon='ERROR')
            self.layout.label(text=controller['st_py_error'][:70])


CLASSES = (BLOSM_OT_UpdatePythonStreets, BLOSM_OT_PythonStreetSection, BLOSM_PT_PythonStreets)


def register():
    for cls in CLASSES:
        if not cls.is_registered:
            bpy.utils.register_class(cls)
    register_runtime()


def unregister():
    for cls in reversed(CLASSES):
        if cls.is_registered:
            bpy.utils.unregister_class(cls)
    unregister_runtime()
