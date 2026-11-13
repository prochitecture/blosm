"""Manual Python Update and authored Section properties; no edit observers."""
import bpy
import bmesh
from .scene import controller_for, NetworkScene, register_runtime, unregister_runtime
from .terrain.blender import MODE_ITEMS


def _import_active():
    from ....progressive.blender import is_busy
    from ....app.blender import app
    session = getattr(app, 'osm_session', None)
    return is_busy() or (session is not None and not session.finished)


def _require_idle():
    if _import_active():
        raise RuntimeError('Finish or cancel the active import before editing streets')


class BLOSM_PG_StreetTerrain(bpy.types.PropertyGroup):
    target: bpy.props.PointerProperty(name='Original terrain', type=bpy.types.Object,
        poll=lambda self,obj: obj.type == 'MESH' and not obj.get('st_py_controller'),
        description='Unmodified terrain surface; empty keeps the prepared street elevation')
    mode: bpy.props.EnumProperty(name='Terrain placement', items=MODE_ITEMS, default='PROFILE')
    clearance: bpy.props.FloatProperty(name='Clearance (m)',default=.02,min=0.)
    maximum_grade: bpy.props.FloatProperty(name='Maximum reference grade',default=.10,min=.001,max=1.,subtype='FACTOR',
        description='Rise/run of the added longitudinal baseline; existing crown relief is preserved')
    control_spacing: bpy.props.FloatProperty(name='Reference spacing (m)',default=4.,min=.1)
    sample_spacing: bpy.props.FloatProperty(name='Terrain sample spacing (m)',default=1.,min=.05)
    maximum_lift: bpy.props.FloatProperty(name='Maximum lift above ground (m)',default=25.,min=.01)


class BLOSM_OT_UpdatePythonStreets(bpy.types.Operator):
    bl_idname = 'blosm.update_python_streets'
    bl_label = 'Update'
    bl_description = 'Regenerate changed streets, dependent surfaces and attached furniture'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return not _import_active() and controller_for(context.active_object) is not None

    def execute(self, context):
        controller = controller_for(context.active_object)
        try:
            _require_idle()
            update = NetworkScene.update(controller)
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        count = update.timings.get('published_objects', 0) + update.timings.get('removed_objects', 0)
        material_count = update.timings.get('material_objects_updated', 0)
        self.report({'INFO'}, f'Updated {count} street objects; {material_count} material assignments'
                    if count or material_count else 'Already up to date')
        return {'FINISHED'}


class BLOSM_OT_LoadStreetMaterials(bpy.types.Operator):
    bl_idname = 'blosm.load_street_materials'
    bl_label = 'Load asphalt materials'
    bl_description = 'Load the three PBR and three diffuse-only asphalt materials from the street material library'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return not _import_active()

    def execute(self, context):
        from pathlib import Path
        from ....app.asset_library import BlenderAssetLoader
        from .. import assets
        settings = getattr(context.scene, 'blosm', None)
        external = getattr(settings, 'streetMaterialFile', '')
        try:
            _require_idle()
            # Explicit user library browsing embeds the chosen materials. Import
            # generation retains its normal asset storage policy.
            loader = BlenderAssetLoader(Path(__file__).resolve().parents[3] / 'assets',
                                        assets.material_catalog(external))
            loader.load_many(assets.ASPHALTS)
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, 'Loaded six asphalt materials')
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

    @classmethod
    def poll(cls, context):
        return not _import_active() and controller_for(context.active_object) is not None

    def invoke(self, context, event):
        try:
            _require_idle()
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
            _require_idle()
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
        if _import_active():
            self.layout.label(text='Finish or cancel the import to edit streets', icon='INFO')
            self.layout.enabled = False
        self.layout.operator('blosm.update_python_streets', icon='FILE_REFRESH')
        self.layout.operator('blosm.python_street_section')
        self.layout.prop(controller, 'blosm_street_pavement')
        self.layout.operator('blosm.load_street_materials')
        for key, label in (('st_py_cell_size', 'Cell size (m)'), ('st_py_resolution', 'Resolution'),
                           ('st_py_exaggeration', 'Exaggeration')):
            self.layout.prop(controller, '["' + key + '"]', text=label)
        terrain = controller.blosm_street_terrain
        box = self.layout.box()
        box.prop(terrain, 'target')
        if terrain.target:
            box.prop(terrain, 'mode')
            box.prop(terrain, 'clearance')
            box.prop(terrain, 'sample_spacing')
            if terrain.mode == 'PROFILE':
                box.prop(terrain, 'maximum_grade')
                box.prop(terrain, 'maximum_lift')
                box.prop(terrain, 'control_spacing')
                box.label(text='Junctions and complete Bundles move together')
            else:
                box.label(text='Terrain cross-slopes affect road profiles')
            box.label(text='Press Update to apply settings and validate')
        self.layout.label(text='Edit source controls, then press Update')
        if controller.get('st_py_error'):
            self.layout.label(text='Output is outdated: see error', icon='ERROR')
            self.layout.label(text=controller['st_py_error'][:70])


CLASSES = (BLOSM_PG_StreetTerrain, BLOSM_OT_UpdatePythonStreets, BLOSM_OT_LoadStreetMaterials,
           BLOSM_OT_PythonStreetSection, BLOSM_PT_PythonStreets)


def register():
    for cls in CLASSES:
        if not cls.is_registered:
            bpy.utils.register_class(cls)
    if not hasattr(bpy.types.Object, 'blosm_street_terrain'):
        bpy.types.Object.blosm_street_terrain = bpy.props.PointerProperty(type=BLOSM_PG_StreetTerrain)
    if not hasattr(bpy.types.Object, 'blosm_street_pavement'):
        bpy.types.Object.blosm_street_pavement = bpy.props.PointerProperty(name='Carriageway material',
            type=bpy.types.Material, description='Optional material override; press Update to apply to this network')
    register_runtime()


def unregister():
    if hasattr(bpy.types.Object, 'blosm_street_pavement'):
        del bpy.types.Object.blosm_street_pavement
    if hasattr(bpy.types.Object, 'blosm_street_terrain'):
        del bpy.types.Object.blosm_street_terrain
    for cls in reversed(CLASSES):
        if cls.is_registered:
            bpy.utils.unregister_class(cls)
    unregister_runtime()
