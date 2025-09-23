"""
Pro-only baking utilities for Blosm: bake multiple materials into a single
material with one diffuse texture. Designed for 3D Tiles now, reusable for OSM later.
Shown in UI only for premium builds when this module is present.
"""

import bpy

_REGISTERED = False


def _active_mesh(context):
    ob = context.object
    return ob if ob and getattr(ob, 'type', None) == 'MESH' else None


class BLOSM_OT_BakeToSingleMaterial(bpy.types.Operator):
    bl_idname = "blosm.bake_to_single_material"
    bl_label = "Bake To Single Diffuse"
    bl_description = "Bake all materials of the active mesh into a single material with one diffuse texture"
    bl_options = {"REGISTER", "UNDO"}

    image_size: bpy.props.IntProperty(
        name="Image size",
        description="Square bake resolution (px)",
        min=128, max=16384, default=4096, step=1,
        subtype='PIXEL'
    )
    bake_margin: bpy.props.IntProperty(
        name="Bake margin",
        description="Margin (pixels) around UV islands during bake",
        min=0, max=64, default=2
    )
    uv_margin: bpy.props.FloatProperty(
        name="UV pack margin",
        description="UV packing margin (relative 0..1)",
        min=0.0, max=0.05, default=0.001, precision=4
    )
    cage_extrusion: bpy.props.FloatProperty(
        name="Cage extrusion",
        description="Distance for selected-to-active bake to catch surface details",
        min=0.0, max=1.0, default=0.01, precision=3
    )
    replace_active: bpy.props.BoolProperty(
        name="Replace active object",
        description="Replace the active object with its baked copy",
        default=False
    )

    @classmethod
    def poll(cls, context):
        return _active_mesh(context) is not None

    def execute(self, context):
        # Capture current selection as potential sources
        pre_selected_sources = [o for o in getattr(context, 'selected_objects', []) if getattr(o, 'type', None) == 'MESH']
        src = _active_mesh(context)
        if not src:
            self.report({'ERROR'}, "Active mesh object required")
            return {'CANCELLED'}

        if not getattr(src.data, 'polygons', None) or len(src.data.polygons) == 0:
            self.report({'ERROR'}, "Active mesh has no faces to bake")
            return {'CANCELLED'}

        # Duplicate as bake target (receiver). If multiple sources are selected,
        # create duplicates for all and join them into a single target object.
        sources = pre_selected_sources or [src]
        dup_objects: list[bpy.types.Object] = []
        for s in sources:
            d = s.copy()
            d.data = s.data.copy()
            d.name = f"{s.name}_BAKE_PART"
            if s.users_collection:
                s.users_collection[0].objects.link(d)
            else:
                context.scene.collection.objects.link(d)
            dup_objects.append(d)

        # If more than one duplicate, join into a single target
        if len(dup_objects) > 1:
            # Join requires selection
            bpy.ops.object.select_all(action='DESELECT')
            for d in dup_objects:
                d.select_set(True)
            context.view_layer.objects.active = dup_objects[0]
            bpy.ops.object.join()
            tgt = dup_objects[0]
            tgt.name = f"Baked_Joined"
        else:
            tgt = dup_objects[0]
            tgt.name = f"{src.name}_BAKE"

        # Read settings from addon properties if present (UI-driven pre-run settings)
        addon = getattr(getattr(context, 'scene', None), 'blosm', None)
        img_size = int(getattr(addon, 'bake_image_size', self.image_size)) if addon else self.image_size
        bake_margin = int(getattr(addon, 'bake_margin_px', self.bake_margin)) if addon else self.bake_margin
        uv_margin = float(getattr(addon, 'bake_uv_margin', self.uv_margin)) if addon else self.uv_margin
        cage_extrusion = float(getattr(addon, 'bake_cage_extrusion', self.cage_extrusion)) if addon else self.cage_extrusion
        replace_active = bool(getattr(addon, 'bake_replace_active', self.replace_active)) if addon else self.replace_active

        # Create a single new material with an image texture
        mat = bpy.data.materials.new(name=f"{src.name}_Baked")
        mat.use_nodes = True
        nt = mat.node_tree
        nodes, links = nt.nodes, nt.links
        for n in list(nodes):
            nodes.remove(n)
        out = nodes.new("ShaderNodeOutputMaterial")
        out.location = (400, 0)
        bsdf = nodes.new("ShaderNodeBsdfPrincipled")
        bsdf.location = (120, 0)
        img_node = nodes.new("ShaderNodeTexImage")
        img_node.location = (-220, 0)
        size = max(128, min(16384, img_size))
        img = bpy.data.images.new(name=f"{src.name}_Baked_diffuse", width=size, height=size, alpha=True)
        img_node.image = img
        # Avoid circular dependency during bake: do NOT connect image before baking
        links.new(bsdf.outputs.get('BSDF'), out.inputs.get('Surface'))
        # Make the image node active for baking
        img_node.select = True
        nt.nodes.active = img_node

        # Assign the single material to the target object
        tgt.data.materials.clear()
        tgt.data.materials.append(mat)

        # Prepare UVs on target
        with _mode(context, 'OBJECT'):
            _select_only(context, tgt)
            bpy.ops.object.editmode_toggle()
            bpy.ops.mesh.select_all(action='SELECT')
            # Blender 4.x smart project arg names
            bpy.ops.uv.smart_project(angle_limit=66, island_margin=0.0, area_weight=0.0, correct_aspect=True, scale_to_bounds=True)
            bpy.ops.uv.pack_islands(rotate=True, margin=uv_margin)
            bpy.ops.object.editmode_toggle()
            if tgt.data.uv_layers:
                tgt.data.uv_layers.active_index = 0
                tgt.data.uv_layers[0].active_render = True

        # Store current engine and bake settings to restore later
        scene = context.scene
        prev_engine = scene.render.engine
        prev_bake = {
            'type': getattr(scene.cycles, 'bake_type', 'DIFFUSE'),
            'use_selected_to_active': scene.render.bake.use_selected_to_active,
            'use_pass_direct': scene.render.bake.use_pass_direct,
            'use_pass_indirect': scene.render.bake.use_pass_indirect,
            'use_pass_color': scene.render.bake.use_pass_color,
        }

        try:
            scene.render.engine = 'CYCLES'
            scene.cycles.bake_type = 'DIFFUSE'
            scene.render.bake.use_selected_to_active = True
            # Only color is needed (no direct/indirect lighting)
            scene.render.bake.use_pass_direct = False
            scene.render.bake.use_pass_indirect = False
            scene.render.bake.use_pass_color = True

            # Perform bake: select sources (all selected meshes) then target active
            bpy.ops.object.select_all(action='DESELECT')
            for o in sources:
                try:
                    o.select_set(True)
                except Exception:
                    pass
            tgt.select_set(True)
            context.view_layer.objects.active = tgt
            # Ensure the correct image node is active on target
            mat.node_tree.nodes.active = img_node

            bpy.ops.object.bake(type='DIFFUSE')

            # After bake: replace original or leave a copy
            # If multiple sources, keep baked copy to avoid replacing a single source only
            replace_effective = replace_active and len(sources) == 1

            if replace_effective:
                src.data = tgt.data
                src.data.name = f"{src.name}_BakedMesh"
                src.data.materials.clear()
                src.data.materials.append(mat)
                # Link baked image for display
                try:
                    mat.node_tree.links.new(img_node.outputs.get('Color'), bsdf.inputs.get('Base Color'))
                except Exception:
                    pass
                _safe_unlink_and_remove(tgt)
            else:
                _select_only(context, tgt)
                try:
                    mat.node_tree.links.new(img_node.outputs.get('Color'), bsdf.inputs.get('Base Color'))
                except Exception:
                    pass

        finally:
            scene.render.engine = prev_engine
            scene.cycles.bake_type = prev_bake['type']
            scene.render.bake.use_selected_to_active = prev_bake['use_selected_to_active']
            scene.render.bake.use_pass_direct = prev_bake['use_pass_direct']
            scene.render.bake.use_pass_indirect = prev_bake['use_pass_indirect']
            scene.render.bake.use_pass_color = prev_bake['use_pass_color']

        self.report({'INFO'}, "Baking finished")
        return {'FINISHED'}


def _safe_unlink_and_remove(obj: bpy.types.Object):
    try:
        for coll in list(obj.users_collection) if obj.users_collection else []:
            coll.objects.unlink(obj)
        bpy.data.objects.remove(obj, do_unlink=True)
    except Exception:
        pass


class _mode:
    def __init__(self, context, mode):
        self.context = context
        self.mode = mode
        self.prev = context.mode
    def __enter__(self):
        if self.context.mode != self.mode:
            bpy.ops.object.mode_set(mode=self.mode)
    def __exit__(self, exc_type, exc, tb):
        if self.context.mode != self.prev:
            bpy.ops.object.mode_set(mode=self.prev)


def _select_only(context, obj):
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    context.view_layer.objects.active = obj


_classes = (
    BLOSM_OT_BakeToSingleMaterial,
)


def ensure_registered():
    global _REGISTERED
    if _REGISTERED:
        return
    for c in _classes:
        if not hasattr(bpy.types, c.__name__):
            bpy.utils.register_class(c)
    _REGISTERED = True


def register():
    ensure_registered()


def unregister():
    global _REGISTERED
    for c in reversed(_classes):
        try:
            bpy.utils.unregister_class(c)
        except Exception:
            pass
    _REGISTERED = False
