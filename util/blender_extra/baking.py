"""
Pro-only baking utilities for Blosm.

Supports baking multiple materials into a single material with a single
texture. Initially designed for 3D Tiles but structured so additional data
types and bake channels can be added later.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Sequence, Set, Tuple

import bpy
from bpy.types import Material


_REGISTERED = False


def _active_mesh(context) -> Optional[bpy.types.Object]:
    obj = getattr(context, "object", None)
    if obj and getattr(obj, "type", None) == 'MESH':
        return obj
    return None


class BakeChannel(Enum):
    DIFFUSE = "DIFFUSE"
    EMIT = "EMIT"


@dataclass
class BakeSettings:
    image_size: int
    bake_margin: int
    uv_margin: float
    cage_extrusion: float
    replace_active: bool


@dataclass
class SceneBakeState:
    engine: str
    bake_type: str
    use_selected_to_active: bool
    use_pass_direct: bool
    use_pass_indirect: bool
    use_pass_color: bool
    margin: int
    cage_extrusion: float


class BakeJob:
    """Base class encapsulating shared bake steps."""

    channel: BakeChannel = BakeChannel.DIFFUSE

    def __init__(self, context, sources: Sequence[bpy.types.Object], target: bpy.types.Object, settings: BakeSettings, base_name: str):
        self.context = context
        self.sources = list(sources)
        self.target = target
        self.settings = settings
        self.base_name = base_name
        self.material: Optional[bpy.types.Material] = None
        self.image_node: Optional[bpy.types.ShaderNodeTexImage] = None

    # -- public API -----------------------------------------------------

    def execute(self) -> bpy.types.Image:
        self.material, self.image_node = self._build_material()
        self._assign_material(self.material)
        _prepare_uv(self.context, self.target, self.settings.uv_margin)
        state = self._capture_scene_state()
        try:
            self._apply_scene_settings()
            _select_sources_for_bake(self.context, self.sources, self.target)
            bpy.ops.object.bake(type=self._bake_type())
        finally:
            self._restore_scene_state(state)
        self._link_image_after_bake()
        return self.image_node.image

    # -- hooks for subclasses ------------------------------------------

    def _bake_type(self) -> str:
        return self.channel.value

    def _configure_nodes(self, nodes: bpy.types.Nodes, links: bpy.types.NodeLinks) -> Tuple[bpy.types.Material, bpy.types.ShaderNodeTexImage]:
        raise NotImplementedError

    def _apply_scene_settings(self) -> None:
        scene = self.context.scene
        scene.render.engine = "CYCLES"
        scene.cycles.bake_type = self._bake_type()
        scene.render.bake.use_selected_to_active = True
        scene.render.bake.use_pass_direct = False
        scene.render.bake.use_pass_indirect = False
        scene.render.bake.use_pass_color = self.channel == BakeChannel.DIFFUSE
        scene.render.bake.margin = int(self.settings.bake_margin)
        scene.render.bake.cage_extrusion = float(self.settings.cage_extrusion)

    def _post_build(self, material: bpy.types.Material, image_node: bpy.types.ShaderNodeTexImage) -> None:
        pass

    def _post_bake(self) -> None:
        raise NotImplementedError

    # -- helpers --------------------------------------------------------

    def _build_material(self) -> Tuple[bpy.types.Material, bpy.types.ShaderNodeTexImage]:
        mat = bpy.data.materials.new(name=f"{self.base_name} Bake")
        mat.use_nodes = True
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        for node in list(nodes):
            nodes.remove(node)
        image_node = self._configure_nodes(nodes, links)
        suffix = "Bake Diffuse" if self.channel is BakeChannel.DIFFUSE else "Bake Emit"
        image = bpy.data.images.new(
            name=f"{self.base_name} {suffix}",
            width=self.settings.image_size,
            height=self.settings.image_size,
            alpha=True,
        )
        image_node.image = image
        self._post_build(mat, image_node)
        image_node.select = True
        mat.node_tree.nodes.active = image_node
        return mat, image_node

    def _assign_material(self, material: bpy.types.Material) -> None:
        self.target.data.materials.clear()
        self.target.data.materials.append(material)

    def _capture_scene_state(self) -> SceneBakeState:
        scene = self.context.scene
        bake = scene.render.bake
        return SceneBakeState(
            engine=scene.render.engine,
            bake_type=getattr(scene.cycles, "bake_type", "DIFFUSE"),
            use_selected_to_active=bake.use_selected_to_active,
            use_pass_direct=bake.use_pass_direct,
            use_pass_indirect=bake.use_pass_indirect,
            use_pass_color=bake.use_pass_color,
            margin=bake.margin,
            cage_extrusion=bake.cage_extrusion,
        )

    def _restore_scene_state(self, state: SceneBakeState) -> None:
        scene = self.context.scene
        scene.render.engine = state.engine
        scene.cycles.bake_type = state.bake_type
        scene.render.bake.use_selected_to_active = state.use_selected_to_active
        scene.render.bake.use_pass_direct = state.use_pass_direct
        scene.render.bake.use_pass_indirect = state.use_pass_indirect
        scene.render.bake.use_pass_color = state.use_pass_color
        scene.render.bake.margin = state.margin
        scene.render.bake.cage_extrusion = state.cage_extrusion

    def _link_image_after_bake(self) -> None:
        self._post_bake()


class DiffuseBakeJob(BakeJob):
    channel = BakeChannel.DIFFUSE

    def _configure_nodes(self, nodes, links):
        output = nodes.new("ShaderNodeOutputMaterial")
        output.location = (400, 0)
        bsdf = nodes.new("ShaderNodeBsdfPrincipled")
        bsdf.location = (120, 0)
        bsdf.inputs[0].default_value = (1.0, 1.0, 1.0, 1.0)
        image = nodes.new("ShaderNodeTexImage")
        image.location = (-220, 0)
        links.new(bsdf.outputs.get("BSDF"), output.inputs.get("Surface"))
        self._bsdf = bsdf
        self._image = image
        return image

    def _post_bake(self) -> None:
        try:
            self.material.node_tree.links.new(self.image_node.outputs.get("Color"), self._bsdf.inputs.get("Base Color"))
        except Exception:
            pass


class EmitBakeJob(BakeJob):
    channel = BakeChannel.EMIT

    def _configure_nodes(self, nodes, links):
        output = nodes.new("ShaderNodeOutputMaterial")
        output.location = (400, 0)
        emission = nodes.new("ShaderNodeEmission")
        emission.location = (120, 0)
        emission.inputs[0].default_value = (1.0, 1.0, 1.0, 1.0)
        image = nodes.new("ShaderNodeTexImage")
        image.location = (-220, 0)
        links.new(emission.outputs.get("Emission"), output.inputs.get("Surface"))
        self._emission = emission
        self._image = image
        return image

    def _post_bake(self) -> None:
        try:
            self.material.node_tree.links.new(self.image_node.outputs.get("Color"), self._emission.inputs.get("Color"))
        except Exception:
            pass


JOB_REGISTRY = {
    BakeChannel.DIFFUSE: DiffuseBakeJob,
    BakeChannel.EMIT: EmitBakeJob,
}


class MaterialSignature:
    __slots__ = ("has_diffuse", "has_emit")

    def __init__(self, has_diffuse: bool, has_emit: bool):
        self.has_diffuse = has_diffuse
        self.has_emit = has_emit


class BakeChannelSelector:
    def choose_channel(self, materials: Sequence[Material]) -> Optional[BakeChannel]:
        raise NotImplementedError


class TilesAutoSelector(BakeChannelSelector):
    SAMPLE_COUNT = 3

    def choose_channel(self, materials: Sequence[Material]) -> Optional[BakeChannel]:
        pool = [m for m in materials if m]
        if not pool:
            return None
        sample = self._sample_materials(pool)
        if not sample:
            return None
        signatures = [_analyze_material(mat) for mat in sample]
        if signatures and all(sig.has_emit and not sig.has_diffuse for sig in signatures):
            return BakeChannel.EMIT
        if signatures and all(sig.has_diffuse and not sig.has_emit for sig in signatures):
            return BakeChannel.DIFFUSE
        return None

    def _sample_materials(self, materials: Sequence[Material]) -> List[Material]:
        if len(materials) <= self.SAMPLE_COUNT:
            return list(materials)
        seed_input = "".join(sorted(m.name for m in materials if m.name))
        seed = sum(ord(ch) for ch in seed_input)
        rand = random.Random(seed)
        return rand.sample(materials, self.SAMPLE_COUNT)


SELECTOR_REGISTRY = {
    "3d-tiles": TilesAutoSelector(),
}


class BLOSM_OT_BakeToSingleMaterial(bpy.types.Operator):
    bl_idname = "blosm.bake_to_single_material"
    bl_label = "Bake To Single Diffuse"
    bl_description = "Bake multiple materials into a single material/texture"
    bl_options = {"REGISTER", "UNDO"}

    image_size: bpy.props.IntProperty(
        name="Image size",
        description="Square bake resolution (px)",
        min=128,
        max=16384,
        default=4096,
        subtype='PIXEL',
    )
    bake_margin: bpy.props.IntProperty(
        name="Bake margin",
        description="Margin (pixels) around UV islands during bake",
        min=0,
        max=64,
        default=2,
    )
    uv_margin: bpy.props.FloatProperty(
        name="UV pack margin",
        description="UV packing margin (relative 0..1)",
        min=0.0,
        max=0.05,
        default=0.001,
        precision=4,
    )
    cage_extrusion: bpy.props.FloatProperty(
        name="Cage extrusion",
        description="Distance for selected-to-active bake",
        min=0.0,
        max=10.0,
        default=0.01,
        precision=3,
    )
    replace_active: bpy.props.BoolProperty(
        name="Replace active object",
        description="Replace the active object with its baked copy",
        default=False,
    )

    @classmethod
    def poll(cls, context):
        return _active_mesh(context) is not None

    def execute(self, context):
        sources = _gather_sources(context)
        if not sources:
            self.report({'WARNING'}, "No mesh objects selected for baking")
            return {'CANCELLED'}

        channel = self._pick_channel(context, sources)
        if not channel:
            self.report({'ERROR'}, "Unable to determine bake channel for the selected configuration")
            return {'CANCELLED'}

        base_name = sources[0].name
        base_mesh_name = sources[0].data.name if sources[0].data else base_name
        target = _duplicate_sources_for_bake(context, sources, base_name, base_mesh_name)
        settings = _gather_settings(context, self)
        job_cls = JOB_REGISTRY[channel]
        job = job_cls(context, sources, target, settings, base_name)
        baked_image = job.execute()

        _handle_replacement(context, settings.replace_active, sources, target, job.material)
        self.report({'INFO'}, f"Baking finished ({channel.name.title()}): {baked_image.name}")
        return {'FINISHED'}

    def _pick_channel(self, context, sources: Sequence[bpy.types.Object]) -> Optional[BakeChannel]:
        addon = getattr(context.scene, "blosm", None)
        data_type = getattr(addon, "dataType", None)
        selector = SELECTOR_REGISTRY.get(data_type)
        materials = _collect_materials_from_objects(sources)
        if selector:
            return selector.choose_channel(materials)
        # Default to diffuse if no selector is registered yet.
        return BakeChannel.DIFFUSE if materials else None


def _gather_settings(context, operator) -> BakeSettings:
    addon = getattr(context.scene, "blosm", None)
    image_size = int(getattr(addon, "bake_image_size", operator.image_size)) if addon else operator.image_size
    bake_margin = int(getattr(addon, "bake_margin_px", operator.bake_margin)) if addon else operator.bake_margin
    uv_margin = float(getattr(addon, "bake_uv_margin", operator.uv_margin)) if addon else operator.uv_margin
    cage_extrusion = float(getattr(addon, "bake_cage_extrusion", operator.cage_extrusion)) if addon else operator.cage_extrusion
    replace_active = bool(getattr(addon, "bake_replace_active", operator.replace_active)) if addon else operator.replace_active
    return BakeSettings(image_size, bake_margin, uv_margin, cage_extrusion, replace_active)


def _gather_sources(context) -> List[bpy.types.Object]:
    selected = [obj for obj in context.selected_objects if getattr(obj, "type", None) == 'MESH']
    if selected:
        return selected
    active = _active_mesh(context)
    return [active] if active else []


def _duplicate_sources_for_bake(context, sources: Sequence[bpy.types.Object], base_name: str, base_mesh_name: str) -> bpy.types.Object:
    duplicates: List[bpy.types.Object] = []
    for src in sources:
        dup = src.copy()
        dup.data = src.data.copy()
        dup.name = f"{src.name} Bake"
        dup.data.name = f"{src.data.name} Bake"
        if src.users_collection:
            src.users_collection[0].objects.link(dup)
        else:
            context.scene.collection.objects.link(dup)
        duplicates.append(dup)

    if len(duplicates) == 1:
        return duplicates[0]

    bpy.ops.object.select_all(action='DESELECT')
    for dup in duplicates:
        dup.select_set(True)
    context.view_layer.objects.active = duplicates[0]
    bpy.ops.object.join()
    target = duplicates[0]
    target.name = f"{base_name} Bake"
    if target.data:
        target.data.name = f"{base_mesh_name} Bake"
    return target


def _collect_materials_from_objects(objects: Sequence[bpy.types.Object]) -> List[Material]:
    seen: Set[int] = set()
    materials: List[Material] = []
    for obj in objects:
        for slot in getattr(obj, "material_slots", []) or []:
            mat = slot.material
            if mat and mat.as_pointer() not in seen:
                seen.add(mat.as_pointer())
                materials.append(mat)
    return materials


def _analyze_material(material: Optional[Material]) -> MaterialSignature:
    if not material or not material.use_nodes:
        return MaterialSignature(False, False)
    outputs = [node for node in material.node_tree.nodes if node.type == 'OUTPUT_MATERIAL']
    if not outputs:
        return MaterialSignature(False, False)
    visited: Set[int] = set()
    stack_nodes: List[bpy.types.Node] = []
    for node in outputs:
        surface = node.inputs.get('Surface')
        if not surface:
            continue
        for link in surface.links:
            stack_nodes.append(link.from_node)
    has_diffuse = False
    has_emit = False
    while stack_nodes:
        node = stack_nodes.pop()
        pointer = node.as_pointer()
        if pointer in visited:
            continue
        visited.add(pointer)
        if node.type in {'BSDF_PRINCIPLED', 'BSDF_DIFFUSE'}:
            has_diffuse = True
        if node.type == 'EMISSION':
            has_emit = True
        for input_socket in node.inputs:
            for link in input_socket.links:
                stack_nodes.append(link.from_node)
    return MaterialSignature(has_diffuse, has_emit)


def _prepare_uv(context, target: bpy.types.Object, uv_margin: float) -> None:
    with _mode(context, 'OBJECT'):
        _select_only(context, target)
        bpy.ops.object.editmode_toggle()
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.uv.smart_project(angle_limit=1.15192, island_margin=0.0, area_weight=0.0, correct_aspect=True, scale_to_bounds=True)
        bpy.ops.uv.pack_islands(rotate=True, margin=uv_margin)
        bpy.ops.object.editmode_toggle()
        if target.data.uv_layers:
            target.data.uv_layers.active_index = 0
            target.data.uv_layers[0].active_render = True


def _select_sources_for_bake(context, sources: Sequence[bpy.types.Object], target: bpy.types.Object) -> None:
    bpy.ops.object.select_all(action='DESELECT')
    for obj in sources:
        try:
            obj.select_set(True)
        except Exception:
            pass
    target.select_set(True)
    context.view_layer.objects.active = target


def _handle_replacement(context, replace_active: bool, sources: Sequence[bpy.types.Object], target: bpy.types.Object, material: bpy.types.Material) -> None:
    if not replace_active or len(sources) != 1:
        _select_only(context, target)
        return
    source = sources[0]
    source.data = target.data
    source.data.name = f"{source.name}_BakedMesh"
    source.data.materials.clear()
    source.data.materials.append(material)
    _safe_unlink_and_remove(target)
    _select_only(context, source)


def _safe_unlink_and_remove(obj: bpy.types.Object) -> None:
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
    for cls in _classes:
        if not hasattr(bpy.types, cls.__name__):
            bpy.utils.register_class(cls)
    _REGISTERED = True


def register():
    ensure_registered()


def unregister():
    global _REGISTERED
    for cls in reversed(_classes):
        try:
            bpy.utils.unregister_class(cls)
        except Exception:
            pass
    _REGISTERED = False
