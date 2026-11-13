"""Reversible view-layer visibility filters for existing and incoming cameras.

    Draft GUI values become active only on Apply. The applied snapshot lives in
    the .blend, alongside object references for visibility owned by this filter.
    No handlers, downloads, name-based ownership or scene-wide draw-time scans.
"""
from functools import lru_cache
import json

import bpy
from bpy.props import BoolProperty, CollectionProperty, PointerProperty, StringProperty

from .filter_data import CameraFilter, CATEGORIES
from .tiles import identifier


def is_imported_camera(obj):
    try:
        if (obj is not None and obj.type == 'CAMERA'
                and 'mapillary_metadata' in obj and 'mapillary_orientation_source' in obj):
            identifier(obj.get('mapillary_image_id'))
            return True
    except (ValueError, ReferenceError):
        pass
    return False


def draft(settings):
    return dict(date_from=settings.date_from.strip(), date_to=settings.date_to.strip(),
                categories=[key for key in CATEGORIES if getattr(settings, key)],
                include_undated=settings.include_undated)


def encode(options):
    return json.dumps(options, sort_keys=True, separators=(',', ':'))


@lru_cache(maxsize=16)
def compiled(snapshot):
    return CameraFilter.build(**json.loads(snapshot)) if snapshot else CameraFilter.build()


def apply_new_camera(obj, view_layer):
    """One predicate check per camera in the shared import commit path.

    hide_set uses the view layer's object base, not the render visibility or the
    global disable-in-viewports flag. Blender synchronizes newly linked bases
    when accessing them; a full dependency-graph update is not required here.
    """
    settings = getattr(view_layer, 'blosm_mapillary_filters', None)
    if settings is None or not settings.applied:
        return
    spec = compiled(settings.applied)
    if not spec.matches(obj.get('mapillary_captured_at'), obj.get('mapillary_camera_model')):
        if not obj.hide_get(view_layer=view_layer):
            obj.hide_set(True, view_layer=view_layer)
            settings.hidden.add().object = obj


def set_visibility(view_layer, changes):
    """Batch large visibility edits with Blender's native operators.

    hide_set tags the entire layer for resynchronization on every object; a
    Python hide_set loop becomes quadratic. Native hide/clear iterate bases
    once. Snapshot all hidden/selected flags first, restore the requested final
    hidden set, then restore selection. Nothing is drawn in the interim and
    nested operators do not add undo entries. Non-selectable/disabled objects
    use the data API fallback. Local View and headless contexts without a 3D
    area also use the conservative fallback.
    """
    if not changes:
        return
    area = next((a for a in bpy.context.screen.areas if a.type == 'VIEW_3D'
                 and not a.spaces.active.local_view), None) if bpy.context.screen else None
    if len(changes) < 32 or area is None or bpy.context.mode != 'OBJECT':
        for obj, hidden in changes.items():
            obj.hide_set(hidden, view_layer=view_layer)
        return
    objects = list(view_layer.objects)
    selected = [obj for obj in objects if obj.select_get(view_layer=view_layer)]
    desired = {obj: obj.hide_get(view_layer=view_layer) for obj in objects}
    desired.update(changes)
    region = next(r for r in area.regions if r.type == 'WINDOW')
    with bpy.context.temp_override(area=area, region=region, view_layer=view_layer):
        try:
            bpy.ops.object.hide_view_clear('EXEC_DEFAULT', False, select=False)
            for obj in selected:
                obj.select_set(False, view_layer=view_layer)
            for obj, hidden in desired.items():
                if hidden:
                    obj.select_set(True, view_layer=view_layer)
            bpy.ops.object.hide_view_set('EXEC_DEFAULT', False, unselected=False)
            # Hidden collections/disabled objects cannot be selected by the
            # native operator. Set their flags too, for when they are revealed.
            for obj, hidden in desired.items():
                if obj.hide_get(view_layer=view_layer) != hidden:
                    obj.hide_set(hidden, view_layer=view_layer)
        finally:
            for obj in objects:
                if obj.select_get(view_layer=view_layer):
                    obj.select_set(False, view_layer=view_layer)
            for obj in selected:
                obj.select_set(True, view_layer=view_layer)


def apply_existing(view_layer, spec):
    """Only restore flags that we set; retain pre-existing manual hiding.

    Retain references temporarily absent from the layer (excluded collections)
    so a later Apply/Reset can restore them after their collection is included.
    Deleted references are dropped; camera renames require no special handling.
    """
    settings = view_layer.blosm_mapillary_filters
    owned = {item.object.as_pointer(): item.object for item in settings.hidden
             if item.object is not None}
    matches = hidden = 0
    changes = {}
    for obj in view_layer.objects:
        pointer = obj.as_pointer()
        if not is_imported_camera(obj):
            # A camera whose identifying metadata was edited is no longer a
            # target, but still release any visibility flag owned by us.
            if pointer in owned:
                changes[obj] = False
                del owned[pointer]
            continue
        matched = spec.matches(obj.get('mapillary_captured_at'), obj.get('mapillary_camera_model'))
        matches += int(matched)
        if matched:
            if pointer in owned:
                changes[obj] = False
                del owned[pointer]
        else:
            hidden += 1
            if not obj.hide_get(view_layer=view_layer):
                changes[obj] = True
                owned[pointer] = obj
    set_visibility(view_layer, changes)
    settings.hidden.clear()
    for obj in owned.values():
        settings.hidden.add().object = obj
    return matches, hidden


class BLOSM_PG_FilterHiddenCamera(bpy.types.PropertyGroup):
    object: PointerProperty(type=bpy.types.Object)


class BLOSM_PG_MapillaryFilters(bpy.types.PropertyGroup):
    expanded: BoolProperty(name='Camera filters', default=False)
    date_from: StringProperty(name='From', description='Capture date, YYYY-MM-DD (UTC, inclusive); blank for no lower bound', maxlen=10)
    date_to: StringProperty(name='To', description='Capture date, YYYY-MM-DD (UTC, inclusive); blank for no upper bound', maxlen=10)
    include_undated: BoolProperty(name='Include undated', default=False,
        description='Include missing or invalid capture dates when a date range is set')
    ordinary: BoolProperty(name='Ordinary', default=True, description='Perspective photographs')
    panoramic: BoolProperty(name='Panoramic', default=True, description='Spherical / equirectangular panoramas')
    fisheye: BoolProperty(name='Fisheye', default=True, description='Photographs with a fisheye camera model')
    other: BoolProperty(name='Other / unknown', default=True, description='Missing or unrecognized camera models')
    applied: StringProperty(options={'HIDDEN'})
    hidden: CollectionProperty(type=BLOSM_PG_FilterHiddenCamera, options={'HIDDEN'})


class BLOSM_OT_ApplyMapillaryFilters(bpy.types.Operator):
    bl_idname = 'blosm.apply_mapillary_filters'
    bl_label = 'Apply filters'
    bl_description = 'Filter existing and incoming Mapillary cameras in this view layer; all cameras are still imported'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT'

    def execute(self, context):
        settings = context.view_layer.blosm_mapillary_filters
        options = draft(settings)
        try:
            spec = CameraFilter.build(**options)
        except ValueError as error:
            self.report({'ERROR'}, str(error))
            return {'CANCELLED'}
        matched, excluded = apply_existing(context.view_layer, spec)
        settings.applied = encode(options)
        self.report({'INFO'}, '%d cameras match; %d excluded by filters' % (matched, excluded))
        return {'FINISHED'}


class BLOSM_OT_ResetMapillaryFilters(bpy.types.Operator):
    bl_idname = 'blosm.reset_mapillary_filters'
    bl_label = 'Reset filters'
    bl_description = 'Remove camera filters in this view layer; keep cameras that were already manually hidden'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT'

    def execute(self, context):
        settings = context.view_layer.blosm_mapillary_filters
        apply_existing(context.view_layer, CameraFilter.build())
        settings.applied = ''
        settings.date_from = settings.date_to = ''
        settings.include_undated = False
        for key in CATEGORIES:
            setattr(settings, key, True)
        return {'FINISHED'}


def draw(layout, context):
    settings = context.view_layer.blosm_mapillary_filters
    if (context.scene.blosm.dataType != 'mapillary' and not settings.applied
            and not is_imported_camera(context.view_layer.objects.active)):
        return
    box = layout.box()
    box.prop(settings, 'expanded', emboss=False,
             icon='TRIA_DOWN' if settings.expanded else 'TRIA_RIGHT')
    if not settings.expanded:
        return
    box.label(text='Capture date · YYYY-MM-DD (UTC)')
    box.prop(settings, 'date_from')
    box.prop(settings, 'date_to')
    row = box.row()
    row.enabled = bool(settings.date_from.strip() or settings.date_to.strip())
    row.prop(settings, 'include_undated')
    box.label(text='Photo type')
    column = box.column(align=True)
    for key in CATEGORIES:
        column.prop(settings, key)
    row = box.row(align=True)
    row.operator('blosm.apply_mapillary_filters', icon='FILTER')
    row.operator('blosm.reset_mapillary_filters', text='Reset', icon='LOOP_BACK')
    if settings.applied:
        box.label(text=('Changes not applied' if encode(draft(settings)) != settings.applied
                        else 'Applied to this view layer'))
    else:
        box.label(text='Existing and incoming cameras')


CLASSES = (BLOSM_PG_FilterHiddenCamera, BLOSM_PG_MapillaryFilters,
           BLOSM_OT_ApplyMapillaryFilters, BLOSM_OT_ResetMapillaryFilters)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.ViewLayer.blosm_mapillary_filters = PointerProperty(type=BLOSM_PG_MapillaryFilters)


def unregister():
    del bpy.types.ViewLayer.blosm_mapillary_filters
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
    compiled.cache_clear()
