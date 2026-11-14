"""Reversible view-layer visibility filters for existing and incoming cameras.

    Draft GUI values become active only on Apply. The applied snapshot lives in
    the .blend, alongside object references for visibility owned by this filter.
    No downloads, name-based ownership or scene-wide draw-time scans. Calendar
    migration uses a load handler and one deferred registration/appended-scene pass.
"""
from functools import lru_cache
from calendar import monthrange
from datetime import date, datetime, timezone
import json

import bpy
from bpy.app.handlers import persistent
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, IntProperty, PointerProperty, StringProperty

from .filter_data import ANY_YEAR, FIRST_YEAR, CameraFilter, CATEGORIES, month_date
from .tiles import identifier
from ..viewport.visibility import hidden_get, underlying_changes

MONTHS=tuple((f'{i:02}',label,'',i) for i,label in enumerate((
    'January','February','March','April','May','June',
    'July','August','September','October','November','December'),1))


def utc_year():
    return datetime.now(timezone.utc).year


@lru_cache(maxsize=None)
def years_for(year):
    # Blender borrows callback strings. Retain every generated list for this
    # process, including when a long-running Blender crosses New Year's Day.
    return [(ANY_YEAR,'No limit','No date boundary',0)]+[
        (str(value),str(value),'',value) for value in range(FIRST_YEAR,year+1)]


def year_items(_settings,_context):
    return years_for(utc_year())


def day_count(settings,side):
    year=getattr(settings,'year_'+side)
    # Hidden day choices retain their value when No limit is selected.
    if year==ANY_YEAR or not year.isdecimal():return 31
    return monthrange(int(year),int(getattr(settings,'month_'+side)))[1]


@lru_cache(maxsize=None)
def days_for(count):
    return tuple((f'{value:02}',str(value),'',value) for value in range(1,count+1))


def from_day_items(settings,_context):return days_for(day_count(settings,'from'))
def to_day_items(settings,_context):return days_for(day_count(settings,'to'))


def calendar_changed(settings,_context):
    settings.calendar_version=2
    if settings.adjusting_calendar:return
    settings.adjusting_calendar=True
    try:
        for side,default in (('from',1),('to',31)):
            if getattr(settings,'year_'+side)==ANY_YEAR:continue
            maximum=day_count(settings,side)
            # Read the stored enum number before resolving its identifier: the
            # previous month may have supplied a day absent from the new list.
            raw=settings.get('day_'+side,default)
            chosen=maximum if side=='to' and settings.to_month_end else max(1,min(int(raw),maximum))
            setattr(settings,'day_'+side,f'{chosen:02}')
    finally:settings.adjusting_calendar=False


def from_day_changed(settings,context):
    calendar_changed(settings,context)


def to_day_changed(settings,context):
    if not settings.adjusting_calendar:settings.to_month_end=False
    calendar_changed(settings,context)


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
    return dict(date_from=month_date(settings.year_from,settings.month_from,day=settings.day_from),
                date_to=month_date(settings.year_to,settings.month_to,day=settings.day_to,end=True),
                categories=[key for key in CATEGORIES if getattr(settings, key)],
                include_undated=settings.include_undated)


def migrate_settings(settings):
    """Populate calendar drafts without changing any applied visibility/filter.

    Applied snapshots keep their exact predicate until the next explicit Apply.
    Legacy draft bounds outside the new year list remain represented by
    No limit; the applied saved range is still retained and labelled in the UI.
    """
    if settings.calendar_version>=2:
        return
    if settings.calendar_version==0:
        settings.adjusting_calendar=True
        try:
            for side,legacy in (('from',settings.date_from),('to',settings.date_to)):
                text=legacy.strip()
                try:bound=date.fromisoformat(text)
                except ValueError:continue
                if bound.isoformat()==text and FIRST_YEAR<=bound.year<=utc_year():
                    setattr(settings,'year_'+side,str(bound.year))
                    setattr(settings,'month_'+side,f'{bound.month:02}')
                    setattr(settings,'day_'+side,f'{bound.day:02}')
                    if side=='to':settings.to_month_end=False
        finally:settings.adjusting_calendar=False
    calendar_changed(settings,None)
    snapshot=json.loads(settings.applied) if settings.applied else {}
    current=draft(settings)
    settings.legacy_applied=bool(snapshot and any(current[key]!=snapshot.get(key,'')
        for key in ('date_from','date_to')))
    settings.calendar_version=2


@persistent
def migrate_all(*_):
    for scene in bpy.data.scenes:
        for layer in scene.view_layers:
            migrate_settings(layer.blosm_mapillary_filters)


def migrate_once():
    migrate_all()
    return None


def schedule_migration():
    if not bpy.app.timers.is_registered(migrate_once):
        bpy.app.timers.register(migrate_once,first_interval=0.)


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
    changes=underlying_changes(view_layer,changes)
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
            if not hidden_get(obj,view_layer):
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
    year_from: EnumProperty(name='From',description='First capture year; No limit leaves the lower boundary open',
        items=year_items,default=0,update=calendar_changed)
    year_to: EnumProperty(name='To',description='Last capture year; No limit leaves the upper boundary open',
        items=year_items,default=0,update=calendar_changed)
    month_from: EnumProperty(name='From month',description='Month of the first inclusive capture date (UTC)',
        items=MONTHS,default='01',update=calendar_changed)
    month_to: EnumProperty(name='To month',description='Month of the last inclusive capture date (UTC)',
        items=MONTHS,default='12',update=calendar_changed)
    day_from: EnumProperty(name='From day',description='First capture day, inclusive (UTC)',
        items=from_day_items,default=1,update=from_day_changed)
    day_to: EnumProperty(name='To day',description='Last capture day, inclusive (UTC)',
        items=to_day_items,default=31,update=to_day_changed)
    to_month_end: BoolProperty(default=True,options={'HIDDEN'})
    adjusting_calendar: BoolProperty(default=False,options={'HIDDEN','SKIP_SAVE'})
    # Retained only for migrating .blend files from the previous text-field UI.
    date_from: StringProperty(options={'HIDDEN'})
    date_to: StringProperty(options={'HIDDEN'})
    calendar_version: IntProperty(default=0,options={'HIDDEN'})
    legacy_applied: BoolProperty(default=False,options={'HIDDEN'})
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
        try:
            options = draft(settings)
            spec = CameraFilter.build(**options)
        except ValueError as error:
            self.report({'ERROR'}, str(error))
            return {'CANCELLED'}
        matched, excluded = apply_existing(context.view_layer, spec)
        settings.applied = encode(options)
        settings.calendar_version=2
        settings.legacy_applied=False
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
        settings.year_from=settings.year_to=ANY_YEAR
        settings.to_month_end=True
        settings.month_from,settings.month_to='01','12'
        settings.adjusting_calendar=True
        try:settings.day_from,settings.day_to='01','31'
        finally:settings.adjusting_calendar=False
        settings.legacy_applied=False
        settings.calendar_version=2
        settings.include_undated = False
        for key in CATEGORIES:
            setattr(settings, key, True)
        return {'FINISHED'}


def draw(layout, context):
    settings = context.view_layer.blosm_mapillary_filters
    if settings.calendar_version<2 and (settings.date_from or settings.date_to or settings.applied):
        # Covers old scenes appended after registration without mutating IDs
        # inside panel drawing. Normal file loads migrate via load_post.
        schedule_migration()
    if (context.scene.blosm.dataType != 'mapillary' and not settings.applied
            and not is_imported_camera(context.view_layer.objects.active)):
        return
    box = layout.box()
    box.prop(settings, 'expanded', emboss=False,
             icon='TRIA_DOWN' if settings.expanded else 'TRIA_RIGHT')
    if not settings.expanded:
        return
    box.label(text='Capture date (UTC)')
    for side,label in (('from','From'),('to','To')):
        row=box.split(factor=.17,align=True)
        row.label(text=label+':')
        if getattr(settings,'year_'+side)!=ANY_YEAR:
            fields=row.split(factor=.31,align=True)
            fields.prop(settings,'year_'+side,text='')
            rest=fields.split(factor=.69,align=True)
            rest.prop(settings,'month_'+side,text='')
            rest.prop(settings,'day_'+side,text='')
        else:
            row.prop(settings,'year_'+side,text='')
    row = box.row()
    row.enabled = settings.year_from!=ANY_YEAR or settings.year_to!=ANY_YEAR
    row.prop(settings, 'include_undated')
    box.label(text='Photo type')
    column = box.column(align=True)
    for key in CATEGORIES:
        column.prop(settings, key)
    row = box.row(align=True)
    row.operator('blosm.apply_mapillary_filters', icon='FILTER')
    row.operator('blosm.reset_mapillary_filters', text='Reset', icon='LOOP_BACK')
    if settings.applied:
        box.label(text=('Saved date range active' if settings.legacy_applied else
                        'Changes not applied' if encode(draft(settings)) != settings.applied
                        else 'Applied to this view layer'))
        if settings.legacy_applied:box.label(text='Apply to use these dates')
    else:
        box.label(text='Existing and incoming cameras')


CLASSES = (BLOSM_PG_FilterHiddenCamera, BLOSM_PG_MapillaryFilters,
           BLOSM_OT_ApplyMapillaryFilters, BLOSM_OT_ResetMapillaryFilters)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.ViewLayer.blosm_mapillary_filters = PointerProperty(type=BLOSM_PG_MapillaryFilters)
    bpy.app.handlers.load_post.append(migrate_all)
    schedule_migration()


def unregister():
    if migrate_all in bpy.app.handlers.load_post:bpy.app.handlers.load_post.remove(migrate_all)
    if bpy.app.timers.is_registered(migrate_once):bpy.app.timers.unregister(migrate_once)
    del bpy.types.ViewLayer.blosm_mapillary_filters
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
    compiled.cache_clear()
