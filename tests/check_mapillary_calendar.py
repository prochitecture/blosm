"""Native calendar dropdowns, leap days, exact filtering and file migration.

Run in an isolated Blender UI process with factory startup, disabled autoexec,
event simulation and the five chat-specific user-resource paths. No HTTP/token.
"""
import argparse
import importlib
import json
from pathlib import Path
import sys
import time
import traceback
from types import SimpleNamespace
from unittest.mock import patch

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
addon=importlib.import_module(ROOT.name);addon.register()
filters=importlib.import_module(ROOT.name+'.mapillary.filters')
data=importlib.import_module(ROOT.name+'.mapillary.filter_data')
checks=[]
state={'phase':0,'ticks':0,'start':time.monotonic()}


def check(label,condition):
    assert condition,label
    checks.append(label)


def settings():return bpy.context.view_layer.blosm_mapillary_filters
def apply():return bpy.ops.blosm.apply_mapillary_filters('EXEC_DEFAULT',True)
def reset():return bpy.ops.blosm.reset_mapillary_filters('EXEC_DEFAULT',True)


class Layout:
    """Record the actual draw routine's conditional properties and enablement."""
    def __init__(self,calls=None):self.calls=[] if calls is None else calls;self.enabled=True
    def box(self,**_):return Layout(self.calls)
    def row(self,**_):return Layout(self.calls)
    def column(self,**_):return Layout(self.calls)
    def split(self,**_):return Layout(self.calls)
    def prop(self,settings,name,**_):self.calls.append((name,self.enabled))
    def label(self,**_):pass
    def operator(self,*_,**__):return SimpleNamespace()


def drawn():
    layout=Layout();filters.draw(layout,bpy.context)
    return dict(layout.calls)


def cameras():
    return sorted((o for o in bpy.context.scene.objects if o.get('calendar_fixture')),key=lambda o:o.name)


def visibility():return [not o.hide_get() for o in cameras()]


def fixtures():
    dates=[data.date_bound('2024-02-01','From')-1,data.date_bound('2024-02-01','From'),
           data.date_bound('2024-03-01','To')-1,data.date_bound('2024-03-01','To'),None]
    for index,stamp in enumerate(dates):
        obj=bpy.data.objects.new('Calendar fixture '+str(index),bpy.data.cameras.new('Calendar'))
        bpy.context.scene.collection.objects.link(obj)
        obj['calendar_fixture']=True
        obj['mapillary_image_id']=str(9007199254740993+index)
        obj['mapillary_metadata']='{}';obj['mapillary_orientation_source']='computed_rotation'
        obj['mapillary_camera_model']='perspective';obj['mapillary_captured_at']='' if stamp is None else str(stamp)


def basic():
    s=settings();s.expanded=True;bpy.context.scene.blosm.dataType='mapillary'
    check('default years are unbounded and months/days are January 1 / December 31',
          (s.year_from,s.year_to,s.month_from,s.month_to,s.day_from,s.day_to)==('ANY','ANY','01','12','01','31'))
    ui=drawn()
    check('No limit hides both month/day selectors and disables undated control',
          not any(key in ui for key in ('month_from','month_to','day_from','day_to')) and not ui['include_undated'])
    current=filters.utc_year();items=filters.year_items(s,None)
    check('year list is No limit then every year from 2014 through current UTC year',
          [item[0] for item in items]==['ANY']+[str(y) for y in range(2014,current+1)]
          and items[0][1]=='No limit' and items[-1][3]==current)
    s.year_from='2014'
    with patch.object(filters,'utc_year',return_value=current+1):
        future=filters.year_items(s,None)
        check('New Year automatically extends list without shifting saved enum identity',future[-1][0]==str(current+1)
              and s.year_from=='2014' and future[1][3]==2014)
    check('dynamic enum strings retain their cached lifetime',filters.year_items(s,None) is items)
    category_layer=bpy.context.scene.view_layers.new('Legacy type-only filter')
    category_settings=category_layer.blosm_mapillary_filters
    category_settings.applied=filters.encode(dict(date_from='',date_to='',categories=['ordinary'],include_undated=False))
    filters.migrate_settings(category_settings)
    check('legacy type-only edits do not produce a saved-date warning',not category_settings.legacy_applied
          and json.loads(category_settings.applied)['categories']==['ordinary'])
    s.year_from='2024'
    ui=drawn()
    check('specific From year shows its month/day only',all(key in ui for key in ('month_from','day_from'))
          and 'month_to' not in ui and 'day_to' not in ui and ui['include_undated'])
    s.month_from='02'
    check('February in leap year offers days 1 through 29',[item[0] for item in filters.from_day_items(s,None)]==[f'{i:02}' for i in range(1,30)])
    s.day_from='29';s.year_from='2023'
    check('changing leap year to non-leap clamps February 29 to 28',s.day_from=='28'
          and len(filters.from_day_items(s,None))==28)
    s.month_from='01';s.day_from='31';s.month_from='04'
    check('changing to a 30-day month clamps day 31',s.day_from=='30' and len(filters.from_day_items(s,None))==30)
    s.year_to='2024';s.month_to='02'
    check('To defaults to leap February month end',s.day_to=='29' and s.to_month_end)
    s.year_to='2023'
    check('default month end follows year change to February 28',s.day_to=='28')
    s.year_to='2024';s.month_to='04'
    check('default month end follows a longer month',s.day_to=='30')
    s.day_to='15';s.month_to='03'
    check('an explicitly chosen To day is retained through month changes',s.day_to=='15' and not s.to_month_end)
    s.year_from=s.year_to='ANY'
    check('returning to No limit hides both lower controls again',not any(key in drawn()
          for key in ('month_from','month_to','day_from','day_to')))
    s.year_from='2023';s.year_to='2024'
    check('choosing a year again retains month/day values',s.month_from=='04' and s.day_from=='30'
          and s.month_to=='03' and s.day_to=='15')
    reset()
    check('Reset restores No limit and January 1 / December 31 defaults',
          (s.year_from,s.year_to,s.month_from,s.month_to,s.day_from,s.day_to)==('ANY','ANY','01','12','01','31') and s.to_month_end)
    fixtures();s.year_from=s.year_to='2024';s.month_from=s.month_to='02';apply()
    check('applied calendar bounds include the complete From/To days',visibility()==[False,True,True,False,False]
          and json.loads(s.applied)['date_to']=='2024-02-29')
    s.include_undated=True;apply()
    check('Include undated remains functional',visibility()==[False,True,True,False,True])
    before=visibility();snapshot=s.applied;s.month_from='03'
    try:result=apply()
    except RuntimeError:result={'CANCELLED'}
    check('reversed calendar range is rejected before visibility or snapshot changes',
          result=={'CANCELLED'} and before==visibility() and snapshot==s.applied)
    reset();s.year_to='2024';s.month_to='02';apply()
    check('No limit lower bound admits earlier captures',visibility()==[True,True,True,False,False])
    reset();s.year_from='2024';s.month_from='02';apply()
    check('No limit upper bound admits later captures',visibility()==[False,True,True,True,False])
    reset()
    check('two No limit bounds also admit undated captures',all(visibility()))
    legacy=filters.encode(dict(date_from='2023-06-15',date_to='2024-02-20',categories=list(data.CATEGORIES),include_undated=False))
    filters.apply_existing(bpy.context.view_layer,filters.compiled(legacy))
    s.applied=legacy;s.date_from='2023-06-15';s.date_to='2024-02-20';s.calendar_version=0
    prior=visibility();filters.migrate_settings(s)
    check('legacy text dates migrate with exact day precision',
          (s.year_from,s.month_from,s.day_from,s.year_to,s.month_to,s.day_to)==('2023','06','15','2024','02','20'))
    check('migration preserves the applied snapshot and visibility',s.applied==legacy and prior==visibility()
          and not s.legacy_applied)
    apply()
    check('reapplying an exact migrated date range preserves its behavior',s.applied==legacy and prior==visibility())
    # Recreate old-file storage, with no new enum keys. load_post must migrate it
    # before the sidebar uses it. This also checks real saved numeric enum IDs.
    s.adjusting_calendar=True
    for name in ('year_from','year_to','month_from','month_to','day_from','day_to'):s.property_unset(name)
    s.adjusting_calendar=False;s.calendar_version=0
    state['legacy']=legacy
    state['path']=str(args.artifacts/'legacy_calendar_fixture.blend')
    bpy.ops.wm.save_as_mainfile(filepath=state['path'])


def finish(error=None):
    addon.unregister()
    result=dict(status='failed' if error else 'passed',error=error,checks=checks,
                blender=bpy.app.version_string,elapsed=time.monotonic()-state['start'])
    (args.artifacts/'calendar.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        assert time.monotonic()-state['start']<90,'Calendar UI check timeout'
        state['ticks']+=1
        if state['ticks']<5:return .05
        state['ticks']=0
        window=bpy.context.window
        area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
        region=next(r for r in area.regions if r.type=='WINDOW')
        with bpy.context.temp_override(window=window,area=area,region=region):
            p=state['phase']
            if p==0:
                window.event_simulate(type='ESC',value='PRESS');window.event_simulate(type='ESC',value='RELEASE')
            elif p==1:
                basic();bpy.ops.wm.open_mainfile(filepath=state['path'])
            elif p==2:
                s=settings()
                check('file-load migration restores year/month/day dropdowns automatically',s.calendar_version==2
                      and filters.draft(s)==json.loads(state['legacy']) and s.applied==state['legacy'])
                bpy.ops.ed.undo_push(message='Before calendar filter edit')
                s.day_to='21';apply();bpy.ops.ed.undo()
                check('Undo restores calendar choices and applied snapshot',settings().day_to=='20' and settings().applied==state['legacy'])
                bpy.ops.ed.redo()
                check('Redo restores the new selected day',settings().day_to=='21'
                      and json.loads(settings().applied)['date_to']=='2024-02-21')
                area.spaces.active.show_region_ui=True
                settings().expanded=True
            elif p==3:
                next(r for r in area.regions if r.type=='UI').active_panel_category='Blosm'
                area.tag_redraw()
            elif p==4:
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'calendar_dates.png'))
                settings().year_from=settings().year_to='ANY';area.tag_redraw()
            elif p==5:
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'calendar_no_limit.png'))
                check('actual sidebar hides month and day after No limit is chosen',not any(name in drawn()
                      for name in ('month_from','month_to','day_from','day_to')))
                finish();return None
        state['phase']+=1
        return .05
    except Exception:
        finish(traceback.format_exc());return None


bpy.app.timers.register(tick,first_interval=.5,persistent=True)
