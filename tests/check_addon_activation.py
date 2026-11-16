"""Exercise Blender's actual restricted addon activation, not direct register().

Run with factory startup, disabled autoexec and an isolated user profile. Works
in background mode; a UI run also captures the registered sidebar tab. All
output belongs in the --artifacts directory, outside the source checkout.
"""
import argparse
import importlib
import json
from pathlib import Path
import sys
import time
import traceback

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import addon_utils
import bpy

parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
checks=[]
errors=[]
state={'start':time.monotonic(),'enabled':False}


def check(label,condition):
    assert condition,label
    checks.append(label)


def failed(exc):
    errors.append(''.join(traceback.format_exception(type(exc),exc,exc.__traceback__)))


def enable():
    module=addon_utils.enable(ROOT.name,default_set=True,handle_error=failed)
    check('Normal addon activation succeeds',module is not None and not errors)
    state['enabled']=True
    check('Preference and runtime enabled flags agree',addon_utils.check(ROOT.name)==(True,True))
    check('Blosm sidebar panel is registered',hasattr(bpy.types,'BLOSM_PT_Extent'))
    panel=bpy.types.BLOSM_PT_Extent
    check('Panel belongs to the Object Mode 3D View sidebar',
          (panel.bl_space_type,panel.bl_region_type,panel.bl_context,panel.bl_category)==
          ('VIEW_3D','UI','objectmode','Blosm'))
    check('Scene and filter RNA properties are registered',
          hasattr(bpy.types.Scene,'blosm') and hasattr(bpy.types.ViewLayer,'blosm_mapillary_filters'))
    return importlib.import_module(ROOT.name+'.viewport.visibility')


def disable(visibility):
    addon_utils.disable(ROOT.name,default_set=True,handle_error=failed)
    state['enabled']=False
    check('Normal disable succeeds',not errors)
    check('No deferred recovery survives disable',not bpy.app.timers.is_registered(visibility._recover_timer))
    check('Sidebar panel and temporary records are unregistered',
          not hasattr(bpy.types,'BLOSM_PT_Extent') and not hasattr(bpy.types.ViewLayer,'blosm_temporary_hidden'))
    check('Visibility handlers are removed',all(callback not in group for group,callback in visibility.handlers()))


def finish(error=None):
    if state['enabled']:
        try:addon_utils.disable(ROOT.name,default_set=True,handle_error=failed)
        except Exception:errors.append(traceback.format_exc())
    result=dict(status='failed' if error or errors else 'passed',error=error,errors=errors,
                blender=bpy.app.version_string,background=bpy.app.background,checks=checks,
                elapsed=time.monotonic()-state['start'])
    (args.artifacts/'activation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)
    if bpy.app.background:
        if result['status']!='passed':raise RuntimeError(error or str(errors))
    else:bpy.ops.wm.quit_blender()


def camera(name):
    obj=bpy.data.objects.new(name,bpy.data.cameras.new(name))
    bpy.context.scene.collection.objects.link(obj)
    bpy.context.view_layer.update()
    return obj


def tick():
    try:
        check('Deferred recovery timer completed',not bpy.app.timers.is_registered(visibility._recover_timer))
        check('Saved stale camera visibility recovered',not state['stale'].hide_get())
        check('Active visibility ownership survived initial recovery',state['live'].hide_get())
        check('Only the live ownership record remains',len(bpy.context.view_layer.blosm_temporary_hidden)==1)
        disable(visibility)
        check('Disable restores the live camera',not state['live'].hide_get())
        restored=enable()
        check('Reactivation schedules fresh recovery',bpy.app.timers.is_registered(restored._recover_timer))
        # A stale record must also recover when disabled before the first timer.
        obj=state['stale'];obj.hide_set(True)
        record=bpy.context.view_layer.blosm_temporary_hidden.add()
        record.name='disable-before-first-tick';record.object=obj;record.hidden=False
        disable(restored)
        check('Immediate disable recovers stale saved visibility',not obj.hide_get())
        enable()
        if not bpy.app.background:
            window=bpy.context.window
            window.event_simulate(type='ESC',value='PRESS')
            window.event_simulate(type='ESC',value='RELEASE')
            area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
            area.spaces.active.show_region_ui=True
            # Allow a redraw so the vertical Blosm tab is present in the image.
            bpy.app.timers.register(capture,first_interval=.5)
        else:finish()
    except Exception:finish(traceback.format_exc())
    return None


def capture():
    try:
        bpy.ops.screen.screenshot(filepath=str(args.artifacts/'sidebar.png'))
        check('Sidebar screenshot captured',(args.artifacts/'sidebar.png').is_file())
        finish()
    except Exception:finish(traceback.format_exc())
    return None


try:
    visibility=enable()
    check('Initial recovery is deferred',bpy.app.timers.is_registered(visibility._recover_timer))
    disable(visibility)
    visibility=enable()
    stale=camera('Stale saved visibility')
    stale.hide_set(True)
    record=bpy.context.view_layer.blosm_temporary_hidden.add()
    record.name='old-viewer';record.object=stale;record.hidden=False
    live=camera('Live temporary visibility')
    state.update(stale=stale,live=live,lease=visibility.TemporaryHide(live,bpy.context.view_layer))
    if bpy.app.background:
        # Background scripts have no timer event loop. Execute the callback once
        # after addon_utils has restored bpy.data, then retire the one-shot timer.
        visibility._recover_timer()
        bpy.app.timers.unregister(visibility._recover_timer)
        tick()
    else:bpy.app.timers.register(tick,first_interval=.2)
except Exception:finish(traceback.format_exc())
