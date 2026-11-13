"""View an already downloaded real capture; optionally refresh through live HTTP.

Supply --scene with the photo validation fixture. Generated evidence only goes
to --artifacts. Credentials, when supplied, are read once and never persisted.
"""
import argparse
import importlib
import json
from pathlib import Path
import statistics
import sys
import time
import traceback
import ctypes

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
import gpu
parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
parser.add_argument('--scene',type=Path,required=True)
parser.add_argument('--token-file',type=Path)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
bpy.context.preferences.filepaths.file_preview_type='NONE'
addon=importlib.import_module(ROOT.name);addon.register()
viewport=importlib.import_module(ROOT.name+'.viewport')
pano=importlib.import_module(ROOT.name+'.mapillary.panorama')
photo=importlib.import_module(ROOT.name+'.mapillary.photo')
from blosm.viewport.ui import InputEvent
with bpy.data.libraries.load(str(args.scene),link=False) as (source,destination):
    destination.objects=[name for name in source.objects if '1079237249490765' in name]
    destination.images=source.images
obj=next(obj for obj in destination.objects if obj is not None and pano.is_panorama(obj))
bpy.context.scene.collection.objects.link(obj)
bpy.context.view_layer.update()
token=None
if args.token_file:
    token=args.token_file.read_text(encoding='utf-8-sig').strip()
    photo.get_token=lambda *_:token
    for image in tuple(bpy.data.images):
        if image.get('mapillary_image_id')==obj['mapillary_image_id']:bpy.data.images.remove(image)
window=bpy.context.window
area=next(a for a in window.screen.areas if a.type=='VIEW_3D')
region=next(r for r in area.regions if r.type=='WINDOW')
state={'phase':0,'ticks':0,'start':time.perf_counter(),'frames':[],'draws':[],'move':[]}


def memory():
    if sys.platform!='win32':return None
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(n,ctypes.c_size_t) for n in
            ('peak_working','working','peak_paged','paged','peak_nonpaged','nonpaged','commit','peak_commit')]
    counters=Counters();counters.cb=ctypes.sizeof(counters)
    if ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.c_void_p(-1),ctypes.byref(counters),counters.cb):
        return {'working_bytes':counters.working,'private_commit_bytes':counters.commit}


def finish(error=None):
    if viewport.current():viewport.current().close('USER')
    addon.unregister()
    result={'status':'failed' if error else 'passed','error':error,'blender':bpy.app.version_string,
            'backend':gpu.platform.backend_type_get(),'renderer':gpu.platform.renderer_get(),
            'image_id':obj['mapillary_image_id'],'fresh_download':bool(args.token_file),
            'viewport_pixels':[region.width,region.height],
            'first_usable_seconds':state.get('first'), 'repeats':state.get('repeats',0)}
    result['memory_before_repeats']=state.get('memory_before')
    result['memory_after_repeats']=state.get('memory_after')
    for key in ('frames','draws','move'):
        values=state[key]
        if values:result[key]={'count':len(values),'median_ms':statistics.median(values)*1000,
                               'p95_ms':sorted(values)[min(len(values)-1,int(len(values)*.95))]*1000}
    text=json.dumps(result,indent=2)
    assert not token or token not in text
    (args.artifacts/'real.json').write_text(text,encoding='utf8')
    print(text,flush=True)
    bpy.ops.wm.quit_blender()


def tick():
    try:
        assert time.perf_counter()-state['start']<150,'Real-photo check timed out'
        state['ticks']+=1
        with bpy.context.temp_override(window=window,area=area,region=region):
            p=state['phase']
            if p==0:
                window.event_simulate(type='ESC',value='PRESS');window.event_simulate(type='ESC',value='RELEASE')
                state['phase']=1
            elif p==1:
                state['open_started']=time.perf_counter()
                state['surface']=surface=pano.Panorama(obj)
                state['session']=viewport.open(bpy.context,surface)
                original_draw=surface.draw
                def measured(frame,painter):
                    started=time.perf_counter()
                    original_draw(frame,painter)
                    state['draws'].append(time.perf_counter()-started)
                surface.draw=measured
                state['phase']=2
            elif p==2:
                surface=state['surface']
                assert not surface.failed,surface.status
                assert not state['session'].closed,'Viewer closed'
                if surface.ready:
                    state['first']=time.perf_counter()-state['open_started']
                    state['phase']=3;state['ticks']=0
            elif p==3 and state['ticks']>8:
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'real_panorama.png'))
                state['surface'].opacity=.35
                state['session'].invalidate()
                state['phase']=4;state['ticks']=0
            elif p==4 and state['ticks']>8:
                bpy.ops.screen.screenshot(filepath=str(args.artifacts/'real_transparent.png'))
                state['surface'].opacity=1.
                state['surface'].handle(InputEvent('LEFTMOUSE','PRESS',500,400))
                state['phase']=5;state['ticks']=0;state['previous']=time.perf_counter()
            elif p==5:
                now=time.perf_counter();state['frames'].append(now-state['previous']);state['previous']=now
                state['surface'].handle(InputEvent('MOUSEMOVE','NOTHING',500+state['ticks']*3,400))
                state['move'].append(time.perf_counter()-now)
                if state['ticks']>=120:
                    state['surface'].handle(InputEvent('LEFTMOUSE','RELEASE',860,400))
                    state['session'].close('USER')
                    count=len(bpy.data.images)
                    state['memory_before']=memory()
                    for _ in range(20):
                        item=viewport.open(bpy.context,pano.Panorama(obj));item.close('USER')
                    assert len(bpy.data.images)==count,'Repeated opening leaked images'
                    state['repeats']=20
                    import gc
                    gc.collect()
                    state['memory_after']=memory()
                    bpy.ops.wm.save_as_mainfile(filepath=str(args.artifacts/'real_panorama.blend'))
                    finish();return None
        return .016 if state['phase']==5 else .1
    except Exception as exc:
        message=photo.redact_text(exc,headers={'Authorization':'OAuth '+token} if token else {})
        finish(message);return None

bpy.app.timers.register(tick,first_interval=1.)
