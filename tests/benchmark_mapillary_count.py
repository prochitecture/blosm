"""Offline timer benchmark against the pre-optimization counting implementation."""
import argparse
import ast
import importlib
import json
from pathlib import Path
import subprocess
import sys
import time
import types

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
addon=importlib.import_module(ROOT.name)
addon.register()
adapter=importlib.import_module(ROOT.name+'.mapillary.progressive')
Projection=importlib.import_module(ROOT.name+'.util.transverse_mercator').TransverseMercator
parser=argparse.ArgumentParser()
parser.add_argument('--mode',choices=('baseline','optimized'),required=True)
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory=str(args.artifacts)
data_root=Path('C:/dev/workspace/artifacts/mapillary_import/count_optimization')
baseline=json.loads((data_root/'baseline.json').read_text())
payload=(data_root/'coverage.mvt').read_bytes()
class Transport:
    max_requests=1
    pending=None
    def submit(self,request):
        if '/images?' in request.url: self.pending=None
        else:self.pending=(1,payload,None)
        return 1
    def poll(self):
        result,self.pending=self.pending,None
        return result
    def close(self):self.pending=None
bounds=baseline['bounds']
origin=((bounds[1]+bounds[3])/2,(bounds[0]+bounds[2])/2)
task=adapter.MapillaryImport(bpy.context,bounds,Projection(lat=origin[0],lon=origin[1]),origin,False,
                             'offline-fixture-token',transport=Transport())
if args.mode=='baseline':
    def source(path):
        return subprocess.check_output(['git','show','9dd2509:'+path],cwd=ROOT,text=True)
    old_tiles=types.ModuleType('baseline_tiles')
    exec(compile(source('mapillary/tiles.py'),'baseline_tiles.py','exec'),old_tiles.__dict__)
    tree=ast.parse(source('mapillary/progressive.py'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='MapillaryImport')
    method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='count_tile')
    namespace=dict(adapter.__dict__)
    namespace['image_points']=old_tiles.image_points
    exec(compile(ast.Module(body=[method],type_ignores=[]),'baseline_count.py','exec'),namespace)
    task.count_tile=types.MethodType(namespace['count_tile'],task)
started=previous=None
gaps=[]
turns=[]
def tick():
    global started,previous
    now=time.monotonic()
    if started is None:started=now
    if previous is not None:gaps.append(now-previous)
    previous=now
    step=time.monotonic()
    task.job.step()
    turns.append(time.monotonic()-step)
    if task.job.error or task.progress.counted_tiles or now-started>240:
        report={'mode':args.mode,'elapsed_seconds':time.monotonic()-started,
                'count_seconds':task.metrics['count_seconds'],'candidates':task.progress.candidates,
                'max_turn_seconds':max(turns),'heartbeats':len(gaps),
                'max_heartbeat_gap':max(gaps,default=0),'error':str(task.job.error) if task.job.error else None,
                'status':task.status,'blender':bpy.app.version_string}
        assert task.progress.candidates==len(baseline['ids']),report
        task.cancel()
        (args.artifacts/'timer.json').write_text(json.dumps(report,indent=2))
        print(json.dumps(report),flush=True)
        bpy.ops.wm.quit_blender()
        return None
    return .03
bpy.app.timers.register(tick,first_interval=.5)
