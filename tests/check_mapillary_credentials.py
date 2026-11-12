"""Isolated Blender operator check using the real local file; no HTTP needed."""
import argparse
import importlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
import bpy
import addon_utils

parser=argparse.ArgumentParser()
parser.add_argument('--artifacts',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
args.artifacts.mkdir(parents=True,exist_ok=True)
addon_utils.enable(ROOT.name,default_set=True)
credentials=importlib.import_module(ROOT.name+'.mapillary.credentials')
adapter=importlib.import_module(ROOT.name+'.mapillary.progressive')
source=Path(r'C:\dev\workspace\blosm\mapillary\mapillary_token.txt')
prefs=bpy.context.preferences.addons[ROOT.name].preferences
bpy.context.scene.blosm.dataType='mapillary'
bpy.context.scene.blosm.progressiveImport=False
seen=[]
reads=[]
real_read=Path.read_text


def read(path,*a,**kw):
    if path==source: reads.append(True)
    return real_read(path,*a,**kw)


class CaptureImport:
    def __init__(self,context,bounds,projection,origin,set_origin,token):
        seen.append(token)
    def run_blocking(self): return True


with patch.object(credentials,'TOKEN_FILE',source), patch.object(Path,'read_text',read), patch.object(adapter,'MapillaryImport',CaptureImport):
    prefs.mapillaryAccessToken='fixture-preference-token'
    assert bpy.ops.blosm.import_data()=={'FINISHED'}
    assert not reads and seen[-1]=='fixture-preference-token'
    prefs.mapillaryAccessToken=''
    assert bpy.ops.blosm.import_data()=={'FINISHED'}
    assert len(reads)==1 and seen[-1] and not any(c.isspace() for c in seen[-1])
    assert bpy.ops.blosm.import_data()=={'FINISHED'}
    assert len(reads)==1 and seen[-1]==seen[-2]
    prefs.mapillaryAccessToken='fixture-changed-preference'
    assert bpy.ops.blosm.import_data()=={'FINISHED'}
    assert seen[-1]=='fixture-changed-preference' and len(reads)==1
    prefs.mapillaryAccessToken=''
    assert bpy.ops.blosm.import_data()=={'FINISHED'}
    assert seen[-1]==seen[1] and len(reads)==1
    report={'status':'passed','file_read_count':len(reads),
            'checks':['preference priority','lazy fallback','one read across repeated imports',
                      'preference override after fallback','cached fallback after clearing preferences']}
    (args.artifacts/'credentials.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report))
seen.clear()
addon_utils.disable(ROOT.name,default_set=True)
