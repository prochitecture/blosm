"""Launch only isolated test processes; never alter Blender's installed files.

python tests/run_mapillary_checks.py [--ui] [--count 5000]
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('--ui',action='store_true')
parser.add_argument('--count',type=int,default=5000)
parser.add_argument('--blender',default=r'C:\Program Files\Blender Foundation\Blender 5.2\blender.exe')
args=parser.parse_args()
root=Path(__file__).resolve().parents[1]
artifacts=Path(r'C:\dev\workspace\artifacts\mapillary_import')/('ui' if args.ui else 'background')
profile=Path(r'C:\dev\workspace\artifacts\blender_profiles\01a0edf2-527a-7e10-be5c-d4e614d29008')
env=os.environ.copy()
for key,part in (('BLENDER_USER_RESOURCES',''),('BLENDER_USER_CONFIG','config'),
                 ('BLENDER_USER_SCRIPTS','scripts'),('BLENDER_USER_EXTENSIONS','extensions'),
                 ('BLENDER_USER_DATAFILES','datafiles')):
    path=profile/part
    path.mkdir(parents=True,exist_ok=True)
    env[key]=str(path)
env['PYTHONDONTWRITEBYTECODE']='1'
artifacts.mkdir(parents=True,exist_ok=True)
cmd=[args.blender,'--factory-startup','--disable-autoexec','--online-mode','--python-exit-code','1']
if not args.ui: cmd+=['--background']
cmd+=['--python',str(root/'tests'/('check_mapillary_ui.py' if args.ui else 'check_mapillary_import.py')),
      '--','--artifacts',str(artifacts),'--count',str(args.count)]
with (artifacts/'blender.log').open('w',encoding='utf-8') as log:
    process=subprocess.Popen(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,
        creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    try:
        code=process.wait(timeout=600)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()
        raise SystemExit('Isolated Blender check exceeded 600 seconds; see '+str(artifacts))
print('Blender exit:',code,'; report directory:',artifacts)
if args.ui and code==0:
    report=artifacts/'ui.json'
    if not report.is_file() or json.loads(report.read_text())['status']!='passed':
        code=1
if code:
    print((artifacts/'blender.log').read_text(encoding='utf-8')[-7000:])
raise SystemExit(code)
