"""Switch a saved real network to direct GN; verify a script-free reopen."""
import argparse
import json
from pathlib import Path
import sys
import time
import bpy
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT.parent),str(ROOT/'tests')]
from blosm.way.street_network.python_backend import ui
from blosm.way.street_network.python_backend.scene import NetworkScene
from blosm.way.street_network.python_backend.terrain.blender import TerrainSampler


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--scene',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:]);args.output.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()));ui.register()
    controllers=[o for o in bpy.data.objects if o.get('st_py_source_data')]
    start=time.perf_counter()
    for c in controllers:
        c.blosm_street_terrain.mode='DIRECT';NetworkScene.update(c)
    bpy.context.view_layer.update()
    groups={m.node_group.as_pointer() for c in controllers for o in c['st_py_outputs'].objects for m in o.modifiers if m.type=='NODES'}
    assert len(groups)==1
    topology={o.name:len(o.data.polygons) for c in controllers for o in c['st_py_outputs'].objects if o.type=='MESH'}
    path=args.output/'python_street_import_direct.blend'
    bpy.context.preferences.filepaths.file_preview_type='NONE'
    bpy.ops.wm.save_as_mainfile(filepath=str(path.resolve()))
    ui.unregister()
    bpy.ops.wm.open_mainfile(filepath=str(path.resolve()),use_scripts=False)
    bpy.context.view_layer.update()
    assert not hasattr(bpy.types.Object,'blosm_street_terrain')
    checked=0;error=0.
    for c in [o for o in bpy.data.objects if o.get('st_py_source_data')]:
        sampler=None;scale=bpy.context.scene.unit_settings.scale_length
        for o in c['st_py_outputs'].objects:
            if o.type!='MESH':continue
            modifier=o.modifiers[0]
            inputs={s.name:getattr(modifier.properties.inputs,s.identifier).value
                for s in modifier.node_group.interface.items_tree
                if s.item_type=='SOCKET' and s.in_out=='INPUT' and s.name in ('Terrain','Prepared Datum','Clearance')}
            if sampler is None:sampler=TerrainSampler(inputs['Terrain'],c,scale)
            raw=np.empty((len(o.data.vertices),3),np.float32);o.data.vertices.foreach_get('co',raw.ravel())
            evaluated=o.evaluated_get(bpy.context.evaluated_depsgraph_get()).data
            assert len(evaluated.polygons)==topology[o.name]
            xyz=np.empty_like(raw);evaluated.vertices.foreach_get('co',xyz.ravel())
            expected=sampler(raw[:,:2]*scale)+(raw[:,2]-inputs['Prepared Datum']+inputs['Clearance'])*scale
            error=max(error,float(np.max(np.abs(xyz[:,2]*scale-expected))))
            assert error < 1e-4,error
            hit=np.empty(len(raw),bool);evaluated.attributes['st_terrain_hit'].data.foreach_get('value',hit)
            assert hit.all()
            np.testing.assert_allclose(raw[:,:2],xyz[:,:2],atol=1e-5/scale,rtol=0)
            checked+=len(raw)
    report=dict(passed=True,shared_groups=len(groups),meshes=len(topology),vertices=checked,
        maximum_height_error_m=error,script_free_reopen=True,seconds=time.perf_counter()-start)
    (args.output/'checks.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('SAVED_TERRAIN_DIRECT_OK',json.dumps(report),flush=True)


if __name__=='__main__':main()
