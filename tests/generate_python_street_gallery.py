"""Regenerate the complete Python street gallery from authored fixture inputs.

Run in background Blender 5.2 with --output <new-directory>. Optionally select
fixtures with repeated --case and verify captured source/mesh/furniture records
with --reference-catalogue. No saved output mesh is imported.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import sys
import bpy
import numpy as np
from mathutils import Quaternion, Vector

parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--reference-catalogue',type=Path)
parser.add_argument('--case',action='append')
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else sys.argv[1:])
assert bpy.app.background,'Run this fixture generator in background Blender'
ROOT=args.root.resolve()
sys.path.insert(0,str(ROOT.parent));sys.path.insert(0,str(ROOT/'tests'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import blosm
assert Path(blosm.__file__).resolve().parent==ROOT
from blosm.way.street_network.python_backend import model,ui
from blosm.way.street_network.python_backend.scene import NetworkScene,SourceObjects,clear_runtime
from check_python_street_suite import source_digest
import _python_street_catalogue as catalogue
SOURCE=source_digest()


def snapshot(controller):
    meshes,furniture={},{}
    for obj in controller['st_py_outputs'].objects:
        key=obj['st_py_key']
        assert not any(modifier.type=='NODES' for modifier in obj.modifiers)
        if obj.type=='MESH':
            assert all(polygon.loop_total==4 for polygon in obj.data.polygons)
            xyz=np.empty(3*len(obj.data.vertices),dtype=np.float32)
            faces=np.empty(len(obj.data.loops),dtype=np.int32)
            obj.data.vertices.foreach_get('co',xyz);obj.data.loops.foreach_get('vertex_index',faces)
            digest=hashlib.sha256(struct.pack('<QQ',len(xyz),len(faces)))
            digest.update(xyz.tobytes());digest.update(faces.tobytes());meshes[key]=digest.hexdigest()
        elif obj.type=='EMPTY':
            furniture[key]=(tuple(obj.location),tuple(obj.rotation_euler),obj.get('st_py_asset'))
    return json.loads(json.dumps(dict(source=SourceObjects.read(controller).to_dict(),mesh=meshes,furniture=furniture)))


def bounds(controller):
    objects=(*controller['st_py_sources'].all_objects,*controller['st_py_outputs'].all_objects)
    points=[obj.matrix_world@Vector(corner) for obj in objects
            if obj.type=='MESH' and len(obj.data.vertices) for corner in obj.bound_box]
    xy=np.asarray([point[:2] for point in points])
    assert len(xy) and np.all(np.isfinite(xy))
    return xy.min(axis=0),xy.max(axis=0)


class GalleryLayout:
    def __init__(self,columns=4,padding=40.):
        self.columns,self.padding=columns,padding
        self.x=self.y=self.row_height=0.
        self.count=0

    def place(self,controller):
        low,high=bounds(controller);size=high-low
        controller.location=(self.x-low[0],self.y-high[1],0.)
        self.x+=size[0]+self.padding
        self.row_height=max(self.row_height,size[1])
        self.count+=1
        if self.count%self.columns==0:
            self.x=0.;self.y-=self.row_height+self.padding;self.row_height=0.


assert not any(obj.get('st_py_source_data') for obj in bpy.data.objects),'Use a fresh background process'
args.output.mkdir(parents=True,exist_ok=False)
for obj in list(bpy.data.objects):bpy.data.objects.remove(obj,do_unlink=True)
if not ui.BLOSM_PT_PythonStreets.is_registered:ui.register()
pavement=bpy.data.materials.new('Python streets / pavement')
pavement.diffuse_color=(.18,.25,.29,1.)
grass=bpy.data.materials.new('Python streets / grass')
grass.diffuse_color=(.22,.4,.17,1.)
materials={'pavement':pavement,'grass':grass}
networks=catalogue.networks(model)
catalogue_ids={network.uuid for network in networks}
assert len(networks)==len(catalogue_ids)
if args.case:
    assert len(set(args.case))==len(args.case)
    assert set(args.case)<={network.uuid for network in networks}
    networks=tuple(network for network in networks if network.uuid in args.case)
reference=None
if args.reference_catalogue:
    reference_report=json.loads(args.reference_catalogue.read_text(encoding='utf8'))
    assert reference_report['passed'] and reference_report['fixtures']==len(catalogue_ids)
    reference={row['case']:row for row in reference_report['records']}
    assert reference.keys()==catalogue_ids
controllers=[];records=[];layout=GalleryLayout()
try:
    for index,network in enumerate(networks):
        controller=NetworkScene.create(network,name=network.uuid,
            grouping='grouped' if index%2 else 'street',materials=materials,generate=False)
        bpy.context.view_layer.update()
        actual_source=json.loads(json.dumps(SourceObjects.read(controller).to_dict()))
        if reference:
            expected=reference[network.uuid]
            if actual_source!=expected['source']:
                (args.output/'source_mismatch.json').write_text(json.dumps(dict(case=network.uuid,
                    expected=expected['source'],actual=actual_source),indent=2),encoding='utf8')
                raise AssertionError('Authored source differs from the canonical gallery: '+network.uuid)
        bpy.context.view_layer.objects.active=controller
        assert bpy.ops.blosm.update_python_streets()=={'FINISHED'}
        assert not controller.get('st_py_error'),controller.get('st_py_error')
        current=snapshot(controller)
        if reference:
            assert current=={key:expected[key] for key in ('source','mesh','furniture')},network.uuid
        bpy.context.view_layer.update();layout.place(controller)
        bpy.context.view_layer.update()
        assert snapshot(controller)==current,'Layout changed source, buffers or furniture'
        controllers.append(controller)
        records.append(dict(case=network.uuid,quads=sum(len(obj.data.polygons)
            for obj in controller['st_py_outputs'].objects if obj.type=='MESH'),
            source_fingerprint=model.fingerprint(actual_source),reference_exact=reference is not None))
        clear_runtime()
        (args.output/'generation_partial.json').write_text(json.dumps(records,indent=2),encoding='utf8')
        print(json.dumps(records[-1]),flush=True)
    boxes=[bounds(controller) for controller in controllers]
    for index,(low,high) in enumerate(boxes):
        assert all(np.any(high<a) or np.any(b<low) for a,b in boxes[index+1:]),'Fixture layout overlaps'
    readme=bpy.data.texts.new('READ ME - Python street gallery')
    readme.write('All gallery meshes are built from the shared authored Python fixture catalogue.\n'
        'Enable this addon, select an input centerline and edit its vertices or Section attributes.\n'
        'Then press Streets > Update. Arrow empties demonstrate source-relative furniture anchors.\n'
        'Actual fixture bounds determine the layout. All sources use their canonical saved settings.\n'
        'The attached/internal examples retain original Street and Section identities.\n')
    bootstrap=bpy.data.texts.new('OPTIONAL - Python street editing.py')
    bootstrap.write('import importlib.util,sys,bpy\nfrom pathlib import Path\np=Path('
        +repr(str(ROOT/'way/street_network/python_backend'))+')\n'
        'name="blosm_python_street_preview"\n'
        'if not (hasattr(bpy.types,"BLOSM_PT_PythonStreets") or hasattr(bpy.types,"BLOSM_PT_python_streets")):\n'
        '    if name not in sys.modules:\n'
        '        spec=importlib.util.spec_from_file_location(name,p/"__init__.py",submodule_search_locations=[str(p)])\n'
        '        module=importlib.util.module_from_spec(spec)\n'
        '        sys.modules[name]=module\n'
        '        spec.loader.exec_module(module)\n'
        '    from blosm_python_street_preview import ui\n'
        '    ui.register()\n')
    low=np.min([box[0] for box in boxes],axis=0);high=np.max([box[1] for box in boxes],axis=0)
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type=='VIEW_3D':
                area.spaces.active.shading.color_type='MATERIAL'
                view=area.spaces.active.region_3d
                view.view_location=(*(.5*(low+high)),0.)
                view.view_rotation=Quaternion((1.,0.,0.,0.));view.view_perspective='ORTHO'
                view.view_distance=float(max(high-low)*1.1)
    bpy.context.preferences.filepaths.file_preview_type='NONE'
    scene=args.output/'python_streets_gallery.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(scene.resolve()))
    assert source_digest()==SOURCE
    report=dict(passed=True,fixtures=len(records),cases=records,scene=str(scene),
        source_sha256=SOURCE,source_unchanged=True,
        catalogue_sha256=hashlib.sha256(Path(catalogue.__file__).read_bytes()).hexdigest(),
        total_quads=sum(row['quads'] for row in records),all_direct_quads=True,
        no_gn_modifiers=True,layout_nonoverlapping=True,reference_checked=reference is not None,
        reusable_runtime_caches_released_between_fixtures=True,
        no_precomputed_output_mesh_imported=True)
    (args.output/'generation_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps({key:value for key,value in report.items() if key!='cases'}),flush=True)
except Exception:
    import traceback
    (args.output/'failure.txt').write_text(traceback.format_exc(),encoding='utf8')
    raise
