"""Native file-load/deleted-owner cache lifetime without a viewport observer."""
import argparse
import gc
import json
from pathlib import Path
import sys
import weakref
import bpy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent));sys.path.insert(0,str(ROOT/'tests'))
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from blosm.way.street_network.python_backend import model,scene,ui


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:]);args.output.mkdir(parents=True,exist_ok=True)
    bpy.context.preferences.filepaths.file_preview_type='NONE'
    target=args.output/'empty.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(target.resolve()))
    observer=tuple(bpy.app.handlers.depsgraph_update_post)
    ui.register();ui.register()
    assert bpy.app.handlers.load_pre.count(scene._clear_on_file_load)==1
    assert tuple(bpy.app.handlers.depsgraph_update_post)==observer
    controllers=[]
    for i in range(2):
        network=model.Network('cache_lifecycle/'+str(i),
            (model.Section(1,1,((0.,i*20.,0.),(10.,i*20.,0.),(20.,i*20.,0.))),))
        controllers.append(scene.NetworkScene.create(network))
    first,removed=controllers
    identity=first['st_py_network_uuid'];gone=removed['st_py_network_uuid']
    before=mesh_buffers(first)
    pointers={o['st_py_key']:o.as_pointer() for o in first['st_py_outputs'].objects}
    first_ref=weakref.ref(scene._ENGINES[identity]);gone_ref=weakref.ref(scene._ENGINES[gone])
    bpy.data.objects.remove(removed,do_unlink=True)
    scene.NetworkScene.update(first);gc.collect()
    print('LIFECYCLE deleted-owner cleanup',flush=True)
    assert gone not in scene._ENGINES and gone_ref() is None
    assert scene._ENGINES[identity] is first_ref()
    assert equal(before,mesh_buffers(first))
    assert pointers=={o['st_py_key']:o.as_pointer() for o in first['st_py_outputs'].objects}
    project=args.output/'lifecycle.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(project.resolve()))
    bpy.ops.wm.open_mainfile(filepath=str(target.resolve()),load_ui=False,use_scripts=False)
    print('LIFECYCLE other file loaded',flush=True)
    gc.collect()
    assert not scene._ENGINES and first_ref() is None
    bpy.ops.wm.open_mainfile(filepath=str(project.resolve()),load_ui=False,use_scripts=False)
    assert not scene._ENGINES
    first=next(o for o in bpy.data.objects if o.get('st_py_network_uuid')==identity)
    pointers={o['st_py_key']:o.as_pointer() for o in first['st_py_outputs'].objects}
    scene.NetworkScene.update(first)
    assert equal(before,mesh_buffers(first))
    assert pointers=={o['st_py_key']:o.as_pointer() for o in first['st_py_outputs'].objects}
    ui.unregister()
    assert scene._clear_on_file_load not in bpy.app.handlers.load_pre and not scene._ENGINES
    ui.register()
    assert bpy.app.handlers.load_pre.count(scene._clear_on_file_load)==1
    assert tuple(bpy.app.handlers.depsgraph_update_post)==observer
    report=dict(passed=True,deleted_controller_cache_released=True,retained_controller_cache_reused=True,
        opening_another_file_releases_cache=True,garbage_collection_verified=True,
        reopen_exact_meshes=True,reopen_retains_output_identity=True,
        repeat_registration_idempotent=True,unregister_releases_cache=True,no_depsgraph_observer=True)
    (args.output/'lifecycle_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
