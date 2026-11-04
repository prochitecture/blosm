"""Native terrain Update/mode/reopen regression and a reproducible gallery."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import bpy
import numpy as np
from mathutils import Vector

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT.parent),str(ROOT/'tests')]
from blosm.way.street_network.python_backend import model,ui
from blosm.way.street_network.python_backend.scene import NetworkScene,clear_runtime
from blosm.way.street_network.python_backend.terrain.blender import TerrainSampler,RAW_Z
from _python_street_fixtures import fixture_networks,transition_fixture_networks,short_junction_pair


def buffers(controller,evaluated=False):
    bpy.context.view_layer.update()
    result={}
    for obj in controller['st_py_outputs'].objects:
        if obj.type!='MESH':continue
        mesh=obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).data if evaluated else obj.data
        xyz=np.empty((len(mesh.vertices),3),np.float32); mesh.vertices.foreach_get('co',xyz.ravel())
        result[obj['st_py_key']]=xyz
    return result


def make_terrain(name,offset=(0.,0.,0.),extent=100.,scale=1.):
    x=np.linspace(-extent,extent,41)
    xx,yy=np.meshgrid(x,x)
    zz=2.+.025*xx+.04*yy+.4*np.exp(-((xx-12.)**2+(yy+3.)**2)/180.)
    xyz=np.column_stack((xx.ravel(),yy.ravel(),zz.ravel()))/scale
    quads=[(j*41+i,j*41+i+1,(j+1)*41+i+1,(j+1)*41+i) for j in range(40) for i in range(40)]
    mesh=bpy.data.meshes.new(name);mesh.from_pydata(xyz,[],quads);mesh.update()
    obj=bpy.data.objects.new(name,mesh);bpy.context.scene.collection.objects.link(obj)
    obj.location=np.array(offset)/scale
    material=bpy.data.materials.get('Terrain demonstration') or bpy.data.materials.new('Terrain demonstration')
    material.diffuse_color=(.24,.32,.16,1.)
    mesh.materials.append(material)
    return obj


def verify(controller):
    before=buffers(controller)
    original_target=controller.blosm_street_terrain.target
    ground_before=np.empty((len(original_target.data.vertices),3),np.float32)
    original_target.data.vertices.foreach_get('co',ground_before.ravel())
    output_ids={o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    check=NetworkScene.update(controller)
    assert check.timings['published_objects']==0,check.timings
    assert output_ids=={o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    controller.blosm_street_terrain.mode='DIRECT'
    NetworkScene.update(controller)
    raw=buffers(controller); direct=buffers(controller,True)
    sampler=TerrainSampler(controller.blosm_street_terrain.target,controller,bpy.context.scene.unit_settings.scale_length)
    scale=bpy.context.scene.unit_settings.scale_length
    groups=set()
    for obj in controller['st_py_outputs'].objects:
        if obj.type!='MESH':continue
        assert len(obj.modifiers)==1
        groups.add(obj.modifiers[0].node_group.as_pointer())
        key=obj['st_py_key'];datum=json.loads(obj['st_py_metadata'])['terrain_datum']
        expected=sampler(raw[key][:,:2]*scale)+(raw[key][:,2]*scale-datum)+controller.blosm_street_terrain.clearance
        np.testing.assert_allclose(direct[key][:,2]*scale,expected,atol=2e-5,rtol=0)
        np.testing.assert_allclose(direct[key][:,:2],raw[key][:,:2],atol=1e-5/scale,rtol=0)
    assert len(groups)==1
    controller.blosm_street_terrain.mode='PROFILE'
    NetworkScene.update(controller)
    for key,xyz in buffers(controller).items():
        np.testing.assert_allclose(xyz,before[key],atol=2e-5/scale,rtol=0)
    assert not any(o.modifiers for o in controller['st_py_outputs'].objects)
    # Failed coverage cannot publish a partial replacement or advance manifest.
    old_manifest=controller['st_py_manifest'];old=buffers(controller)
    terrain=controller.blosm_street_terrain.target
    terrain.location.x+=10000.;bpy.context.view_layer.update()
    try:
        NetworkScene.update(controller)
    except ValueError as exc:
        assert 'no surface below' in str(exc)
    else: raise AssertionError('Missing terrain should reject Update')
    assert controller['st_py_manifest']==old_manifest
    for key,xyz in buffers(controller).items():np.testing.assert_array_equal(xyz,old[key])
    terrain.location.x-=10000.;bpy.context.view_layer.update()
    NetworkScene.update(controller)
    invalid=terrain.modifiers.new('Invalid feedback test','SHRINKWRAP')
    invalid.target=next(o for o in controller['st_py_outputs'].objects if o.type=='MESH')
    try:
        NetworkScene.update(controller)
    except ValueError as exc:
        assert 'must not depend' in str(exc)
    else:raise AssertionError('Terrain feedback dependency was accepted')
    terrain.modifiers.remove(invalid);bpy.context.view_layer.update();NetworkScene.update(controller)
    # A terrain edit affects published streets and attached furniture only on Update.
    initial=buffers(controller)
    anchors={o['st_py_key']:tuple(o.location) for o in controller['st_py_outputs'].objects if o.type=='EMPTY'}
    terrain.location.z+=.4/scale;bpy.context.view_layer.update()
    for key,xyz in buffers(controller).items():np.testing.assert_array_equal(xyz,initial[key])
    NetworkScene.update(controller)
    for key,xyz in buffers(controller).items():np.testing.assert_allclose((xyz[:,2]-initial[key][:,2])*scale,.4,atol=2e-5)
    for obj in controller['st_py_outputs'].objects:
        if obj.type=='EMPTY':
            np.testing.assert_allclose((obj.location.z-anchors[obj['st_py_key']][2])*scale,.4,atol=2e-5)
    terrain.location.z-=.4/scale;bpy.context.view_layer.update();NetworkScene.update(controller)
    original=buffers(controller)
    controller.blosm_street_terrain.target=None
    NetworkScene.update(controller)
    for obj in controller['st_py_outputs'].objects:
        if obj.type!='MESH':continue
        raw=np.empty(len(obj.data.vertices),np.float32);obj.data.attributes[RAW_Z].data.foreach_get('value',raw)
        np.testing.assert_allclose(buffers(controller)[obj['st_py_key']][:,2]*scale,raw,atol=1e-6)
    controller.blosm_street_terrain.target=terrain
    NetworkScene.update(controller)
    for key,xyz in buffers(controller).items():np.testing.assert_allclose(xyz,original[key],atol=2e-5/scale,rtol=0)
    ground_after=np.empty_like(ground_before);original_target.data.vertices.foreach_get('co',ground_after.ravel())
    np.testing.assert_array_equal(ground_after,ground_before)
    return dict(mode_switch=True,no_op=True,missing_hit_rollback=True,terrain_update=True,
                shared_projection_group=True,remove_terrain_restores_prepared_mesh=True,
                furniture_follows_terrain=True,feedback_rejected=True,original_terrain_unchanged=True)


def verify_edits(controller):
    from check_python_street_scene import section_control
    original=buffers(controller)
    obj,vertex=section_control(controller,1)
    position=vertex.co.copy();vertex.co.y+=.15;obj.data.update()
    NetworkScene.update(controller)
    changed=buffers(controller)
    assert not np.array_equal(changed['road:1'],original['road:1'])
    vertex.co=position;obj.data.update();NetworkScene.update(controller)
    for key,xyz in buffers(controller).items():np.testing.assert_allclose(xyz,original[key],atol=3e-5,rtol=0)
    layer=obj.data.attributes['st_width']
    selected=[edge.index for edge in obj.data.edges if obj.data.attributes['st_section_id'].data[edge.index].value==1]
    for i in selected:layer.data[i].value+=.2
    obj.data.update();NetworkScene.update(controller)
    assert not np.array_equal(buffers(controller)['road:1'],original['road:1'])
    for i in selected:layer.data[i].value-=.2
    obj.data.update();NetworkScene.update(controller)
    return dict(centerline_edit=True,width_edit=True)


def verify_units_and_transforms():
    original_scene=bpy.context.scene
    scene=bpy.data.scenes.new('Terrain units test');bpy.context.window.scene=scene
    scene.unit_settings.scale_length=.01
    source=model.Network('centimetres', (model.Section(1,1,((-25.,0.,0.),(0.,1.,0.),(25.,0.,0.)),datum=7.),),
                         settings=model.Settings(exaggeration=1.))
    controller=NetworkScene.create(source,generate=False)
    controller.location=(200.,300.,400.)
    controller.rotation_euler.z=.2
    target=make_terrain('Transformed original terrain',extent=100.,scale=.01)
    target.location=(200.,300.,400.);target.rotation_euler.z=.1;target.scale=(1.1,.9,1.3)
    controller.blosm_street_terrain.target=target;bpy.context.view_layer.update()
    NetworkScene.update(controller)
    result=verify(controller)
    for obj in list(scene.objects):bpy.data.objects.remove(obj,do_unlink=True)
    bpy.context.window.scene=original_scene;bpy.data.scenes.remove(scene)
    clear_runtime()
    return dict(scene_scale=.01,prepared_datum=7.,transformed_terrain=True,rotated_controller=True,**result)


def present_gallery(controllers, output):
    scene=bpy.context.scene
    for controller in controllers:
        font=bpy.data.curves.new('Fixture label','FONT')
        font.body=controller.name.replace('terrain_','').replace('_',' ')
        font.size=3.;font.align_x='LEFT'
        text=bpy.data.objects.new(font.body,font);scene.collection.objects.link(text)
        text.location=controller.location+Vector((-65.,70.,8.))
    readme=bpy.data.texts.new('README - Street terrain placement')
    readme.write('Select a source Street, generated mesh or network controller.\n'
        'Open the 3D View sidebar: Streets > Python streets.\n'
        'Choose Preserve profiles (Python) or Direct projection (GN), then press Update.\n'
        'Edit source XY controls or Section properties and press Update.\n'
        'Original terrain stays unchanged. Clear its target and Update to restore construction elevations.\n'
        'Complete Bundle regions currently move at one shared elevation.\n'
        'Use the updated Blosm addon from C:/dev/workspace/blosm.\n')
    camera_data=bpy.data.cameras.new('Terrain inspection camera')
    camera=bpy.data.objects.new('Terrain inspection camera',camera_data);scene.collection.objects.link(camera)
    target=Vector((20.,0.,2.));camera.location=target+Vector((65.,-105.,85.))
    camera.rotation_euler=(target-camera.location).to_track_quat('-Z','Y').to_euler()
    camera_data.type='ORTHO';camera_data.ortho_scale=115.;scene.camera=camera
    scene.render.engine='BLENDER_WORKBENCH';scene.render.resolution_x=1400;scene.render.resolution_y=950
    scene.render.resolution_percentage=100
    shading=scene.display.shading
    shading.light='STUDIO';shading.color_type='MATERIAL';shading.show_shadows=True;shading.show_cavity=True
    shading.cavity_type='BOTH';shading.show_specular_highlight=True
    scene.render.image_settings.file_format='PNG';scene.render.filepath=str((output/'profile_placement.png').resolve())
    bpy.ops.render.render(write_still=True)
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type=='VIEW_3D':
                area.spaces.active.region_3d.view_location=target
                area.spaces.active.region_3d.view_rotation=camera.rotation_euler.to_quaternion()
                area.spaces.active.region_3d.view_distance=130.
                area.spaces.active.shading.color_type='MATERIAL'
    for obj in bpy.context.selected_objects:obj.select_set(False)
    bpy.context.view_layer.objects.active=controllers[0];controllers[0].select_set(True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,
        default=ROOT.parent/'artifacts/streets/terrain_placement/gallery')
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:]);args.output.mkdir(parents=True,exist_ok=True)
    for obj in list(bpy.data.objects):bpy.data.objects.remove(obj,do_unlink=True)
    ui.register()
    from blosm.app.asset_library import BlenderAssetLoader
    from blosm.way.street_network import assets
    loader=BlenderAssetLoader(ROOT/'assets',assets.material_catalog())
    materials={'pavement':loader.load(assets.PAVEMENT),'grass':loader.load(assets.GRASS)}
    networks=(fixture_networks(model)[-1],short_junction_pair(model),*transition_fixture_networks(model)[:2],*fixture_networks(model)[:6])
    controllers=[];report=[]
    for i,network in enumerate(networks):
        network=replace(network,uuid='terrain_'+network.uuid,settings=replace(network.settings,exaggeration=1.))
        offset=(i%3*240.,i//3*240.,0.)
        controller=NetworkScene.create(network,name=network.uuid,generate=False,materials=materials)
        controller.location=offset
        terrain=make_terrain('Original terrain / '+network.uuid,offset)
        controller.blosm_street_terrain.target=terrain
        bpy.context.view_layer.update()
        update=NetworkScene.update(controller)
        report.append(dict(name=network.uuid,placement=json.loads(controller['st_py_terrain_report'])))
        controllers.append(controller)
        print('TERRAIN_FIXTURE',network.uuid,update.timings['published_objects'],flush=True)
    report[0]['lifecycle']=verify(controllers[0])
    report[0]['source_edits']=verify_edits(controllers[0])
    report.append(verify_units_and_transforms())
    source=fixture_networks(model)[-1]
    source=replace(source,uuid='terrain_direct_comparison',settings=replace(source.settings,exaggeration=1.))
    comparison=NetworkScene.create(source,name=source.uuid,generate=False,materials=materials)
    comparison.location=(240.,720.,0.)
    comparison.blosm_street_terrain.target=make_terrain('Original terrain / direct comparison',(240.,720.,0.))
    comparison.blosm_street_terrain.mode='DIRECT';bpy.context.view_layer.update();NetworkScene.update(comparison)
    controllers.append(comparison)
    present_gallery(controllers,args.output)
    path=args.output/'street_terrain_modes.blend'
    bpy.context.view_layer.objects.active=controllers[0];controllers[0].select_set(True)
    bpy.context.preferences.filepaths.file_preview_type='NONE'
    bpy.ops.wm.save_as_mainfile(filepath=str(path.resolve()))
    signatures={c.name:buffers(c) for c in controllers}
    names=[c.name for c in controllers]
    bpy.ops.wm.open_mainfile(filepath=str(path.resolve()));clear_runtime()
    for name in names:
        controller=bpy.data.objects[name]
        update=NetworkScene.update(controller)
        assert update.timings['published_objects']==0,(name,update.timings)
        for key,xyz in buffers(controller).items():np.testing.assert_array_equal(xyz,signatures[name][key])
    report.append(dict(reopen=True,original_terrain_retained=True,fixture_count=len(names)))
    (args.output/'checks.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('TERRAIN_CHECKS_OK',str(path),flush=True)


if __name__=='__main__':main()
