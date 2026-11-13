"""Native UV publication/lifecycle checks and a checker-textured street gallery.

Run in factory-startup Blender with --python-exit-code 1. All generated outputs
default to workspace/artifacts/streets/uv_mapping. No material library changes.
"""
import argparse
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys
import time

import bpy
import numpy as np
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT.parent), str(ROOT / 'tests')]
# A worktree need not be named blosm. Resolve this checkout explicitly rather
# than accidentally importing the sibling release checkout.
if 'blosm' not in sys.modules:
    spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
                                                submodule_search_locations=[str(ROOT)])
    package = importlib.util.module_from_spec(spec)
    sys.modules['blosm'] = package
    spec.loader.exec_module(package)

from blosm.way.street_network.python_backend import model, ui
from blosm.way.street_network.python_backend.scene import NetworkScene, SourceObjects, clear_runtime
from _python_street_fixtures import fixture_networks, transition_fixture_networks
from check_python_street_scene import section_control
from check_street_terrain import make_terrain


def stage(name):
    print('UV_STAGE', name, flush=True)


def mesh_snapshot(controller, evaluated=False):
    bpy.context.view_layer.update()
    result = {}
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH':
            continue
        mesh = obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).data if evaluated else obj.data
        xyz = np.empty((len(mesh.vertices), 3), dtype=np.float32)
        mesh.vertices.foreach_get('co', xyz.ravel())
        ids = np.empty(len(mesh.loops), dtype=np.int32)
        mesh.loops.foreach_get('vertex_index', ids)
        assert 'UVMap' in mesh.uv_layers, (controller.name, obj.name)
        layers = {}
        for layer in mesh.uv_layers:
            values = np.empty((len(mesh.loops), 2), dtype=np.float32)
            layer.data.foreach_get('uv', values.ravel())
            assert np.isfinite(values).all(), (obj.name, layer.name)
            layers[layer.name] = values
        result[obj['st_py_key']] = dict(xyz=xyz, ids=ids, layers=layers)
    return result


def compare(a, b, *, geometry=True):
    assert a.keys() == b.keys()
    for key in a:
        if geometry:
            np.testing.assert_array_equal(a[key]['xyz'], b[key]['xyz'], err_msg=key)
            np.testing.assert_array_equal(a[key]['ids'], b[key]['ids'], err_msg=key)
        assert a[key]['layers'].keys() == b[key]['layers'].keys()
        for name in a[key]['layers']:
            np.testing.assert_array_equal(a[key]['layers'][name], b[key]['layers'][name], err_msg=key + '/' + name)


def verify_mapping(controller):
    network = SourceObjects.read(controller)
    source = network.section_index
    snapshot = mesh_snapshot(controller)
    roads, planar, quads = 0, 0, 0
    scale = bpy.context.scene.unit_settings.scale_length
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH':
            continue
        key = obj['st_py_key']
        entry = snapshot[key]
        metadata = json.loads(obj['st_py_metadata'])
        uv = entry['layers']['UVMap']
        quads += len(obj.data.polygons)
        if key.startswith('road:') and 'terrain_stations' in metadata:
            section = source[int(key.split(':')[1])]
            stations = np.asarray(metadata['terrain_stations'])
            across = metadata['transverse_intervals']
            transverse = section.offset + .5 * section.width * np.linspace(-1., 1., across + 1)
            points = np.column_stack((np.repeat(stations, across + 1), np.tile(transverse, len(stations))))
            expected = points[entry['ids']]
            np.testing.assert_allclose(uv, expected, atol=4e-6, rtol=1e-7, err_msg=key)
            roads += 1
        else:
            # Planar charts use controller-local physical XY, so relocating a
            # fixture or changing the scene unit does not move its pattern.
            expected = entry['xyz'][entry['ids'], :2].astype(float) * scale
            np.testing.assert_allclose(uv, expected, atol=1e-5, rtol=1e-7, err_msg=key)
            planar += 1
    return dict(meshes=len(snapshot), roads=roads, planar=planar, quads=quads)


def verify_edits(controller):
    stage('edits / initial snapshot')
    before = mesh_snapshot(controller)
    ids = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    result = NetworkScene.update(controller)
    stage('edits / no-op Update complete')
    assert result.timings['published_objects'] == 0
    assert ids == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    stage('edits / runtime cleared')
    result = NetworkScene.update(controller)
    stage('edits / cold Update complete')
    assert result.timings['published_objects'] == 0
    compare(before, mesh_snapshot(controller))

    source, vertex = section_control(controller, 1)
    stage('edits / control located')
    original = vertex.co.copy()
    vertex.co.y += .15
    source.data.update()
    NetworkScene.update(controller)
    stage('edits / moved-control Update complete')
    edited = mesh_snapshot(controller)
    assert not np.array_equal(edited['road:1']['xyz'], before['road:1']['xyz'])
    verify_mapping(controller)
    NetworkScene.update(controller, force=True)
    stage('edits / moved-control forced rebuild complete')
    compare(edited, mesh_snapshot(controller))
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    stage('edits / restored-control Update complete')
    compare(before, mesh_snapshot(controller))

    attr = source.data.attributes
    widths = [(edge.index, attr['st_width'].data[edge.index].value) for edge in source.data.edges
              if attr['st_section_id'].data[edge.index].value == 1]
    for i, value in widths:
        attr['st_width'].data[i].value = value + .4
    source.data.update()
    NetworkScene.update(controller)
    stage('edits / width Update complete')
    wider = mesh_snapshot(controller)
    verify_mapping(controller)
    assert np.ptp(wider['road:1']['layers']['UVMap'][:, 1]) > np.ptp(before['road:1']['layers']['UVMap'][:, 1])
    NetworkScene.update(controller, force=True)
    stage('edits / width forced rebuild complete')
    compare(wider, mesh_snapshot(controller))
    for i, value in widths:
        attr['st_width'].data[i].value = value
    source.data.update()
    NetworkScene.update(controller)
    stage('edits / restored-width Update complete')
    compare(before, mesh_snapshot(controller))

    exaggeration = controller['st_py_exaggeration']
    controller['st_py_exaggeration'] = 10.
    NetworkScene.update(controller)
    stage('edits / exaggeration Update complete')
    compare(before, mesh_snapshot(controller), geometry=False)
    controller['st_py_exaggeration'] = exaggeration
    NetworkScene.update(controller)
    stage('edits / restored-exaggeration Update complete')
    compare(before, mesh_snapshot(controller))

    resolution = controller['st_py_resolution']
    controller['st_py_resolution'] = 2
    NetworkScene.update(controller)
    stage('edits / R2 Update complete')
    dense = mesh_snapshot(controller)
    verify_mapping(controller)
    for operation in (np.min, np.max):
        np.testing.assert_allclose(operation(before['road:1']['layers']['UVMap'], axis=0),
                                   operation(dense['road:1']['layers']['UVMap'], axis=0), atol=1e-5, rtol=0)
    controller['st_py_resolution'] = resolution
    NetworkScene.update(controller)
    stage('edits / restored-resolution Update complete')
    compare(before, mesh_snapshot(controller))

    return dict(no_op=True, cold_update=True, moved_control=True, width_edit=True,
                full_rebuild_equal=True, restoration_exact=True, exaggeration_independent=True,
                resolution_independent_extent=True)


def verify_terrain(controller):
    stage('terrain / initial snapshot')
    # Exercise recovery of all maps, including discontinuous corner data,
    # rather than merely regenerating UVMap from vertex positions.
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH':
            continue
        primary = np.empty((len(obj.data.loops), 2), np.float32)
        obj.data.uv_layers['UVMap'].uv.foreach_get('vector', primary.ravel())
        secondary = primary * .37 + .11
        secondary[::4, 0] += .2
        layer = obj.data.uv_layers.new(name='FeaturePreview', do_init=False)
        layer.uv.foreach_set('vector', secondary.ravel())
    prepared = mesh_snapshot(controller)
    target = make_terrain('UV test / Original terrain', tuple(controller.location))
    controller.blosm_street_terrain.target = target
    controller.blosm_street_terrain.mode = 'PROFILE'
    bpy.context.view_layer.update()
    NetworkScene.update(controller)
    stage('terrain / PROFILE Update complete')
    placed = mesh_snapshot(controller)
    compare(prepared, placed, geometry=False)
    assert any(not np.array_equal(prepared[k]['xyz'], placed[k]['xyz']) for k in prepared)
    clear_runtime()
    stage('terrain / PROFILE runtime cleared')
    target.location.z += .25
    bpy.context.view_layer.update()
    NetworkScene.update(controller)
    stage('terrain / moved target cold PROFILE Update complete')
    compare(prepared, mesh_snapshot(controller), geometry=False)
    target.location.z -= .25
    bpy.context.view_layer.update()
    NetworkScene.update(controller)
    stage('terrain / restored target Update complete')
    compare(placed, mesh_snapshot(controller))

    controller.blosm_street_terrain.mode = 'DIRECT'
    NetworkScene.update(controller)
    stage('terrain / DIRECT Update complete')
    compare(prepared, mesh_snapshot(controller), geometry=False)
    compare(prepared, mesh_snapshot(controller, evaluated=True), geometry=False)
    stage('terrain / DIRECT evaluated snapshot complete')
    clear_runtime()
    NetworkScene.update(controller)
    stage('terrain / cold DIRECT Update complete')
    compare(prepared, mesh_snapshot(controller, evaluated=True), geometry=False)
    controller.blosm_street_terrain.target = None
    NetworkScene.update(controller)
    stage('terrain / cleared-target Update complete')
    compare(prepared, mesh_snapshot(controller))
    # Save the gallery's bent example using the profile-preserving mode.
    controller.blosm_street_terrain.target = target
    controller.blosm_street_terrain.mode = 'PROFILE'
    NetworkScene.update(controller)
    stage('terrain / final PROFILE Update complete')
    return dict(profile_preserves_uv=True, direct_evaluated_preserves_uv=True,
                named_corner_maps_survive_terrain=True,
                cold_terrain_hydration=True, clearing_restores_prepared_mesh=True)


def verify_units(materials):
    stage('units / create scene')
    previous = bpy.context.scene
    scene = bpy.data.scenes.new('UV centimetre validation')
    bpy.context.window.scene = scene
    scene.unit_settings.scale_length = .01
    network = model.Network('uv_centimetres', (model.Section(1, 1, ((-20., 0., 0.), (0., 0., 0.), (20., 0., 0.)), offset=1.5),),
                            settings=model.Settings(exaggeration=1.))
    controller = NetworkScene.create(network, name=network.uuid, materials=materials)
    stage('units / generated centimetre network')
    controller.location = (200., 300., 0.)
    controller.rotation_euler.z = .3
    result = verify_mapping(controller)
    stage('units / mapping verified')
    for obj in list(scene.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.context.window.scene = previous
    bpy.data.scenes.remove(scene)
    stage('units / temporary scene removed')
    clear_runtime()
    return dict(physical_metres=True, translated_rotated_controller=True, **result)


def verify_named_corner_layers():
    stage('named layers / begin')
    from blosm.way.street_network.python_backend.mesh import Mesh
    from blosm.way.street_network.python_backend.scene import SceneWriter
    xyz = np.array(((0., 0., 0.), (1., 0., 0.), (2., 0., 0.),
                    (0., 1., 0.), (1., 1., 0.), (2., 1., 0.)))
    quads = np.array(((0, 1, 4, 3), (1, 2, 5, 4)), dtype=np.int32)
    primary = xyz[quads, :2].copy()
    primary[1, :, 0] += 10.  # Two UVs at geometric vertices 1 and 4.
    secondary = primary * .25 + .5
    mesh = Mesh(xyz, quads, uv_layers={'UVMap': primary, 'CrossingPreview': secondary})
    published = SceneWriter._mesh('UV corner payload validation', mesh, 1.)
    stage('named layers / published standalone mesh')
    try:
        assert len(published.vertices) == 6 and len(published.polygons) == 2
        for name, expected in mesh.uv_layers.items():
            actual = np.empty((8, 2), dtype=np.float32)
            published.uv_layers[name].data.foreach_get('uv', actual.ravel())
            np.testing.assert_array_equal(actual, expected.reshape(-1, 2).astype(np.float32))
        assert primary[0, 1, 0] != primary[1, 0, 0]
    finally:
        bpy.data.meshes.remove(published)
    stage('named layers / mesh verified and removed')
    return dict(named_layers=True, shared_vertex_chart_seam=True, unchanged_connectivity=True)


def checker_material():
    material = bpy.data.materials.new('UV inspection / metre checker')
    material.use_nodes = True
    material.diffuse_color = (.32, .38, .44, 1.)
    nodes, links = material.node_tree.nodes, material.node_tree.links
    nodes.clear()
    out = nodes.new('ShaderNodeOutputMaterial')
    shader = nodes.new('ShaderNodeBsdfPrincipled')
    shader.inputs['Roughness'].default_value = .8
    uv = nodes.new('ShaderNodeUVMap')
    uv.uv_map = 'UVMap'
    image = bpy.data.images.new('UV inspection / 1 m tile', width=128, height=128)
    y, x = np.mgrid[:128, :128]
    cells = ((x // 64 + y // 64) % 2).astype(bool)
    pixels = np.empty((128, 128, 4), dtype=np.float32)
    pixels[:] = (.66, .73, .77, 1.)
    pixels[cells] = (.12, .21, .27, 1.)
    pixels[(x % 64 < 2) | (y % 64 < 2)] = (.04, .07, .09, 1.)
    image.pixels.foreach_set(pixels.ravel())
    image.update()
    image.pack()
    texture = nodes.new('ShaderNodeTexImage')
    texture.image = image
    texture.interpolation = 'Linear'
    texture.extension = 'REPEAT'
    nodes.active = texture
    uv.location = (-600, 0)
    texture.location = (-200, 0)
    shader.location = (100, 0)
    out.location = (400, 0)
    links.new(uv.outputs['UV'], texture.inputs['Vector'])
    links.new(texture.outputs['Color'], shader.inputs['Base Color'])
    links.new(shader.outputs['BSDF'], out.inputs['Surface'])
    return material


def render_cpu(scene, output, *, samples=12, width=1000):
    """Render the packed UV inspection material without a headless GPU context."""
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'CPU'
    scene.cycles.samples = samples
    scene.cycles.use_denoising = False
    scene.cycles.max_bounces = 2
    scene.cycles.diffuse_bounces = 1
    scene.cycles.glossy_bounces = 1
    scene.render.resolution_x = width
    scene.render.resolution_y = round(width * .8)
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = 'PNG'
    scene.render.filepath = str((output / 'uv_mapping.png').resolve())
    if scene.world is None:
        scene.world = bpy.data.worlds.new('UV inspection world')
    scene.world.use_nodes = True
    background = scene.world.node_tree.nodes.get('Background')
    background.inputs['Color'].default_value = (.16, .2, .25, 1.)
    background.inputs['Strength'].default_value = .65
    sun = bpy.data.objects.get('UV inspection sun')
    if sun is None:
        light = bpy.data.lights.new('UV inspection sun', 'SUN')
        sun = bpy.data.objects.new(light.name, light)
        scene.collection.objects.link(sun)
    sun.data.energy = 2.5
    sun.data.angle = .35
    sun.rotation_euler = (.4, -.5, -.5)
    stage('CPU render / begin')
    bpy.ops.render.render(write_still=True)
    stage('CPU render / complete')


def present(controllers, output, render):
    scene = bpy.context.scene
    for controller in controllers:
        curve = bpy.data.curves.new('UV fixture label', 'FONT')
        curve.body = controller.name.replace('UV / ', '')
        curve.size = 3.
        obj = bpy.data.objects.new(curve.body, curve)
        scene.collection.objects.link(obj)
        obj.location = controller.location + Vector((-48., 53., 1.))
    readme = bpy.data.texts.new('README - Python street UV mapping')
    readme.write('Select a street and inspect its UVMap in the UV Editor.\n'
                 'UV units are physical metres: a 1 m tile contains 0.5 m checker squares.\n'
                 'Roads use station/transverse charts; other surfaces use local XY charts.\n'
                 'Chart seams at intersections and Section boundaries are intentional.\n'
                 'Enable Material Preview to see the packed checker image.\n'
                 'Streets > Python streets > Update regenerates edited input Sections.\n'
                 'The bent example also demonstrates profile-preserving terrain placement.\n')
    camera_data = bpy.data.cameras.new('UV inspection camera')
    camera = bpy.data.objects.new('UV inspection camera', camera_data)
    scene.collection.objects.link(camera)
    target = Vector((17., 2., 2.))
    camera.location = target + Vector((45., -70., 100.))
    camera.rotation_euler = (target - camera.location).to_track_quat('-Z', 'Y').to_euler()
    camera_data.type = 'ORTHO'
    camera_data.ortho_scale = 100.
    scene.camera = camera
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == 'VIEW_3D':
                area.spaces.active.region_3d.view_location = target
                area.spaces.active.region_3d.view_rotation = camera.rotation_euler.to_quaternion()
                area.spaces.active.region_3d.view_distance = 115.
                area.spaces.active.shading.type = 'MATERIAL'
    for obj in bpy.context.selected_objects:
        obj.select_set(False)
    bpy.context.view_layer.objects.active = controllers[0]
    controllers[0].select_set(True)
    if render:
        render_cpu(scene, output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT.parent / 'artifacts/streets/uv_mapping')
    parser.add_argument('--render', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    started = time.perf_counter()
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    ui.register()
    material = checker_material()
    grass = bpy.data.materials.new('UV inspection / grass')
    grass.diffuse_color = (.16, .25, .1, 1.)
    materials = dict(pavement=material, grass=grass)
    basic = fixture_networks(model)
    straight = model.Network('straight_offset',
        (model.Section(1, 1, ((-35., 0., 0.), (0., 0., 0.), (35., 0., 0.)), offset=1.5),))
    networks = (basic[-1], straight, *basic[:6], basic[6], transition_fixture_networks(model)[0])
    controllers, reports = [], []
    for i, network in enumerate(networks):
        network = replace(network, uuid='uv_' + network.uuid,
                          settings=replace(network.settings, exaggeration=1., resolution=1))
        controller = NetworkScene.create(network, name='UV / ' + network.uuid[3:], generate=False,
                                         materials=materials)
        controller.location = (i % 3 * 150., i // 3 * 140., 0.)
        bpy.context.view_layer.update()
        NetworkScene.update(controller)
        reports.append(dict(name=controller.name, **verify_mapping(controller)))
        controllers.append(controller)
        print('UV_FIXTURE', controller.name, reports[-1], flush=True)
    reports[0]['edits'] = verify_edits(controllers[0])
    stage('all edit checks complete')
    reports[0]['terrain'] = verify_terrain(controllers[0])
    stage('all terrain checks complete')
    reports.append(dict(units=verify_units(materials)))
    stage('unit checks complete')
    reports.append(dict(corner_layers=verify_named_corner_layers()))
    stage('named layer checks complete')
    present(controllers, args.output, args.render)
    stage('presentation/render complete')
    expected = {c.name: mesh_snapshot(c) for c in controllers}
    names = [c.name for c in controllers]
    path = (args.output / 'street_uv_mapping.blend').resolve()
    bpy.ops.wm.save_as_mainfile(filepath=str(path))
    stage('scene saved')
    clear_runtime()
    bpy.ops.wm.open_mainfile(filepath=str(path), use_scripts=False)
    stage('scene reopened')
    clear_runtime()
    for name in names:
        controller = bpy.data.objects[name]
        compare(expected[name], mesh_snapshot(controller))
        stage('reopen verified / '+name)
        result = NetworkScene.update(controller)
        assert result.timings['published_objects'] == 0, (name, result.timings)
        compare(expected[name], mesh_snapshot(controller))
    report = dict(passed=True, fixtures=reports, reopened=True,
                  fixture_count=len(names), seconds=time.perf_counter() - started,
                  scene=str(path), rendered=args.render)
    (args.output / 'checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print('UV_CHECKS_OK', json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
