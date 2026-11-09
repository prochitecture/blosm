"""Qualify the six asphalt materials and build a 160 m street comparison.

Blender 5.2 background; outputs and reports belong in artifacts/streets.
--reopen checks native data in a separate process, with no authoring scripts.
--render renders one variant at a time to keep texture memory bounded.
"""
import argparse
import json
import runpy
from pathlib import Path
import sys
from types import SimpleNamespace

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_street_uv import ROOT, model, ui, verify_mapping
from blosm.app.asset_library import BlenderAssetLoader
from blosm.way.street_network import assets
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime
from blosm.way.street_network.python_backend.renderer import PythonStreetRenderer


def snapshot(controller):
    return {o['st_py_key']: (o.as_pointer(), o.data.as_pointer(),
            np.array([v.co[:] for v in o.data.vertices]),
            np.array([loop.uv[:] for loop in o.data.uv_layers['UVMap'].data]))
            for o in controller['st_py_outputs'].objects if o.type == 'MESH'}


def assert_same(a, b):
    assert a.keys() == b.keys()
    for key in a:
        assert a[key][:2] == b[key][:2], key
        for x, y in zip(a[key][2:], b[key][2:]):
            np.testing.assert_array_equal(x, y)


def audit(library):
    with bpy.data.libraries.load(str(library)) as (data, _):
        assert len(data.materials) == 10, data.materials
        assert len(data.images) == 9, data.images
        assert not data.node_groups and not data.objects and not data.meshes
    loader = BlenderAssetLoader(ROOT / 'assets', assets.material_catalog(str(library)))
    loaded = loader.load_many(assets.ASPHALTS)
    for reference in assets.ASPHALTS:
        material = loaded[reference]
        assert material.asset_data and material['blosm_license'] == 'CC0-1.0'
        nodes = material.node_tree.nodes
        texture_nodes = [n for n in nodes if n.type == 'TEX_IMAGE']
        variant = material['blosm_material_variant']
        assert len(texture_nodes) == (3 if variant == 'PBR' else 1)
        assert len([n for n in nodes if n.type == 'BSDF_PRINCIPLED']) == (variant == 'PBR')
        assert len([n for n in nodes if n.type == 'BSDF_DIFFUSE']) == (variant == 'Diffuse')
        assert nodes['Base color'].image.colorspace_settings.name == 'sRGB'
        tile = tuple(material['blosm_tile_m'])
        scale = nodes['Physical tile size'].inputs['Scale'].default_value
        np.testing.assert_allclose(scale[:2], [1 / tile[0], 1 / tile[1]], rtol=1e-6)
        assert next(n for n in nodes if n.type == 'UVMAP').uv_map == 'UVMap'
        for n in texture_nodes:
            assert n.image.packed_file and tuple(n.image.size) == (2048, 2048)
            assert n.extension == 'REPEAT'
            if n.name != 'Base color':
                assert n.image.colorspace_settings.name == 'Non-Color'
    for i in range(0, 6, 2):
        assert loaded[assets.ASPHALTS[i]].node_tree.nodes['Base color'].image is \
               loaded[assets.ASPHALTS[i+1]].node_tree.nodes['Base color'].image
    assert len({n.image.as_pointer() for m in loaded.values() for n in m.node_tree.nodes if n.type == 'TEX_IMAGE'}) == 9
    return loaded


def camera(target, distance, scale):
    scene = bpy.context.scene
    obj = scene.camera
    if obj is None:
        obj = bpy.data.objects.new('Asphalt comparison camera', bpy.data.cameras.new('Asphalt comparison camera'))
        scene.collection.objects.link(obj)
        scene.camera = obj
    target = Vector(target)
    obj.location = target + Vector(distance)
    obj.rotation_euler = (target - obj.location).to_track_quat('-Z', 'Y').to_euler()
    obj.data.type = 'ORTHO'
    obj.data.ortho_scale = scale
    return obj


def render(args):
    bpy.ops.wm.open_mainfile(filepath=str(args.output / 'asphalt_long_streets.blend'), use_scripts=False)
    controllers = sorted([o for o in bpy.data.objects if o.get('asphalt_fixture')], key=lambda o: o['asphalt_fixture'])
    selected = controllers[args.render]
    active = set(selected['st_py_outputs'].objects)
    for obj in bpy.data.objects:
        if obj.type in ('MESH', 'FONT') and obj not in active:
            obj.hide_render = True
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'CPU'
    scene.cycles.samples = args.samples
    scene.cycles.use_denoising = False
    scene.cycles.max_bounces = 2
    scene.cycles.diffuse_bounces = 1
    scene.cycles.glossy_bounces = 1
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = 'PNG'
    y = selected.location.y
    for view in ('long', 'detail'):
        if view == 'long':
            camera((80, y + 2, 0), (0, -70, 125), 165)
            scene.render.resolution_x, scene.render.resolution_y = 1280, 340
        else:
            camera((47, y + .4, .04), (0, -11, 18), 6)
            scene.render.resolution_x, scene.render.resolution_y = 720, 480
        scene.render.filepath = str(args.output / f'{args.render:02d}_{view}.png')
        bpy.ops.render.render(write_still=True)
        # Release decoded pixels; packed bytes remain in the saved fixture.
        for image in bpy.data.images:
            if image.packed_file:
                image.buffers_free()
    (args.output / f'render_{args.render:02d}.json').write_text(json.dumps(dict(
        passed=True, material=selected['asphalt_material'], samples=args.samples, device='CPU')), encoding='utf8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reopen', action='store_true')
    parser.add_argument('--render', type=int, choices=range(6))
    parser.add_argument('--samples', type=int, default=16)
    parser.add_argument('--regression', choices=('materials', 'registration'))
    args = parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    if args.regression:
        name = 'test_street_material_library.py' if args.regression == 'materials' else 'test_street_backend_registration.py'
        sys.argv = [str(ROOT / 'tests' / name)]
        if args.regression == 'materials':
            sys.argv += ['--', '--asset-root', str(args.library.parent), '--output',
                         str(args.output / 'shared_library_regression.blend')]
        return runpy.run_path(str(ROOT / 'tests' / name), run_name='__main__')
    if args.render is not None:
        return render(args)
    if args.reopen:
        bpy.ops.wm.open_mainfile(filepath=str(args.output / 'asphalt_long_streets.blend'), use_scripts=False)
        meshes = [o for o in bpy.data.objects if o.get('st_py_key') and o.type == 'MESH']
        assert len(meshes) == 6
        assert all(o.data.uv_layers.get('UVMap') and o.data.materials for o in meshes)
        assert not bpy.data.libraries and not bpy.data.node_groups
        assert len({o.data.materials[0].name for o in meshes}) == 6
        images = {n.image for obj in meshes for material in obj.data.materials
                  for n in material.node_tree.nodes if n.type == 'TEX_IMAGE'}
        assert len(images) == 9 and all(i.packed_file for i in images)
        ui.register()
        controller = next(o for o in bpy.data.objects if o.get('asphalt_fixture') == '00')
        before = snapshot(controller)
        controller.blosm_street_pavement = bpy.data.materials[assets.ASPHALTS[1].name]
        result = NetworkScene.update(controller)
        assert result.timings['material_objects_updated'] == 1
        assert_same(before, snapshot(controller))
        (args.output / 'reopen_checks.json').write_text(json.dumps(dict(passed=True, meshes=6)), encoding='utf8')
        print('ASPHALT REOPEN PASSED', flush=True)
        return
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    ui.register()
    loaded = audit(args.library)
    controllers, report = [], []
    for index, reference in enumerate(assets.ASPHALTS):
        network = model.Network(f'asphalt_{index}', (model.Section(1, 1,
            ((0., 0., 0.), (40., 0., 0.), (80., 4., 0.), (120., 4., 0.), (160., 0., 0.)),
            width=7., rise=.08),))
        controller = NetworkScene.create(network, name=reference.name, materials={'pavement': loaded[reference]})
        controller.location.y = index * 14.
        controller['asphalt_fixture'] = f'{index:02d}'
        controller['asphalt_material'] = reference.name
        controllers.append(controller)
        report.append(dict(material=reference.name, length_m=160, **verify_mapping(controller)))
    # Material-only Update must retain exact mesh/UV buffers and object identity.
    controller = controllers[0]
    before = snapshot(controller)
    controller.blosm_street_pavement = loaded[assets.ASPHALTS[1]]
    result = NetworkScene.update(controller)
    assert result.timings['material_objects_updated'] == 1
    assert_same(before, snapshot(controller))
    assert next(iter(controller['st_py_outputs'].objects)).data.materials[0] is loaded[assets.ASPHALTS[1]]
    # Source edits also retain the explicit material override.
    source = next(iter(controller['st_py_sources'].objects))
    source.data.vertices[2].co.y += .05
    source.data.update()
    NetworkScene.update(controller)
    assert all(o.data.materials[0] is loaded[assets.ASPHALTS[1]] for o in controller['st_py_outputs'].objects)
    source.data.vertices[2].co.y -= .05
    source.data.update()
    NetworkScene.update(controller)
    controller.blosm_street_pavement = None
    NetworkScene.update(controller)
    assert next(iter(controller['st_py_outputs'].objects)).data.materials[0] is loaded[assets.ASPHALTS[0]]
    # A fresh importer uses the same dedicated override field.
    app = SimpleNamespace(basePath=str(ROOT), assetLibraryMode='APPEND', streetMaterialFile=str(args.library),
                          streetCarriagewayMaterial=loaded[assets.ASPHALTS[4]])
    renderer = PythonStreetRenderer(app)
    assert renderer.load_materials()['pavement'] is loaded[assets.ASPHALTS[4]]
    # Mark every asset/variant explicitly, without modifying the source textures.
    for index, reference in enumerate(assets.ASPHALTS):
        text = bpy.data.curves.new(reference.name + ' label', 'FONT')
        text.body = reference.name.replace('Prochitecture Streets / Asphalt / ', '') + ' — 160 m'
        text.size = 1.2
        obj = bpy.data.objects.new(text.body, text)
        bpy.context.scene.collection.objects.link(obj)
        obj.location = (0, index * 14 + 5.5, .03)
    scene = bpy.context.scene
    scene.world = bpy.data.worlds.new('Asphalt neutral environment')
    scene.world.use_nodes = True
    background = scene.world.node_tree.nodes.get('Background')
    background.inputs['Color'].default_value = (.6, .65, .72, 1.)
    background.inputs['Strength'].default_value = .6
    light = bpy.data.objects.new('Asphalt comparison sun', bpy.data.lights.new('Asphalt comparison sun', 'SUN'))
    scene.collection.objects.link(light)
    light.rotation_euler = (.5, -.4, -.5)
    light.data.energy = 3.
    light.data.angle = .15
    camera((80, 35, 0), (0, -100, 180), 185)
    scene.view_settings.view_transform = 'AgX'
    scene.unit_settings.system = 'METRIC'
    scene.render.engine = 'CYCLES'
    note = bpy.data.texts.new('START HERE — asphalt material comparison')
    note.write('Six generated 160 m carriageways: New, Worn, Cracked; PBR and Diffuse for each.\n'
        'All source maps are packed. Physical tile sizes: 2.10 m / 2.08 m / 3.00 m.\n'
        'Material Preview or Rendered view shows textures. Source controls remain editable.\n'
        'Select a controller or generated street: Streets > Python streets > Carriageway material, then Update.\n'
        'Blosm import panel: Load asphalt materials, choose Carriageway material, import.\n'
        'Diffuse uses only base color in a Diffuse BSDF; PBR adds roughness and OpenGL tangent normals.\n'
        'Crack repetition is intentionally visible; no stochastic texture blending applied.\n')
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == 'VIEW_3D':
                area.spaces.active.region_3d.view_location = (80, 35, 0)
                area.spaces.active.region_3d.view_distance = 190
                area.spaces.active.shading.type = 'MATERIAL'
    for obj in bpy.context.selected_objects:
        obj.select_set(False)
    controllers[0].select_set(True)
    bpy.context.view_layer.objects.active = controllers[0]
    for image in bpy.data.images:
        if image.packed_file:
            image.buffers_free()
    bpy.ops.wm.save_as_mainfile(filepath=str(args.output / 'asphalt_long_streets.blend'))
    (args.output / 'checks.json').write_text(json.dumps(dict(passed=True, fixtures=report,
        material_only_update=True, source_edit=True, import_override=True), indent=2), encoding='utf8')
    clear_runtime()
    print('ASPHALT MATERIAL CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()
