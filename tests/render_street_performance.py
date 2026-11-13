"""Controlled saved-output renders; never run a generator or execute scene code.

Use the same script/settings for the reference and candidate in separate factory
processes. Copies are used for inspection exaggeration and shader substitutions;
the input scenes, their materials and runtime asset library are read-only.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector


def material(name, color, roughness=.42):
    result = bpy.data.materials.new(name)
    result.use_nodes = True
    shader = result.node_tree.nodes.get('Principled BSDF')
    shader.inputs['Base Color'].default_value = (*color, 1.)
    shader.inputs['Roughness'].default_value = roughness
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--material-library', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=32)
    parser.add_argument('--shot', action='append')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    old_scene = bpy.context.scene
    sources = sorted((o for o in old_scene.objects if o.type == 'MESH' and o.get('st_py_key')),
                     key=lambda o: o['st_py_key'])
    assert sources
    scene = bpy.data.scenes.new('Controlled street comparison')
    bpy.context.window.scene = scene
    scene.unit_settings.scale_length = old_scene.unit_settings.scale_length
    objects, original, relief = {}, {}, {}
    neutral = material('Comparison / neutral pavement', (.24, .29, .32))
    grass = material('Comparison / grass', (.20, .32, .12), .85)
    wire = material('Comparison / wire', (.005, .008, .01), .8)
    for source in sources:
        key = source['st_py_key']
        obj = bpy.data.objects.new('Compared / ' + key, source.data.copy())
        obj.matrix_world = source.matrix_world.copy()
        scene.collection.objects.link(obj)
        objects[key] = obj
        xyz = np.empty((len(obj.data.vertices), 3), np.float32)
        obj.data.vertices.foreach_get('co', xyz.ravel())
        original[key] = xyz
        attr = obj.data.attributes.get('st_terrain_prepared_z')
        dz = np.zeros(len(xyz), np.float32)
        if attr:
            attr.data.foreach_get('value', dz)
            metadata = json.loads(source.get('st_py_metadata', '{}'))
            dz -= float(metadata.get('datum', 0.))
        relief[key] = dz
        obj.data.materials.clear()
        obj.data.materials.append(grass if key.startswith(('separator:', 'junction_grass:')) else neutral)
        obj.data.materials.append(wire)
    # The benchmark's original terrain provides the same visible contact and
    # occlusion in both runs; unrelated buildings and source lines are omitted.
    terrain = old_scene.objects.get('Terrain')
    if terrain is not None and terrain.type == 'MESH':
        ground = bpy.data.objects.new('Comparison / original terrain', terrain.data.copy())
        ground.matrix_world = terrain.matrix_world.copy()
        scene.collection.objects.link(ground)
        ground.data.materials.clear()
        ground.data.materials.append(material('Comparison / terrain', (.065, .078, .085), .95))
    # Every copied object is ordered by its stable source/output key, avoiding
    # random UUID/object insertion differences in deterministic sampling.
    bpy.context.view_layer.update()
    names = [f'Prochitecture Streets / Asphalt / Worn / {v}' for v in ('PBR', 'Diffuse')]
    with bpy.data.libraries.load(str(args.material_library.resolve()), link=False) as (data, target):
        assert set(names) <= set(data.materials)
        target.materials = names
    asphalt = dict(zip(('pbr', 'diffuse'), target.materials))
    assert all(n.image.packed_file for n in asphalt['pbr'].node_tree.nodes if n.type == 'TEX_IMAGE')
    scene.world = bpy.data.worlds.new('Comparison / world')
    scene.world.use_nodes = True
    scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (.15, .17, .20, 1.)
    scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value = .6
    sun = bpy.data.objects.new('Comparison / sun', bpy.data.lights.new('Comparison / sun', 'SUN'))
    scene.collection.objects.link(sun)
    sun.data.energy = 2.5
    sun.data.angle = .08
    sun.rotation_euler = Vector((-.45, -.6, -1.)).to_track_quat('-Z', 'Y').to_euler()
    camera = bpy.data.objects.new('Comparison / camera', bpy.data.cameras.new('Comparison / camera'))
    scene.collection.objects.link(camera)
    scene.camera = camera
    camera.data.type = 'ORTHO'
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'CPU'
    scene.cycles.samples = args.samples
    scene.cycles.seed = 1701
    scene.cycles.use_animated_seed = False
    scene.cycles.use_denoising = False
    scene.cycles.use_adaptive_sampling = False
    scene.cycles.max_bounces = 2
    scene.cycles.diffuse_bounces = 1
    scene.cycles.glossy_bounces = 1
    scene.render.resolution_x, scene.render.resolution_y = 1200, 800
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_mode = 'RGBA'
    scene.render.image_settings.color_depth = '8'
    scene.view_settings.view_transform = 'AgX'
    scene.view_settings.look = 'AgX - Medium High Contrast'
    scene.view_settings.exposure = 0.
    scene.view_settings.gamma = 1.
    shots = (
        ('overview', None, 1., 'neutral', (0., -.3, 1.), 1.1),
        ('junction_oblique', 'junction_mesh:155', 1., 'neutral', (.45, -.75, .8), 1.3),
        ('junction_relief', 'junction_mesh:155', 10., 'neutral', (.45, -.75, .6), 1.3),
        ('junction_side', 'junction_mesh:155', 10., 'neutral', (.2, -1., .08), 1.3),
        ('separator', 'separator:14', 10., 'neutral', (0., -.2, 1.), 1.2),
        ('junction_wire', 'junction_mesh:19', 10., 'wire', (0., 0., 1.), 1.4),
        ('asphalt_pbr', 'junction_mesh:19', 1., 'pbr', (.3, -.7, .7), .35),
        ('asphalt_diffuse', 'junction_mesh:19', 1., 'diffuse', (.3, -.7, .7), .35),
    )
    all_points = np.array([obj.matrix_world @ Vector(corner)
                          for obj in objects.values() for corner in obj.bound_box])
    records = []
    for label, key, exaggeration, shader, direction, margin in shots:
        if args.shot and label not in args.shot:
            continue
        if key is not None:
            assert key in objects, key
        points = all_points if key is None else np.array([
            objects[key].matrix_world @ Vector(corner) for corner in objects[key].bound_box])
        target = Vector((points.min(axis=0) + points.max(axis=0)) * .5)
        width = float(np.ptp(points[:, :2], axis=0).max()) * margin
        for token, obj in objects.items():
            xyz = original[token].copy()
            xyz[:, 2] += (exaggeration - 1.) * relief[token]
            obj.data.vertices.foreach_set('co', xyz.ravel())
            obj.data.update()
            is_grass = token.startswith(('separator:', 'junction_grass:'))
            obj.data.materials[0] = grass if is_grass else asphalt.get(shader, neutral)
        overlay = None
        if shader == 'wire':
            overlay = objects[key].modifiers.new('Inspection edges only', 'WIREFRAME')
            overlay.thickness = .025
            overlay.use_replace = False
            overlay.use_even_offset = True
            overlay.material_offset = 1
        camera.location = target + Vector(direction).normalized() * max(width, 100.)
        camera.rotation_euler = (target - camera.location).to_track_quat('-Z', 'Y').to_euler()
        camera.data.ortho_scale = width
        camera.data.clip_end = max(10000., width * 10.)
        scene.render.filepath = str(args.output / (label + '.png'))
        bpy.context.view_layer.update()
        bpy.ops.render.render(write_still=True)
        if overlay is not None:
            objects[key].modifiers.remove(overlay)
        records.append(dict(name=label, key=key, exaggeration=exaggeration,
            shader=shader, target=list(target), scale=width,
            camera_location=list(camera.location), image=label+'.png'))
        print('STREET_VISUAL', label, flush=True)
    (args.output / 'renders.json').write_text(json.dumps(dict(scene=str(args.scene.resolve()),
        samples=args.samples, seed=1701, mesh_count=len(objects), shots=records), indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
