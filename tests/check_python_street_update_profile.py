"""Warm Update timing and cProfile, preserving a real scene's exact inputs."""
import argparse
import cProfile
import json
import importlib.util
import hashlib
from pathlib import Path
import pstats
import sys
import time
import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
                                             submodule_search_locations=[str(ROOT)])
package = importlib.util.module_from_spec(spec)
sys.modules['blosm'] = package
spec.loader.exec_module(package)
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from _python_street_edit import SectionProbeEdit
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime
from blosm.way.street_network.python_backend import backend, scene
from blosm.way.street_network.python_backend import intersections, updates
from blosm.way.street_network.python_backend.terrain import blender as terrain_backend


def visual_signature(controller):
    result = {}
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH':
            continue
        digest = hashlib.sha256()
        for layer in obj.data.uv_layers:
            value = np.empty((len(obj.data.loops), 2), np.float32)
            layer.uv.foreach_get('vector', value.ravel())
            digest.update(layer.name.encode('utf8') + value.tobytes())
        value = np.empty((len(obj.data.corner_normals), 3), np.float32)
        obj.data.corner_normals.foreach_get('vector', value.ravel())
        digest.update(value.tobytes())
        digest.update(repr([m.name if m else None for m in obj.data.materials]).encode('utf8'))
        result[obj['st_py_key']] = digest.hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--section', type=int)
    parser.add_argument('--edit-kind', choices=('both', 'vertex', 'width'), default='both')
    parser.add_argument('--skip-profile', action='store_true')
    comparison = parser.add_mutually_exclusive_group()
    comparison.add_argument('--compare-validation', action='store_true')
    comparison.add_argument('--compare-work', action='store_true')
    comparison.add_argument('--compare-terrain-cache', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    package.register()
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    controller = max((o for o in bpy.data.objects if o.get('st_py_source_data')),
                     key=lambda c: len(json.loads(c['st_py_source_data'])['sections']))
    before = mesh_buffers(controller)
    visual_before = visual_signature(controller)
    clear_runtime()
    start = time.perf_counter()
    NetworkScene.update(controller)
    cold_seconds = time.perf_counter() - start
    assert equal(before, mesh_buffers(controller))
    assert visual_before == visual_signature(controller)
    print(json.dumps(dict(phase='warmed', cold_seconds=cold_seconds)), flush=True)
    terrain_edit_seconds = []
    if args.compare_terrain_cache:
        # Reopening intentionally has no persistent numerical cache. Exercise
        # a real terrain translation and its restoration to populate the same
        # constraint cache that a fresh Import creates, before paired timings.
        target = controller.blosm_street_terrain.target
        assert target is not None
        # The saved benchmark hides Terrain at the depsgraph level. Blender
        # does not evaluate object transforms while hide_viewport is set.
        # Make this disposable test's terrain evaluable before editing it.
        was_hidden = target.hide_viewport
        target.hide_viewport = False
        bpy.context.view_layer.update()
        position = target.location.copy()
        try:
            for raised in (True, False):
                target.location = position
                if raised:
                    target.location.z += .01
                bpy.context.view_layer.update()
                start = time.perf_counter()
                NetworkScene.update(controller)
                terrain_edit_seconds.append(time.perf_counter() - start)
                print(json.dumps(dict(phase='terrain-cache-warmup', raised=raised,
                    seconds=terrain_edit_seconds[-1])), flush=True)
                if raised:
                    assert not equal(before, mesh_buffers(controller))
                else:
                    assert equal(before, mesh_buffers(controller))
                    assert visual_before == visual_signature(controller)
        finally:
            target.location = position
            bpy.context.view_layer.update()
            target.hide_viewport = was_hidden
            bpy.context.view_layer.update()
    noop_times = []
    for _ in range(3):
        start = time.perf_counter()
        NetworkScene.update(controller)
        noop_times.append(time.perf_counter() - start)
        assert equal(before, mesh_buffers(controller))
        assert visual_before == visual_signature(controller)
    edit = SectionProbeEdit(controller, args.section,
        distance=0. if args.edit_kind == 'width' else .05,
        width=0. if args.edit_kind == 'vertex' else .05)
    records = []
    batched = backend.validate_overlaps
    frame = intersections.Approach.frame
    order = updates.DependencyIndex._topological_order
    terrain_fit = terrain_backend.fit
    def uncached_fit(*args, **kwargs):
        kwargs['cache'] = None
        return terrain_fit(*args, **kwargs)
    if args.compare_validation:
        from _python_overlap_reference import validate_reference
        modes = [('reference', False), ('batched', False), ('batched', False), ('reference', False), ('batched', True)]
    elif args.compare_work:
        from _python_work_reference import uncached_approach_frame, sorted_list_topological_order
        modes = [('old-work', False), ('reused', False), ('reused', False), ('old-work', False), ('reused', True)]
    elif args.compare_terrain_cache:
        modes = [('uncached-terrain', False), ('cached-terrain', False),
                 ('cached-terrain', False), ('uncached-terrain', False)]
    else:
        modes = [('batched', False), ('batched', True)]
    if args.skip_profile:
        modes = [mode for mode in modes if not mode[1]]
    reference_edited = None
    reference_visual = None
    for method, profiled in modes:
        terrain_backend.fit = uncached_fit if method == 'uncached-terrain' else terrain_fit
        backend.validate_overlaps = scene.validate_overlaps = validate_reference if method == 'reference' else batched
        if args.compare_work:
            intersections.Approach.frame = uncached_approach_frame if method == 'old-work' else frame
            updates.DependencyIndex._topological_order = (sorted_list_topological_order if method == 'old-work' else order)
        edit.apply()
        profiler = cProfile.Profile() if profiled else None
        cpu_start = time.process_time()
        start = time.perf_counter()
        if profiler:
            profiler.enable()
        NetworkScene.update(controller)
        if profiler:
            profiler.disable()
        seconds = time.perf_counter() - start
        cpu_seconds = time.process_time() - cpu_start
        current = mesh_buffers(controller)
        current_visual = visual_signature(controller)
        assert not equal(before, current)
        if reference_edited is None:
            reference_edited = current
            reference_visual = current_visual
        else:
            assert equal(reference_edited, current), 'Measured methods changed published mesh buffers'
            assert reference_visual == current_visual, 'Measured methods changed UVs, normals or materials'
        records.append(dict(method=method, profiled=profiled, seconds=seconds, cpu_seconds=cpu_seconds,
                            changed_geometry=[key for key in current if key not in before or
                                not np.array_equal(current[key][0], before[key][0]) or
                                not np.array_equal(current[key][1], before[key][1])],
                            update=json.loads(controller['st_py_last_update'])))
        print(json.dumps(dict(phase='measured', method=method, profiled=profiled,
                              seconds=seconds, cpu_seconds=cpu_seconds)), flush=True)
        if profiler:
            profiler.dump_stats(str(args.output / 'update.prof'))
            with (args.output / 'profile.txt').open('w') as stream:
                pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats('cumulative').print_stats(70)
        edit.restore()
        NetworkScene.update(controller)
        assert equal(before, mesh_buffers(controller))
        assert visual_before == visual_signature(controller)
    backend.validate_overlaps = scene.validate_overlaps = batched
    intersections.Approach.frame = frame
    updates.DependencyIndex._topological_order = order
    terrain_backend.fit = terrain_fit
    counts = [(r['update']['timings']['overlap_check']['candidate_pairs'],
               r['update']['timings']['publication_overlap_check']['candidate_pairs']) for r in records]
    assert len(set(counts)) == 1, 'Reference/batched candidate sets differ'
    assert all(r['changed_geometry'] == records[0]['changed_geometry'] for r in records)
    report = dict(passed=True, scene=str(args.scene.resolve()), section=edit.section.id,
                  checkout=str(ROOT), edit_kind=args.edit_kind, noop_seconds=noop_times,
                  scoped_import_performance=True,
                  cold_seconds=cold_seconds, measurements=records, exact_restoration=True,
                  compared_work_reuse=args.compare_work, exact_edited_buffers_across_methods=True,
                  terrain_cache_abba=args.compare_terrain_cache, exact_uvs_normals_materials=True,
                  terrain_edit_and_restore_seconds=terrain_edit_seconds,
                  publication_bounds_reused_in_both_work_methods=args.compare_work)
    (args.output / 'profile_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(dict(passed=True, section=edit.section.id, cold_seconds=cold_seconds,
                         seconds=[r['seconds'] for r in records])), flush=True)


if __name__ == '__main__':
    from blosm.app.import_performance import ImportPerformance
    with ImportPerformance():
        main()
