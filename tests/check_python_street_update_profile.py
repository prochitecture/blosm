"""Warm Update timing and cProfile, preserving a real scene's exact inputs."""
import argparse
import cProfile
import json
from pathlib import Path
import pstats
import sys
import time
import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from _python_street_edit import SectionProbeEdit
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime
from blosm.way.street_network.python_backend import backend, scene
from blosm.way.street_network.python_backend import intersections, updates


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--section', type=int)
    comparison = parser.add_mutually_exclusive_group()
    comparison.add_argument('--compare-validation', action='store_true')
    comparison.add_argument('--compare-work', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    controller = max((o for o in bpy.data.objects if o.get('st_py_source_data')),
                     key=lambda c: len(json.loads(c['st_py_source_data'])['sections']))
    before = mesh_buffers(controller)
    clear_runtime()
    start = time.perf_counter()
    NetworkScene.update(controller)
    cold_seconds = time.perf_counter() - start
    assert equal(before, mesh_buffers(controller))
    print(json.dumps(dict(phase='warmed', cold_seconds=cold_seconds)), flush=True)
    edit = SectionProbeEdit(controller, args.section)
    records = []
    batched = backend.validate_overlaps
    frame = intersections.Approach.frame
    order = updates.DependencyIndex._topological_order
    if args.compare_validation:
        from _python_overlap_reference import validate_reference
        modes = [('reference', False), ('batched', False), ('batched', False), ('reference', False), ('batched', True)]
    elif args.compare_work:
        from _python_work_reference import uncached_approach_frame, sorted_list_topological_order
        modes = [('old-work', False), ('reused', False), ('reused', False), ('old-work', False), ('reused', True)]
    else:
        modes = [('batched', False), ('batched', True)]
    reference_edited = None
    for method, profiled in modes:
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
        assert not equal(before, current)
        if reference_edited is None:
            reference_edited = current
        else:
            assert equal(reference_edited, current), 'Measured methods changed published mesh buffers'
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
    backend.validate_overlaps = scene.validate_overlaps = batched
    intersections.Approach.frame = frame
    updates.DependencyIndex._topological_order = order
    counts = [(r['update']['timings']['overlap_check']['candidate_pairs'],
               r['update']['timings']['publication_overlap_check']['candidate_pairs']) for r in records]
    assert len(set(counts)) == 1, 'Reference/batched candidate sets differ'
    assert all(r['changed_geometry'] == records[0]['changed_geometry'] for r in records)
    report = dict(passed=True, scene=str(args.scene.resolve()), section=edit.section.id,
                  cold_seconds=cold_seconds, measurements=records, exact_restoration=True,
                  compared_work_reuse=args.compare_work, exact_edited_buffers_across_methods=True,
                  publication_bounds_reused_in_both_work_methods=args.compare_work)
    (args.output / 'profile_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(dict(passed=True, section=edit.section.id, cold_seconds=cold_seconds,
                         seconds=[r['seconds'] for r in records])), flush=True)


if __name__ == '__main__':
    main()
