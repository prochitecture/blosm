"""Benchmark local tile import; network and scene cleanup are excluded.

Blender --background --factory-startup --python-exit-code 1 --python this_file --
    --artifacts OUTPUT_DIR --helsinki TILE_DIR --repeats 20

Run one Blender process at a time. Only benchmark fixtures and the historical
operator baseline write tile files. The new production renderer never does.
"""
import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import random
import statistics
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gltf_test_support import (
    bpy, BASELINE_REF, BlenderRenderer, fixture, convert_b3dm, reset_scene,
    operator_import, baseline_renderer, manager_for, install_scene_settings, import_gltf,
)

parser = argparse.ArgumentParser()
parser.add_argument("--artifacts", type=Path, required=True)
parser.add_argument("--helsinki", type=Path, required=True)
parser.add_argument("--repeats", type=int, default=20)
parser.add_argument("--warmups", type=int, default=3)
args = parser.parse_args(sys.argv[sys.argv.index("--")+1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
install_scene_settings()
ReleaseRenderer = baseline_renderer()
rng = random.Random(241104)
rows, samples, summaries = [], [], []


def run_group(label, operations):
    # Warm every path; then shuffle order within each round, retaining pairing.
    for _ in range(args.warmups):
        for operation in operations.values():
            reset_scene()
            with redirect_stdout(io.StringIO()):
                operation()
    for iteration in range(args.repeats):
        order = list(operations)
        rng.shuffle(order)
        for name in order:
            reset_scene()
            with redirect_stdout(io.StringIO()):
                started = time.perf_counter()
                operations[name]()
                elapsed = time.perf_counter() - started
            rows.append({"sample": label, "iteration": iteration, "path": name, "seconds": elapsed})
    group = {}
    for name in operations:
        values = [r["seconds"] * 1000 for r in rows if r["sample"] == label and r["path"] == name]
        quartiles = statistics.quantiles(values, n=4, method="inclusive")
        group[name] = {"median_ms": statistics.median(values),
                       "q1_ms": quartiles[0], "q3_ms": quartiles[2]}
    summaries.append({"sample": label, "paths": group})
    print("BENCHMARK_GROUP", json.dumps(summaries[-1]), flush=True)


def glb_operations(content, source_file):
    temp = args.artifacts / "operator_roundtrip.glb"
    def roundtrip():
        try:
            temp.write_bytes(content)
            operator_import(temp)
        finally:
            temp.unlink(missing_ok=True)
    return {
        "operator_file": lambda: operator_import(source_file),
        "direct_file": lambda: import_gltf(source_file),
        "operator_roundtrip": roundtrip,
        "direct_memory": lambda: import_gltf(source_file, content=content),
    }


def renderer_operation(renderer_type, content):
    def operation():
        renderer = renderer_type("Benchmark", False, "")
        manager = manager_for(renderer, args.artifacts, content)
        renderer.prepare(manager)
        try:
            if renderer_type is BlenderRenderer:
                renderer.renderB3dm(manager, "fixture", "tile.b3dm", None)
            else:
                renderer.renderB3dm(manager, "fixture", "tile.b3dm", None, False)
            assert renderer.num_imported_tiles == 1
        finally:
            renderer.finalize(manager)
    return operation


doc, binary, synthetic = fixture()
synthetic_file = args.artifacts / "synthetic.glb"
synthetic_file.write_bytes(synthetic)
run_group("synthetic_textured_glb", glb_operations(synthetic, synthetic_file))

for tile in sorted(args.helsinki.glob("*.b3dm")):
    content = tile.read_bytes()
    extracted, rtc = convert_b3dm(content)
    glb = bytes(extracted)
    source_file = args.artifacts / (tile.stem + ".glb")
    source_file.write_bytes(glb)
    samples.append({
        "tile": tile.name, "b3dm_bytes": len(content), "glb_bytes": len(glb),
        "sha256": hashlib.sha256(content).hexdigest(), "rtc_center": rtc,
        "source": "https://kartta.hel.fi/3d/mesh/Helsinki_2024/1/" + tile.name,
    })
    run_group(tile.stem + "_glb", glb_operations(glb, source_file))
    run_group(tile.stem + "_b3dm_renderer", {
        "release_operator": renderer_operation(ReleaseRenderer, content),
        "direct_memory": renderer_operation(BlenderRenderer, content),
    })

report = {
    "blender": bpy.app.version_string,
    "blender_build_hash": bpy.app.build_hash.decode(),
    "baseline_ref": BASELINE_REF,
    "repeats": args.repeats, "warmups_per_path": args.warmups,
    "timing": "perf_counter wall time; warmed paths; deterministic shuffled order; cleanup and network excluded",
    "scope": "GLB pipeline isolation plus full B3DM renderer, fresh scene per trial, background Blender",
    "dataset": "City of Helsinki, Helsinki_2024; CC BY-SA 4.0 as supplied by the user",
    "samples": samples, "summary": summaries, "raw": rows,
}
(args.artifacts / "benchmark.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print("BENCHMARK_COMPLETE", bpy.app.version_string, flush=True)
