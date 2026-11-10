# Street generation performance baseline

## Karl-Marx-Allee — 10 November 2026

Measured the current Python backend in the dedicated `streets_development`
checkout. Before adding benchmark tools, `buildings_gn` and development both
pointed to `c2981e754acfdfb109c2666d29341455bb945f21`; no merge was necessary.
No generation algorithm was modified or optimized for this measurement.

Input: `osm_extracts/facade_visibility/berlin_karl_marx_allee.osm`.
SHA-256: `d22cf78ef0bd8daa2473dea0a85687f5483f1bd8567c2d7d8deebf9e6c2583dc`.

| Fresh process | OSM preparation | Street generation | Dependency graph | Total ready |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.587 s | **66.947 s** | 0.024 s | 67.557 s |
| 2 | 0.438 s | **110.637 s** | 0.020 s | 111.096 s |
| 3 | 0.474 s | **62.235 s** | 0.022 s | 62.732 s |

The median street generation time is **66.947 s**, with an observed range of
**62.235–110.637 s**. These are three samples, not a statistical guarantee.
Median total ready time is **67.557 s**.

Every run reports COMPLETE coverage, two connected networks, 328 generated
meshes and **361,698 quads**. Native checks verify quad faces, UVMap presence and
material assignments after the timed interval. Matching counts do not establish
bit-identical geometry between runs.

### What was measured

Street generation encloses normal renderer construction, preparation, render,
finalization and cleanup. It includes physical-network adaptation, material
loading, input/controller creation, numerical mesh construction, UV generation,
validation, material assignment and Blender mesh publication.

OSM preparation encloses layer creation, local XML parsing and normal manager
processing. An explicit view-layer update is then timed separately. Total ready
is their sum. Blender startup, Python/module/app setup, benchmark assertions,
file saving, rendering and interactive viewport redraw are excluded. This is
not an end-to-end GUI import measurement. The final process saves an inspection
scene after measuring; that save took 0.055 s.

### Settings and machine

- Blender 5.2.0 LTS; backend dev77.
- Intel Core i5-12500H, 16 logical processors, approximately 15.7 GiB visible RAM.
- Resolution 1; cell size 0.7 m; smooth cubic interpolation; 0.08 m crown rise;
  exaggeration 1. No terrain, buildings, render or network downloads.
- Existing library-default pavement/grass materials. Asphalt texture rendering
  and material-preview performance are outside this baseline.
- Blender thread limit 16; OpenBLAS/OMP/MKL child limits 1. Python generation is
  largely serial; the Blender limit is not a claim that it uses 16 cores.
- Fresh process and factory startup per run, using this chat's isolated config,
  user resources, extensions, scripts, datafiles and temporary directories.
- Safety monitoring: 3 GiB allocation/resident cap, 1 GiB system reserve. No
  working-set throttle, CPU affinity, cProfile or guarded-allocation tracing.
  All runs finished without reaching the caps. Maximum resident high-water was
  **544.3 MiB**; maximum sampled private allocation **600.0 MiB**. The latter is
  a sampled value, not a guaranteed process high-water mark.

OS file caches were not flushed. Other user activity and CPU power/clock behavior
were not controlled. The resource monitor labels its own checks as resource
qualification, not performance profiling; wall-clock phase timings come from
the independent benchmark timers.

### Where the time goes

| Measurement | Median |
| --- | ---: |
| Numerical backend, including its UV/global validation | 58.624 s |
| Blender publication | 4.866 s |
| Adaptation, materials, inputs and other renderer work | 3.019 s |
| Junction-surface tasks, inside the numerical backend | 30.899 s |
| Separator tasks, inside the numerical backend | 16.776 s |

These are independently calculated medians and nested measurements: do not add
the task rows to backend/publication totals. Junction surfaces and separators
account for about 71% of generation time in the median-generation run. Simple
carriageway mesh tasks total approximately 0.371 s at their median.

The slow second run merits investigation before selecting an optimization.
For example, task `separator:8` took **4.305 / 17.618 / 3.876 s**, whereas
`separator:14` took **7.275 / 7.002 / 6.729 s**. The cause has not been established.
Do not explain this solely as global machine slowdown or claim that the input
geometry was bit-identical. Compare per-task inputs/iteration behavior and then
profile the expensive junction and separator work in separate diagnostic runs.
Keep unprofiled timing runs distinct from profiler-instrumented results.

## Reproduction

The native entry point is
[`benchmark_python_street_import.py`](../tests/benchmark_python_street_import.py).
It bootstraps this exact checkout, times the unchanged importer and wraps only
two component-level calls for phase separation. Existing per-task timers are
aggregated without adding per-vertex callbacks. Its wrappers are restored after
the measured import.

The host entry point is
[`benchmark_streets.py`](../tools/benchmark_streets.py). Use a **fresh output
directory** and your own chat's profile. Do not run this command in the reserved
main Blosm checkout.

```powershell
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' `
  C:/dev/workspace/blosm_streets_dev/tools/benchmark_streets.py `
  --osm C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm `
  --output C:/dev/workspace/artifacts/streets/performance/karl_marx_allee_new_run `
  --profile-root C:/dev/workspace/artifacts/blender_profiles/01a0807c-abda-7822-8855-02127742633c `
  --runs 3 --threads 16 --resolution 1
```

The host launcher sets process-local profile environment variables, invokes the
existing bounded Blender runner, starts workers serially and saves the last
scene outside the timed interval. No user preferences, installed addons or
Blender installation files are modified. Its `--guard` option can override the
default sibling `prochitecture_gn` runner path. `--summarize-only` reads existing
reports without launching Blender.

Current evidence is under
`C:/dev/workspace/artifacts/streets/performance/karl_marx_allee_20261110/`:

- `summary.json`: aggregate timings, counts and memory observations.
- `run_01` through `run_03`: phase/task timings and native process resource logs.
- `run_03/street_benchmark.blend`: generated inspection scene.

JSON, logs and scenes are development outputs and must not be staged. Preserve
this source guide and the benchmark tools in the repository for future comparisons.
