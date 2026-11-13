# Street generation performance baseline

For current dev78 measurements and the implementation plan toward the requested
10x speedup, see [STREETS_PERFORMANCE_PLAN.md](STREETS_PERFORMANCE_PLAN.md).
The sections below retain the earlier baselines and scheduling investigation.

The controlled follow-up is in
[STREETS_TIMING_DIAGNOSIS.md](STREETS_TIMING_DIAGNOSIS.md). It reproduces large
timing differences through Windows CPU placement and execution policy while
independently checking identical geometry, UVs and construction results. The
measurements below preserve the earlier investigation; their initially
unexplained variability should be read alongside that follow-up.

## Saved-scene comparison — 10–11 November 2026

The original 66.947 s median below is a flat-ground, programmatic renderer
baseline. It does **not** reproduce the selected terrain in
`C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend` and must not be
presented as the expected duration of that scene's GUI Import action. The user
confirmed that the reported import of over 3 minutes 30 seconds was streets only.

The revised benchmark can load that scene with `load_ui=False` and
`use_scripts=False`, register this checkout's property definitions and copy the
effective street settings to the prepared Python street renderer. Other saved
feature switches and the saved general UI mode do not enable buildings, forests
or aviation in this benchmark. The source scene is never saved or overwritten.
The same OSM file and SHA-256 as the original baseline were used.

All comparison runs used revision `3a777d98f99d82b92f7072fb9992644ff91b1092`,
Blender 5.2, resolution 1, effective RNA cell size 0.699999988079071 m, crown rise
0.07999999821186066 m, default materials and the same isolated profile. The
selected original terrain has 1,298 vertices and no modifiers; placement is
PROFILE with 1 m sampling, 0.02 m clearance and a 10% reference-grade limit.
Blender/numerical thread limits remained 16/1, with no working-set throttle.

| Run, in execution order | Street generation | Numerical construction | Publication, including placement | Terrain placement, inside publication |
| --- | ---: | ---: | ---: | ---: |
| Saved scene, terrain enabled, first run | **194.884 s** | 176.285 s | 14.937 s | 9.803 s |
| Saved scene, terrain disabled, control | **66.786 s** | 58.557 s | 5.227 s | 0 s |
| Saved scene, terrain enabled, repeat | **75.027 s** | 59.255 s | 12.672 s | 7.319 s |

OSM preparation took 0.48–0.57 s and the explicit dependency-graph update
0.014–0.020 s. All three completed with COMPLETE coverage, two networks,
328 output meshes and 361,422 quads. The small difference from the old baseline's
361,698 quads accompanies the effective Blender float properties instead of
literal Python `.7` and `.08`; it is not evidence of missing roads.

The slow run's junction-surface tasks totalled **92.039 s**, versus **31.626 s**
in its terrain-enabled repeat. Separator tasks totalled **60.407 s**, versus
**16.812 s**. Terrain fitting accounts for a modest additional cost, **not** the
whole discrepancy. The unprofiled terrain run itself varied from 75 to 195 s.
The slow run consumed 187.625 process CPU seconds for preparation/generation,
versus 74.047 s in the repeat; this was not simply a long idle wait for a UI
redraw. These observations do not identify a particular hardware, OS or
numerical-library cause.

The saved controller input dictionaries match exactly after removing only
their per-import UUIDs. An independent native inspection also compared all
328 outputs' float32 XY, pre-placement height (the stored prepared-Z attribute
where terrain was enabled), face-index arrays and construction metadata,
excluding terrain-reference/datum/bounds metadata. They match exactly across
these three runs. Thus changed footprints, density or accepted construction
results do not account for the large timing difference. This does not by itself
prove that every intermediate numerical operation followed an identical path.

The memory guard completed normally: peak resident memory was approximately
597–598 MiB with terrain and 517 MiB without. No safety-cap termination or
deliberate paging throttle occurred. NumPy was the same bundled 2.3.4 build in
all runs. Its thread environment was also identical; differing thread limits
cannot explain the variation between these particular runs.

### Reproducing and continuing the investigation

Use the command in the Reproduction section below with these additional
arguments, a fresh output directory and `--runs 1`:

```powershell
--scene C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend
```

`--scene` uses its effective street resolution, cell size and other street
properties instead of the hard-coded defaults, including `--resolution`.
Add `--without-terrain` for the corresponding control. `--numeric-threads N`
allows an explicit, recorded numerical-library thread limit (default 1).
This remains a background prepared-renderer measurement, excluding Blender
startup, module setup, assertions, saving, rendering and interactive UI redraw.

For a **separate diagnostic**, add `--diagnostic-profile`. It writes
`numerical.prof` and `numerical_profile.txt` for numerical backend calls only.
Its settings and summary explicitly label profiler instrumentation. Profiler
overhead is included in its durations, so do not mix it with the unprofiled
baseline. The wrappers are restored when the measured import finishes.

Evidence is under
`C:/dev/workspace/artifacts/streets/performance/user_scene_comparison/`:

- `saved_settings.json`, `effective_settings.json`: input scene inspection.
- `terrain_rna_01`, `no_terrain_rna_01`, `terrain_rna_02`: independent unprofiled
  runs, reports, guarded resource logs and generated inspection scenes.
- `mesh_comparison.json`, `mesh_comparison_summary.json`: independent native
  comparison of sources, prepared geometry and construction metadata.
- `*_sources.json`: UUID-normalized source dictionaries extracted from saved
  controller payloads.
- `profile_01`: the separate profiler-instrumented diagnostic described below.

### Numerical diagnostic — 11 November 2026

A separate run with `--diagnostic-profile` completed with COMPLETE coverage,
the same 328 meshes and 361,422 quads. It recorded **68,061,366 function calls**
(68,050,217 primitive calls). The numerical wrapper took 174.981 s under
instrumentation; total generation took 199.567 s. These durations include
substantial profiler overhead and are **not** a fourth unprofiled timing sample.

| Profiled routine | Calls | Cumulative time, including children |
| --- | ---: | ---: |
| Compound quad-chart construction, `partition.py:CompoundQuadChart.build` | 44 | 98.661 s |
| Polygon visibility kernel, `quad_charts.py:visibility_center` | 1,443 | 54.581 s |
| Visibility decomposition, `decomposition.py:VisibilitySplitter.split` | 126 | 52.340 s |
| Fixed-rim quad quality, `quad_charts.py:FixedRimQuality.improve` | 165 | 29.121 s |
| NumPy `linalg.norm`, across the numerical backend | 2,695,572 | 22.010 s |

These are nested call-tree costs and must not be added. They identify repeated
polygon-kernel clipping, partition/decomposition work and quad-quality
relaxation as concrete optimization targets. The profile was collected only
inside numerical preparation, so its costs do not come from terrain raycasts
or viewport redraw. It does not establish the cause of the unprofiled
75-versus-195-second execution-rate variation. Follow-up should time/cache or
reduce repeated visibility/kernel work with fixed inputs, retaining footprint,
quad-quality and surface checks; measure changes using unprofiled repeated runs.

At this checkpoint, the remaining investigation was the large variability in
numerical construction time. The controlled follow-up above addresses that
mechanism. Profile the expensive junction and separator routines, compare operation
counts with their accepted construction metadata, and separate numerical
workload changes from execution-rate changes before choosing an optimization.
Do not attribute the user's entire 3.5-minute duration to terrain, materials,
NumPy threading, memory shortage or viewport work without further evidence.
No runtime generation algorithm has been changed for this investigation, and
the user's main checkout and Blender installation remain untouched.

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
