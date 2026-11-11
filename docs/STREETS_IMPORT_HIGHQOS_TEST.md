# Temporary performance request during Import

This is the pre-integration experiment. The resulting production preference,
implementation and validation are documented in [Import performance](IMPORT_PERFORMANCE.md).

## Test scope — 11 November 2026

This experiment tests the proposed Windows HighQoS request during the **normal
Blosm Import operator**. It is separate from the earlier diagnostic comparison
of CPU classes and explicit EcoQoS. At the time of this experiment, no
production renderer or GUI default had been changed; the request was opt-in
test code in the isolated checkout.

The development branch is `streets_development`, synchronized to reference
`buildings_gn` at `3a777d98f99d82b92f7072fb9992644ff91b1092`. The reference had
no additional commits before this task. The main checkout and normally
installed Blender/addon are not used for the experiment.

Input scene: `C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend`.
OSM: `C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm`.
Blender 5.2.0 LTS, 3D realistic mode, Python streets only, resolution 1,
0.7 m cells, smooth interpolation, 0.08 m crown, exaggeration 1. The selected
Terrain uses PROFILE placement, 1 m sampling and 0.02 m clearance. Default
bundled pavement/grass materials are used. These are effective Blender RNA
settings, including their float32 representation.

## Matched GUI background protocol

The planned sequence is **system-managed / HighQoS / HighQoS / system-managed**.
Each Import runs in a fresh factory Blender process. One idle factory Blender
with a cube scene stays in front throughout all four measurements; it has a
separate sub-profile and no addon. Only this owned idle window is activated.
The user's existing Blender windows are not used or changed.

Each importing process waits until the idle target's known PID is foreground
before starting. Light visibility-kernel entry/exit observations record the
foreground PID and the importing thread's core class. The host verifier rejects
a run if those observations do not identify the same idle target throughout.
This is stricter than inferring background state from a window title or intent.

Both modes have the same observations. There is no cProfile or per-vertex
callback. Whole Import time includes the scoped request and restoration costs.
Blender startup, source-scene preparation, result saving and independent mesh
verification occur outside the timed interval. Blender's thread limit is 16;
OpenBLAS/OMP/MKL child-process environment limits are 1. No CPU affinity,
process-priority or working-set override is applied. OS/shader caches are not
flushed, and unrelated user activity is not controlled.

The serial host monitor checks each owned importing process against a 3 GiB
resident/private cap and a 1 GiB system reserve. It records memory observations
without forcing paging. The idle target exits at the end of the session; each
importing process saves its result and exits before the next one starts.

## Request and restoration

The test-only `ProcessQoS` context manager reads the existing process policy,
sets the execution-speed control bit and clears its throttling state for
HighQoS. It preserves unrelated policy bits. The original masks are restored
when the Import call returns or raises. It does not change thread/process CPU
priority or the system power plan.

The native Blender cleanup qualification has passed for:

- An exception inside the performance scope and repeated cleanup.
- A HighQoS scope entered from an existing explicit EcoQoS policy, restoring
  that EcoQoS policy before the outer scope restores the original state.
- A real failed `bpy.ops.blosm.import_data()` using an intentionally missing
  OSM file, with exact policy restoration after the raised error.

The timed Import receipts independently record before/during/after policies.
The host verifier requires the original policy after every successful Import,
and unchanged process/thread priority. The separate native mesh checker must
then verify equivalent sources, prepared and terrain-deformed positions,
quad topology, UVs, material names and construction metadata.

## Reproduction

The host entry point is
[`test_street_import_qos.py`](../tools/test_street_import_qos.py). Launch it on
the desktop with a fresh output directory and this chat's isolated profile:

```powershell
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' `
  C:/dev/workspace/blosm_streets_dev/tools/test_street_import_qos.py `
  --scene C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend `
  --osm C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm `
  --output C:/dev/workspace/artifacts/streets/performance/import_highqos_new_run `
  --profile-root C:/dev/workspace/artifacts/blender_profiles/01a0807c-abda-7822-8855-02127742633c
```

If a newly opened importing window takes focus, activate the owned
`focus_target.blend` window. `progress.json` identifies the current case and
whether it is waiting for the target or importing. Do not activate a normal
user window as part of this protocol. The GUI launcher needs desktop access;
a temporary sandboxed command may not leave a visible window.

Relevant source:

- [`open_street_performance_scene.py`](../tests/open_street_performance_scene.py):
  prepares the fixture and scopes the optional request around the real Import.
- [`street_execution_diagnostics.py`](../tests/street_execution_diagnostics.py):
  Windows policy reader and reversible `ProcessQoS` context manager.
- [`check_street_import_qos_cleanup.py`](../tests/check_street_import_qos_cleanup.py):
  real Blender failure/nesting/restoration qualification.
- [`check_street_timing_equivalence.py`](../tests/check_street_timing_equivalence.py):
  independent saved-scene comparison, without executing the generator.

Evidence is under
`C:/dev/workspace/artifacts/streets/performance/import_highqos_20261111/`.
`cleanup/cleanup.json` contains the failure-path qualification. The matched
session is `matched_gui_01`, with per-case native logs, timing receipts,
resource observations and saved results. `summary.json` is written only after
all four cases pass the host checks. Mesh equivalence is a separate gate.
Generated JSON, logs and Blender scenes must not be staged.

## Recorded GUI results

The first matched pair passed every foreground and policy gate:

| Case | Execution policy during Import | Whole Import |
| --- | --- | ---: |
| 1 | System-managed | 118.623 s |
| 2 | Temporary HighQoS | 74.563 s |

This pair shows a **37.1% reduction** in elapsed Import time. Both runs kept
the same owned idle Blender in front at every recorded visibility-kernel
observation; the importing thread was on class-1 cores in both cases.
The read-back process masks changed from `(0, 0)` to `(1, 0)` during HighQoS
and returned to `(0, 0)` afterwards. Process and thread priorities were unchanged.

The third Import completed in 78.352 s and restored its original policy, but
the foreground PID changed during that run. The strict host gate rejected it;
**its timing is excluded from the matched comparison**. The fourth GUI case
was not started. `matched_gui_01/progress.json` therefore remains `FAILED`,
and there is no passing four-run GUI summary. The launcher also requires that
host qualification before allowing `--summarize-only` to publish a summary.

Computer Use was then ended at the user's request. Further repetitions run
without UI control, as separate headless controls; they do not retroactively
complete the rejected GUI sequence. The first pair is one observation per mode,
so the percentage is a local paired result, not a repeatability estimate.

## Headless repeat controls

The noninteractive controls run the same complete scene and OSM through
`bpy.ops.blosm.import_data()` in four fresh background Blender processes,
in system/high/high/system order. Use
[`benchmark_streets.py`](../tools/benchmark_streets.py), with `--operator`,
`--scene`, `--runs 1`, `--threads 16`, `--numeric-threads 1`, `--resolution 1`
and `--qos system` or `--qos high`. Each case has a fresh output directory and
the same isolated profile, 3 GiB process cap and 1 GiB system reserve.
No profiler, coarse execution-diagnostic wrappers or CPU affinity is enabled.

These controls have no visible importing window and no fixed foreground
target. They test the normal operator's output, restoration and execution
times in that different context. The headless helper applies the request before
the timed operator and restores it after the dependency-graph update; its
operator timing excludes the policy setup and cleanup. The GUI pair above
includes both policy setup and cleanup in its timed Import wrapper.

Evidence directories are `headless_01_system`, `headless_02_high`,
`headless_03_high` and `headless_04_system`, each with `run_01/timings.json`,
`run_01/process/resources.json` and a saved scene. Read these as a separate
comparison, without pooling their samples with the GUI pair.

| Order | Policy | Whole Import | Street renderer |
| --- | --- | ---: | ---: |
| 1 | System-managed | 170.862 s | 169.050 s |
| 2 | Temporary HighQoS | 80.913 s | 80.323 s |
| 3 | Temporary HighQoS | 78.651 s | 78.107 s |
| 4 | System-managed | 79.598 s | 79.039 s |

Both HighQoS receipts show the exact original process policy restored after
the measured operation. Every run reports COMPLETE coverage, two ready
networks, 328 meshes and 361,422 quads. No resource guard aborted or forced a
working-set limit. The measured headless peaks were below 0.7 GiB resident
and private memory; the GUI cases peaked below 1.5 GiB private memory.

The HighQoS samples cluster around 79–81 s. The automatic samples range from
79.6 to 170.9 s, and the final automatic sample matches the HighQoS times.
The request therefore has little measured benefit when the automatic policy
already produces fast execution. The headless sequence has two observations
per mode and uncontrolled foreground activity, temperature and cache history;
it does not isolate the cause of every difference between those four runs.

## Output equivalence and completed checks

The independent native Blender checker passed for **all seven completed
Imports**: the two qualified GUI cases, the third GUI case used only as
geometry evidence, and all four headless controls. Each saved result has two
source networks, 328 generated meshes and 361,422 quads. After normalizing the
per-import source UUID, the following compare exactly:

- Source settings and topology.
- Published float32 vertex positions before and after terrain placement.
- Quad face indices, loop UVs and material names.
- Construction metadata, excluding per-import terrain references and bounds.

No generation is rerun during this comparison. The full evidence is
`equivalence/report.json`; its separate resource-guard receipt is
`equivalence/process/resources.json`. `evaluation.json` consolidates the
qualified GUI pair, excluded timing, headless controls and cleanup/equivalence
gates. It explicitly records that the planned four-run GUI matrix did not pass.

The new helpers compile, `git diff --check` passes, and the host summarizer
rejects the failed GUI session rather than publishing a successful matrix
summary. Test outputs remain in artifacts; source code and this report stay
in the isolated development checkout. No production Import code, main-checkout
files, user Blender sessions, normal Blender profile or persistent power
settings were changed.

## Interpretation

Small sample counts can show a local effect and check restoration; they are
not a statistical guarantee across machines, temperatures, workloads or
OS power settings. Read raw times and ranges alongside any medians.
The experiment tests whether the request benefits this Import in the background,
without claiming that window focus forces a particular processor class or
retrospectively identifying the exact policy of the original 3.5-minute run.

The qualified GUI pair supports the request as a promising way to reduce
scheduling-related slowdowns during a user-triggered Import on this machine.
The separate repetitions support checking variability as well as average time.
They do not establish a universal speedup or a guaranteed maximum Import time.
Any production integration should retain an optional, scoped request, exact
restoration and an unsupported-platform fallback. The code in this experiment
was test-only when these measurements were recorded. The later production
integration is described in [Import performance](IMPORT_PERFORMANCE.md).

See [the earlier diagnosis](STREETS_TIMING_DIAGNOSIS.md) for the initial
foreground/background and CPU-policy measurements.
