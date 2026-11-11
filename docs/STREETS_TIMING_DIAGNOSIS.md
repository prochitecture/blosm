# Karl-Marx-Allee timing investigation — 11 November 2026

## Finding and scope

The large timing difference can be reproduced by changing the Windows execution
policy and the processor class running the numerical work. It does not require
different street geometry, sampling or missing output. The Python backend has a
substantial serial workload, so the speed of its executing core matters greatly.

The user confirmed that the import used **3D realistic mode, without buildings**
and that another Blender instance was started while the importing instance was
working. That focus/visibility change is relevant to Windows power management.
It is consistent with the mechanism demonstrated by the tests. The user's exact
historical GUI run was not instrumented. The requested interactive check has
now completed in 69.367 s; a subsequent direct two-instance background test
completed in 119.923 s. Both are described below.

The established finding is **execution-rate sensitivity to CPU placement and
Windows execution policy**, with unchanged calculation counts and geometry.
The fixed-core test reproduces a 232-second import on an efficient core; the
explicit HighQoS/EcoQoS tests reproduce a larger 71-to-372-second range. These
are sufficient causes of large runtimes on this input. They do not retrospectively
identify the exact core, clock or automatic QoS tag in the user's original
210-second run. In particular, do not claim that losing focus invariably moves
Blender to efficient cores: the direct background GUI case remained on
performance-class cores.

Microsoft documents that thread quality of service (QoS) can affect core
selection and processor power management. Focused, visible-but-unfocused and
fully occluded/minimized windows can receive different automatic QoS levels.
An explicit EcoQoS request favors efficient cores and frequencies; an explicit
HighQoS request marks performance-critical work. See
[Windows Quality of Service](https://learn.microsoft.com/en-us/windows/win32/procthread/quality-of-service)
and the official
[SetProcessInformation examples](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-setprocessinformation).

An automatic policy is **not** the same as explicit EcoQoS. A returned
ControlMask/StateMask of zero means system-managed policy, not proof of either
fast or slow execution. Indeed, system-managed runs in this investigation
included both fast and slow imports.

## Controlled tests

Checkout `streets_development`, synchronized to `buildings_gn` commit
`3a777d98f99d82b92f7072fb9992644ff91b1092`. No runtime geometry algorithm was
changed. The tests use Blender 5.2.0 LTS and its bundled NumPy 2.3.4 on the
Intel Core i5-12500H laptop. Windows reports logical CPUs 0–7 as efficiency
class 1 (performance cores, including their hardware threads), and 8–15 as
class 0 (efficiency cores).

Input scene: `C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend`.
Input OSM: `C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm`,
SHA-256 `d22cf78ef0bd8daa2473dea0a85687f5483f1bd8567c2d7d8deebf9e6c2583dc`.
Resolution 1, 0.7 m cell size, smooth cubic interpolation, 0.08 m crown,
PROFILE placement on the scene's original terrain. Effective Blender float
property values, rather than Python decimal literals, are retained.

The normal `bpy.ops.blosm.import_data()` is invoked in 3D realistic mode with
streets enabled and buildings, forests, other vegetation, aviation, railways
and water disabled. The bundled street material library is used. Missing OSM
members and network downloads are disabled. The source scene is loaded with
embedded scripts disabled and is never overwritten. Setup/preferences changes
are confined to the isolated test process and profile.

| Experiment | Street generation | Numerical preparation | Terrain, nested in publication |
| --- | ---: | ---: | ---: |
| Normal Import, system-managed, no instrumentation | 255.581 s | 204.830 s | 24.929 s |
| Fixed CPU 0, coarse diagnostics | 102.990 s | 81.245 s | 10.331 s |
| Fixed CPU 8, coarse diagnostics | 232.058 s | 192.505 s | 17.935 s |
| HighQoS, no affinity, coarse diagnostics | 70.742 s | 56.000 s | 6.801 s |
| System-managed repeat, no affinity, coarse diagnostics | 72.515 s | 57.313 s | 7.000 s |
| EcoQoS, no affinity, coarse diagnostics | 372.155 s | 299.732 s | 36.604 s |
| HighQoS after EcoQoS, no affinity, no instrumentation | 73.706 s | 57.981 s | 7.310 s |

The fixed-core pair uses one Blender thread and one numerical-library thread.
The no-affinity experiments use 16 Blender threads and one numerical-library
thread. Thus the valid fixed-core comparison is CPU 0 versus CPU 8; the valid
policy comparison uses the same no-affinity/thread settings. Diagnostics are
identical within each pair. They count coarse calls, take processor snapshots
and run small calibration kernels; their elapsed times include that overhead.
Do not mix them into the old uninstrumented baseline's median.

All completed runs have COMPLETE coverage, two networks, **328 meshes and
361,422 quads**, with UVs and material assignments. The guard uses a 3 GiB cap
and 1 GiB system reserve; no deliberate working-set throttle or CPU priority
change is applied. CPU affinity experiments affect only the launched test
child. QoS experiments affect only the current test process and restore the
exact original process policy after the timed operation, including on failure.

An independent factory-Blender check loaded all **seven** completed scenes
without running the generator or embedded scripts. The UUID-normalized source
records, published float32 prepared and terrain-deformed positions, face-index
arrays, UV layers, material names and construction-metadata digests match
exactly. Its receipt is `equivalence/report.json` under the evidence root.
Geometry simplification or dropped streets cannot account for these results.

### Same operations, different execution rate

The CPU 0/8 and HighQoS/system-managed comparisons have matching call counts
and array-shape distributions for the measured routines:

| Routine | Calls | CPU 0 wall time | CPU 8 wall time |
| --- | ---: | ---: | ---: |
| Polygon visibility kernel | 1,443 | 15.561 s | 47.380 s |
| Visibility splitting | 126 | 14.177 s | 42.880 s |
| Quad energy evaluation | 32,928 | 6.230 s | 11.870 s |
| Quad gradient evaluation | 6,349 | 6.766 s | 12.060 s |
| Fixed-rim quad quality improvement | 165 | 19.090 s | 35.790 s |
| Compound quad-chart construction | 44 | 30.398 s | 75.280 s |
| Boundary height field | 44 | 3.143 s | 5.620 s |

These are nested measurements and must not be added. The numerical workload
has not gained extra iterations in the slower fixed-core run: the measured
energy/gradient call counts and array-shape distributions match. Its pure
Python and NumPy calibration kernels also execute more slowly.

At the visibility-kernel call boundaries, CPU 0 is class 1 in all 2,886 checks;
CPU 8 is class 0 in all checks. HighQoS has class 1 in all those checks, and its
fast system-managed repeat has class 1 in 2,885/2,886 checks. These are processor
snapshots at entry/exit, not continuous CPU-residency percentages.

Process priority remains NORMAL_PRIORITY_CLASS (32), and thread priority
remains 0. HighQoS is not a high/realtime priority setting and does not bind the
process to a fixed core. The recorded API policy changes from (ControlMask 0,
StateMask 0) to (1, 0), then restores to (0, 0).

The explicit EcoQoS counterexample uses (ControlMask 1, StateMask 1), then
restores (0, 0). It has the same measured call counts and shape distributions as
the HighQoS run. At visibility-kernel entry/exit, only 29/2,886 processor checks
fall on class 1; the remainder fall on efficient cores. The energy-saving policy
can therefore slow even unchanged CPU work substantially. It is a stronger
explicit policy than some automatic background classifications, so its
6-minute duration is not claimed to be the user's exact automatic QoS state.

The subsequent uninstrumented HighQoS import returns to 73.706 s in a fresh
process. This reversible policy comparison, the fixed-core comparison and the
matching operation counts together distinguish execution-rate changes from
additional geometry work. Normal Import overhead outside preparation and the
street renderer is under one second in these examples.

## Reproduction and evidence

Use this chat's checkout and isolated profile, with a fresh output directory:

```powershell
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' `
  C:/dev/workspace/blosm_streets_dev/tools/benchmark_streets.py `
  --osm C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm `
  --scene C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend `
  --operator --execution-diagnostics --qos high `
  --output C:/dev/workspace/artifacts/streets/performance/new_policy_test `
  --profile-root C:/dev/workspace/artifacts/blender_profiles/01a0807c-abda-7822-8855-02127742633c `
  --runs 1 --threads 16
```

- `--qos system` uses automatic policy; `--qos high` and `--qos eco` are explicit,
  reversible **test-process** overrides. They are not addon preferences.
- Fixed-core comparison: omit `--qos`, use `--threads 1 --affinity-cpu 0`, then
  CPU 8 in a separate fresh output directory.
- Omit `--execution-diagnostics` for uninstrumented timing checks. Keep
  `--diagnostic-profile` runs separate; it adds substantial cProfile overhead.
- `--operator` uses the real import operator. Without it, the older prepared
  renderer entry point is measured. Both are identified in the report.

Code is in
[`benchmark_python_street_import.py`](../tests/benchmark_python_street_import.py),
[`benchmark_streets.py`](../tools/benchmark_streets.py) and
[`street_execution_diagnostics.py`](../tests/street_execution_diagnostics.py).
Independent scene comparison is provided by
[`check_street_timing_equivalence.py`](../tests/check_street_timing_equivalence.py).

Evidence root:
`C:/dev/workspace/artifacts/streets/performance/timing_cause_20261111/`.
Each successful experiment contains phase/task timings, native stdout/stderr,
resource-monitor logs and an inspection scene. Instrumented ones additionally
contain `execution_diagnostics.json` with core topology, API policy reads,
calibrations and call/shape/processor distributions. Generated JSON, scenes,
logs and profiler outputs must not be staged.

The successful experiment folders are `operator_02`, `cpu0_diagnostic_02`,
`cpu8_diagnostic_01`, `high_qos_diagnostic_01`, `system_qos_diagnostic_01`,
`eco_qos_diagnostic_01` and `high_qos_plain_01`. `controlled_comparison.json`
records the matching call/shape distributions. `equivalence` contains the
independent seven-scene comparison and resource logs. The failed `operator_01`
was an early benchmark-bootstrap error, before geometry generation; it is not
included in the timing table or comparisons.

## Interactive confirmation

The isolated interactive fixture was opened for the user from this checkout,
independently of the main Blosm folder, using
[`open_street_performance_scene.py`](../tests/open_street_performance_scene.py).
It uses the same factory profile, Blender 5.2, source scene, street settings and
numeric-thread environment as the controlled tests. Import is **not** started
automatically. No QoS, affinity or CPU-priority override is applied.

The prepared copy and launch receipt are in
`C:/dev/workspace/artifacts/streets/performance/timing_cause_20261111/manual_foreground_03/`.
`ready.json` confirms the exact checkout/revision, input path and effective
settings. The visible window title identifies `street_foreground_test.blend`.
The temporary **Street Test** sidebar has **Run street import** and
**Reset test scene** buttons. Reset reloads only this dedicated prepared copy.

Run street import invokes the normal `bpy.ops.blosm.import_data()` and reports
elapsed operator time. Lightweight observations at visibility-kernel entry/exit
record the current core class and the foreground process. Counter keys are
`(is_this_blender_foreground, foreground_pid, logical_cpu, efficiency_class)`.
There are no per-vertex hooks, cProfile, persistent observers or dependency-graph
handlers added by the fixture. Function aliases are restored after the call,
including aliases imported lazily during generation. These observations add
some overhead and are intended to verify foreground/core behavior, not to
replace the uninstrumented baseline.

Each user-triggered run writes `foreground_test_001.json`, then a new numbered
report for subsequent runs. It includes elapsed/CPU time, the normal import
coverage/task report, process power-policy reads and foreground/core counters.
The panel also displays the last elapsed time. The actual user-run result must
be inspected before attributing its timing to a focus change; the prepared file
and a visible window alone are not evidence that Import has run.

The user's completed run is recorded in `manual_foreground_03/foreground_test_001.json`
and its compact `summary.json`. Total Import time was **69.367 s**, with
67.484 process CPU seconds. The renderer reports 67.562 s; numerical preparation
totals 53.526 s and terrain placement 6.486 s. Coverage is COMPLETE, with two
READY networks and 361,422 quads. All 2,886 visibility-kernel processor checks
report performance-class cores. The process policy remains system managed,
ControlMask/StateMask (0, 0), before and after the test. There was no performance
policy override in this run.

The raw foreground PID was Blender's own PID before Import, but PID 1964
throughout the kernel observations. A read-only process query verified PID 1964
as `dwm`. These flags must **not** be interpreted as proof that the user switched
away. Windows can replace an unresponsive top-level window with a ghost window
while its main thread is busy; see
[Microsoft's message-queue documentation](https://learn.microsoft.com/en-us/windows/win32/winmsg/about-messages-and-message-queues).
The DWM observations are consistent with that proxy behavior; this test did not
record the proxy's window class. The importing thread's core-class readings are
independent of the foreground-window proxy and remain usable. Raw receipts are
preserved rather than relabelled as a continuous focus measurement.

### Direct two-instance background comparison

`interactive_background_01` is a separate GUI experiment, not an explicit
HighQoS/EcoQoS or affinity experiment. The normal importing instance uses this
checkout and the same prepared settings. An idle factory Blender instance uses
its own sub-profile and a default cube scene, with no addon or generation.
The idle instance was activated before the delayed Import began and remained
in front. Its known PID, 33500, appears in **all 2,886** recorded foreground
checks, independently confirming which process was in front.

| GUI measurement | User-triggered case | Behind the idle Blender instance |
| --- | ---: | ---: |
| Whole Import operator | 69.367 s | 119.923 s |
| Process CPU time | 67.484 s | 119.078 s |
| Renderer report | 67.562 s | 118.641 s |
| Numerical preparation | 53.526 s | 94.032 s |
| Terrain placement | 6.486 s | 10.901 s |
| Coverage / quads | COMPLETE / 361,422 | COMPLETE / 361,422 |

The main importing thread's sampled cores are class 1 in both cases. No
continuous effective-frequency measurement was taken, so no particular clock
speed can be inferred from these receipts. The background case changed both
focus and the presence of an idle instance; it is not proof that focus alone
must produce exactly this multiplier. Its CPU time grows with its wall time,
which is consistent with slower execution rather than a large idle wait.

The background result was saved as `interactive_background_01/main/street_result.blend`.
A further independent native comparison against `high_qos_plain_01` passes for
sources, all prepared and terrain-deformed vertex positions, quad topology,
UVs, material names and construction metadata. By comparison with the earlier
seven-scene result, **all eight saved scenes** are equivalent. The user's
unsaved foreground scene is supported by its complete import receipt and quad
counts, not an independent saved-scene hash.

The importing instance and idle target closed themselves after recording and
saving the result. Their native logs confirm normal exits, and the returned
window inventory confirmed that those temporary windows were gone. The user's
other applications were not closed or reconfigured. Source automation uses the
explicit `--auto-run-after 30 --save-result --quit-after-run` flags on the
interactive harness, plus
[`street_performance_focus_target.py`](../tests/street_performance_focus_target.py).
The default interactive harness still waits for the user to click Import.

Two launch issues were resolved before handoff. Changing sidebar visibility
on an inactive screen during CLI startup caused a native editor-context crash;
setup now runs once after UI initialization and only with a live window/area
context. The sandboxed command did not leave a visible desktop window, so the
authorized interactive launch was performed outside that temporary command
sandbox, retaining the same isolated environment. The resulting window and
button were visually verified. Neither issue occurred during street generation
or involved the user's existing Blender processes or installation.

## Consequence for future work

The previous 67-second flat-ground figure describes a fast execution case;
it does not promise the same duration under other Windows execution conditions.
The interactive check confirms that approximately 69 seconds is achievable
with the selected terrain too. The measured background GUI case is slower,
while the fixed-core and power-policy controls demonstrate substantially larger
slowdowns without extra geometry work. Report foreground/background context and
execution-policy constraints alongside future benchmark numbers.

A possible addon improvement is a scoped Windows HighQoS request during an
explicit user-triggered generation, restoring the original policy in `finally`.
It should be optional/policy-aware, preserve all unrelated policy bits and have
a no-op fallback on unsupported systems. Persistently pinning Blender to one
processor is not appropriate: core numbering and available resources vary.
No such runtime or GUI behavior change has been made as part of this diagnosis.
The user's main checkout, Blender installation, existing Blender processes,
power plan and persistent OS settings remain unchanged.

The follow-up [temporary Import performance-request experiment](STREETS_IMPORT_HIGHQOS_TEST.md)
records the scoped HighQoS test, policy cleanup checks, a qualified GUI pair,
the rejected GUI repeat and separate noninteractive controls.

Algorithmic optimization remains useful independently of scheduling. Prioritize
the repeated polygon visibility/kernel clipping and quad-quality work identified
in the separate numerical profile, retaining the existing footprint, profile,
quality, UV and material checks. See [the original measurements](STREETS_PERFORMANCE.md).
