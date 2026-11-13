# Cooperative Python street import

Implementation record, **13 November 2026**, Python street runtime
`0.1.0.dev82`. Development is in the chat's `streets_development` checkout,
synchronized with committed `buildings_gn` at
`a2eab58743a762e6cc04321824e87b4bf8d16d0f`. This document does not record a
transfer into the user's main checkout. **The initial complete-component
integration is qualified by the targeted checks below.** The record keeps its
responsiveness limits explicit; further latency improvements are separate work.

The [shared OSM session](OSM_PROGRESSIVE_IMPORT.md) now drives Python street
preparation, numerical work, terrain fitting, validation and detached object
staging through cooperative producers. A complete connected component remains
the publication unit. The scheduler can return to Blender during construction;
the component's generated street meshes become visible together when that
component is ready. There are no provisional street mesh previews.

The [12 November requirements audit](OSM_PROGRESSIVE_STREETS_REQUIREMENTS.md)
remains the historical design and merge contract. Its descriptions of the
then-current synchronous implementation should not be read as the status of
this implementation. The generic scheduling and ownership protocol is described
in [progressive/README.md](../progressive/README.md).

## Using progressive street import

Enable the existing **Progressive import** option for OpenStreetMap and import
with the Python street backend. The import starts a background timer in the
current Blender window; progress labels describe preparation and numerical work
before a complete component's generated meshes appear. Import options, including
material selection, are captured when the operation starts.

Use **Cancel import** to stop at the next safe checkpoint. Completed components
remain, and the incomplete component's pending allocations are removed. An
individual numerical builder or the final component swap can delay reaching
that checkpoint. After the import finishes or is cancelled, edit the source
controls and press **Update** as before. The Update button remains blocking and
is disabled while an import is active. Scripts that need an immediate complete
result should use blocking import.

## Preparation and dependency order

`WayManager.process_steps(session)` initializes ways in batches and drives each
preparation action in its existing order. `StreetGenerator.do_steps` places
checkpoints between its established network preparation phases. The blocking
methods drain these same producers. Way preparation still runs before consumers
that require it, including building facade analysis when road output is disabled.

`WayManager.render_steps` and the `StreetRenderer` facade delegate to a backend's
optional producer. An older backend without that producer retains the shared
synchronous unit fallback. Backend selection, feature registration, shared
materials and terrain assets remain outside the generic scheduler.

[`way.progressive.producer_for`](../way/progressive.py) also checks the method
resolution order. If a custom subclass overrides only the blocking action or
renderer method, that override runs once as an atomic fallback instead of an
older inherited producer bypassing it. The subclass can opt back into
cooperative dispatch by providing its own matching `*_steps` method.

The Python renderer reads the complete prepared network, resolves physical
ownership, converts to the numerical model, partitions declared connected
components and prepares their retained boundary support. It does not interpret
one OSM way or one Blender object as an independent street construction unit.
Sections connected through intersections, connections or Bundle membership
remain in the same component.

## Producer API

All producers below run on Blender's main thread when used by the scene adapter.
Numerical modules remain usable without Blender. Checkpoints yield `None` and
results are returned through the generator's return value. Callers must close
unfinished producers on cancellation. `on_progress`, where supported, receives
a readable phase string; it does not mark an object as committed.

| Producer | Result and boundary |
| --- | --- |
| `UpdateCoordinator.prepare_steps(network, force=False, on_progress=None)` | Returns a `PreparedUpdate`; checkpoints in source/dependency preparation and between numerical tasks; accepted coordinator state remains unchanged |
| `PythonMeshBackend.prepare_steps(network, force=False, on_progress=None)` | Returns the validated numerical update, including UV preparation and overlap checking |
| `ClippedApproachSupport.prepare_steps(network, on_progress=None)` | Returns the supported source network and report; checkpoints between original-context preparation, junction checks, source extensions and numerical retries |
| `validate_overlaps_steps(...)` | Returns overlap statistics; checkpoints while collecting geometry, indexing and checking candidates |
| `ScenePlacement.create_steps(..., fit_cache=None, ownership=None, on_progress=None)` | Returns raw and placed updates; checkpoints through hydration, terrain fitting, projection preflight and product placement |
| `SceneWriter.publish_steps(controller, network, engine, update, validate=None, on_progress=None, parent_ownership=None)` | Returns the placed update after complete component publication; stages one detached output at a time |
| `NetworkScene.update_steps(controller, force=False, materials=None, validate=None, on_progress=None, parent_ownership=None)` | Reads authored inputs, prepares the engine update, validates revisions and publishes it |
| `PythonStreetRenderer.render_steps(manager, data, session)` | Drives component work and integrates accepted results, diagnostics and visible revision with the OSM session |

Optional arguments in this table are keyword-only where declared by the source.
The existing `prepare`, `publish` and `NetworkScene.update` methods drain the
same algorithms synchronously. The explicit Update button remains a blocking
operator; adding its own timer-driven UI is separate work. Its synchronous
publication hook is preserved for existing extensions and regression tests.

`None` is a soft checkpoint: the shared job may process several before its
time/step budget ends. After accepting a component, the renderer yields
`YieldToUI()` so Blender can display the completed result before starting another
component. Numerical task completion, detached staging and visible publication
are distinct events.

## Publication and cancellation

For each component the scene adapter follows this sequence:

1. Read authored sources and capture the controller/input revision.
2. Prepare dependent numerical products without advancing the accepted engine.
3. Fit and validate the complete terrain corridor, expanding the rebuilt set
   when fitted neighbors or furniture must move together.
4. Validate the actual float32 publication coordinates, quad edge ratios and
   overlap against relevant unchanged output.
5. Allocate detached replacement meshes and objects, including UVs, materials,
   ownership metadata and terrain attributes, with checkpoints between outputs.
6. Revalidate inputs, then link replacements, detach retired objects, apply
   controller/material state and commit the raw numerical update without yielding.
7. Accept the allocation ledger and retire detached old IDs before returning.
   The renderer accepts the component and advances the OSM revision before its
   next yield.

The engine commits `placement.raw`, never the fitted coordinates. This preserves
the prepared elevation and prevents terrain offsets accumulating on later
Updates. Existing source identity, dependencies, UV conventions, material roles,
DIRECT projection modifiers and publication tolerances remain in use.

Each new mesh and object is registered immediately in a private allocation
ledger and, during import, the parent session's ledger. DIRECT terrain assets
are captured only around their synchronous load. No global scene snapshot
spans a yield. Closing a pending producer removes its owned unaccepted IDs;
objects created by the user between checkpoints are not inferred to belong to
the import. Pending tags also participate in the shared history cleanup.

The final swap contains no yield. Old objects remain allocated until metadata,
material assignments and the engine commit succeed, so a precommit exception
can restore their collection links, material slots, controller state and engine
revision. Cancellation before the swap retains the previous accepted output.
On initial import, the parent session removes a cancelled component's controller
and sources while earlier accepted components remain.

Disposal failures after a successful swap do not roll back the accepted result.
Retired objects have already been detached, and failures are recorded in
`update.timings['publication_cleanup_errors']`. This keeps a cleanup error from
turning into deletion of a coherent accepted component.

Expected component construction failures retain the renderer's diagnostic
source/controller state and FAILED report entry. The renderer can continue to
later components. Interrupted imports retain accepted component reports with
partial coverage. Neither a failed diagnostic component nor staged geometry is
counted as a successful generated component.

`StreetInputChanged` distinguishes invalidated inputs or destinations from a
numerical construction failure. Deleted/unlinked component, output or import
collections, changed sources and other stale publication targets abort the
active import. They are not converted into retained FAILED source diagnostics
followed by more construction. Pending allocations are removed and previously
accepted components remain.

## Input guards and retained behavior

The OSM session validates the scene, units, destination and original terrain.
The scene adapter additionally checks controller state, source/output collection
identity, output ownership, material selection and pending ID liveness. Before
commit it rereads authored source geometry and compares the captured network
fingerprint. Every final publication also compares the evaluated terrain
signature, including work inside an OSM session. Blender ID inputs such as the
street material override are captured with the import request, alongside scalar
settings. A competing Update of the same network is rejected.

The street Update, Section properties and material-loading operators are
disabled during an active import, and their execute paths enforce the same
guard. Source controls changed through other Blender operations invalidate
publication when checked; the importer does not overwrite the user's source
edits. The shared driver retains its scene/Edit Mode pause behavior, cancellation
and undo/file-load/disable/window-close handling.

No-op and material-only Updates use the same publication path and retain
unchanged geometry objects. Reopened scenes restore manifest identity and small
product metadata; process-local caches are not required to rebuild correctly.
Failed or cancelled initial preparation releases its new engine cache, and the
renderer releases caches belonging to pending component controllers on teardown.

## Responsiveness limits

The implementation is cooperative. It does not guarantee a maximum callback or
cancellation latency inside an individual operation. In particular:

- An individual numerical task builder still runs synchronously. A difficult
  intersection or Bundle can dominate a callback even when neighboring tasks
  yield frequently.
- Each shared street preparation phase, prepared-network conversion and component
  partition remains a coarse operation. Boundary-support preparation now yields
  between its phases and numerical tasks; one footprint fit or source-extension
  calculation still runs synchronously.
- Large NumPy concatenation/sort calls, evaluated-terrain extraction/BVH setup,
  an individual terrain sample batch and one Blender mesh upload are not
  preempted by the scheduler.
- Source/controller creation and the final multi-object swap/retirement are
  synchronous. Final source validation also reads the component's authored data.
- A dominant connected component has no generated mesh preview while its
  numerical/terrain work is pending. Component-level publication is the current
  milestone; provisional previews and smaller commit units are outside its scope.

Checkpoints preserve existing construction choices and validation. They do not
reduce density, quality checks or terrain samples, parallelize Blender access,
or establish a solver throughput speedup. More checkpoints can add overhead;
measure active callback duration, cancellation latency, first visible component,
elapsed throughput and memory separately. A background run can establish
geometry and transaction behavior but does not establish GUI responsiveness,
navigation or native undo correctness.

## Validation results

The following records were produced on 13 November 2026 with Blender 5.2.0 LTS
in the isolated development profile. They compare blocking and cooperative
execution of the same dev82 runtime. They do not replace the historical dev81
performance qualification or establish a speedup against it.

### Native fixture and transaction checks

[`native_02/checks.json`](../../artifacts/streets/progressive_import_20261113/native_02/checks.json)
records **11 passing tests**. These cover the normal OSM adapter, cancellation
before starting, hidden numerical work and detached staging, preservation of
earlier accepted components and unrelated user data, precommit mesh failure,
scene-switch pausing, unit changes and target collection deletion.

Exact blocking/cooperative comparison covers T, unequal-width four-arm and
bent-two-intersection fixtures in both Exact polyline and Smooth interpolation.
Terrain comparison covers PROFILE and DIRECT at scene scales 1.0 and 0.1,
including unchanged original terrain. Comparison includes stored coordinates,
quad indices, corner normals, smooth flags, face materials, all UV maps, prepared
relief attributes, transforms, source data, output metadata and manifests. Only
the newly allocated network UUID is normalized in source/manifest identity;
timing data is not used as a geometry equality criterion.

[`extended_02/checks.json`](../../artifacts/streets/progressive_import_20261113/extended_02/checks.json)
extends this to **54 fixture/interpolation combinations**: **53 exact mesh
comparisons and one exact rejection comparison**. It includes the junction
families, divided carriageways, Bundle overlaps/contacts, dimensional seams,
attached overlap, clipped source support and compound internal routes, in both
interpolation modes. The `seamed_approach` fixture in Exact polyline mode rejects
its strict 2 m curb-radius fit in both drivers with the same exception and
diagnostic; it is not counted as successful mesh generation. The unchanged
committed reference at `a2eab58` rejects the same input with the same message,
as recorded in
[`committed_rejection/checks.json`](../../artifacts/streets/progressive_import_20261113/committed_rejection/checks.json).
No radius tolerance or geometry rule was relaxed to pass the matrix.

### GUI timer and native history checks

[`gui_01/ui_checks.json`](../../artifacts/streets/progressive_import_20261113/gui_01/ui_checks.json)
records a passing factory-profile GUI run. The unmodified production timer
presented three fixture components incrementally while viewport navigation
continued. Repeated native Undo/Redo preserved exact completed geometry. A
separate staged-cancellation scenario exercised scene-switch pause/resume,
preserved the earlier component and an unrelated user object/text, and repeated
Undo/Redo without restoring incomplete pending output.

For the normal fixture run, the first accepted component took **1.559 s**, the
largest import callback **0.441 s**, and the largest observed heartbeat gap
**0.456 s**. Cancellation at the deliberately held detached-staging checkpoint
took **0.004 s**. That cancellation measurement starts at an available safe
checkpoint; it does not measure interruption of a running numerical builder.
These fixture observations are not city-scale GUI latency guarantees.

### Real GUI timing without simulated navigation

Real-extract GUI timing uses `--no-navigation` to retain the observer heartbeat
without continuously changing the viewport. This is separated from the
navigation stress tests below. The scheduler and runtime policy are unchanged;
the test harness additionally records nested wall-time scopes and source hashes.

[`gui_berlin_04/ui_checks.json`](../../artifacts/streets/progressive_import_20261113/gui_berlin_04/ui_checks.json)
and [`gui_canal_01/ui_checks.json`](../../artifacts/streets/progressive_import_20261113/gui_canal_01/ui_checks.json)
record passing production-timer imports with the window shown, the scoped
performance request enabled, and simulated navigation disabled. Both accepted
two components, observed their intermediate presentation, and matched their
exact background scene hashes.

| Extract/run | First component commit (s) | First scene-visible heartbeat (s) | Import elapsed (s) | Largest callback (s) | Largest heartbeat gap (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Berlin Karl-Marx-Allee / GUI 04 | 77.140 | 77.315 | 77.657 | 3.784 | 3.805 |
| Canal Street / GUI 01 | 9.157 | 9.237 | 9.411 | 4.126 | 4.147 |

The [Berlin resource report](../../artifacts/streets/progressive_import_20261113/gui_berlin_04/process/resources.json)
records 1,125,871,616 resident high-water bytes (1.049 GiB) and 1,578,500,096
maximum sampled private bytes (1.470 GiB). The
[Canal resource report](../../artifacts/streets/progressive_import_20261113/gui_canal_01/process/resources.json)
records 817,848,320 resident high-water bytes (0.762 GiB) and 1,275,846,656
maximum sampled private bytes (1.188 GiB).

| Normal-import wall scope | Berlin GUI 04 (s) | Canal GUI 01 (s) |
| --- | ---: | ---: |
| Cumulative core job steps | 46.828 | 7.056 |
| Cumulative controller ticks | 47.628 | 7.115 |
| Derived interval outside controller ticks | 30.047 | 2.310 |

Core and controller times are nested scopes and must not be summed. The
outside-controller interval includes UI idle/redraw and harness work. These
are wall-time observations, not CPU profiling or a universal latency bound.
Berlin's ordinary GUI run took about **77.7 s**, while its continuous viewport
automation stress runs took **367–387 s**. The production runtime and scheduler
policy were unchanged between them. This distinction is required when reading
the timing record; it is not a claim of a solver optimization or an explanation
of a particular OS scheduling mechanism. The largest indivisible callbacks
still lasted several seconds.

### Berlin continuous-navigation stress checks

[`gui_berlin_01/ui_checks.json`](../../artifacts/streets/progressive_import_20261113/gui_berlin_01/ui_checks.json)
[`gui_berlin_02/ui_checks.json`](../../artifacts/streets/progressive_import_20261113/gui_berlin_02/ui_checks.json)
and [`gui_berlin_03/ui_checks.json`](../../artifacts/streets/progressive_import_20261113/gui_berlin_03/ui_checks.json)
record passing Berlin Karl-Marx-Allee imports using the unchanged production
timer while the test heartbeat continuously increased viewport `view_distance`
by 0.001 and requested redraw, with a 0.01 s heartbeat interval requested. Actual
callback intervals depend on Blender's event loop. These are continuous viewport
automation stress runs, not ordinary GUI throughput runs. All three accepted
two components, observed the intermediate one-component state, and matched the full exact scene hash
from the background parity run. This establishes timer integration and
final-result parity for the real extract, without claiming foreground user
responsiveness.

| Observation | Run 01 | Run 02 | Run 03 |
| --- | ---: | ---: | ---: |
| Window condition | Hidden | Hidden | Shown without activation |
| Import performance request | Disabled | Enabled | Enabled |
| First component commit | 384.351 s | 367.093 s | 366.590 s |
| First visible scene component sampled by heartbeat | 384.972 s | 367.262 s | 366.749 s |
| Import elapsed | 386.775 s | 368.998 s | 367.067 s |
| Largest import callback | 3.179 s | 3.377 s | 3.461 s |
| Largest heartbeat gap | 3.193 s | 3.391 s | 3.476 s |
| Peak resident memory (GiB) | 0.979 | 1.054 | 1.050 |
| Maximum sampled private memory (GiB) | 1.479 | 1.477 | 1.470 |

The visibility observation samples scene objects eligible for display, not
pixels presented to a user. Runs 01 and 02 used hidden windows; run 03 requested
`SW_SHOWNOACTIVATE`, and the child window was observed as visible during live
monitoring. The finalized resource file retains its last sample after exit,
when no child window remains. The process guard
labels its measurements as qualification rather than a performance benchmark.
Enabling the scoped Import performance request did not remove the long elapsed
time: run 02 still took about 369 seconds. Showing the window in run 03 also
left elapsed time near 367 seconds. Window state or the preference alone does
not explain the long stress-run times. The same-extract no-navigation result
above separates normal timer qualification from continuous viewport automation
stress; it does not isolate every rendering, observer or event-loop cost.
No OS scheduling cause, universal foreground latency or solver speedup is inferred.

The [run 01 resource report](../../artifacts/streets/progressive_import_20261113/gui_berlin_01/process/resources.json)
[run 02 resource report](../../artifacts/streets/progressive_import_20261113/gui_berlin_02/process/resources.json)
and [run 03 resource report](../../artifacts/streets/progressive_import_20261113/gui_berlin_03/process/resources.json)
record the 3 GiB process cap, 1 GiB reserve, resident process high-water and
sampled private-memory maximum. Run 01's
[stdout](../../artifacts/streets/progressive_import_20261113/gui_berlin_01/process/stdout.log)
and [stderr](../../artifacts/streets/progressive_import_20261113/gui_berlin_01/process/stderr.log)
and run 02's [stdout](../../artifacts/streets/progressive_import_20261113/gui_berlin_02/process/stdout.log)
and [stderr](../../artifacts/streets/progressive_import_20261113/gui_berlin_02/process/stderr.log),
and run 03's [stdout](../../artifacts/streets/progressive_import_20261113/gui_berlin_03/process/stdout.log)
and [stderr](../../artifacts/streets/progressive_import_20261113/gui_berlin_03/process/stderr.log)
are retained beside their reports.

### Real extract parity and background timing

[`real_01/checks.json`](../../artifacts/streets/progressive_import_20261113/real_01/checks.json)
records exact blocking/cooperative parity for Berlin Karl-Marx-Allee, Canal
Street and Nashville Wyoming Avenue. All three imports succeeded with two
accepted components and no errors. The report retains input and exact-output
SHA-256 values for reproducibility.

| Extract | First component commit (s) | Largest callback (s) | Cooperative elapsed (s) |
| --- | ---: | ---: | ---: |
| Berlin Karl-Marx-Allee | 38.869 | 3.291 | 38.961 |
| Canal Street | 6.273 | 3.958 | 6.303 |
| Nashville Wyoming Avenue | 16.383 | 1.266 | 16.725 |

These are **background** measurements that advance the real session with
`job.step(max_steps=1)`. The first-component value marks when generated streets
are committed and eligible for display; no viewport presentation was timed in
this run. The largest callbacks demonstrate the remaining indivisible work.
They do not establish GUI timer gaps, a cancellation bound or production frame
time on these extracts.

The parity harness disables the Import performance preference, uses street
resolution 1 and 0.7 m cells, and keeps the other recorded fixture/import settings
fixed between the two drivers. Its elapsed times are not directly comparable
with earlier performance runs using different preferences, settings or test
drivers. No solver speedup or memory improvement is inferred from this table.

### Follow-up guards, rollback and existing lifecycle

[`native_03/retry/checks.json`](../../artifacts/streets/progressive_import_20261113/native_03/retry/checks.json)
records **six additional passing methods** covering captured carriageway
material selection, deleted material pointers, source edits after staging,
mixed native building/street parity for both PADL GN and Python buildings,
and invalidation after unlinking each of the still-live component, output and
street-root collections. A failure injected after the real engine commit
restored the last accepted geometry, manifest and engine state; retry succeeded.
The three collection variants run within one method.

[`lifecycle_01/checks.json`](../../artifacts/streets/progressive_import_20261113/lifecycle_01/checks.json)
records a passing rerun of the existing native lifecycle harness with terrain
enabled. It covers native file/cache lifetime, material-only Update retaining
the same mesh IDs and UV data, saved-scene Undo/Redo, publication failure
rollback and retry. These are additional regression checks, not replacements
for the progressive cancellation/history cases above.

The shared `test_osm_scheduler` suite also passed all **three** tests. Its first
attempt could not connect to its loopback fixture under the restricted
execution policy (`WinError 10013`);
the authorized loopback rerun passed. That initial environment restriction was
not a scheduler-code failure.

The focused pure-Python suites passed **22 tests**: five optional-producer
dispatch checks, seven backend registry/lifecycle checks and ten numerical
progressive checks. They cover custom synchronous override fallback, explicit
producer opt-in, unchanged backend dispatch, numerical equivalence, cancellation
and accepted-revision behavior. Their sources are
[`test_python_street_progressive_dispatch.py`](../tests/test_python_street_progressive_dispatch.py),
[`test_street_backends.py`](../tests/test_street_backends.py) and
[`test_python_street_progressive.py`](../tests/test_python_street_progressive.py).

The Berlin GUI observations remain separate from background parity timings and
the smaller GUI fixture. These completed checks qualify the initial integration
milestone: complete-component publication with cooperative construction,
same-runtime geometry parity, cancellation/ownership and native history. They
do not establish a throughput optimization or turn the soft scheduling budget
into a latency bound. Provisional previews and finer commit units are outside
this milestone.

The implementation's targeted checks are
[`tests/test_python_street_progressive.py`](../tests/test_python_street_progressive.py),
[`tests/check_street_progressive.py`](../tests/check_street_progressive.py) and
[`tests/check_street_progressive_ui.py`](../tests/check_street_progressive_ui.py),
alongside the existing publication, lifecycle, UV, terrain and native import
tests listed in the [requirements audit](OSM_PROGRESSIVE_STREETS_REQUIREMENTS.md#112-street-tests-after-the-merge).
Further runs should keep the same-runtime exact comparison and preserve the
separation between geometry/transaction tests and responsiveness measurements.

## Reproducing the checks

The following PowerShell command launches only a new Blender process from this
chat's development checkout. It supplies a child-only environment, factory
startup and disabled automatic script execution. It writes all output under
`artifacts`, enables only the test's checkout inside that process, and leaves
the normal Blender installation and profile unchanged. Other chats must replace
the checkout and profile paths together with their own isolated paths.

Set `suite` to `native`, `extended`, `real`, `lifecycle` or `gui` before running.
The variants use the existing Windows resource guard with a 3 GiB process cap
and 1 GiB system reserve. The `gui` variant intentionally opens a
separate visible factory-profile Blender window for the timer/history test; it
does not attach to an existing Blender session. This visible rerun is a different
window condition from the recorded hidden-window Berlin qualification.

```powershell
@'
import os
from datetime import datetime
from pathlib import Path
import subprocess
import sys

suite = 'native'
gui_osm = None  # Optional absolute Path for a real-extract GUI run.
gui_expected_hash = None  # exact_geometry_sha256 from the background report.
gui_performance_request = False
gui_no_navigation = True  # Ordinary real-extract timing; fixture still navigates.
repo = Path('C:/dev/workspace/blosm_streets_dev')
blender = Path('C:/Program Files/Blender Foundation/Blender 5.2/blender.exe')
guard = Path('C:/dev/workspace/prochitecture_gn/experiments/streets/blosm_integration/bounded_blender.py')
profile = Path('C:/dev/workspace/artifacts/blender_profiles/01a0807c-abda-7822-8855-02127742633c/documentation_reruns')
output = Path('C:/dev/workspace/artifacts/streets/progressive_import_20261113') / (
    'rerun_' + suite + '_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
output.mkdir(parents=True, exist_ok=True)
env = os.environ.copy()
for key, folder in {
    'BLENDER_USER_RESOURCES': 'resources', 'BLENDER_USER_CONFIG': 'config',
    'BLENDER_USER_SCRIPTS': 'scripts', 'BLENDER_USER_EXTENSIONS': 'extensions',
    'BLENDER_USER_DATAFILES': 'datafiles', 'TEMP': 'temp', 'TMP': 'temp',
}.items():
    path = profile / folder
    path.mkdir(parents=True, exist_ok=True)
    env[key] = str(path)
env['PYTHONDONTWRITEBYTECODE'] = '1'
script = repo / 'tests/check_street_progressive.py'
args = ['--artifacts', str(output)]
if suite == 'extended':
    args += ['--only', 'fixture_matrix', '--extended']
elif suite == 'real':
    args += ['--only', 'real_extracts']
    for source in (
        'facade_visibility/berlin_karl_marx_allee.osm',
        'streets/canal_street.osm', 'streets/nashville_wyoming_ave.osm',
    ):
        args += ['--osm', str(Path('C:/dev/workspace/osm_extracts') / source)]
elif suite == 'lifecycle':
    script = repo / 'tests/check_street_performance_lifecycle.py'
    args = ['--output', str(output)]
elif suite == 'gui':
    script = repo / 'tests/check_street_progressive_ui.py'
    if gui_osm is not None:
        args += ['--osm', str(gui_osm)]
        if gui_no_navigation:
            args += ['--no-navigation']
    if gui_expected_hash is not None:
        args += ['--expected-hash', gui_expected_hash]
    if gui_performance_request:
        args += ['--performance-request']
elif suite != 'native':
    raise ValueError('Unknown suite: ' + suite)
guard_args = ['--blender', str(blender), '--script', str(script),
    '--output', str(output / 'process'), '--memory-gib', '3',
    '--reserve-mib', '1024', '--threads', '4', '--numeric-threads', '1', '--', *args]
command = [sys.executable, str(guard), *guard_args]
if suite == 'gui':
    # Use the same guard in a new helper, omitting only --background from its
    # Blender command. The guard's source file and installed Blender stay unchanged.
    gui_guard = '''
import importlib.util, sys
path = sys.argv.pop(1)
spec = importlib.util.spec_from_file_location('street_gui_guard', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
original = module.worker_command
module.worker_command = lambda options, arguments: [
    value for value in original(options, arguments) if value != '--background']
sys.exit(module.main())
'''
    command = [sys.executable, '-c', gui_guard, str(guard), *guard_args]
subprocess.run(command, cwd=repo, env=env, check=True,
               creationflags=subprocess.CREATE_NO_WINDOW)
'@ | & 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' -
```

Read `checks.json` or `ui_checks.json`, not only the process exit code. Resource
evidence is in `process/resources.json`, with `stdout.log` and `stderr.log`
beside it. `--only` accepts a test-method
substring and may be repeated for focused reruns. The real-extract command uses
the existing local OSM files and does not download them.

For a Berlin GUI parity rerun, set `suite = 'gui'`,
`gui_osm = Path('C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm')`
and `gui_expected_hash = '5a9b3400076a78fb1d27a204d86a579920b0ba603f0e55ab07ec60f57098c4aa'`.
These supply the test's `--osm` and `--expected-hash` flags. Real-extract GUI mode
runs the normal production-timer import only; the synthetic staging/Undo scenario
remains the no-`--osm` GUI test. `gui_performance_request = True` adds
`--performance-request` to enable the normal scoped preference inside that test
process. Record the setting and window condition when comparing timings.
The default `gui_no_navigation = True` adds `--no-navigation` for a selected
real extract. Set it to `False` only when intentionally repeating the continuous
viewport automation stress test. The observer heartbeat and scene-visibility
checks still run when navigation is disabled.

The GUI harness writes `wall_timing` and `source_hashes_at_start` into
`ui_checks.json` and a compact `ui_progress.json` while running. Its
`core_job_step`, `osm_job_step` and `controller_tick` scopes measure nested
elapsed wall time, record callback counts/gaps and the longest callbacks, and
count yielded scheduler commands. Derived differences separate wrapper time
from nested job time and time outside controller callbacks. Do not add the
nested totals or label outside-callback time as OS scheduling or CPU cost:
that interval includes UI idle/redraw and test-harness work. Instrumentation
wraps existing calls; it does not change production scheduler budgets or policy.

The numerical and shared scheduler tests run in ordinary Python with the
checkout as the working directory. NumPy is required for the numerical suite;
the scheduler test requires permission to connect to its loopback fixture.

```powershell
Set-Location -LiteralPath C:/dev/workspace/blosm_streets_dev
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' -B -m unittest discover -s tests -p test_python_street_progressive.py
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' -B -m unittest discover -s tests -p test_python_street_progressive_dispatch.py
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' -B -m unittest discover -s tests -p test_street_backends.py
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' -B -m unittest discover -s tests -p test_osm_scheduler.py
```

Generated scenes, logs and reports belong under
`C:/dev/workspace/artifacts/streets/progressive_import_20261113/`. Blender
qualification uses this chat's isolated factory-startup profile; it does not
modify the default Blender installation, normal user profile or main Blosm
checkout.
