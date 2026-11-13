# Progressive OpenStreetMap import

Implementation on `progressive_import`, based on committed `release`
`1bf4bb2a68d4f5c9546b282cec59e6857eb7d163`. The design audit and larger integration
plan remain in [OSM_PROGRESSIVE_IMPORT_PLAN.md](OSM_PROGRESSIVE_IMPORT_PLAN.md).

## Behavior

Progressive import is enabled by default. Its checkbox appears above the data
type selector for both OSM and 3D Tiles. With OSM progressive import selected,
the **Import as a single object** checkbox is hidden and its effective setting
is false before managers, actions or layers are constructed. The saved preference
is preserved for the next blocking import.

Objects are completed in source/manager order. A building includes its parts
and local dome/onion roof joins before its commit. Road features can contain
several splines. 2D Geometry Nodes creates a modifier and face attributes for
each separate footprint object, sharing the loaded node group.

The scene remains usable between callbacks. Individual Blender operations,
native roof solvers, asset loading and GN/dependency-graph evaluation are
synchronous and can exceed the scheduler's soft time budget. This is cooperative
importing, not preemption of Blender or guaranteed frame time.

Realistic forests use separate projected area brushes and a shared scattering
system per imported layer. Each accepted patch updates cumulative area/density;
rounding is applied to total area, and particle settings are copied before
changing a prior/shared asset. Terrain changes made by a failed patch are rolled
back without removing earlier accepted patches. Source brushes are hidden, as
in the existing forest representation; the scattered forest is the visible result.

Uncached Import for export textures render in a private Blender process. A
building remains staged/hidden until its textures are complete and loaded.
Compositor scene snapshots, temporary images/files, process/profile lifetime and
texture publication belong to the import. Completed texture cache files remain
available after cancellation/Undo; Blender history does not undo external files.
Snapshot serialization and material preparation still run synchronously.

Cancel stops pending work and retains accepted results. Renaming/deleting an
accepted object does not cause it to be recreated. Deleting/unlinking the import
collection stops with an error. Terrain geometry/transform changes invalidate
pending placement; already completed objects and the user's terrain edits remain.
Importing does not temporarily move the user's terrain origin.

One Undo removes uninterrupted importing. Scene edits during the import keep
their own undo steps and divide importing into consecutive groups. Undo/Redo
cancels an active job; Redo restores stored scene data rather than downloading
again. Pending staged IDs are explicitly tagged and removed from restored
undo/save snapshots, so incomplete staged objects do not become canonical output.

Edit Mode and scene/view-layer switching pause commits. File loading, disabling
the addon and closing the owning window stop the job. There is one terminal
result report, including cancellation and retained partial results. Failure
messages include the original Blosm error explanation in both the progressive
popup/operator report and the blocking report. Missing required folders, input
files or setup scripts show only the original message, without an import-results
summary. These prerequisites are checked before enabling the progress panel,
registering the import timer or adding an undo checkpoint. Initialization that
allocates Blender data remains scheduled. Failures after preflight (including
parsing failures before the first object) retain the summary and original
explanation. Another
Blosm import or legacy modal imagery import cannot replace the active app state.

## Acquisition and settings

Main OSM and missing-member queries use the shared HTTP helper. Missing members
use bounded read-only POST batches and an atomically published supplementary OSM
file, preserving reuse of local companion files. XML is read in chunks and
decoded on the main thread; it no longer constructs the whole XML tree first.
The necessary OSM node/way/relation graph remains in memory.

Local imports without uncached missing members start no HTTP process. Initial
configuration, parsing, relation repair, preparation, rendering and finalization
all participate in the session. Existing data-folder and asset validation apply
to both modes.

Progressive operator `FINISHED` means scheduled. Use blocking import for scripts
that require an immediate final result. Interactive progressive requests in
background Blender are rejected. `commandLineMode=True` selects the blocking
driver, while a selected progressive option still forces separate objects.
To request combined meshes, set `progressiveImport=False` explicitly.

Panel options and default-level distributions are captured at start. Changing
them affects the next import. Random choices use session state, isolated from
other addons' random calls between callbacks.

## Provider contract and buildings_gn integration

`app/osm_import.py:OsmImport` orchestrates the complete operation. Existing public
manager methods remain blocking wrappers over the same producers. Optional
`process_steps(session)`, `render_steps(session)` and
`render_extra_steps(session)` hooks may yield:

- `None`: soft processing checkpoint;
- `YieldToUI()`: immediate end of the current callback;
- shared `Request`, `Download` or `Drain`: existing network protocol.

`session.unit(callback)` wraps a synchronous coherent publication. It captures
only allocations from that synchronous call, owns them explicitly afterward,
finishes local roof work, waits for any export textures, validates inputs,
accepts the unit and advances visible revision/count before yielding to Blender.
A unit can contain several objects. Work progress need not correspond to an
object, and a producer can consume a complete prepared manager/network.

For work that itself spans turns, register each allocation immediately through
`session.ownership.own(kind, datablock)`. Never leave the legacy synchronous
`ownership.capture()` context open across a yield. Do not clean a yielding
transaction by comparing all current Blender data with a job-start snapshot.
The provider owns its numerical dependencies, private buffers and success-only
finalization. The session owns cancellation and disposal of pending IDs.

An older/custom backend without cooperative hooks runs once as an atomic fallback.
If it overrides a legacy method but only inherits a producer from an older base
class, its actual override takes precedence. That preserves the newer algorithm;
the fallback's long call does not acquire fine-grained responsiveness automatically.

The session establishes the renderer destination before `prepare`, initializes
projected bounds, preserves manager/action order, supports layer/renderer
finalization and manager cleanup, and scopes an available `import_performance`
policy to each active callback. These match the newer building/street lifecycle.
Both drivers share this session. `is_busy()` exposes the active-import guard for
later provider Update/editor operators.

The street handoff is at
`C:/dev/workspace/blosm_streets_dev/docs/OSM_PROGRESSIVE_STREETS_REQUIREMENTS.md`.
Keep `buildings_gn` setup/backend registration and its revised building algorithms
when resolving conflicts. In particular:

- PADL Python/GN producers must keep complete outline/part groups, editable owners,
  generated geometry, persistence Texts, statuses and shared asset dependencies
  coherent. GN status evaluation still needs its local synchronous update.
- Realistic/PML on that branch has two passes. Preserve preprocessing and
  `renderExtrudedVolumes`; do not inherit an old per-building producer that
  bypasses a newer render override. Its separate-object footprint-GN finalization
  and base-layer finalization hook require the baseline checks from the plan.
- Way preparation remains necessary for facade analysis even if road meshes are
  disabled. Street producers must retain the actual dependency graph, component
  publication, raw-versus-terrain-fitted state and stored-precision checks.
- Replace the street backend's synchronous global-difference transaction before
  adding inner yields. A connected component can be an atomic commit while its
  numerical construction yields many times.

The release implementation supplies the extension contract without copying the
branch algorithms. The native PADL integration below implements that contract;
the street integration remains separate work.

### buildings_gn merge qualification (13 November 2026)

The merge of `release` at `e5fe4ac` into `buildings_gn` at `bcc2fa6`
preserves native PADL Python/GN generation, complete outline/part groups, the
revised realistic building-part algorithm, renderer-owned preprocessing and
two-pass realistic rendering. At merge time PADL simple 3D, two-pass realistic
buildings and native streets used the atomic manager fallback. The PADL
integration below replaces that fallback for simple 3D. Two-pass realistic
buildings and native streets still run atomically; cancellation cannot interrupt
those native calls.

The merge adapts the branch's BMesh finalization for 2D GN buildings, separate
curves and forest areas, and passes the renderer to realistic footprint
finalization. Package-relative scheduler imports preserve Blender extension
isolation; standalone tests load a package shell instead of executing startup.

Validation uses Blender 5.2 in this chat's isolated profiles. Native PADL adapter
and building-part tests cover both backends, status handling, terrain, persistence
and reopening. Import performance-policy checks exercise the actual operator and
controller. Shared OSM tests cover real extracts, acquisition, grouping, errors,
terrain, forests, trees, 2D GN, lifecycle/undo and rectangular realistic
Python-volume import/export.
The realistic test explicitly selects `Result.VolumeMesh` and a compatible asset
bundle. The tracked `assets/prochitecture_buildings.blend` lacks the legacy
`Blosm Building Single Object` and `Blosm Building Separate Object` node groups;
the realistic `FootprintWithGn` output path remains unqualified with this binary.
This is distinct from the PADL simple-3D GN backend tested above.
The legacy realistic `VolumeMesh` route also fails on some Karl Marx Allee
roofs because its item renderer calls the obsolete `BuildingRendererNew.createFace`
API. Its bake/undo ownership check therefore uses the rectangular fixture;
the merge does not claim to repair the legacy complex-roof renderer.

Reports are under `C:/dev/workspace/artifacts/osm_progressive_import/merge_buildings_gn/`.
The broad host-Python discovery also includes older Blender-only and numerical
street suites that require `bpy`/`numpy`; it is not a passing full-suite claim.
Shared scheduler, transport, extension and performance checks are run separately
with their required environment, including localhost permission for HTTP fixtures.

### Native PADL progressive buildings (13 November 2026)

`building/padl_integration/manager.py` implements `render_steps(session)` for
both `geometry_nodes` and `python`. After resolving outline/part ownership, each
complete building group is passed to `session.unit()`. An outline and its owned
parts are never split between commits. GN owners, Python editable inputs and
generated volumes, persistence Texts, generation statuses and dependencies are
published together. The existing progressive controller provides cancellation,
viewport redraw, Edit Mode pausing and grouped undo. Progressive mode always
uses separate buildings; the stored blocking single-object preference is retained.

Blocking imports drain the same generation code in one ownership scope and
retain batch status-reader reuse. Progressive GN imports verify the shared graph
for the batch and reuse a status reader in a private temporary scene. The private
scene contains one status owner, its empty mesh and the shared terrain dependencies
when needed. Each synchronous building call borrows that building's input mesh,
reads diagnostics and releases the mesh before yielding. It does not add temporary
owners to the growing import scene. Updates to the original root or nested generator NodeTree IDs
trigger revalidation before the next building. Deferred scene/object/mesh updates
from publication and normal viewport work do not invalidate the generator.
NodeTree updates matter even when their geometry/shading/transform flags are all
false. Deleted GN dependencies and changed embedded program Texts stop the import
with a partial-result error. The private scene, owner, empty mesh and reader group
are tagged transient: they persist between import turns, are removed on batch
teardown and are stripped from restored undo snapshots. Shared generator and
terrain dependencies remain owned by completed buildings. Deleting a private
helper stops further generation with an error and cleans up the remaining helpers.
Terrain status uses the import scene's frame/subframe, FPS, time remapping and
viewport subdivision simplification settings; the private scene follows changes
without changing the active scene or its time.
Only the private graph's notifications are ignored; edits to shared generators
in another user scene still invalidate the asset. The first progressive building
updates the main view layer once to consume asset-initialization notifications;
subsequent status reads update only the private layer.
The report Text and app result/status attributes are updated after each building,
so cancellation and failures retain the report for completed work. Renamed or
deleted accepted owners are handled by RNA identity/liveness rather than names;
the cumulative report still records buildings the user later deletes.

Failures before the first accepted building reclaim pending nested GN groups and
their persistence Texts. Shared ownership cleanup revisits dependencies released
by removing their consumers, without purging unrelated orphan or user data.

Footprint ownership preparation, initial asset loading, a complete building's
Python generation or GN evaluation, and its status read remain synchronous.
The viewport becomes interactive between building groups; this is not a hard
latency bound for a single complex building. Terrain projection uses the existing
backend and shared session validity checks. Native streets and legacy two-pass
realistic generation are outside this integration.

Qualification outputs are under
`C:/dev/workspace/artifacts/osm_progressive_import/padl_progressive/`.
The Karl Marx Allee fixture produces 138 complete groups from 188 accepted
footprints. Both backends match their blocking outputs and reports. In the
isolated dataset run, first progressive results appeared at about 2.27 seconds
for GN and 0.43 seconds for Python; total generation was about 38.11 and 1.22
seconds respectively. These headless scheduler measurements include initial
asset preparation and are observations from one run, not viewport latency
guarantees or a backend speed comparison.

The GUI performance investigation found that the original broad graph-update
guard revalidated the entire shared asset after each building's deferred scene
updates. On the first 20 Karl Marx Allee groups it took 48.09 seconds, including
44.49 seconds in 20 asset validations. Restricting invalidation to the actual
generator IDs reduced that same GUI run to 5.31 seconds with one validation;
the final registered benchmark repeated it in 5.48 seconds (about 9x faster).
Before the private status scene optimization, the full 138-group GUI import
completed in 72.88 seconds, including 32.57 seconds of status evaluation and
34.53 seconds in undo checkpoints. The speedup above describes the measured
20-group comparison, not every dataset.
Per-building status evaluation and undo checkpoints remain enabled. Report links
are only assigned when needed, avoiding writes to unchanged owner metadata.
The GUI suite now asserts that uninterrupted import, ordinary scene edits and
accepted-object renaming/deletion do not revalidate an unchanged generator;
root/nested graph, interface, document and deletion checks protect invalidation.
Profiles and event probes are under
`C:/dev/workspace/artifacts/osm_progressive_import/padl_progressive_performance/`.

The subsequent private-scene optimization avoids rebuilding the growing main
scene's dependency relations merely to add/remove a diagnostic owner. The final
138-group GUI run in Blender 5.2 took 43.59 seconds (an earlier repeat took 45.70),
versus 72.88 seconds before this change. Status work fell from 32.57 to 3.35
seconds. Average observed publication intervals fell from 510 to 294 ms, and
the 95th percentile from 886 to 494 ms. A Python run of the same extract took
7.20 seconds with 46 ms average intervals. GN therefore improves by about 40%
in total time but does not yet match Python responsiveness. Initial GN results
still take about 2.8 seconds. Main-scene work during undo checkpoints dominates
the remaining cost (34.60 seconds); individual building evaluation and initial
asset loading also remain synchronous. These are observations on this machine,
not guaranteed input-event latency. Editable GN modifiers, per-building reports
and the existing undo boundaries remain enabled. All four native building suites
pass, including animated/remapped-time and simplified terrain, helper deletion,
failure cleanup, dataset geometry parity and actual GUI history operations.
Measurements and validation provenance are in
`C:/dev/workspace/artifacts/osm_progressive_import/gn_responsiveness/qualification.json`.

The separate `blosm_build` repository uses explicit runtime file inventories.
Its common inventory includes `app/osm_import.py` and `progressive/ownership.py`;
the Pro inventory also includes `building2/bake_job.py` and
`building2/bake_worker.py`. These four entries were synchronized after the OSM
implementation at the user's request. Build from an OSM-capable source revision
such as `1af9075`; earlier revisions do not contain these required files.
The Base qualification payload for the implementation task
adds the two common files explicitly and is an evaluation artifact, not an
official release package.

## Tests and measurements

The release catalog now contains:

- `padl_progressive`: both native backends, evaluated blocking/progressive
  geometry parity, complete outline/part commits, retained editable inputs and
  reports, cancellation and allocation rollback including cold GN dependencies;
- `padl_progressive_ui`: actual operator/timer execution for both backends,
  uninterrupted and user-edited Undo/Redo, cancellation, accepted-object renaming
  and deletion, destination deletion and report persistence through history;
- `padl_progressive_datasets`: full Karl Marx Allee evaluated-output and report
  parity for both backends, translated/sloped terrain, shared GN terrain helpers
  and partial cancellation; records observed callback and first-result timings;
- `padl_status_scene`: private evaluator reuse, diagnostic parity, borrowed input
  release, exception/deletion cleanup, terrain dependencies and asset invalidation;
- `padl_progressive_performance` (opt-in GUI benchmark): times actual timer
  execution, asset validation, status evaluation and undo checkpoints. The
  default stops after 20 Karl Marx Allee groups; `--extended` completes all 138.
  Set `BLOSM_PADL_BENCHMARK_MODE=blocking` in the launch environment for the
  synchronous separate-object comparison, or `BLOSM_PADL_BENCHMARK_BACKEND=python`
  for the Python backend. It records `padl_performance.json`, including the
  intervals between observed building publications (not a hard UI latency bound);
- `osm_progressive`: real extract parity, local parsing, missing-member POST and
  companion reuse, captured settings, per-object GN, resource/transaction failures,
  terrain validity, forests, realistic materials and real uncached export;
- `osm_progressive_ui`: explicit GUI process testing actual timers/operator entry,
  visibility, grouping/user edits, pause, cancellation, deletion, staged bake
  history and disabling;
- `test_osm_scheduler.py` in the standalone suite: local scheduling, immediate
  boundaries, unknown commands, POST retry policy and actual helper roundtrip.

Existing OSM/roof/terrain/GN/Buildify/export/validation and 3D Tiles suites remain
in the catalog. Older blocking suites explicitly disable progressive mode.
Pro-only cases in the common new suite require `--assets`; Base runs the common
cases without importing Pro modules.

Example from this chat's checkout:

```powershell
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python313/python.exe' -B tests/run_release_tests.py `
  --blender 'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe' `
  --suite osm_progressive --suite osm_progressive_ui `
  --assets C:/prochitecture/assets/assets_release `
  --artifacts C:/dev/workspace/artifacts/osm_progressive_import/my_run `
  --profile-root C:/dev/workspace/artifacts/blender_profiles/01a1a828-899a-7200-8fed-ec81d14ab15a/my_run
```

Use a new output directory per run. All five Blender resource overrides and
factory/autoexec flags also apply to bake workers. Generated outputs stay in
artifacts; source/tests/docs stay in this checkout. No normal Blender profile or
main checkout is changed by these tests.

`osm_metrics.json` records first-result time, maximum callback, peak parent RSS on
Windows and total duration
for successful session cases. The GUI report separately measures actual event
loop execution. These measurements include native operations but are not hard
latency guarantees or a constant-memory claim for the whole OSM graph.
