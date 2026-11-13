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

Before merging, recheck the current street handoff at
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

This release implementation does not copy or change those branch algorithms.
The synthetic bulk/override/dependency tests qualify the extension contract;
actual progressive PADL and street qualification belongs to their integration
after evaluation. A clean source merge is not that qualification.

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
