# Progressive OSM import: implementation and merge plan

Prepared on **12 November 2026**. Status: **source-audited implementation plan;
runtime implementation and runtime qualification have not been performed**.

This is the original design/audit snapshot. Subsequent implementation and
qualification are recorded in [OSM_PROGRESSIVE_IMPORT.md](OSM_PROGRESSIVE_IMPORT.md).

## 1. Objective and authoritative inputs

Extend Blosm's existing progressive-import infrastructure to OpenStreetMap while
preserving the blocking path and preparing compatible integration into
`buildings_gn`. Downloading, parsing, preparation, generation, publication and
finalization all need explicit treatment. A timer around the final render loop
does not make the entire import responsive.

The user's option rule is mandatory: **when Progressive import is selected for
OSM, hide Import as a single object and treat that option as false in code**.
Keep the saved preference intact so switching progressive mode off restores it.
This applies before setup and layer construction, not just to the final join.

| Source | Audited state / role |
| --- | --- |
| Local `release` | `1bf4bb2a68d4f5c9546b282cec59e6857eb7d163` |
| Local `buildings_gn` | `b9517de05889a3a7f0d7bc8107d923f7af871f41` |
| This chat's development checkout | `C:/dev/workspace/blosm_progressive`, branch `progressive_import` |
| Main checkout | `C:/dev/workspace/blosm`; read-only, currently on another branch, so working files were not used as the release baseline |
| Street requirements | `C:/dev/workspace/blosm_streets_dev/docs/OSM_PROGRESSIVE_STREETS_REQUIREMENTS.md` |
| Street source audit | Matching committed branch plus relevant current uncommitted street work identified by that handoff; inspected read-only |

The street handoff includes qualified optimizations not all represented by its
parent commit. Recheck their committed status at integration time. Do not copy
neighboring uncommitted files into this work. Neither Git reference was changed
for this audit, and no new worktree was created.

Shared documentation on the audited release is `progressive/README.md` and
`downloads/README.md`; principal implementations are `progressive/core.py`,
`progressive/blender.py`, `downloads/transport.py` and
`threed_tiles/progressive.py`. This development checkout is older than the
audited reference: before implementation, preserve its unfinished work and
merge committed `release` changes according to its local `AGENTS.md`.

This plan does not implement Mapillary, change street geometry algorithms,
replace the building backends, or introduce 3D Tiles caching. It includes the
integration contract needed to avoid a scheduler redesign after the merge.

## 2. What the current code actually does

### 2.1 Release entry, acquisition and parsing

`BLOSM_OT_ImportData.execute` in `__init__.py` runs setup, resets projection,
copies scene settings into the singleton `app`, and calls `importOsm`.
`BlenderApp.initOsm` currently validates the data/assets directories, downloads
the main OSM file synchronously, initializes caches, and can load the 2D GN
library. Consequently a progressive timer started after `initOsm` would still
block during acquisition and asset loading.

The remainder of `importOsm` is approximately:

```text
choose setup -> change mode/selection -> create Osm/managers/layers
-> establish or derive projection -> parse main XML
-> acquire/parse missing relation members -> repair relations
-> finalize bounds -> set up terrain/BVH -> initialize layers
-> process managers -> render -> write scene origin -> clean
```

`parse/osm/__init__.py:Osm.parse` builds an entire ElementTree before iterating
nodes, ways and relations. Parsing also executes registered conditions and
manager callbacks, assigns layers, calculates bounds and stores incomplete
relations. It cannot simply be moved wholesale into a worker: these callbacks
belong to the configured import and may touch Blender or mutable application
state. Building part association and street preparation need a complete,
repaired source graph before publication.

`app/__init__.py:downloadOsmWays` sends an Overpass **POST** for missing ways;
the shared `downloads.Request` and worker currently support GET requests only.
`loadMissingWays` also reuses a supplementary local file. Preserve that OSM
storage behavior deliberately; the no-cache rule for 3D Tiles is not an OSM
storage redesign. Current release requires an existing configured data folder
for every OSM import, including local files.

### 2.2 Release lifecycle and feature paths

`BaseApp.process` executes managers in registration order. Release
`BlenderApp.render` executes renderer `prepare`, `Renderer.begin`, managers'
`render`, `Renderer.end`, managers' `renderExtra`, and renderer `cleanup`.
Branch-specific differences are listed in section 10; do not replace the
target branch's ordering with this release sequence during merge.

| Feature | Actual release path | Initial coherent publication unit |
| --- | --- | --- |
| 2D polygons/areas | `manager.BaseManager.render`, `renderer.Renderer2d`, `renderer.layer.MeshLayer` | Complete way/relation with all rings and modifiers |
| Roads/railways | `manager.WayManager`, `renderer/curve_renderer.py`, `renderer/curve_layer.py` | Complete source feature, possibly multiple splines |
| Simple buildings | `setup/base.py`, `building/manager.py`, `building/renderer.py` | Logical building including outline, parts and roof objects |
| Legacy realistic buildings | `setup/premium_default.py`, `setup/premium.py`, `realistic/building/renderer.py` | Same complete-building boundary |
| Experimental realistic/PML | `setup/realistic_dev.py`, `building2/manager.py`, `building2/renderer.py`, `action/volume` | Building with all style-generated footprints, materials and required auxiliary geometry |
| Separate trees | `renderer/node_renderer.py:SingleTreeRenderer` | Tree instance using shared asset meshes |
| Realistic forests | `setup/premium.py:setup_forests`, `realistic/manager.py:AreaManager`, `realistic/renderer.py:ForestRenderer` | Separate area geometry plus a coherent shared scattering-state update |
| 2D building GN | `building/gn_2d.py`, `util/geometry_nodes.py` | One object's geometry, face attributes and configured modifier |
| Import for export | Export layer/item renderers and `building2/texture_exporter.py` | Building with its required completed textures/materials |

The base manager traverses eligible relations, then ways, then rendered nodes,
honoring element-specific renderers. Preserve that order and action ordering.
Current realistic renderers share item pools, reusable roof generators, asset
selection counters and class-level action lists; they are not reentrant.

### 2.3 Concrete prerequisites and hazards

These are source findings, **not newly reproduced Blender failures**:

- Both OSM UI branches expose `singleObject`; 2D GN controls additionally depend
  on its saved value. `Layer` copies the option, simple building code reads
  `app.singleObject`, and realistic setup installs `Offset` conditionally.
  Normalize effective settings before all those decisions.
- `GnBldg2dLayer.attributeValues` accumulates a layer-wide list but
  `finalizeBlenderObject` writes it into each current object's FACE domain.
  Multiple separate buildings require object-local, face-aligned values.
- `AreaManager.createLayer` forcibly enables `singleObject`.
  `renderExtra` assumes one `layer.obj`; removing the flag alone breaks the
  forest pipeline. The inspected file also imports the 3D Tiles `BaseManager`
  despite calling the OSM constructor signature. Qualify and, if confirmed,
  fix that blocking-baseline defect explicitly before forest parity tests.
- `Terrain.init(createBvhTree=True)` temporarily transforms the existing
  terrain mesh/origin until `cleanup`. A user must not interact with that
  temporary state between callbacks. `initProjectionProxy` also collects all
  building vertices and applies Shrinkwrap in one call.
- Dome/onion `RoofMesh` objects are queued for end-of-import joining through
  `Renderer.toJoin`, keyed by mutable target names. A progressive building
  must finish its local joins before becoming a committed result.
- `TextureExporter.renderTexture` synchronously calls `bpy.ops.render.render`.
  Yielding between buildings does not interrupt that render.
- `Renderer.cleanup` may clean shared terrain; lifetime must cover all remaining
  consumers rather than whichever renderer happens to clean up first.
- Shared app/class state, global RNG use, rotating asset selection and mutable
  parser elements require session ownership or strictly scoped use. A generator
  retaining a context manager across `yield` does not restore it while idle.

## 3. User-visible behavior and effective settings

### 3.1 Option rules

Expose the existing Progressive import checkbox for OSM in the import box above
the data-type row. Use data-type-neutral tooltip/status/cancellation text.
Keep its default enabled. Hide both `singleObject` controls when the OSM
progressive option is selected; do not merely disable them.

Capture two distinct concepts:

```text
progressive_requested = saved progressiveImport option
run_progressively = progressive_requested and not commandLineMode
effective_single_object = saved singleObject and not progressive_requested
```

The separate-object rule follows the selected option even when command-line
execution uses the blocking driver. A script requesting combined objects must
explicitly disable progressive mode. Keep the saved RNA `singleObject` value
unchanged, including on failure/cancellation; never toggle it temporarily to
communicate the policy. Custom setup/layer creation must receive and respect
the effective policy, and validation must catch a progressive aggregate layer
being re-enabled internally.

| Progressive option | Saved singleObject | Execution | Effective singleObject |
| --- | --- | --- | --- |
| Off | False | Blocking | False |
| Off | True | Blocking | True |
| On | Either | Interactive progressive | False |
| On + commandLineMode | Either | Blocking driver | False |

Match existing 3D Tiles entry policy: an interactive progressive request without
an interactive window is rejected with instructions to select blocking mode;
`commandLineMode` explicitly selects blocking execution. Progressive entry
requires Object Mode; do not force a user's mode change behind a timer.

Keep 2D GN controls available when progressive mode is selected even though
effective aggregation is false. Their availability must describe mode/backend
capability, not depend solely on `addon.singleObject`. Preserve blocking
combined behavior and its existing UI unless a separately tested change is needed.

### 3.2 Meaning of progress

Display phase and useful counts: downloading, parsing, resolving relations,
preparing buildings/ways, generating, publishing, scattering/baking, finalizing.
Totals may be unknown. Do not promise a percentage based only on object count.

Separate `work_done`, `visible_revision` and `committed_units`. A complete
building may own several objects; numerical street tasks may own none.
Already committed objects are editable and stay on explicit cancellation.
The returned operator `FINISHED` means scheduled in progressive mode; expose
the final outcome/result independently of the original operator instance.

Report success, partial completion, cancellation or failure once, with retained
feature counts and aggregated errors. Per-feature detail goes to diagnostics.
Normal cancellation is not a geometry failure to catch and continue past.

## 4. Recommended architecture

Names below are proposed APIs, not claims that these files/classes already exist.

| Location | Responsibility |
| --- | --- |
| `app/osm_import.py` | OSM session, phase orchestration, effective settings, acquired files/parser/managers, result aggregation; same producer in both modes |
| `progressive/core.py` | Generic bounded scheduling, local checkpoints, HTTP commands and explicit return-after-publication command |
| `progressive/blender.py` | Existing timer/controller, context/lifecycle, undo and notifications; only small generic extensions |
| `progressive/ownership.py` | Explicit owned-datablock transaction utility, if shared by several adapters; no OSM/street algorithms |
| `parse/osm/__init__.py` | Reusable element decoder and incremental parsing producer; blocking wrapper drains it |
| Existing manager/renderer modules | Optional cooperative hooks for their own work and publication units |
| `downloads` | Backward-compatible method/body support and request handling; no OSM interpretation |

Do not add a second timer, retry implementation, thread pool inside Blender or
provider-specific scheduler. Retain import compatibility under Blender extension
package namespaces through relative imports.

### 4.1 Session lifetime

The session owns immutable copied settings, target scene/view-layer/window,
projection/units, stable input identities/revisions, progress/report sink,
transport, resource ledger, RNG state, producers and terminal outcome.

Create it before expensive or mutating initialization. Keep constructor and
operator preflight small, then schedule initialization as work. Register
resources immediately as created; do not rely exclusively on a generator's
`finally`, because closing an unstarted generator never enters its body.
Use idempotent terminal cleanup that also covers constructor, initial undo
checkpoint, timer-registration and first-step failures.

Audit calls reading `bpy.context.scene.blosm` after settings capture, including
terrain setup. Pass captured values explicitly. Replace reliance on `app.op`
with an independent report sink. Retain the singleton compatibility surface
where necessary, but prevent competing Blosm operations from reinitializing it.
Guard older modal overlay activity and later street/PADL updates as well as
ordinary Import. Normal unrelated scene editing remains allowed.

### 4.2 Producer and scheduler contract

Provide optional `process_steps(session)`, `render_steps(session)` and
finalization/disposal hooks, or equivalent registered adapters. A default
adapter calls an existing synchronous method once; do not introduce required
abstract methods into every old/custom renderer.

The producer chooses whether it consumes one element, a building group or a
complete prepared manager/network. Shared code never assumes a signature
`render(element)` or enumerates street product kinds.

Preserve current commands `None`, URL/`Request`, `Download` and `Drain`. Add an
explicit local immediate-return boundary, for example `YieldToUI`, so a commit
can display before further work. Dispatch it explicitly; today's `Job` treats
any unknown non-None value as a download request. Reject unsupported commands
clearly. Progress remains adapter state rather than arbitrary yielded objects.

Support CPU/local-file work without launching a network helper: use a lazy
transport facade or optional transport with an explicit local path. Existing
3D Tiles producers and polling/backpressure semantics must stay compatible.

Keep the current approximately 8 ms soft callback budget and maximum generator
advances initially; tune from measurements. `None` does not guarantee immediate
return. No scheduler can preempt an ongoing Blender operator, native roof solver,
GN evaluation or library append. Give long operations named progress stages
and measure them; add internal checkpoints only where invariants permit.

### 4.3 Publication and resource ownership

A transaction may stage several objects and publish them as one coherent unit.
Track every owned mesh/object/curve/text/node group/image/temporary scene and
file, plus edits to existing state. Distinguish pending allocations, shared
dependencies, committed allocations and replaced user data.

Never delete everything created since a global `bpy.data` snapshot across
callbacks. User-created objects between ticks must survive failures. Resolve
objects/collections by validated identity rather than delayed name lookups.

Recommended sequence:

```text
capture/validate inputs -> prepare (many checkpoints)
-> generate and validate -> stage detached resources (bounded batches)
-> revalidate inputs -> coherent synchronous commit -> record revision/count
-> return to Blender -> release temporary state
```

Cancellation rolls back pending owned resources only. Commit includes metadata
and all related geometry; after commit, cleanup errors are retained-result
errors, not justification for undoing earlier accepted units. Do not allow a
time budget to weaken geometry validation or dependency closure.

Stage first without visible provisional previews. Even detached Blender IDs can
enter a user-created undo snapshot: validate cleanup/Undo/Redo of staged data,
or keep preparation in Python and constrain RNA staging to a commit callback
until ownership reconciliation is proven. A preview feature is separate work
with its own save/history/cleanup policy, not an incidental side effect.

Persist scene geographic origin with the first retained publication (including
owned terrain where applicable). Partial/cancelled imports must remain correctly
georeferenced. If nothing was committed, remove only this session's temporary
metadata and leave the original scene settings intact.

### 4.4 Context, undo, invalidation and cleanup

Apply context preservation and any performance request around each synchronous
callback, restoring before idle. Do not restore a selection snapshot captured
at the very beginning of the import; preserve user edits between callbacks.

Reuse grouped checkpoints: one Undo for uninterrupted importing, separate steps
for intervening user edits. Increment the visible revision for every relevant
scene mutation, not only object creation. Include changes during finalization
or explicit cancellation in the history policy; current controller finalization
does not automatically checkpoint them. Silent undo/file-load teardown must
not push a conflicting undo step or open a result popup.

Retain pause on scene/view-layer mismatch or Edit Mode, and cancellation on
Undo/Redo, file load, addon disable or owning-window close. Validate identities
before dereferencing Blender data after each pause.

Renaming/deleting committed outputs does not recreate them. Deleting/unlinking
the destination collection stops the import with an error. Renaming it or moving
it within the same scene is valid. Detect terrain/source/node-group deletion or
relevant changes before further publication. Prefer conservative stop-and-report
for invalidated dependencies in the first version, retaining prior results;
do not silently combine geometry prepared from incompatible terrain revisions.
Use scoped dirty/revision tracking, not a full-scene hash every timer tick.

## 5. Acquisition, parsing and preparation

1. Split configuration/storage validation from network acquisition in `initOsm`.
   Preserve release's required data directory, file naming, supplementary file
   reuse, setup selection, licensing/assets and extension loading behavior.
2. Extend the shared request protocol with explicit method/body support for the
   read-only Overpass query. Preserve GET defaults. Existing IPC JSON has a
   64 KiB request limit: specify an encoding/size limit and split large missing-ID
   queries deterministically. Never implicitly enable retries for arbitrary
   non-idempotent POST requests; opt this read-only query into the retry policy.
3. Download through the existing helper. For the first implementation, handle
   bounded responses using its byte budget and write the OSM file in bounded
   chunks to an owned partial file before atomic publication. Avoid one large
   main-thread `write`. Oversized downloads must fail clearly, not silently
   truncate or silently bypass the limit. Streaming-to-file is a follow-up
   shared transport capability if supported dataset sizes require it.
4. Replace the full XML-tree load with a chunk-fed `XMLPullParser`/equivalent
   incremental reader. Decode a completed top-level element before clearing its
   children; remove processed XML nodes from the root so empty nodes do not
   accumulate. Keep the actual OSM node/way/relation graph needed by managers.
5. Use bounded file-read/feed/record-processing chunks and checkpoints inside
   expensive member/bounds work. Share the same element decoder and semantics
   with blocking parsing. Preserve deletion records, tags, condition order,
   custom manager callbacks, missing members, relation identity, bounds and
   projection decisions. Handle parse errors and closure/cancellation explicitly.
6. Repair missing members before dependent manager processing. Preserve the
   existing ability to report supplementary-download errors and continue with
   valid data, but make partial coverage explicit. Worker/protocol failure is
   terminal. Do not treat a cancelled request as an ordinary missing member.
7. Keep source graph and required manager data alive through all consumers.
   Do not publish buildings or streets directly as XML records arrive.

Recommended first parser implementation is cooperative main-thread decoding,
which avoids serializing setup callbacks and the full domain graph. It is not
a claim of background parsing or guaranteed latency. Measure large files and
single huge relations. If decoding/file I/O cannot meet the measured budget,
add a distinct process protocol returning bounded plain-record batches, with
source/task/revision IDs, pause backpressure, stale-result rejection, limits
and cancellation/reaping. Domain callbacks and Blender work remain on the main
thread. Do not send the singleton app or Blender references to a worker.

## 6. Release feature implementation

### 6.1 Base polygons, curves and ordinary buildings

Refactor manager loops into producers while retaining the synchronous wrappers.
Preserve element-specific renderers and `preRender -> render* -> postRender`.
The committed unit includes geometry, OSM tags, materials, parenting and required
modifiers. Large multipolygons, curve subdivision and `MeshLayer.slice` require
inner checkpoint opportunities; a half-created object is not a complete feature.

Checkpoint complete-graph building-part association, including explicit relations,
shared nodes and fallback BVH work. Generate one whole logical building at a
time with all parts and a common terrain/offset policy. Finish dome/onion joins
inside that unit, with selection scoped locally and stable target identity.
Do not defer work on a committed object until after the user can rename/delete it.

Keep per-renderer mutable roof/item state exclusive to its current unit. Preserve
construction order and style/asset decisions. Use job-owned RNG streams or save,
bind and restore RNG state per callback without holding the binding across a
yield; prove fixed-seed parity and independence from unrelated random use.

### 6.2 Realistic/PML buildings

Retain style resolution, `Building.actions`, `Footprint.actions`, preprocessing,
terrain, offset, volume and roof-generation order. `Offset` creates the separate
object, so it must be installed after effective settings are normalized.

Add checkpoints to manager/action preparation, then initially publish complete
buildings. Audit `Volume`'s temporary parser-element `t`/`ls` changes before adding
inner yields; use local representations or guaranteed restoration. Bound item
pool/generated-buffer retention. Release per-building work only after its
consumers finish; preserve shared materials and assets.

Record native polygonization, triangulation, asset append and BMesh conversion
as possible indivisible calls. They remain explicit limitations until measured
and refactored. Run action/exporter/cache cleanup once on every terminal path.

### 6.3 Terrain used by all features

Replace the import-long terrain-origin mutation with a nonmutating coordinate
adapter/BVH snapshot reproducing existing placement. Verify translated, rotated,
scaled and shared-mesh terrain; do not introduce an unrelated placement change
under the scheduler refactor. Required Blender evaluation is main-thread work.

Prepare projection-proxy inputs incrementally, own the proxy, and preserve its
shrinkwrap semantics. If chunking proxy evaluation is needed, prove placement
equivalence rather than assuming independent chunks are identical. Validate
terrain revision before publishing later units. Keep envelopes, proxies, BVHs
and assets alive until their final consumer; centralize disposal responsibility.
Do not adopt or delete an unrelated envelope just because its name matches.

### 6.4 Geometry Nodes in 2D mode

Support progressive separate objects as a first-class path, retaining the loaded
node group and configuring a modifier on every object. Use per-object face rows
for `building:levels`, including holes/multipolygons and any topology change
before modifier attachment. No layer-wide attribute list may spill into later
objects. Verify missing/height-derived levels against current behavior.

Share verified node groups, not mutable per-object modifier input state. Retain
Blender-version handling in `util.geometry_nodes`; do not impose PADL's Blender
5.2 requirement on this existing release feature. Force evaluation only for
required geometry reads, not an unconditional update after every object. Measure
viewport/GN evaluation outside the callback as well as Python callback duration.

### 6.5 Forests and individual trees

First qualify the suspected `AreaManager` base-class import defect and repair
it in an isolated prerequisite if confirmed. Ordinary base-mode vegetation
polygons are distinct from realistic Dynamic Paint/particle forests.

Replace realistic forest `layer.obj` assumptions with a layer-owned collection
of explicitly owned **separate area objects**. Prepare/intersect/slice/project
each whole area, including holes, progressively. Keep layer-level surface-area
accounting and shared scattering state independent of geometry aggregation.
Do not turn `singleObject` back on internally.

Retain the existing density equation `ceil(treeDensity / 10000 * total_area)`;
rounding each patch separately changes the result. Stage a new brush and publish
its membership, cumulative area, density weights and particle-count change
coherently. Preserve required modifiers/surfaces/vertex groups and record every
new slot on an existing terrain. Copy settings when mutating a shared particle
settings datablock would affect earlier imports. Rollback must restore exactly
the state changed by the pending unit, preserving prior committed patches.

Prototype partial Dynamic Paint/scattering evaluation early. If its constraints
require a whole-layer update, expose preparation separately, then perform one
coherent scattering commit; do not describe hidden brushes as visible trees or
claim fully incremental scattering. Acceptance for responsive forests requires
measured usable partial scattering, or an explicitly agreed limitation before
shipping this mode. Replacing the vegetation system with a new GN scatterer is
not an implicit part of this task.

Own automatically created flat terrain and derive its extent only after local
file bounds are complete. Keep it if committed forest results depend on it;
otherwise clean it on cancellation. Individual trees can publish independently
using shared meshes, with stable asset cycling, random rotation and scale.

### 6.6 Import for export

Keep export renderer selection, facade UV layout, material semantics and texture
cache keys. Treat texture creation as an explicit dependency before publishing
a building that needs it. Cached textures can use the ordinary building path.

For an uncached texture, first separate preparation, rendering, validation and
file publication from the nested material call. Refactor consumers to request
or suspend for the texture rather than calling a blocking render deep inside
the feature generator. Limit outstanding bakes and reuse completed results.

Prototype native asynchronous render-job integration in an isolated Blender:
prepare the template under a short context scope, start/poll a job, and resume
only after its completion signal. Verify owning-window behavior, cancellation,
render errors, input edits, teardown and missing-window cases. Do not hold
selection/performance overrides across that asynchronous wait.

If native render integration cannot preserve scene interaction/lifecycle,
use a separate isolated Blender bake process with a bounded job protocol and
explicitly serialized/copied inputs. This is not the ordinary HTTP helper.
Keep the parent render-template state immutable, handle cancellation and stale
results, and qualify cross-platform packaging before adopting it. Choose one
qualified implementation; do not maintain two untested paths.

Until that gate passes, the existing synchronous bake is only a compatibility
fallback with a clearly reported pause, not completed nonblocking export support.
Blocking export continues to drain the same dependency workflow synchronously.

Use owned temporary files/scenes and atomic cache-file publication. Preserve
files referenced by committed materials; remove only incomplete owned outputs.
Do not purge previously existing templates/caches or delete valid files merely
because scene objects were undone. Document that Blender Undo does not revert
the external texture cache. Test real uncached rendering, not just mocks/cache hits.

## 7. Blocking mode, compatibility and failures

Blocking and progressive drivers must run the same producers, algorithms,
validation and finalization. Retain existing public synchronous methods as
wrappers. The non-Blender command-line path must not import Blender timers or
unconditionally construct Blender session objects.

Classify failures:

- Invalid global configuration, malformed primary data, invalidated required
  inputs or helper/protocol failure: stop, dispose pending work, retain and
  report previously committed results where valid.
- Recoverable feature geometry/asset failure: roll back that whole pending unit,
  record diagnostic identity and continue only where remaining dependencies are
  valid. Do not skip a failed prerequisite and publish its dependents.
  Preserve any explicitly supported provider-level partial-success policy, such
  as PADL's per-footprint diagnostics within a coherent building-group result;
  do not mistake a reported recoverable subfeature failure for a transaction
  failure or silently count it as fully successful.
- Postcommit disposal/report failure: retain accepted geometry and report the
  cleanup problem without pretending the whole unit never committed.
- Explicit cancellation: cease scheduling, discard deferred/stale results,
  dispose pending work, retain completed units and emit one terminal result.

Success-only `renderExtra`/joins/finalization must not run blindly on cancellation.
Some providers need partial-result metadata finalization; make it a separate
explicit hook. Always run initialized-resource disposal exactly once.

## 8. Integration contract for buildings_gn buildings

### 8.1 PADL simple 3D: both Python and GN backends

The branch's simple 3D path is `setup/base.py` ->
`building/padl_integration/manager.py:PadlBuildingManager` ->
`BuildingGenerationService`; this is not the release simple renderer or 2D GN.

`OsmFootprintSource.collect` and `BuildingFootprintAssembler.assemble` collect
and validate complete ownership, explicit relations, outer/inner rings and
parts before `group_buildings`. Keep this complete-source preparation barrier,
with inner cooperative checkpoints. Source IDs/order must not depend on ticks.

With separate-object mode, one group means a **whole building**, not one object
per footprint. A group's outline and parts share terrain datum and editable
source identity. Preserve deterministic seeds, PADL attributes/presence flags,
program identity and source-key namespaces.

Publish a group transaction containing its editable owner/mesh, Python generated
output/mesh when applicable, persistence Texts (`padl_source_keys`,
`padl_document`), report references and newly required helpers. Current service
failure cleanup covers only part of these resources; extend explicit ownership
before exposing internal yields. Retain dependencies of accepted groups.

Atomic publication concerns the consistency of owner, generated output and
statuses, not mandatory all-or-nothing geometry success. Preserve PADL's existing
per-footprint failure policy: valid siblings may be published with explicit
partial-result diagnostics when one footprint fails. An unexpected transaction
failure still rolls back the pending group.

| Backend | Audit finding | Adapter requirement |
| --- | --- | --- |
| Python | `PythonBuildingBackend` reads editable footprints, invokes PADL roof adapters, combines MeshData, creates a generated object and hides the source owner | Keep generated world-space/identity-transform semantics, status and roof failure policy; commit owner/output/metadata together |
| Python terrain | `BlenderTerrainSampler` builds evaluated terrain BVH per group | Consider session snapshot reuse only with unchanged coordinate/coverage semantics and revision validation; preserve hole/part boundary checks |
| GN | `GeometryNodesBackend` uses verified compiled trees and a batch `StatusReader` wrapper | Keep the wrapper as an explicitly owned session resource; dispose on cancellation/start failure without duplicating libraries per tick |
| GN status | `gn_runtime.py` creates a temporary object, updates/evaluates it, reads status attributes and removes it | Treat this local evaluation as an atomic main-thread stage; do not remove the required update or yield with the inspection object/context half configured |

Native roof solvers and mathutils-dependent code cannot simply run in the
ordinary Python download helper. Measure long calls and refactor only at valid
backend boundaries. Preserve `BuildingProgramLibrary`,
`VerifiedNodeGroupLibrary` asset identity/validation and extension-safe
`padl.services.configure_blosm` resolution. Blender 5.2 gating remains specific
to the GN backend rather than the shared scheduler.

### 8.2 Realistic/PML and footprint GN on buildings_gn

This is a third building pipeline. It prepares ways for facade classification,
including when street output is disabled; building managers perform global
outline/part/shared-edge analysis and action preprocessing.

With `renderAfterExtrude=False`, generation is two-pass: `render(building)`
prepares style/volume state and `building.renderInfo`; later
`renderExtrudedVolumes(building)` publishes results. Preserve this distinction
and bounded checkpoints. Do not flip the flag or fuse passes merely to display
the first building earlier. Streaming requires proof that all preprocessing,
facade, material and geometry dependencies remain equivalent. Record this path's
existing retained-buffer peak rather than claiming a new constant-memory bound.

The `preferableResult=FootprintWithGn` route uses `item_renderer/footprint.py`,
`item_renderer/gn`, `building2/content_generator.py`, and Python fallback volumes.
It is separate from PADL. Publish all footprint/holes/parts/fallback geometry and
attributes of a building coherently.

Two merge prerequisites need targeted baseline fixtures:

1. Separate-object `building2/renderer.py:postRender` calls `layer.finalize()`
   without its renderer, while `RealisticBuildingLayer.finalize(None)` returns;
   the renderer-supplied call is currently in the combined-object branch.
   Establish and correct the per-object GN finalization contract as required.
2. The branch lacks release's `finalizeBlenderObject` call path in its base layer
   code. Reconcile this deliberately so the GN2D per-object attribute/modifier
   fix survives the merge; do not accept either side wholesale.

These are source-level risks, not measured failures. Qualify both Python/GN
PADL and realistic footprint-GN separately in the compatibility milestone.

## 9. Street/intersection integration contract

Read the complete external street handoff before implementing or merging. Its
P01-P16 requirements map to this plan as follows:

| Requirements | Implementation obligation |
| --- | --- |
| P01 phases, P02 bulk producer, P03 dependencies, P04 checkpoints | Sections 4-5: optional whole-manager producers, preserved preparation order, CPU-only inner steps |
| P05 work/visible/committed distinction | Sections 3.2 and 4.3 |
| P06 main thread, P07 shared blocking algorithms | Sections 4-5 and 7 |
| P08 ownership, P09 multiobject commit | Section 4.3; synthetic bulk transaction test |
| P10 input validity/exclusion, P11 lifecycle | Sections 4.1 and 4.4 |
| P12 determinism, P13 memory | Sections 5-6 and per-backend retention rules |
| P14 backend separation | Section 4 and provider-specific integration |
| P15 performance scope | Section 10: per-callback/finalization scope |
| P16 measured responsiveness | Section 12 |

The latest user option rule supersedes any handoff suggestion to progressively
flush aggregate single-object layers. Atomic multiobject publication is still
required: a component can commit many separate objects together.

`WayManager.process` and `StreetGenerator.do` prepare the complete network;
`WayManager.render` passes the manager to the `StreetRenderer` facade. Preserve
that facade/registry, feature providers, physical source inventory and boundary
support. No street solver/task enum belongs in the shared OSM scheduler.

The initial merge may use a synchronous street fallback for compatibility. Full
progressive street support additionally requires localized producer refactoring:

1. Expose resumable component adaptation and dependency execution in
   `PythonStreetRenderer._render`, `NetworkScene.update` and `UpdateCoordinator`.
   Preserve actual task dependencies/topological order. A partially executed
   task must not advance its completion/digest state.
2. Replace the synchronous global-difference `ImportTransaction` before any
   yield enters its scope. Register controller, sources, outputs, texts and
   helpers through explicit ownership.
3. Checkpoint numerical tasks, UV/overlap validation, terrain fitting,
   stored-float32 publication checks and staging. A connected component is an
   atomic publication boundary, not one giant callback.
4. Commit the complete validated component. Preserve `placement.raw` as the
   numerical accepted state, distinct from terrain-fitted published geometry.
   Keep old accepted outputs/manifest until the coherent replacement is ready.
5. Release published mesh buffers at the existing safe boundary. Treat failure
   during later retirement/material refresh separately from precommit rollback.
6. Add active-import exclusion and revision checks to Street Update/editor
   operators. Preserve reconstruction from saved sources/manifest after reopen
   and Undo/Redo; process-local caches cannot become correctness requirements.

A dominant connected component can delay the first visible mesh while numerical
work remains responsive. Smaller permanent commits need closure over seams,
Bundles, overlap checks, terrain references and furniture; they are not required
for the first safe integration. Visible provisional previews are deferred until
their cancellation/undo/save behavior is independently designed and tested.

## 10. Merge strategy and responsibilities

Do release runtime work only in this chat's synchronized development checkout.
Transfer to `release` only on request. Later merge into `buildings_gn` only in
the appropriate authorized development checkout; do not work in another chat's
folder or change this chat's reference branch implicitly.

Prefer small commits by responsibility: prerequisites, scheduler/session,
acquisition/parser, base adapters, realistic/terrain, GN2D, forests, export,
tests/docs. Keep algorithm changes separate from scheduling changes where
possible so parity failures can be attributed.

| Merge-sensitive area | Preserve/reconcile explicitly |
| --- | --- |
| `__init__.py`, `app/blender.py` | OSM entry/session extraction, release data-dir validation, GN/backend checks and correct operator return semantics |
| Render lifecycle | Branch runs `Renderer.begin -> renderer.prepare -> manager.render -> Renderer.end -> manager.renderExtra -> optional layer.finalize -> renderer.finalize -> renderer.cleanup -> manager.cleanup`; release prepares before begin and lacks some hooks |
| `BaseApp.process` | Branch projected extent initialization and manager/action dependency order |
| Setup/building managers | Branch PADL selection and realistic way-before-building setup; do not replace with older release setup |
| Layer/finalization hooks | Separate-object GN2D and realistic footprint-GN attribute/modifier finalization |
| `progressive/blender.py` | Branch Import performance scope around each actual callback/finalization, restored while idle/paused |
| Backend/assets | Street facade/registry, verified building libraries, shared asset loading and source identity |
| Street optimization work | Reconcile the commits containing audited changes, not merely the older parent commit or a copy of local dirty files |

Introduce optional lifecycle hooks/default no-ops without accidentally invoking
existing side effects twice. Preserve success order and dispose-once semantics.
The source-audited render-order difference needs explicit compatibility tests,
not a generic reorder made on assumption.

Merge readiness means the common contract works and existing branch backends
still function. Fully responsive PADL/PML/street generation requires their
adapter milestones too; a clean Git merge alone cannot establish that claim.

## 11. Implementation sequence and exit criteria

| Stage | Work | Exit criterion |
| --- | --- | --- |
| 0: synchronize and establish baselines | Preserve unfinished work, record/merge current committed release, identify qualified fixtures/assets, reproduce suspected forest/GN issues | Blocking behavior recorded; prerequisite fixes distinguished from scheduler changes |
| 1: effective settings and UI | OSM progressive checkbox, hide/override singleObject, maintain GN2D availability, capture entry policy | UI and actual operator truth-table tests pass; saved preference unchanged |
| 2: generic session/lifecycle | Optional producer hooks, local/lazy transport, explicit immediate return, ownership, independent result sink, context/undo/cleanup | Individual + synthetic bulk producers pass cancellation/startup/ownership/history tests |
| 3: OSM acquisition/parser | GET/Overpass POST, bounded file publication, incremental XML/relations, progress | Large/local/server/missing-member fixtures responsive and equivalent; no helper for local-only work |
| 4: base features/terrain | Nonmutating terrain setup, complete polygons/curves/buildings, feature-local roof joins | Blocking separate vs progressive parity, visible incremental results and native undo verified |
| 5: realistic/PML and GN2D | Checkpoint preprocessing/actions, per-object GN face data and modifiers | All release building modes pass geometry/material/attribute/cleanup tests |
| 6: forests | Baseline repair if needed, separate brushes, transactional shared scattering | Density/terrain/cancel/undo parity and measured partial-forest behavior qualified |
| 7: export | Explicit texture dependency flow and qualified asynchronous bake strategy | Uncached real export render keeps import lifecycle correct; cache/material/roundtrip tests pass |
| 8: release qualification | Cross-mode real datasets, failure injection, memory/latency, 3D Tiles regressions, docs | Acceptance matrix satisfied; long indivisible operations/limits reported accurately |
| 9: buildings_gn integration | Merge reconciliation; PADL Python/GN and two-pass PML adapters; component street producer | Both building backends and mixed streets/buildings pass branch tests and progressive GUI qualification |

Run short forest-scattering and export-render feasibility spikes during stages
0-2, not after all easy features are finished. Their results determine the
remaining implementation work. Do not call the whole OSM feature complete while
required modes silently fall back to long blocking work.

## 12. Validation plan

### 12.1 Shared and entry-policy tests

- Exercise all option rows in section 3 through the actual Import operator;
  check both OSM UI branches and programmatic calls, saved-true settings,
  commandLine/background policy and 2D GN controls.
- Test local commands without HTTP, typed dispatch/unknown commands, soft budget
  and immediate return. Existing 3D Tiles commands retain behavior.
- Synthetic bulk producer: complete prepared input, several dependent hidden
  tasks, staged multiobject commit, failure after an earlier accepted group.
  Test a building-like consumer waiting for way preparation with way output off.
- Cancel before first generator advance, at every phase, during staging,
  immediately before/after commit, and with deferred replies. Inject timer/start,
  parser/helper, allocation, metadata, commit and postcommit-cleanup failures.
- Create unrelated sentinel objects/materials/texts between callbacks; prove
  rollback does not remove them. Exercise user changes to shared terrain/assets.
- Actual GUI: navigation/selection changes, Edit Mode and scene/view-layer pause,
  destination rename/delete/unlink, input invalidation, Undo/Redo, file loading,
  addon disable, window close and conflicting Blosm operation exclusion.
- Verify exactly-once terminal notification/cleanup, retry after failure and
  repeated import, no leaked callbacks/helpers/resources. Repeated Undo/Redo
  after cancellation must not resurrect incomplete staged canonical results.
- Compare deterministic fixed-seed geometry/attributes/materials/tags with the
  same backend in blocking **separate-object** mode. Test blocking combined mode
  separately for regression. Do not widen numerical tolerances merely because
  scheduling changed; justify exceptions with measured underlying behavior.

### 12.2 Feature matrix

| Feature | Required fixtures and observations |
| --- | --- |
| Parsing | Ordered OSM, deletion records, missing/invalid members, nested building relations, malformed XML, huge relation, local bounds absent/present, supplementary reuse/POST, cancellation before any output |
| Simple buildings | All supported roof families, holes and multiparts, complete dome/onion joins, terrain and flat scene, renamed/deleted committed buildings |
| Realistic/PML | Legacy and experimental setups, generated footprints, mixed material/style choices, hipped solver failure, action/item/cache cleanup, shared-asset lifetime |
| Curves/areas | Roads/rails/multiline relations, long subdivisions, large multipolygons with holes, terrain clipping/slicing, bevel and material parity |
| GN2D | At least two buildings with different levels, missing/height-derived values, holes, shared node group/independent modifiers, correct FACE counts, evaluated geometry and terrain |
| Forests/trees | Several patches and holes, aggregate density rounding, existing/generated terrain, previous forest import, shared particle settings, cancel after one patch, partial scattering and undo |
| Export | Cached and uncached textures, actual bake failure/cancel, template/file lifetime, export UVs/materials, texture cache reuse and GLB roundtrip |
| PADL Python/GN | Same source/seed/program, outline/parts datum and holes, statuses, native roofs, owner/output/Texts, helper/group reuse, missing/corrupt assets, terrain revisions, reopen/rebuild |
| Branch realistic GN | Two-pass prerequisites, footprint GN plus Python fallback, separate-object finalization, mixed roads/buildings with road output disabled |
| Streets | Individual and dominant components, intersections/Bundles/seams/boundary support, terrain PROFILE/DIRECT, float32/UV/material parity, partial publication, update/save/reopen/history |

Reuse and extend current release tests: `check_release_osm.py`,
`check_release_roofs.py`, `check_release_terrain.py`, `check_release_pml.py`,
`check_release_geometry_nodes.py`, `check_release_export.py`, plus
`test_progressive_transport.py`, `test_download_resilience.py`,
`check_progressive_import.py`, `check_progressive_ui.py` and
`check_concurrent_tiles.py`. Existing generic GN tests do not establish GN2D
importer face-attribute correctness; add a dedicated importer fixture.

On `buildings_gn`, extend `test_padl_source.py`, `test_padl_osm_mapping.py`,
`test_padl_program_cache.py`, `blender_padl_adapters.py`, `blender_padl_parts.py`,
`blender_padl_berlin.py`, `blender_padl_reopen.py`, asset-library and
`check_import_performance.py` tests. Follow the street handoff's publication,
UV, terrain, lifecycle, source-update and numerical suites. Use the established
Berlin, Canal Street and Nashville street inputs, not hard-coded runtime cases.

### 12.3 Performance evidence and realistic limits

Record time to scheduled job, first visible/first committed result, callback
p50/p95/max, longest atomic operation, total duration, cancellation latency,
peak memory including undo, queued bytes/records and residual resources after
cleanup. Profile parsing, association, terrain BVH/proxy, native roofs, GN
evaluation, forest evaluation and texture rendering separately.

Retain the scheduler's current soft budget initially. Set numeric acceptance
thresholds after stage-0 baseline measurements on named hardware/datasets;
do not claim an unmeasured hard frame-time guarantee. GUI heartbeat/interaction
measurements must include depsgraph/GN/scattering work after callbacks and undo
checkpoint cost, not only time inside Python generators. Any remaining long
atomic call must be identified by phase and reported with its measured bound.

The full OSM graph, preparation structures and any required branch two-pass
buffers consume dataset-dependent memory. Bound extra queues/results and avoid
whole-scene copies per tick; do not claim constant overall memory. Qualify one
large building/relation/component as well as many small independent features.

Use this chat's isolated factory Blender profile under
`C:/dev/workspace/artifacts/blender_profiles/01a1a828-899a-7200-8fed-ec81d14ab15a/`,
all five required `BLENDER_USER_*` process-local overrides, and
`--factory-startup --disable-autoexec`. Do not attach to normal Blender sessions
or modify the installed addon/Python/profile. Apply the same isolation to child
bake processes/resource-guard launchers, using task-owned subdirectories where
concurrent processes need separate profiles.

Store development outputs under `C:/dev/workspace/artifacts/osm_progressive_import`
and street qualification under `C:/dev/workspace/artifacts/streets`. Keep source,
test code and documentation in the repository. Do not stage generated scenes,
test binaries, logs or JSON reports.

## 13. Delivery checklist

- [ ] UI and runtime enforce the separate-object rule without changing the saved preference.
- [ ] Acquisition and parsing participate in the session rather than blocking before it starts.
- [ ] Every required release feature has a qualified producer and explicit finalization/cleanup.
- [ ] Terrain and shared assets remain safe while users interact.
- [ ] Complete buildings include parts/roof joins; GN attributes/modifiers are valid per object.
- [ ] Forest and uncached-export behavior is measured and accurately described.
- [ ] Ownership, stale-input validation, undo, cancellation and terminal reporting pass actual GUI tests.
- [ ] Blocking/background and 3D Tiles regressions pass.
- [ ] Synthetic bulk producer proves the street contract before branch integration.
- [ ] Merge preserves both PADL backends, realistic/PML phases, lifecycle additions and performance scope.
- [ ] Branch-specific responsive generation is qualified before claiming seamless progressive operation there.
- [ ] Documentation records source commits, commands, dataset/asset versions, results and remaining limitations.

The planning task produces this document only. None of these unchecked items is
claimed implemented or tested by the source audit.
