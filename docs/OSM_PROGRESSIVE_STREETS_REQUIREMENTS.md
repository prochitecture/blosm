# Progressive OSM import: street and intersection integration requirements

**Implementation update, 13 November 2026:** the dev82 street integration is
documented in [Cooperative Python street import](STREETS_PROGRESSIVE_IMPORT.md).
This document retains the original 12 November source audit and proposed merge
contract; its descriptions of the then-current implementation are historical.

Technical handoff for the chat implementing progressive OSM import on `release`.
Prepared on **12 November 2026** from the street-development checkout based on
`buildings_gn`. This is a design and source audit, **not an implementation of
progressive streets**.

The shared importer should support both individual OSM features and a renderer
that consumes a complete prepared network, performs many dependent numerical
tasks, and publishes several related objects together. Do not make one OSM way,
one Blender object, or one call to `render()` the mandatory scheduling unit.

For the release implementation, start with the
[essential requirements](#2-essential-requirements-for-the-release-implementation),
[shared scheduling contract](#5-proposed-shared-scheduling-contract),
[implementation sequence](#10-suggested-implementation-and-merge-sequence), and
[release-side tests](#111-tests-the-release-chat-can-implement-without-the-street-backend).
The [pipeline audit](#3-what-the-current-street-pipeline-actually-does),
[publication contract](#6-publication-the-principal-street-specific-constraint),
and [source map](#12-source-map-for-the-implementing-chat) support the subsequent
street integration. Finish with the [merge checklist](#13-merge-readiness-checklist).

## 1. Scope, source version and handoff

The immediate work in the neighboring chat is progressive import for all OSM
feature types supported by `release`, including buildings. This document defines
the integration seams needed so that merging that work into `buildings_gn` does
not require redesigning the scheduler for streets.

| Item | Audited state |
| --- | --- |
| Street checkout | `C:/dev/workspace/blosm_streets_dev` |
| Development branch | `streets_development` |
| Reference branch | `buildings_gn` |
| Development HEAD and reference commit | `b9517de05889a3a7f0d7bc8107d923f7af871f41` |
| Local `release` ref inspected | `1bf4bb2a68d4f5c9546b282cec59e6857eb7d163` |
| Audited Python street runtime | `0.1.0.dev81`, including uncommitted performance work |
| Recorded qualified runtime digest | `c6cc49fd129e3e6b08f2cbc9b67b5640b2d2ebac4658b571ef7bfcf93e55309d` |
| Street Blender target | 5.2 or newer; do not impose this requirement on the shared release importer |

The runtime digest and measurements are from the completed
[performance qualification](STREETS_PERFORMANCE_IMPLEMENTATION.md), not a new
progressive-import test. The checkout contains uncommitted qualified runtime,
test and documentation changes. Reading only the commit above will miss them.
The two branch refs were equal at the audit; no synchronization merge was needed.
No runtime code is changed by this handoff.

The neighboring chat must use its own checkout, development branch and isolated
Blender profile. It should read this document, not develop in this chat's folder.
Source links below are relative to the repository and refer to symbols rather
than relying on line numbers that will change during the merge.

Terminology in this document:

- **Current** describes code that exists in the audited checkout.
- **Required** defines a compatibility or correctness constraint for the new
  importer or its eventual street adapter.
- **Proposed** describes an implementation choice; equivalent designs are fine.
- **Release owner** means the shared progressive OSM implementation.
- **Street owner** means the later adapter/refactoring in `buildings_gn`.

The release implementation does **not** need to import the Python street solver,
copy its task graph into shared code, restore a GN backend, or implement new
street geometry. It does need to provide a tested extension point for it.

## 2. Essential requirements for the release implementation

These requirements are the merge contract. Later sections explain the reasons
and the street-specific work that remains.

| ID | Requirement | Owner |
| --- | --- | --- |
| P01 | Schedule parsing/preparation, generation, validation, scene publication and finalization as distinguishable phases. A progressive render loop alone does not make OSM acquisition or processing progressive. | Release |
| P02 | Support a bulk/network producer whose input is complete prepared manager data, alongside ordinary per-element producers. | Release |
| P03 | Preserve manager/action dependencies and branch-specific lifecycle ordering. Do not interleave building processing ahead of required way preparation. | Release and merge |
| P04 | Offer cooperative checkpoints inside a producer, with cancellation and progress independent of visible objects. Do not restrict yields to the outer object loop. | Release API; street implementation later |
| P05 | Separate work progress, preview changes, and successful publication/commit. A completed numerical task is not necessarily a usable Blender object. | Release |
| P06 | Keep Blender data mutations on the main thread. Any worker receives only explicitly serialized/copied data and returns results for main-thread validation/publication. | Release and street |
| P07 | Retain the blocking/background path using the same work producer and algorithms. A fallback synchronous producer remains valid, with its limitations reported honestly. | Release |
| P08 | Own new datablocks explicitly across yields. Never roll back a yielding operation by deleting everything absent from a global initial `bpy.data` snapshot. | Release and street |
| P09 | Allow a producer-defined atomic commit unit containing many objects; cancellation must not expose a half-replaced network or an advanced manifest with old geometry. | Release API; street commit logic |
| P10 | Protect job inputs, target scene, terrain and pending outputs from stale-reference publication. Provide revision validation and an active-import guard usable by later street operators. | Release and street |
| P11 | Preserve undo/redo, file-load, disable, window-close and pause behavior, with exactly-once terminal cleanup and explicit partial-result reporting. | Release and street |
| P12 | Preserve deterministic geometry, materials, profiles, UVs and source identities. Scheduling must not change construction choices or discard checks. | Both test suites |
| P13 | Keep memory bounded; do not prebuild all future mesh results or copy the whole scene per scheduling step. | Release and street |
| P14 | Keep backend selection, feature registration and asset loading outside the generic scheduler. | Release and street |
| P15 | Preserve the scoped Import performance request around active callbacks and finalization, restoring the previous policy while idle/paused. | Release/merge |
| P16 | Measure time to first visible result, worst callback time, cancellation latency, throughput and memory. Object counts alone do not establish responsiveness. | Both test suites |

## 3. What the current street pipeline actually does

### 3.1 Import preparation happens before geometry generation

The current operator is synchronous. The significant call chain is:

```text
BLOSM_OT_ImportData.execute / importOsm
  initialize app and selected setup; acquire OSM if necessary
  create parser/managers/layers; establish projection
  parse OSM; load and process missing relation members
  establish terrain and initialize layers
  App.process
    managers in registered order
      WayManager.process
        initialize ways
        actions, including StreetGenerator.do
  App.render
    Renderer.begin
    all renderer.prepare
    all manager.render
      WayManager.render
        StreetRenderer.render
          selected backend: PythonStreetRenderer
    Renderer.end
    all manager.renderExtra
    single-object layer finalization when enabled
    renderer.finalize; renderer.cleanup; manager.cleanup
  finalize scene geographic origin; App.clean
```

Sources: [OSM operator](../__init__.py),
[base App processing](../app/__init__.py),
[Blender App rendering](../app/blender.py),
[WayManager](../way/manager.py),
[StreetGenerator](../action/generate_streets.py).

`StreetGenerator.do` operates on a network, not on one final road at a time.
It finds intersections, splits and joins source segments, creates Sections and
Streets, resolves connectors, dimensions and Bundles, and retains physical
source information. Parsing a way does not mean its final carriageway footprint
or junction ownership is known.

`setup/realistic_dev.py` explicitly registers way preparation before building
processing: building facade analysis needs initialized ways. This applies even
when roads are needed for facade analysis but road mesh import is disabled.
Preserve this dependency; the scheduler must not equate “do not render streets”
with “skip all way preparation”. See
[realistic setup](../setup/realistic_dev.py) and
[facade setup](../setup/__init__.py).

### 3.2 Source identity and geometry ownership are different

The input vocabulary is:

| Entity | Meaning relevant to scheduling |
| --- | --- |
| OSM way/relation | Source data and provenance; not necessarily a final output unit |
| Street | Ordered logical street structure, potentially containing several Sections and connectors |
| Section | Authored directed centerline and dimensions/profile/interpolation parameters |
| Endpoint incidence | `End(section_id, end)`, with `end` equal to 0 or 1; both ends may belong to junctions |
| Intersection | Set of endpoint incidences and corner parameters |
| Connection | Explicit two-end Section seam, including dimension/frame transitions |
| Bundle | Ordered member traversals and possible interior Sections; coupled pavement/separator construction |
| Component | Connected group through intersections, connections **and Bundle membership** |
| Task | Dependency graph work, identified by `Key(kind, id)`; often no Blender output |
| Part / AnchorTransform | Publishable mesh product / furniture placement result |

See [model](../way/street_network/python_backend/model.py),
[providers](../way/street_network/python_backend/providers.py), and
[`components()`](../way/street_network/python_backend/adapter.py).

A grouping choice in Blender is not an ownership boundary. Source Sections can
be grouped in one source object, per Street, or in mixed arrangements. Generated
ownership is recorded using network/controller identity and product keys, not
object names or mutable vertex indices. Initial imports use per-Street source
grouping. The scheduler must not allocate new Section or junction IDs according
to timer order or derive dependencies from object names.

### 3.3 Backend adaptation and connected components

`PythonStreetRenderer._render` currently performs:

1. Read `PreparedNetworkAdapter.read(manager)` into a snapshot with original
   centerlines, provenance, logical/physical ownership and boundary context.
2. Select the physical inventory using `PhysicalStreetInventory.view`.
3. Convert with `PreparedPythonAdapter.read` to the numerical `Network` model.
4. Resolve the shared material library and create the import collection.
5. Iterate components determined by declared topology and Bundle membership.
6. Run `ClippedApproachSupport.prepare` for each component, create its controller
   and editable sources, configure terrain, then call `NetworkScene.update`.
7. Record per-component READY/FAILED results and aggregate COMPLETE/PARTIAL
   coverage in a Blender text report and collection properties.

Sources: [renderer](../way/street_network/python_backend/renderer.py),
[prepared snapshot](../way/street_network/adapter.py),
[physical inventory](../way/street_network/physical_sources.py),
[Python adapter](../way/street_network/python_backend/adapter.py),
[source support](../way/street_network/python_backend/source_support.py).

The boundary support phase can retry construction and extend an approach using
retained source context. It is not a trivial metadata conversion. Preserve the
policy of completing supported real junctions near the import boundary, with
the existing bounded extension/error behavior. Progressive batching must not
clip every piece independently or turn a pending boundary junction into a
permanent dead end.

Logical aggregation may have absorbed a Street while its physical centerline
remains necessary. Retained physical Streets, junction origins and Bundle
members must survive preparation and all consumers. Do not free a manager just
because its first visible output appeared.

### 3.4 Numerical work already has an explicit dependency graph

`UpdateCoordinator.prepare` constructs tasks from registered feature providers,
checks missing dependencies, sorts them deterministically, fingerprints source
and dependency signatures, and executes changed task builders in order. It
returns a `PreparedUpdate`; `commit` advances accepted state separately.

The graph below is a conceptual dependency map, not an alternative hard-coded
execution order. **The actual `Task.dependencies` remain authoritative.**

```text
authored Sections / intersections / connections / Bundles / settings
          |
          +--> paths, profiles, routed approaches, physical footprints
                            |
                            +--> compound plan (shared ownership)
                                      |
                       final layouts and dependent surface tasks
                         (different dependencies; see table below)
                                      |
                            complete prepared products
                                      |
                         UV fallback + overlap validation
                                      |
                    terrain corridor and reference fitting
                                      |
                    stored-precision publication validation
                                      |
                    stage objects --> coherent commit
```

Sources: [updates](../way/street_network/python_backend/updates.py),
[backend](../way/street_network/python_backend/backend.py),
[intersection tasks](../way/street_network/python_backend/intersections.py),
[road and furniture tasks](../way/street_network/python_backend/providers.py),
[Bundle tasks](../way/street_network/python_backend/bundles.py).

The main parallel/dependent branches after footprint preparation are:

| Product | Important inputs, in addition to its authored/local data |
| --- | --- |
| Junction surfaces | Compound plan and final layouts |
| `junction_surfaces_ready` | All multi-ended junction surfaces |
| Separator port index | Layouts, and where applicable compound plan and the surface barrier |
| Bundle plans | Member paths/profiles and relevant layouts/compound plan; angular/internal cases also need the port index; not every Bundle waits for all surfaces |
| Roads | Local path/profile, endpoint layouts, compound plan where applicable, and relevant Bundle plan; no direct global surface barrier |
| Separators | Bundle plan, port index, relevant roads/layouts/transitions and applicable compound/surface dependencies |
| Furniture | Its actual road/transition/Bundle/compound surface dependencies, including the declared surface barrier where applicable |

This table intentionally avoids imposing a new all-surfaces-before-all-roads
barrier. Respect each provider's declared graph.

Important barriers include all physical footprints required by the compound
plan, and the `junction_surfaces_ready` barrier used by downstream consumers.
Other tasks need only local dependencies. Let the provider specify which case
applies; the shared OSM scheduler should not know these task names.

The products may include invisible planning data, a suppressed/absorbed road,
several visible surfaces for one junction, or one shared pavement for overlapping
Bundle members. Therefore “processed 100 ways” cannot mean “100 completed
carriageways”, and failure of one required task cannot be ignored while its
dependent road is marked complete.

Broad ordering dependencies and narrow invalidation signatures are deliberately
different. A consumer may need the global surface barrier to finish but use only
its actual owner's surface in its cache signature. Preserve both contracts.
The current synchronous loop calculates a digest before calling the builder;
in a resumable version, that must not make the task appear finished. Publish
the completed product/digest together to downstream tasks only after success.

### 3.5 Validation and terrain are part of construction

`PythonMeshBackend.prepare` applies the fallback UV strategy where a provider did
not supply one and validates overlaps across the prepared products. Publication
then constructs `ScenePlacement`, which may enlarge the changed-output set.

The current terrain modes are:

- No terrain: keep the prepared construction geometry.
- **PROFILE:** derive elevation references, solve contacts and grade/clearance
  constraints, preserve full carriageway profiles and coordinated junctions.
  A Bundle and junctions sharing its protected seams can require one rigid
  elevation. Connected road references are fitted consistently.
- **DIRECT:** publish the Python mesh with a small shared GN projection modifier,
  preserve prepared relief, and validate initial terrain coverage. Subsequent
  modifier evaluation is Blender work. This is not the retired GN street solver.

See [terrain specification](STREETS_TERRAIN.md),
[placement adapter](../way/street_network/python_backend/terrain/blender.py),
[corridor](../way/street_network/python_backend/terrain/corridor.py), and
[alignment](../way/street_network/python_backend/terrain/alignment.py).

Terrain fitting cannot be run independently on arbitrary object batches. A new
constraint can raise an adjoining reference or require replacing an unchanged
neighbor. Conversely, lowering terrain must allow the fit to lower; no cached
fitted surface may override that result. Existing caching stores bounded local
constraints, not a permanent final fit.

UVs describe the prepared surface in physical plan distances and survive vertical
placement. Do not unwrap each visible batch again or alter its origin/scale.
Furniture is attached through Section identity, station and side/offset, not
through generated face indices; it follows the same final reference.

The future raised-strip/curb extension can supply `terrain_corridor` metadata
that constrains an existing elevation reference. Such features must contribute
their constraints **before** that reference is fitted; they are not necessarily
independent decorations that can be added after terrain placement.

## 4. Existing progressive infrastructure and merge differences

The repository already has a general scheduler used by 3D Tiles:

| Existing component | Current responsibility | OSM integration requirement |
| --- | --- | --- |
| `progressive.core.Job` | Drives a generator; handles transport results, completion and cancellation | Accept local CPU work without requiring a network request; retain existing transport behavior |
| `progressive.blender.Controller` | Timer lifecycle, active job, context, pause, undo checkpoints, termination | Reuse for OSM through an adapter rather than another independent global timer system |
| `threed_tiles.progressive.TilesImport` | Adapter with `job`, `scene`, `view_layer`, `revision`, `status`, `cancel`, `finish` | Example of adapter shape, not the street unit-of-work model |
| `downloads` | Out-of-process transport and request protocol | OSM acquisition/missing-member requirements need their own qualification |
| `app.import_performance` | Scoped temporary performance request | Keep around each active callback and finalization, not across idle time |

Sources: [scheduler](../progressive/core.py),
[Blender driver](../progressive/blender.py),
[existing adapter](../threed_tiles/progressive.py),
[progressive documentation](../progressive/README.md),
[performance request](IMPORT_PERFORMANCE.md).

### 4.1 Do not confuse generator yields with UI turns

Currently `Job.step` has an 8 ms default soft budget and at most 32 generator
advances. `yield None` provides a cooperative checkpoint to that loop; it does
**not** force the Blender callback to return immediately. Deferred download
consumers currently return from the step after publication, but local CPU steps
do not have an equivalent explicit commit boundary.

For OSM, provide a documented way to end a callback after a visible scene change
when required. Either extend typed scheduler commands compatibly or add an
equivalent runner rule. Do not yield an arbitrary progress/event object into the
current `Job`: any non-`None`, non-`Drain` value is currently treated as a transport
request. Progress can initially be adapter state updated before `yield None`.

One expensive `task.build` or Blender operation can exceed the soft budget. A
timer cannot interrupt it. The API must allow a producer to expose finer safe
steps without changing its mathematical termination conditions.

### 4.2 OSM acquisition and parsing need explicit handling

Currently `App.initOsm` can download the main OSM file before rendering starts.
The operator subsequently parses it and optionally downloads missing relation
members. Starting the progressive scheduler only at `App.render` leaves those
earlier pauses untouched.

The current generic download request supports GET-style transport; missing OSM
ways use a POST body through `App.downloadOsmWays`. Account for this difference
when extending shared transport. Preserve online-access checks, cancellation,
error reporting and request-body semantics. Do not silently drop missing-member
completion when using progressive mode.

Large parsing/preparation calls need incremental parsing or a defined helper
protocol, followed by main-thread installation of results. The chosen parser
strategy belongs to the release implementation. Required for streets: no
publication based on incomplete topology unless a separate provisional-preview
contract explicitly permits and subsequently replaces it. That preview is not
needed for the first street integration.

### 4.3 Lifecycle differences that must be reconciled during the merge

At the audited refs, `release` calls renderer preparation before `Renderer.begin`
and lacks some later finalization/manager-cleanup stages present in
`buildings_gn`. The current street-development ordering is shown in section 3.1.
Do not overwrite it with the older release order during conflict resolution.

`Renderer` also keeps shared class state such as its collection and join list;
the app is a singleton. Preserve registered manager order initially and serialize
imports. Arbitrary round-robin execution of different managers is not yet proven
safe. Explicit inter-manager dependencies can be introduced later.

Legacy single-object layers hold a BMesh until finalization. Yielding between
their OSM elements does not necessarily show new geometry. The release
implementation needs a repeat-safe partial flush/preview strategy if it promises
visibility in that mode. Do not repeatedly call full `Layer.finalize`: it can
repeat slicing, modifiers or other completion work. Likewise, current join
targets include name-based lookup; capture/validate datablock identity across
yields so rename/deletion cannot redirect a join to an unrelated object.

Required lifecycle rules:

1. Initialization and preparation happen once per session, not once per tick.
2. Do not repeat `StreetGenerator.do` or rebuild topology at each timer callback.
3. `renderExtra`, layer finalization, renderer finalization, joining and cleanup
   remain separate hooks, called in the correct order only when their inputs
   are ready.
4. Finalization may itself be long; it must be reportable and optionally
   cooperative. Existing final joins remain non-preemptible operations until
   changed explicitly.
5. On cancellation/error, dispose initialized resources exactly once without
   blindly invoking success-only finalization.
6. Shared terrain/managers/assets remain alive until all consumers finish.
7. Geographic origin must be established consistently if partial imported data
   is retained. Do not wait until successful whole-import completion to record
   the origin needed by already committed geometry.

Also preserve the projected extent initialization in `BaseApp.process` on the
merge target. The non-Blender `app/command_line.py` rendering path must remain
usable without Blender timers or an interactive window.

## 5. Proposed shared scheduling contract

The names below are illustrative. They are **not existing APIs** and should not
be copied into runtime without aligning them with the release implementation.

### 5.1 Session and producer responsibilities

Keep three responsibilities separate:

| Responsibility | Suggested contents |
| --- | --- |
| Import session | Captured settings, projection/units, target scene/view layer/window, cancellation state, progress/report sink, ownership ledger, lifecycle phase |
| Feature producer | Manager-specific preparation and generation; exposes safe checkpoints and declared publication units; owns its numerical dependencies |
| Scene transaction | Staged datablocks, validation state, commit/rollback, manifest and ownership changes |

The session owns the complete import, not a single building or road. A producer
may yield no visible objects for many steps. It may return several objects at
one commit. The shared layer must not require a per-way `render(element)`
signature from bulk renderers.

An optional `process_steps(session)` / `render_steps(session)` convention, or an
equivalent registered adapter, would fit the existing managers and street facade.
The first implementation can provide a synchronous fallback that executes the
existing method once. Do not make a new abstract method mandatory for every
existing backend or custom renderer in the initial merge.

The street facade must delegate an optional progressive capability to its
registered backend. A future GN backend may implement a different producer using
the same lifecycle. Import setup must remain independent of backend technology.

See [backend registry](../way/street_network/backends.py),
[facade](../way/renderer_streets.py), and
[backend extension documentation](STREET_BACKENDS.md).

### 5.2 Information a producer must be able to expose

| Information | Purpose |
| --- | --- |
| Stable session/producer identity | Ownership and diagnostics; not derived from mutable names |
| Phase and current work label | Distinguish “building junction surface” from “publishing objects” |
| Completed work / optional total | Totals can be unknown until preparation or task discovery finishes |
| Visible revision | Trigger redraw/undo handling when scene state changes, including provisional previews |
| Committed-unit count | Track usable completed results independently of preview changes |
| Errors and warnings | Aggregate reporting without a popup for each failed task |
| Input revision validation | Reject stale preparation before continuing/publishing |
| Cooperative cancellation point | Stop only at a consistent checkpoint; clean pending work |
| Terminal outcome | Success, partial result, cancellation, failure; include retained results |

If only the existing single `revision` field is retained, it must still change
for every undo-relevant visible mutation. Keep committed counts as separate
adapter/report metadata. A preview is undo-relevant but is not a completed unit.

Use a deterministic producer in both blocking and progressive mode. Blocking
mode drains it without scheduling UI turns. It must not bypass validation,
terrain placement, rollback or success finalization. Existing scripts must be
able to request blocking behavior explicitly.

If another feature provider uses random style/material choices, preserve its
seed and decision order in job-owned state. User actions or another addon between
callbacks must not perturb those choices through shared global random state.

In progressive mode, the operator returning `FINISHED` means “job scheduled”,
not “all meshes ready”, as already documented for 3D Tiles. Expose reliable
completion/outcome state for tests and callers; do not make callers infer
completion from a collection appearing or a callback returning once.

### 5.3 Scope context changes to one callback

Selection restoration, active object/collection, cursor restoration, Blender
context overrides and temporary performance requests must cover one synchronous
chunk only. A context manager surrounding a generator `yield` does not exit at
that yield and therefore does not restore state while the user interacts.

Never keep such a scope open across UI turns. Re-establish the captured import
context for the next active chunk and restore it before returning. Preserve user
selection changes made between chunks; do not repeatedly restore a stale
selection snapshot captured at the beginning of the entire import.

Capture actual values and validated datablock identities, not just the current
`App` object. Some current functions, including terrain setup, read live RNA
settings later. Audit those reads when extracting phases. Do not rely on the
original Import operator (`app.op`) remaining valid after it returns; give the
session an independent progress/report sink.

### 5.4 UI, background mode and optional workers

Expose progressive mode for OSM with data-type-neutral status and cancellation
wording. The current checkbox and tooltip are tile-specific. Define whether a
request without an interactive window is rejected or explicitly selects blocking
mode; preserve `commandLineMode` behavior and test the actual Import operator.
Define Object Mode requirements without resetting an unrelated user's context.

Local-file/CPU-only work must not require a network helper. A no-op/composite
transport or a compatible scheduler extension can satisfy that requirement.
If parsing or numerical work uses a helper, specify bounded input/output queues,
task/revision IDs, pause backpressure, cancellation acknowledgement/termination,
stale-result rejection, malformed/oversized reply handling and child reaping.
No Blender references may cross that protocol. Worker use is optional; it is
not implied by the existence of an HTTP downloader.

Retain each branch's intended initialization/storage validation. In particular,
the inspected release requires an existing configured OSM data directory; do not
inadvertently import a different fallback policy while moving acquisition out of
`initOsm`. Preserve extension namespace/relative imports and the release's own
Blender compatibility floor.

## 6. Publication: the principal street-specific constraint

### 6.1 Current publication and its guarantees

`SceneWriter.publish` currently:

1. Resolves outputs owned by the controller and rejects invalid/duplicate keys.
2. Applies terrain placement and obtains the actual changed-output closure.
3. Validates stored-coordinate geometry, including quad edge ratio and overlap
   against relevant unchanged outputs at Blender's float32 representation.
4. Creates replacement meshes/objects detached from the destination collection,
   assigns materials and ownership/terrain metadata, and records staged objects.
5. Links replacements, writes source data, manifest, timings and terrain state,
   then commits the engine's prepared numerical state.
6. Removes retired outputs, releases published numerical mesh buffers, and
   refreshes material assignments.

Publication uses the fitted/placed products, but the engine commits
`placement.raw`, the unplaced numerical update. Preserve this distinction in
resumable state: committing fitted heights as raw geometry would accumulate
terrain offsets on later Updates.

Source: [SceneWriter and NetworkScene](../way/street_network/python_backend/scene.py).

The staging/link/controller/engine block has rollback. Retirement, mesh-buffer
release and material refresh occur after that block. **The whole function is
not currently a resumable transaction.** Do not insert a yield in its middle and
assume existing exception handling covers cancellation at that point.

In particular, qualify failures in retirement/material refresh separately from
precommit rollback. After a coherent commit has succeeded, a disposal error must
not trigger a global rollback that deletes previously accepted units. Define
whether remaining cleanup is retried or reported as a retained-result error.

`NetworkScene.create` creates controller, source and output collections/sources;
it relies on the surrounding import transaction for failure cleanup. Those
objects also need explicit ownership when converted to progressive operation.

### 6.2 Global difference rollback must not span UI turns

[`ImportTransaction`](../way/street_network/transaction.py) currently snapshots
all objects, collections, meshes, curves and texts in `bpy.data`. On an uncaught
exception it deletes newly appearing datablocks. Its docstring explicitly says
the transaction is synchronous. `PythonStreetRenderer.render` encloses the
entire `_render` call in it.

With a yield inside that scope, the user could create an unrelated object while
the import is paused. A subsequent rollback would mistakenly delete it.

**Required:** before exposing yields from this scope, replace its cross-turn
ownership mechanism with an explicit session/transaction ledger. Record each
datablock immediately when created, its commit status and parent transaction.
Track replacements separately from existing user data. Rollback removes only
pending data owned by that transaction and restores only metadata it changed.
Avoid nested global-difference snapshots spanning different producers.

Shared materials, images and node groups loaded by the asset service may be
reused by other outputs. Do not purge them by name or because their creation time
falls within the session. Preserve the current asset-service ownership policy;
clean only demonstrably private, unused allocations when appropriate.

### 6.3 Conservative first commit unit: a connected component

The existing whole-component generation/publication boundary is the safest
first integration unit. This is an **atomicity boundary**, not permission to
run its entire calculation in one UI callback.

Proposed state sequence:

```text
prepare component inputs
  -> build numerical tasks (many checkpoints, no accepted scene change)
  -> validate prepared surfaces
  -> fit/validate complete terrain corridor
  -> validate stored representation
  -> stage objects (bounded batches)
  -> commit complete component
  -> release temporary buffers
```

For new imports, earlier committed components can remain after cancellation or
a later component failure. Pending component outputs must not be reported READY.
For updating an existing component, retain its last good outputs and manifest
until replacements are validated and ready for the coherent swap.

Current per-component `ValueError`/`RuntimeError` failures are recorded and the
renderer continues, retaining inspection sources. Uncaught exceptions trigger
the outer synchronous rollback. Progressive mode needs an explicit revised
policy: retain committed units; roll back pending units; mark the import
PARTIAL/FAILED/CANCELLED as appropriate. Do not accidentally classify a user
cancellation as an ordinary geometry failure and continue with the next network.

Preserve structured diagnostics and source inspection for an expected geometry
failure if that diagnostic state is deliberately committed. A cancelled pending
component need not leave orphan controllers or source collections.

### 6.4 Visibility within a large component

A city extract may have one dominant component. Waiting to finish all numerical
work before any street appears is responsive only if calculation itself yields;
it is not early visible mesh generation. The inspected Berlin case has just two
components, so component-level publication alone is not a complete visual-progress
solution.

A safe next step is **batched provisional display after complete component
validation and terrain fitting**:

- Create/link validated objects in batches, clearly associated with an importing
  component. These are previews until the complete component commit.
- Restrict the first preview implementation to **new** component imports.
  Existing-component Update replacements stay detached until the coherent swap,
  avoiding an overlapping display of old and new pavement.
- Keep previews in a separate ownership namespace/collection, excluded from
  canonical output enumeration, manifest recovery and Update. Do not install
  duplicate or uncommitted `st_py_key` objects in the normal outputs collection.
- Do not make its editable sources/Update commands active during this state.
- Track preview revisions for undo and redraw separately from successful commits.
- On cancellation, remove that component's previews and pending data; keep prior
  committed components. On an existing-component Update, keep the old result
  until a coherent replacement is possible.
- Complete metadata/manifest/engine commit only when all required objects exist.

This displays geometry progressively during publication but does not make early
numerical tasks independently publishable. An alternative detached staging
implementation is valid for the first milestone and simpler to reason about.

Publishing smaller *permanent* groups during construction requires a later
street-specific proof of closure: stable mouths and shared seams, compound/Bundle
ownership, neighbor overlap checks, final terrain references, source provenance
and furniture dependencies. Topological task order alone is not this proof.
The release scheduler should support such groups without trying to calculate
them itself.

### 6.5 Keep the final swap coherent

Once cancellation is accepted, do not advance the engine revision or retire old
outputs. Conversely, do not let a callback stop after advancing the manifest but
before replacing its geometry. Define a bounded, non-yielding final commit with
rollback/recovery for its operations, or design a recoverable multi-step commit.

If final linking, deletion or serialization is too large for the callback
budget, measure it and refactor ownership/publication explicitly. Do not resolve
the pause by weakening the transaction. Postcommit resource disposal can be
scheduled separately only when those resources are no longer visible/referenced
and errors cannot invalidate the accepted result.

## 7. Cancellation, edits, undo and saved files

Required behavior for the merged street path:

| Event | Required behavior |
| --- | --- |
| User navigates viewport | Continue at safe checkpoints; preserve context between callbacks |
| Switch scene/view layer or enter Edit Mode | Preserve the established pause policy; do not mutate the inactive target |
| User changes import panel settings | Current job uses captured settings; next import uses new values |
| User edits pending source, controller transform, units or terrain | Invalidate/restart the affected preparation explicitly, or cancel with a clear reason; never commit stale results |
| User deletes target controller/collection/terrain | Validate identity and stop safely; do not adopt an unrelated object with the same name |
| User invokes Street Update/topology commands | Guard pending/import-owned networks; no reentrant update against an in-flight engine. A conservative first implementation can block street mutation commands while an import is active |
| Cancel | Stop at a safe point; preserve committed units and unrelated edits; clean pending data; report what remains |
| Undo/Redo | Cancel before Blender invalidates references; restore accepted scene state; do not resume stale generators or restart imports automatically |
| Load another file / disable addon / close owning window | Stop timers/helpers, release private resources, clear session caches; no result popup during silent teardown |
| Save during import | Define preview persistence explicitly. Saved accepted components must reopen as valid sources/outputs; pending previews must not masquerade as complete networks |
| New Import request | Enforce one active session across OSM, terrain, tiles and other operations sharing importer state |

The current progressive driver already handles several of these lifecycle events
for tiles. Extend its adapter checks, not a separate conflicting lifecycle.

The existing import exclusion is incomplete across older modal operators: the
imagery path has its own timer and does not claim `progressive.active`. A modal
imagery operation started first can therefore overlap a later progressive job.
The requirement above is to establish shared ownership or explicit mutual guards
for operations using the app singleton, not a claim that every path is already
serialized.

For previews, define an explicit history policy. Currently `Controller.tick`
checkpoints a changed revision, but `Controller.close` does not checkpoint
mutations performed by `cancel`/`finish`. Preview removal on explicit cancellation
or failure must leave an undo-consistent terminal state. Repeated Undo/Redo after
the job terminates must not resurrect an orphan preview-only canonical network.
Either prevent such history states or persist a safe incomplete-state cleanup/
recovery protocol. Silent Undo/load teardown must not push new history entries.

The current street Update operator only checks for a controller; it does **not**
already guard against a future progressive import. That is street integration
work, facilitated by the release API's active-job/ownership query.

The conservative saved-file policy is to persist committed components normally,
tag any previews as incomplete and ensure they are cleaned or explicitly
recoverable on reload. No automatic execution from saved files is required.
Do not serialize live generator objects, Python callables or executable module
paths. The existing trusted feature registry and saved provider IDs remain the
authority for rebuilding.

Current native street Update reconstructs its engine from persisted source data,
manifest and output metadata when needed. It detects a manifest/revision mismatch
after undo/redo. File loading clears process-local engines and terrain-query
caches. Progressive work must preserve these semantics rather than depend on
unsaved in-memory tasks surviving reload.

Resources allocated during session construction need cleanup even if execution
never starts. Closing an unstarted generator does not run its body or `finally`.
Test cancellation/start failure before the first producer advance, including
timer-registration failure; independently owned session cleanup must handle it.

## 8. Responsiveness and performance requirements

### 8.1 Known long work is inside tasks, not only between objects

The optimized street benchmark reduced generation from 70.278 to 41.768 seconds
in alternating controls while preserving the output. That measures total work,
not maximum callback latency. The scene contains 328 meshes / 350,571 quads.

Compound chart partitioning, mesh refinement, slope/quality work, separator
candidate construction and terrain fitting can each contain long synchronous
loops. One optimized separator still took about a second in isolated replay;
ten investigated compound builds each took at least 0.25 s. These are examples,
not bounds for arbitrary networks.

Terrain sampling already has internal batches, and some functions are Python
generators. They do not currently yield to Blender's scheduler. NumPy vectorized
operations and Blender mesh/BVH operations remain indivisible from the UI's
perspective until their caller returns.

Required street checkpoints should eventually include:

- preparation/adaptation and component boundary-support attempts;
- dependency task completion, plus inner iterations of expensive builders;
- UV batches and complete-coverage overlap validation batches;
- terrain sampling, constraint assembly/propagation and final placement;
- detached mesh/object creation, materials and metadata installation;
- publication/finalization and safe disposal.

Numerical iteration counts, tie-breaks, candidate order and stopping tolerances
must remain independent of elapsed time. A deadline may pause an iteration; it
must not accept a lower-quality surface because time expired.

### 8.2 Instrument the scheduling layer

Report separately:

1. Input acquisition/parsing/preparation time.
2. Numerical generation and terrain/validation time.
3. Publication and finalization time.
4. Active processing time versus wall time including pauses/UI scheduling.
5. Time to first visible preview and to first committed result.
6. Callback p50/p95/maximum, largest uninterrupted operation and its label.
7. Cancellation request-to-stop latency.
8. Peak memory, including temporary arrays, Blender meshes and undo snapshots.

An 8–16 ms scheduling target is a reasonable initial soft budget, consistent with
the current runner. It is **not** an asserted maximum callback duration. Establish
and publish measured UI/cancellation targets during implementation; long atomic
operations must be named and measured rather than hidden behind an average.
Qualification must use actual timer-driven Blender interaction, not only a loop
that drains the producer in background mode.

Progressive mode is primarily a responsiveness feature. Extra redraws, dependency
graph evaluation, undo snapshots and artificial waits can increase total time.
Batch small tasks and checkpoint visible transactions sensibly; do not create an
undo snapshot for every invisible numerical task. Preserve the existing Import
performance preference in both blocking and timer paths, and restore previous
settings between callbacks and when paused/cancelled/failed.

Compare blocking/progressive throughput with the same input, runtime, materials,
terrain, quality, cache state and performance-request setting. Use alternating
runs because the earlier qualification observed timing drift. Agree any allowed
overhead from evidence rather than assuming progressive execution is faster.

### 8.3 Bound working memory

Do not convert all task iterators to lists of precomputed geometry, retain every
published mesh twice, or prepare all connected components concurrently by default.
Preserve existing bounded pattern/path/terrain caches and mesh-buffer release.
Cache invalidation must include all inputs that currently affect signatures.

Keeping a complete current component's raw mesh set can be necessary for overlap
checks, corridor construction and terrain fitting. The restriction is against
eager whole-import accumulation, not that required barrier. Prepared products
can share/alias nested mesh payloads with query surfaces and cached products;
do not free them merely because their own task finished. Pin data required by
downstream Bundle/separator/furniture consumers, and release only at a proven
lifetime boundary. Raw mesh hydration for terrain must retain its UV and relief
attributes.

The qualified example used approximately 605 MiB resident / 679 MiB sampled
private memory in background import. A GUI import with undo and redraw can use
more; this is not a product memory guarantee. Development qualification retained
a 3 GiB process cap and 1 GiB system reserve. Do not assume users can provide
7–8 GiB of spare memory or remove checks to fit more pending work.

The shared scheduler must permit bounded queues/backpressure and disposal on
cancel. Parallel numerical execution is optional future work, not a prerequisite
for progressive rendering. Do not run `bpy` or mutable scene/manager access in
background threads to work around latency.

## 9. Geometry, assets and extension contracts to preserve

Changing scheduling must preserve:

- Original source centerlines/provenance and bounded local interpolation behavior;
  both Exact polyline and Smooth interpolation modes.
- Rounded junction corners, incoming profile rises, common edge minima and the
  current crown/grade/quality constraints.
- Directly constructed quads and current resolution/shape checks. No reduction
  of density, validation, solver iterations or terrain samples to meet a timer.
- T-junctions, general three-arm, oblique, four-arm, five-or-more-arm junctions,
  unequal widths and divided carriageways with separators; complex compounds,
  overlapping Bundle pavements, internal connections and boundary junctions.
- Existing overlap ownership: avoid duplicate pavement, cracks or unassigned
  gaps while reconciling several Section/Bundle sources.
- UV layer names, per-face-corner ordering, physical scale, normals/shading,
  material slots/indices, object transforms and scene units.
- Terrain mode, prepared relief attributes and modifier assignments, including
  the DIRECT projection modifier on Python meshes.
- Editable source grouping, stable IDs, dependency manifests and the Update
  button's selective rebuilding/last-good-result behavior.
- Furniture anchors and future attachment dependencies. Do not bake their
  identity into temporary face/vertex indices.

The current material file is separate from generator assets, loaded through the
feature-independent asset service. Respect explicit carriageway material
overrides, external libraries, edition availability and shared resource reuse.
Do not add a material library dependency to the generic scheduler. See
[asset libraries](ASSET_LIBRARIES.md), [asphalt materials](STREETS_ASPHALT.md),
and [UV contract](STREETS_UV.md).

Future sidewalks, crossings, lane markings, cycleways, tram/rail surfaces and
terrain adaptation may add tasks and publications. The generic scheduler should
carry opaque feature/product identities and dependencies owned by the producer;
it must not enumerate “road, intersection, grass” as the only supported kinds.
None of those future geometry features is a requirement to implement in the
release progressive-import task.

## 10. Suggested implementation and merge sequence

### R1 — Shared OSM session and lifecycle in the release development branch

- Move the OSM operation behind a session that owns settings, input lifetime,
  reports and completion. Include acquisition/missing members and parsing.
- Preserve a blocking driver and add a progressive driver through the existing
  infrastructure. Add optional manager/renderer producer hooks or adapters.
- Preserve lifecycle ordering and state; document branch differences for merge.
- Add owned-resource transactions, explicit safe cancellation, active-session
  guards and visible/committed progress distinction.
- Keep synchronous fallback support for custom/older renderers. Label large
  fallback calls accurately; they are compatibility, not fine-grained progress.

### R2 — Release feature adapters and contract tests

- Convert release buildings and other supported OSM feature renderers in
  manageable increments, including their postprocessing/joining behavior.
- Qualify actual visible progress, cancellation and undo.
- Add a **synthetic bulk network producer** without importing streets. It must
  consume a fully prepared manager, build dependent hidden work, validate, and
  commit multiple related outputs together. Exercise cancellation before/after
  commit and a failure after earlier commits.
- Merge only after the general contract tests below pass.

### G1 — Merge into buildings_gn without changing street geometry

- Preserve `buildings_gn` manager/render lifecycle additions and the scoped Import
  performance request. Keep street registry/facade and shared asset service.
- Keep the audited street optimization changes; they are not all in the recorded
  parent commit. Reconcile with whatever commits contain them at merge time.
- Initially allow the street renderer to use the synchronous fallback, if needed.
  This is a compatibility milestone only and must not be advertised as responsive
  street generation.
- Run a normal street import and mixed buildings/streets import for regression.

### G2 — Cooperative street construction with component commits

- Split `PythonStreetRenderer._render` and `NetworkScene.update` into resumable
  phases while retaining blocking wrappers.
- Expose checkpointed dependency execution from `UpdateCoordinator`; do not
  change task signatures/order or committed-state semantics.
- Replace the cross-turn `ImportTransaction` snapshot with explicit ownership.
- Add inner checkpoints to measured expensive tasks and terrain/validation work;
  maintain atomic final publication and raw/placed geometry separation.
- Add active-import guards to street editing operators, revision checks and
  complete/partial/cancelled reporting.

### G3 — Progressive display and finer commits where justified

- Batch validated staging/publication and, if useful, offer clearly provisional
  component previews with cancellation/undo handling.
- Measure latency and memory in the actual GUI and resolve remaining long calls.
- Consider smaller permanent publication groups only after dependency, terrain
  and ownership closure is implemented and tested. This is not necessary to merge
  the release scheduler or to obtain responsive construction.

These stages make later street integration localized and predictable. They do
not eliminate the need to refactor current long synchronous numerical calls.
Merely merging an OSM per-object timer loop will not satisfy G2/G3.

## 11. Acceptance tests

### 11.1 Tests the release chat can implement without the street backend

| Test | Required observation |
| --- | --- |
| Individual producer plus bulk producer | Both use one scheduler; bulk work has many tasks and one multi-object commit |
| Scheduler dispatch | Repeated `None` checkpoints obey soft-budget/max-step semantics without transport submission; any new immediate-return/commit command has explicit dispatch |
| Manager dependency fixture | Building-like consumer cannot run before way-like preparation, even if way output is disabled |
| Atomic fallback | Old synchronous renderer runs once with unchanged success/error behavior |
| Blocking/progressive parity | Same deterministic geometry/attributes and final result structure |
| Cancellation checkpoints | Cancel during preparation, hidden work, staging and after one committed group; only pending owned data is removed |
| Unstarted session | Cancel or fail startup/timer registration before first generator advance; already allocated resources are still cleaned |
| Unrelated user edits | Create a sentinel object between steps, force failure, and prove that object/mesh survives |
| Commit boundary | No stale/half manifest or duplicated old/new accepted output; revision signals visible mutation |
| Long task | Recorded callback duration reveals an intentionally long step; soft budget is not presented as preemption |
| Context scope | User selection/navigation between callbacks survives; context/performance overrides are restored before idle |
| Lifecycle | Pause/resume, scene/view-layer switch, target deletion, Undo/Redo, file load, disable and window close clean up once |
| Preview history | Show provisional output, Cancel, then repeatedly Undo/Redo after the job ends; no orphan previews or incomplete canonical network reappear |
| Input changes | Captured settings remain fixed; edited/deleted producer dependencies invalidate safely |
| HTTP/parsing | Main OSM and missing-member paths preserve data/error/cancel semantics, including POST where used |
| Deferred results | Cancel with pending consumers; paused jobs do not mutate the scene; queues remain bounded and stale/helper-error results are discarded safely |
| Success finalization | `renderExtra`, joins, layer and renderer finalization occur once and in order; cancellation does not run success-only work |
| Reporting | Distinguish scheduled, running, successful, partial, cancelled and failed; counts do not treat pending tasks as completed objects |
| Memory/backpressure | Large synthetic producer does not queue all completed meshes or leak pending allocations across cancel/retry |
| Extension namespace | Relative imports and addon loading work under extension package names, not only a top-level `blosm` import |
| Entry policy | Actual OSM Import operator honors progressive/blocking, background/command-line, Object Mode and initialization/storage rules |
| Existing tiles | Shared core/driver changes preserve tile visibility, download concurrency, cancellation, joining and undo tests |

The synthetic producer is especially important: it verifies the architectural
contract in `release` before real street code is available there.

### 11.2 Street tests after the merge

1. Compare blocking and progressive runs of **the same current runtime**, not
   just against an old baseline. Preserve coordinates, quad connectivity,
   normals, shading, UV maps and material assignments. Ignore only explicitly
   transient fields such as new import UUIDs/timings, not geometry metadata
   indiscriminately. No new numerical tolerance is justified by scheduling alone.
2. Exercise the six junction families and four-arm intersections, both
   interpolation modes, bent connected streets, Section seams, Bundle overlaps,
   separators, internal routes and clipped boundary-junction support.
3. Use Berlin Karl-Marx-Allee, Canal Street and Nashville Wyoming Avenue as
   regression inputs, not as hard-coded runtime cases. Include one dominant
   component and several disconnected components. Actual input locations are
   recorded in the existing performance tools/docs.
4. Exercise no terrain, PROFILE and DIRECT placement, source/terrain transforms,
   scene units, terrain raise/lower, missing hits and invalid grade/lift cases.
   Confirm fitted neighbors and furniture commit consistently.
5. Cancel after every phase and selected internal checkpoint; inject a mesh,
   metadata/material, terrain and commit failure. Confirm prior accepted
   components and unrelated user objects remain unchanged, then retry.
   Include cancellation before the first step and after visible previews,
   followed by repeated Undo/Redo with no active job.
6. Save/reopen complete and partial imports; verify no automatic code execution,
   valid manifests, source editing, material-only Update, geometry Update,
   no-op Update and native Undo/Redo. Caches must not be required for correctness.
7. In a factory-profile GUI session, navigate while importing, cancel during the
   slowest observed geometry phase, change scene/Edit Mode, and measure callback
   gaps. Repeat with foreground and background activity under the same Import
   performance setting.
8. Compare controlled overview/junction/separator/wireframe/asphalt renders.
   Scheduling-only changes should preserve the existing visual qualification;
   do not substitute viewport impressions for numerical comparison.

Existing test entry points to extend/reuse:

| File | Relevant coverage |
| --- | --- |
| [`tests/test_progressive_transport.py`](../tests/test_progressive_transport.py) | Shared core/transport, cancellation and helper cleanup |
| [`tests/check_progressive_import.py`](../tests/check_progressive_import.py) | Existing progressive adapter integration |
| [`tests/check_progressive_ui.py`](../tests/check_progressive_ui.py) | Timer-driven GUI, context, grouped Undo/Redo and lifecycle |
| [`tests/check_import_performance.py`](../tests/check_import_performance.py) | Scoped performance policy, including timer steps and failure |
| [`tests/benchmark_python_street_import.py`](../tests/benchmark_python_street_import.py) | Real Import timings and saved scene reports |
| [`tests/check_street_timing_equivalence.py`](../tests/check_street_timing_equivalence.py) | Saved mesh, UV, material and normal comparison |
| [`tests/check_python_street_publication.py`](../tests/check_python_street_publication.py) | Publication failure and rollback |
| [`tests/check_python_street_undo.py`](../tests/check_python_street_undo.py) | Native source/output Undo/Redo |
| [`tests/check_street_performance_lifecycle.py`](../tests/check_street_performance_lifecycle.py) | Terrain-enabled lifecycle, material-only changes, cache disposal and rollback |
| [`tests/check_street_uv.py`](../tests/check_street_uv.py) | Junction families, terrain modes, units, extra UV maps and reopen |
| [`tests/check_python_street_suite.py`](../tests/check_python_street_suite.py) | Numerical fixture/regression suite |

The current passing suites do not validate a future progressive implementation.
Extend them with real cancellation and timer-driven tests rather than claiming
their current results as evidence for the new scheduling behavior.

All generated scenes, renders, timing data and validation reports belong under
`C:/dev/workspace/artifacts`, with street integration results under
`artifacts/streets`. Keep source, this documentation and documentation illustrations
in the repository. Do not stage generated scenes, test binaries or JSON reports.
Use each chat's own factory-startup Blender profile; do not modify or attach to
the user's normal Blender installation/session.

## 12. Source map for the implementing chat

| Concern | File / symbol |
| --- | --- |
| Import entry and captured settings | [`__init__.py`](../__init__.py): `BLOSM_OT_ImportData.execute`, `_execute`, `importOsm` |
| Acquisition, processing order | [`app/__init__.py`](../app/__init__.py): `downloadOsmWays`, `process`; [`app/blender.py`](../app/blender.py): `initOsm`, `render`, `clean` |
| Building dependence on way preparation | [`setup/realistic_dev.py`](../setup/realistic_dev.py), [`setup/__init__.py`](../setup/__init__.py) |
| Way preparation and dispatch | [`way/manager.py`](../way/manager.py): `process`, `render`; [`action/generate_streets.py`](../action/generate_streets.py): `StreetGenerator.do` |
| Backend-neutral boundary | [`way/renderer_streets.py`](../way/renderer_streets.py), [`way/street_network/backends.py`](../way/street_network/backends.py) |
| Snapshot and physical provenance | [`way/street_network/adapter.py`](../way/street_network/adapter.py), [`physical_sources.py`](../way/street_network/physical_sources.py) |
| Numerical input/component partition | [`python_backend/model.py`](../way/street_network/python_backend/model.py), [`adapter.py`](../way/street_network/python_backend/adapter.py) |
| Import component loop and coverage report | [`python_backend/renderer.py`](../way/street_network/python_backend/renderer.py): `PythonStreetRenderer._render` |
| Boundary support attempts | [`python_backend/source_support.py`](../way/street_network/python_backend/source_support.py): `ClippedApproachSupport.prepare` |
| Task ordering and accepted revisions | [`python_backend/updates.py`](../way/street_network/python_backend/updates.py): `DependencyIndex`, `UpdateCoordinator`, `PreparedUpdate` |
| Feature registration and post-task checks | [`python_backend/backend.py`](../way/street_network/python_backend/backend.py): `PythonMeshBackend.prepare`, `release_published_meshes` |
| Junction/road/Bundle dependencies | [`intersections.py`](../way/street_network/python_backend/intersections.py), [`providers.py`](../way/street_network/python_backend/providers.py), [`bundles.py`](../way/street_network/python_backend/bundles.py) |
| Publication, source grouping, rebuild | [`python_backend/scene.py`](../way/street_network/python_backend/scene.py): `SourceObjects`, `SceneWriter.publish`, `NetworkScene`, runtime cache hooks |
| Synchronous-only global rollback | [`way/street_network/transaction.py`](../way/street_network/transaction.py): `ImportTransaction` |
| User editing operators | [`python_backend/ui.py`](../way/street_network/python_backend/ui.py) |
| Terrain placement and revision expansion | [`terrain/blender.py`](../way/street_network/python_backend/terrain/blender.py): `ScenePlacement` |
| Shared terrain constraints and cache | [`terrain/alignment.py`](../way/street_network/python_backend/terrain/alignment.py), [`terrain/cache.py`](../way/street_network/python_backend/terrain/cache.py) |
| Mesh/UV/validation contracts | [`mesh.py`](../way/street_network/python_backend/mesh.py), [`uv.py`](../way/street_network/python_backend/uv.py), [`validation.py`](../way/street_network/python_backend/validation.py) |
| Scheduler and UI lifecycle | [`progressive/core.py`](../progressive/core.py), [`progressive/blender.py`](../progressive/blender.py) |
| Import performance scope | [`app/import_performance.py`](../app/import_performance.py) |
| Generic asset loading | [`app/asset_library.py`](../app/asset_library.py), [`way/street_network/assets.py`](../way/street_network/assets.py) |

## 13. Merge-readiness checklist

Before handing the release implementation back for merging, confirm:

- [ ] A bulk producer can perform many non-visible steps and then atomically
  publish multiple objects without shared code knowing street geometry.
- [ ] CPU-only steps do not need fake download requests; local yield/commit
  semantics are documented and tested with the current `Job` protocol.
- [ ] Blocking and progressive mode share algorithms and lifecycle.
- [ ] Manager ordering and optional renderer hooks preserve the `buildings_gn`
  preparation/finalization/cleanup requirements.
- [ ] Resource ownership across yields is explicit and unrelated user changes
  survive cancellation/failure.
- [ ] Partial results, previews, success, cancellation and failures are distinct.
- [ ] Input invalidation and active-import guards are available to street Update
  and future backend/feature adapters.
- [ ] Undo/load/disable/window-close and context restoration work in actual GUI
  tests, including a pending multi-object transaction.
- [ ] The scoped Import performance request and existing tiles behavior survive.
- [ ] Known non-preemptible phases and measured latency/memory are documented.
- [ ] Remaining G2/G3 street work is explicitly listed; a successful release
  import is not reported as proof of progressive street integration.

With these seams in place, the later street work stays within its adapter,
dependency execution and publication layers. The intersection algorithms,
source contracts, backend registry and shared OSM scheduler can retain their
respective responsibilities.
