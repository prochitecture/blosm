# Proposed Blosm release tests

Prepared on 2026-11-09 for local `release` commit
`a89dc20e97143d7315fd7a99cb42636238c23ae9`.
The clean `blosm_progressive` worktree was at that same commit when inspected.
Local `release` is 16 commits ahead of the locally stored `origin/release` ref;
the remote was not fetched for this review.

GeoJSON is excluded at the owner's request because its code will be removed.
This proposal covers the features present on the inspected release branch.
Tests for Pro features apply to the Pro distribution if it is being published.
New acceptance criteria below are proposals, not claims that the current code
already satisfies them. The initial observations below describe the November 9
review; implemented automation and current commands are maintained in
[tests/README.md](../tests/README.md).

The implemented runner now has `general`, `base` and `pro` profiles. General tests
exercise the checkout. Base and Pro require the actual extension ZIP, install it
in isolated Blender profiles and prohibit checkout runtime imports. Both start
with mandatory package/edition checks. Pro reuses all general functionality;
Base runs the common functionality and dedicated checks that Pro-only modules,
operators and controls are unavailable. Only baking, import-for-export and PML
are inapplicable to Base. Missing common features are failures. GeoJSON remains
excluded from all profiles.

Filesystem paths in this document are relative to the Blosm repository root
unless a different base is explicitly stated. Sibling repositories and artifacts
use `../`; original machine-specific source locations are intentionally omitted.

## What exists today

Most release automation covers the new direct glTF, progressive import and
concurrent download paths. These suites should be reused.

| Existing script | Coverage to retain |
| --- | --- |
| `tests/test_progressive_transport.py` | Seven standalone tests for scheduling, response identity, bounded concurrency, paused backpressure, cancellation, gzip, retries and helper exit. |
| `tests/test_download_resilience.py` | Sixteen standalone tests for deadlines, stalled DNS, incomplete transfers, provider cooldowns, scoped authorization failures, redaction and disconnect cleanup. |
| `tests/check_direct_gltf_import.py` | Direct bytes/file/native operator equivalence, geometry, UVs, colors, textures, attribution, RTC, extensions, invalid input and patch/logger cleanup. |
| `tests/check_direct_gltf_operator.py` | Registered Blosm operator, nested manifests, mixed tile formats, repeated downloads, partial failures and absence of dataset cache access. |
| `tests/check_progressive_import.py` | Incremental geometry, context preservation, joining, cancellation, destination changes, failed-tile cleanup and lifecycle handling. |
| `tests/check_concurrent_tiles.py` | Sequential versus concurrent and blocking versus progressive parity; transformed branches, session queries, out-of-order completion and truncated-response recovery. |
| `tests/check_progressive_ui.py` | Real Blender timers, grouped Undo/Redo, intervening edits, cancellation, Edit Mode pause, scene switching and disable cleanup. |
| `tests/check_live_tiles.py`, `tests/check_live_tiles_ui.py` | Bounded live-provider imports and UI scenarios. |
| `tests/benchmark_direct_gltf_import.py` | Controlled direct-import benchmark against a historical renderer. |

There are also legacy PML translation tests in `pml/tests`. No comparable
automated release suites were found for the OSM, terrain, imagery or GPX import
pipelines. A useful standalone Geometry Nodes test exists outside this repository
at `../blosm_release_gn_compat_test.py`; it has a stale hardcoded
checkout path and should be adapted before reuse.

The 23 standalone tests were run during this review using Blender 4.5's bundled
Python 3.11.11. All passed in 29.866 seconds. An initial sandboxed run could not
establish the helper's loopback connection; an independent subprocess probe
confirmed Windows socket error 10013. The successful rerun allowed those local
connections and used local HTTP fixtures, with no live-provider requests.
Blender integration, GUI, live-provider and benchmark suites were not rerun
during this review. Earlier results in `downloads/LiveValidation.md` and
`threed_tiles/doc/DirectImport.md` are historical evidence, not a new release pass.

Current logs:

- `../artifacts/blosm_release_tests/2026-11-09_a89dc20_audit/standalone-tests-loopback.log` — successful run.
- `../artifacts/blosm_release_tests/2026-11-09_a89dc20_audit/standalone-tests.log` — sandbox-limited attempt.

The owner supplied a real GPX fixture, copied into
`../artifacts/blosm_release_tests/fixtures/gpx/2026-04-16_13-16_Thu.gpx`.
The synthetic GPX and its validation report were deleted at the owner's request.
The real GPX copy has been verified by SHA-256; its Blender import has not yet
been tested. The owner's Buildify setup was also moved to the fixture library;
see the reusable input locations below.

## Proposed test families

P0 means a release-blocking smoke test or regression suite. P1 adds edge cases
and broader feature coverage before publication. P2 establishes longer-running
quality and performance checks. Promote a P1 case to P0 when it exercises a
feature changed in the final release candidate.

| ID | Priority | Test family and concrete acceptance criteria | Existing coverage |
| --- | --- | --- | --- |
| 01 | P0 | **Distribution, installation and registration.** Install each actual release ZIP into an isolated Blender profile with the source checkout absent from the import path. Enable, disable and re-enable twice. Require the expected operators, properties and panels, loadable bundled assets, no duplicate handlers/timers, and no traceback. Check that the ZIP's version and minimum Blender version match the intended release. Include paths with spaces and Unicode. | New; current scripts mainly import a checkout. |
| 02 | P0 | **Direct GLB/B3DM correctness.** Reuse the existing byte/file/operator comparisons. Require correct topology, world transforms, UVs, colors, material assignments, packed texture content and attribution, including absent/zero/nonzero RTC centers. Invalid tiles must release patches and logging handlers. Production 3D Tiles imports must not read/write dataset caches or invoke the glTF operator. | Strong existing suite. |
| 03 | P0 | **3D Tiles execution modes.** Import one deterministic mixed GLB/B3DM hierarchy with progressive on/off, concurrency off/one/four workers, and join on/off; add targeted eight-worker cases. Compare expected tile identities, final geometry and material/texture content. Completion must wait for all submitted content. Only joining should change the expected object grouping. | Strong existing suite; extend the boundary configurations. |
| 04 | P0 | **Progressive UI, cancellation and Undo/Redo.** Require incremental tiles and timer activity, preserved selection/active object/cursor/collection, cancellation that retains committed tiles, and one Undo/Redo for an uninterrupted import including its final join. User edits and consecutive imports keep separate undo groups. Verify pause/resume for mode, scene and view-layer changes; destination deletion/unlink stops cleanly, while rename/move within the scene preserves it. Check actual file loading, window closure and addon disable cleanup in isolated GUI processes. | Strong existing core/UI tests; add missing real lifecycle transitions and Global Undo disabled. |
| 05 | P0 | **Download and tile failure recovery.** Reuse the 23 standalone tests and Blender adapter failure checks: 429/503, deadlines, hung DNS, truncated bodies, gzip/size errors, provider authorization, helper startup/crash and malformed tiles. Assert bounded retries/queues, no late scene writes after cancellation, one aggregate tile-failure result, redacted diagnostics, preserved committed geometry and exited helper processes. | Strong existing suite; keep transport and Blender adapter checks together. |
| 06 | P1 | **Tile selection, URLs and independent coordinate checks.** Use synthetic nested manifests with `uri`/legacy `url`, relative/root-relative/absolute URLs, session changes, box/sphere/region bounds and geometric-error thresholds. Assert the exact requested tile set at inside/outside/boundary cases. Check nested translation/rotation/scale against independently calculated world positions; include scaled spheres and touching boxes. Cover empty and malformed manifests. Document or reject geographic extents outside supported behavior. | Partial fixtures exist; independent selection/geodetic assertions need adding. |
| 07 | P0 | **OSM parsing and relations.** Supply small local OSM fixtures containing nodes, closed ways, split/reversed multipolygon members, inner courtyards, several outer rings, building relations and building parts. Require the correct element grouping and holes. Supply missing relation members through a local fixture server; test recovery enabled/disabled and failure. Empty, duplicate and invalid geometry must have a defined skip/error result and must not corrupt the next import. | New. |
| 08 | P0 | **Simple buildings and roofs.** Import known footprints in 3D simple mode with explicit heights, levels, minimum heights/levels and roof parameters. Test every roof type registered by `BuildingRenderer.initRoofs`, using rectangular, concave and courtyard fixtures where appropriate. Require finite coordinates, valid indices, expected footprint and vertical bounds, correct roof/material assignment and the documented fallback for courtyard/unsupported roofs. Compare single-object and separate-object geometry. Fix randomness for missing-tag cases. | New. |
| 09 | P0 | **2D OSM and other layers.** Import building footprints, water, forest/vegetation areas, highways, paths and railways. Exercise each feature toggle and single/separate-object settings. Require expected meshes/curves, preserved polygon holes, correct layer assignment, no geometry from disabled layers, and expected terrain projection/subdivision. Test the release branch's existing road renderer. | New. |
| 10 | P0 | **Geometry Nodes integration.** Adapt the external helper tests for socket lookup, literal/attribute switching, false/zero/empty values, per-modifier isolation and evaluated geometry updates. Add an actual 2D OSM import with a small GN fixture: `building:levels` must be an INT attribute on the FACE domain with the right value per footprint, and the setup's mode input must be set correctly. Test missing library/group and incompatible sockets. | External helper tests available; end-to-end import is new. |
| 11 | P0 | **Terrain construction.** Generate synthetic gzip-compressed HGT data with known elevations, a slope and void samples. Test 1/3 arc-second settings, representative reduction ratios including extremes, quads/triangles, and an extent spanning four source tiles. Assert dimensions/counts, height values, `height_offset`, consistent face winding, shared seam vertices and no seam gaps. Missing/corrupt HGT must yield a defined error without corrupting existing objects. | New. |
| 12 | P0 | **Shared geographic origin and terrain alignment.** Import terrain, OSM, GPX and 3D Tiles around known control points in several orders. With `relativeToInitialImport` enabled, require stable scene latitude/longitude, horizontal alignment and correct terrain height offsets; with it disabled, verify each import's documented new-origin behavior. Include nonzero terrain transforms, negative coordinates and objects outside terrain coverage. Use independent expected coordinates and per-fixture tolerances. | Partial 3D Tiles placement tests; cross-import coverage is new. |
| 13 | P0 if Pro ships | **Pro assets, textures and vegetation.** Use the exact asset package intended for customers. Check realistic buildings, roof variants, facade textures/UVs, lit-window settings and tree/forest placement with fixed seeds. Require valid material slots, resolvable images and meaningful evaluated geometry. Cover missing/outdated assets and missing material/vegetation libraries. Exercise the shipped style options. Import-for-export and 3D Tiles baking have dedicated suites 22 and 23 below. | PML translator tests cover only a small component; pipeline tests are new. |
| 14 | P1 | **Imagery overlays.** Serve colored/numbered image tiles locally. Require correct tile selection, mosaic order/orientation, pixel values, terrain UV placement and material assignment. Test image cache reuse, 404-driven zoom fallback, exhausted fallback, corrupt images, packed images and save-to-file with a fresh output directory. Verify the image remains usable after reopening the blend file. | New. |
| 15 | P1 | **GPX tracks.** Test namespaced/non-namespaced files, multiple track segments, optional elevation, and curve/mesh modes. Require correct point/edge/spline counts, no connections between separate segments, expected geographic/elevation placement, curve bevel/material and documented projection behavior on terrain. Empty/one-point/malformed tracks must not crash subsequent imports. | New. |
| 16 | P1 | **Input validation and repeatable configuration.** Exercise missing/unwritable data directories, invalid file/URL/bounds, missing glTF addon, absent provider credentials, missing setup script or `setup` function, empty default-level settings and stale terrain selection. Require an actionable error, appropriate operator result, no helper leaks, and a successful valid import afterward. Test progressive/concurrent defaults and settings captured at import start. | Some 3D Tiles preconditions are covered; other paths need tests. |
| 17 | P0 | **Save/reopen and repeated use.** Import representative datasets, save under artifacts and reopen in a fresh Blender process. Require unchanged geometry, materials, packed textures, attribution, GN attributes and geographic origin. Repeat success, failure, cancel, cleanup and another import; require no stale global state, accumulating temporary collections/handlers or orphaned helpers. Also test a previously released scene with saved Blosm settings. | Existing lifecycle checks are partial; persistence and mixed-operation sequences are new. |
| 18 | P1 | **Background and scripting contract.** Invoke the public import operator in factory-startup background Blender. Require blocking completion for command-line mode and explicit rejection of interactive progressive mode in background execution. Run OSM, terrain, GPX and local 3D Tiles examples; immediately inspect final results. Ensure invalid input produces a detectable script/test failure and that sequential jobs can share one Blender process. | Partial operator tests; broaden to other supported import types. |
| 19 | P2 | **Visual regression.** Produce small fixed-camera renders for roof shapes, courtyard buildings, roads on terrain, realistic facades/vegetation, imagery and textured 3D Tiles. Fix seeds, lighting, renderer and resolution. Compare against reviewed, Blender-version-specific baselines with a calibrated tolerance; save image differences for failures. Numeric geometry/material checks remain the primary oracle. | New harness; existing imports supply reusable scenes. |
| 20 | P2 | **Performance and resource use.** Reuse the direct-import benchmark and add small/medium/large OSM and tile fixtures. Record wall time, first-tile time, UI timer gaps, peak parent/helper RSS, datablock counts and cleanup after repeated cycles. Measure one/four/eight workers and joining costs serially. Establish baselines before choosing thresholds; fail hangs/leaked helpers immediately and review sustained memory or timing regressions. | Direct-import benchmark and some integration timings exist; resource/soak checks are new. |
| 21 | P1 | **Live provider smoke tests.** Run bounded imports against the existing Helsinki/Aarhus/AGI/San Francisco fixtures and one small Google region using an available test credential. Include OSM/terrain/imagery providers that the release advertises. Verify actual imports, attribution, valid content and cleanup; classify outages/rate limits separately. Keep network runs separate from deterministic tests and cap the region, request count and duration. | Live 3D Tiles scripts exist; Google and other import providers need coverage. |
| 22 | P0 if Pro ships | **Import for export.** Exercise the separate realistic-building export layer/renderers, generated texture compositing, shared facade/cladding UV layout, material reuse and portable export/reimport. Test cold/warm texture caches, different colors/assets/dimensions, missing templates, write/render failures and scene cleanup. Include the Blender 4.5 and 5.x compositor API paths. Detailed cases follow. | No dedicated release tests found; previously covered too broadly by family 13. |
| 23 | P0 if Pro ships | **3D Tiles baking to one texture.** Exercise the public bake operator with diffuse and emission materials, one/multiple tile objects, copy/replacement settings, UVs, image resolution and packed-image persistence. Compare texture appearance, preserve source geometry and unrelated objects, verify Undo/Redo and restore render settings after success/failure. Detailed cases follow. | No dedicated release tests found. |

## Mandatory single/separate object matrix

Single/separate output is a P0 test dimension across the relevant families,
including import-for-export and baking. It is not satisfied by one generic
"joining works" test. Run each paired fixture from the same clean starting scene,
with the same inputs, settings and random seeds, changing only the grouping option.

OSM uses `scene.blosm.singleObject` ("Import as a single object"). The base
renderers aggregate objects **per layer** when enabled; an import containing
buildings, water and different road classes need not become one scene object.
When disabled, ordinary building/feature objects have their own origin,
layer collection/parent and, where the renderer assigns them, OSM tags.
Realistic area layers can explicitly force aggregation even when the global
option is disabled. Expected counts must therefore be defined per fixture/layer.

3D Tiles uses the independent `scene.blosm.join3dTilesObjects` option. A successful
mesh-only import should produce one joined mesh when enabled and retain its
imported mesh-object grouping when disabled. One tile can contain several mesh
objects, so object count is not synonymous with imported tile count. Existing
hierarchy flattening still applies. Joining currently does not weld duplicate
vertices; that operation is commented out in `BlenderRenderer.joinObjects`.

| Workflow | Single/join enabled | Single/join disabled | Additional dimensions |
| --- | --- | --- | --- |
| OSM 2D footprints and other supported layers | Required | Required | Buildings, water, vegetation, road/rail curves; terrain absent/present. |
| OSM 3D simple buildings | Required | Required | All roof fixtures, parts/courtyards, terrain absent/present. |
| OSM 3D realistic with ordinary materials | Required | Required | Both the standard and experimental asset/style paths if shipped; terrain absent/present. |
| OSM Import for export | Required | Required | Texture generation and reuse, multiple building colors/sizes, terrain absent/present, export/reimport. |
| OSM 2D Geometry Nodes | Required | UI-supported configuration must be respected | The GN selector is shown only with Single object enabled. Switch to disabled after setting it and test retained hidden settings; require a defined safe result, not corrupt attributes or missing buildings. Do not claim separate-object GN support solely from a stale stored value. |
| 3D Tiles | Required | Required | Blocking/progressive x one/four workers; add the concurrency-disabled baseline and targeted eight-worker cases. GLB, B3DM and mixed multi-object tile fixtures. |
| 3D Tiles material replacement and baking | Required as input | Required as input | Bake one joined object versus all separate imported meshes; copy/replacement options, save/reopen and export. |

For every applicable pair, assert:

1. **Expected grouping and ownership.** Correct mesh/curve counts, layer
   collections, parents and separate-object names/tags; no accidental joining of
   unrelated scene objects or a previous import. Include empty and one-feature
   imports as well as multiple features, shared names and multiple materials.
2. **Equivalent world-space results.** Compare the union of evaluated surfaces
   or curve paths, bounds, geographic placement, height/terrain projection and
   face orientation. Separate origins/parents and joined transforms may differ.
   Do not compare raw local coordinates, object iteration order or just vertex
   counts. Use full small-fixture comparisons with explicit float tolerances.
3. **Equivalent appearance and attributes.** Verify semantic face-to-material
   assignment, UVs, colors, textures and relevant named attributes. Joining may
   remap material-slot indices and names. For roofs and realistic components,
   ensure auxiliary objects join the intended building/layer and remain textured.
   Per-object OSM tags in separate mode need their own assertions; do not demand
   identical per-object metadata from an aggregated layer.
4. **Independent editing and state.** Moving/editing one separate building or
   tile must not alter another's geometry unexpectedly. Toggle the option between
   successive imports in the same session; each import must honor its own setting
   without merging with prior results or retaining stale renderer/attribute state.
5. **Lifecycle and persistence in both states.** Save/reopen and GUI Undo/Redo
   preserve the expected grouping and appearance. In progressive 3D Tiles, objects
   remain separate while loading; successful completion joins once when requested,
   and cancellation retains completed objects without joining. Exercise partial
   failures, destination edits and errors around finalization in both settings.
6. **Export and bake equivalence.** For OSM export mode, both groupings must
   yield portable materials/textures after export/reimport. For tile baking,
   require equivalent scene appearance and one output material/image in both
   cases, allowing different atlas layouts. Verify the documented source-retention
   difference: replacement affects one source, while a multiple-source bake keeps
   originals and creates a new result.

Existing tile integration/UI scripts already run important joined/unjoined
cases. Extend them with explicit object-grouping, appearance and multi-object-tile
assertions where needed. The OSM pairs require new automated fixtures. This matrix
adds coverage requirements; it does not report new test executions.

## Dedicated export and baking cases

These suites were made explicit following the owner's review of the initial
proposal. They add proposed tests, not new implementation or completed test runs.
They cover separate paths: building import prepares export materials through
compositing, while the 3D Tiles tool performs a selected-to-active Cycles bake.

For **Import for export**, use a small deterministic OSM/style/asset fixture
containing colored cladding, window overlays, doors of different dimensions,
curtain walls and representative roofs. Run the public import operator with
`importForExport` enabled and disabled, with the same seeds and asset package.

| Case | Proposed assertion |
| --- | --- |
| E1 — Routing and geometry | The shipped UI/settings select the export path when enabled. The export layer uses its `facade` UV map for both facade and cladding; materials must not depend on the normal path's separate `cladding` UV map or `cladding_color` attribute. Compare world-space surfaces/bounds with ordinary import, allowing intentional topology differences. Exercise single/separate objects and ensure a later ordinary import still works. |
| E2 — Texture compositing | Check pure cladding, facade plus cladding/color, direct-use diffuse textures, doors, curtain walls and roof cladding. Require readable nonempty images with the intended format/dimensions, expected fixture colors/features, correct texture scale and door placement, and assigned export materials. Compare reviewed appearance with ordinary mode within a documented tolerance; specify intended losses for shader effects not carried by the export material. |
| E3 — Cache and material identity | Run with an empty texture directory and then with existing baked textures/materials. Identical inputs should reuse the correct result. Different colors, asset packages containing same-named textures, door sizes and relevant texture scales must not silently reuse an incompatible result. Include a corrupt cached texture and an explicit expected recovery/error policy. |
| E4 — Compositor compatibility | Run actual texture generation with the distributed `building_material_templates.blend` on Blender 4.5, 5.0 and 5.2. Verify both legacy `scene.node_tree`/file slots and newer `compositing_node_group`/file-output items through the production code. Check PNG/JPEG output, temporary filenames and output relocation without extra stale files. |
| E5 — Portable result | Export a small building fixture to GLB and at least one additional format the release promises, such as FBX. Reimport in a fresh Blender process with only the exported package and no access to the original asset/generated-texture directories. Require correct world geometry, UV-driven appearance, material assignment and resolvable/embedded images. Check ordinary blend save/reopen separately. |
| E6 — Failure and cleanup | Missing template scene/nodes, missing source image, unwritable output directory and injected render failure must produce an actionable error. User scenes/settings/assets must survive; only operation-owned temporary scenes/files may be removed. Include a user scene sharing a template name, repeated operations, and a valid import following failure. |

Source anchors: `setup/realistic_dev.py`, `building2/layer.py`,
`building2/texture_exporter.py`, `item_renderer/texture/export/`.
The current compositor code writes generated textures beneath the configured
data directory and loads scenes from the external asset package. Tests must
therefore point both assets and data at isolated test locations, with generated
outputs beneath artifacts. It is insufficient to check only that mesh import
finishes or that a material has an image node.

For **3D Tiles baking**, generate tiny meshes with known colored/checker textures
and use actual imported GLB/B3DM fixtures for integration. Use 128/256-pixel CPU
bakes for fast correctness checks, plus a representative larger bake in the
performance tier. Exercise `bpy.ops.blosm.bake_to_single_material()` with scene
bake settings explicitly set; the implementation prefers those settings over
operator property values when they are present.

| Case | Proposed assertion |
| --- | --- |
| B1 — Channel detection | Pure diffuse materials choose DIFFUSE and pure emission materials choose EMIT. Empty, mixed and unsupported shader configurations have a defined error policy. Include Principled emission, node groups, disconnected nodes and multiple outputs to establish which real imported materials are supported. Test more than three materials with one incompatible outlier and changed ordering/names so sampling cannot silently produce an incorrect bake. |
| B2 — Baked appearance | One tile with several materials and several selected tile objects produce one target material and one image at the configured resolution. Assert expected interior pixel colors/features and compare unlit rendered appearance before/after baking, including seams. Diffuse color baking must exclude scene lighting; emission output must preserve source color. Test margins, close surfaces and cage extrusion to detect missed faces and projection bleed. |
| B3 — UVs and transforms | Require finite, usable target UVs with the intended active/render layer, adequate coverage and no unintended island overlap. Test absent/multiple source UV maps, rotated/scaled/parented tiles and large geographic offsets. World-space geometry must be preserved, and original meshes/UVs must remain unchanged when producing a copy. Test the UV-margin setting explicitly rather than assuming it is effective. |
| B4 — Copy/replacement behavior | With replacement off, all sources remain unchanged and the new baked result becomes active. With replacement on and one source, preserve object identity/transform and change only its intended mesh/material data; shared-mesh users must remain intact. For multiple sources, verify the current behavior of keeping originals and making a joined result even when replacement is requested. Non-mesh and unselected objects must remain untouched. |
| B5 — Settings and lifecycle | Verify scene settings and operator fallback, resolution, bake margin and cage extrusion. Restore render engine, bake type and all modified bake flags/margins on success and injected failure. In a real GUI process, one Undo restores the prior scene and Redo restores the baked result. Check registration/disable/re-enable and valid selection/mode handling. |
| B6 — Failure cleanup | Inject failure during duplication/join, UV preparation, bake and image packing; also test missing source textures and unavailable bake capability. No success report may mask an incomplete result. Preserve source geometry and user settings, and account for temporary objects, meshes, images and materials, not just the target object. A valid bake must succeed afterward. |
| B7 — Persistence and export | Require the final image to be packed, save/reopen in a fresh process, and verify pixels/material links after the original tile textures are unavailable. Export/reimport the baked object using the formats the release promises and compare appearance/placement. Repeat baking to check for accidental data accumulation. |
| B8 — Full user workflow | Import mixed GLB/B3DM fixtures, optionally join or use the existing "Replace materials with export-ready" tool, then bake and save/export. Verify that texture assignments survive material replacement and that baking works after both blocking and completed progressive imports. During an active progressive job, require either safe handling or explicit prevention; there must be no lost tiles or corrupted undo history. |

Source anchors: `util/blender_extra/baking.py`, the 3D Tiles Tools panel and
material-replacement operator in `gui/__init__.py`.
Two code observations motivate specific tests: `TilesAutoSelector` samples at
most three materials, and `_prepare_uv` passes zero margin to Smart UV Project
while its margin-aware packing call is commented out. These are source-level
observations; visual correctness, cleanup and version compatibility have not
been validated by running this new suite.

## Proposed hosting and storage

Keep test source in the same Blosm Git repository as the addon so a commit or
release tag identifies both implementation and matching tests. The existing
`tests/` directory is a suitable starting point. Store the common runner,
standalone tests, Blender background tests, GUI tests and fixture generators
there, with the plan and runner instructions in `docs/`. Keep CI workflow
definitions in `.github/workflows/`. The inspected local release commit has
test scripts but no checked-in CI configuration. This is a hosting proposal;
no runner, workflow, remote repository or external storage was configured.

| Content or workload | Recommended location |
| --- | --- |
| Test code, fixture generators, expected-value definitions and fixture hashes | Blosm repository, alongside the code they qualify. Use the owner's existing OSM repository and artifact fixtures at the locations below. Exclude development tests and CI files from the customer distribution ZIP through the actual packaging route. |
| Fast standalone/local-HTTP tests and small background Blender smoke tests | GitHub Actions with GitHub-hosted runners; execute on relevant pushes and pull requests. Pin Blender versions/builds rather than depend on whatever is preinstalled. |
| Full Blender version matrix, real UI/Undo tests, compositing, baking and visual/performance tests | A dedicated Windows test machine, initially invoked locally and optionally connected as a self-hosted GitHub Actions runner. Use a dedicated desktop session for UI tests and isolated Blender profiles. Run heavy jobs serially for reproducible timings. Restrict this runner to trusted commits. |
| Pro asset packages and larger licensed or downloaded fixtures | Provisioned storage available to the test machine, with recorded versions/hashes and controlled access. Avoid embedding customer asset packages or credentials in the test repository or uploaded reports. |
| Reusable GPX and Geometry Nodes `.blend` inputs | `../artifacts/blosm_release_tests/fixtures/`, separate from disposable run outputs. Use the supplied real GPX and Buildify file listed below, and record their hashes. |
| OSM input files | Use `../osm_extracts` directly; keep the OSM dataset in its dedicated repository. |
| Per-run generated fixtures, imported scenes, renders, logs and JSON reports | `../artifacts/blosm_release_tests/<commit>/<run>/`. Keep these outputs out of Git. |
| CI results for review | GitHub Actions summaries and uploaded workflow artifacts, especially failure logs and image differences. Choose a retention period for routine runs and retain final release qualification evidence in durable project storage. Do not use expiring CI artifacts as the sole fixture source or release archive. |

GitHub supports both hosted execution environments and self-hosted machines
with custom software/hardware, as documented in its
[hosted runner](https://docs.github.com/en/actions/concepts/runners/github-hosted-runners)
and [self-hosted runner](https://docs.github.com/en/actions/concepts/runners/self-hosted-runners)
guides. Its [workflow artifacts](https://docs.github.com/en/actions/tutorials/store-and-share-data)
support test-output uploads and configurable retention. The workload split above
is a recommendation for this addon, not a claim that hosted runners cannot run
Blender or perform CPU baking.

Start with a common runner working locally on the available Blender installations.
Then connect the fast deterministic suite to GitHub Actions, followed by the
dedicated Blender runner if unattended execution is wanted. Run full release
qualification on demand against a named candidate commit/ZIP. Live-provider
checks should have a separate explicit trigger and bounded request budget.

## Locations for reusable test inputs

The owner selected these actual input locations on 2026-11-09:

| Input | Location and handling |
| --- | --- |
| OSM files | `../osm_extracts` — use the existing `prochitecture/osm_extracts` repository directly, without copying it into Blosm. The supplied path used singular `osm_extract`; the verified existing checkout uses plural `osm_extracts`. |
| Real GPX | `../artifacts/blosm_release_tests/fixtures/gpx/2026-04-16_13-16_Thu.gpx` — copied from the owner-supplied file; original retained. |
| GN 2D setup | `../artifacts/blosm_release_tests/fixtures/gn2d/buildify_1.0.blend` — moved from the owner-supplied file into the fixture library. |

Verified SHA-256 values after transfer:

- GPX: `4c933280057fecf74ab26a077cecdf1d57e4a98e89c56b36570fa3cc6b2c903d`.
- Buildify: `d3438b0a02ca0ab8f7b790cd10fc87dbff68c668c7142044d79fade75ab5147a`.

The synthetic GPX and its synthetic-only validation report have been deleted.
Their earlier expected point counts and coordinates do not apply to the real file.
These transfers establish fixture locations, not successful Blender validation.
Inspect the Buildify node-group interface and saved Blender version when adding
the GN import test. Keep this supplied reference intact and use working copies
for tests that modify scene data.

Record the OSM repository revision and hashes of the selected extracts in each
test run. The checkout currently also has untracked `dev/` and `multipolygon/`
directories, so a commit ID alone does not identify every available input.
Do not modify or clean the owner's OSM checkout as part of test setup.

The selected OSM inputs are now specified in `tests/fixtures/osm_catalog.toml`,
with relative paths, SHA-256 values, companion files, tiers and import policy.
Use this catalog rather than discovering arbitrary OSM files at runtime.
Use fixture catalogs to record each
fixture's ID, version, relative file path, SHA-256, origin/license where applicable,
required Blender/node-group interface, and expected counts or other assertions.
Resolve the OSM repository and artifact fixture library through separate
configurable root arguments rather than hardcoding this workstation's paths in
test code.

Keep the reusable fixture library available between runs. Generated run scenes,
image differences and reports belong in per-run output directories. A future CI
machine can reproduce generated inputs or obtain versioned external fixtures from
durable controlled storage, verify the catalog hashes, and run the same tests.
If a fixture is required but unavailable, record an explicit incomplete/failed
qualification rather than silently dropping the case. Keep all development
fixtures out of the customer addon ZIP unless a file is independently a required
distributable runtime asset.

## Selected OSM files

The initial inventory parsed the XML inputs from `osm_extracts` revision
`0fec21466717dde28297c5a35b9ca07ae8920520`. The catalog selects 12 required core
fixtures and two extended fixtures, totaling 17 files including companions.
All catalog hashes and paths were verified locally. This establishes the input
set; Blender imports and rendered-output baselines have not yet been validated.
The inventory report is under
`../artifacts/blosm_release_tests/fixture_inventory/osm_inventory.json`.

All paths below are relative to `../osm_extracts`.

| Tier | Exact main input | Assigned coverage |
| --- | --- | --- |
| Required core | `facade_visibility/berlin_karl_marx_allee.osm` | Main integration fixture explicitly required by the owner. 2D, 2D with Buildify where supported, 3D simple, 3D realistic and Import for export. Test Single object enabled/disabled, plus terrain absent/present and save/reopen. |
| Required core | `fit_rectangles/rectangle_01.osm` | Minimal building baseline and fast GN/material smoke tests. |
| Required core | `fit_rectangles/quadrangle_01.osm` | Four-sided footprint geometry and roof placement. |
| Required core | `feature_detection/complex_01.osm` | Complex single-building footprint. |
| Required core | `feature_detection/curved_01.osm` | Curved footprint with many outline vertices. |
| Required core | `feature_detection/shared_edges_03.osm` | Three building outlines; adjacent geometry, grouping and independent editing. |
| Required core | `fit_rectangles/washington_pentagon.osm` | Building multipolygon with 1 outer and 41 inner way members; preserve holes and use mode-specific roof expectations. Member counts are not assumed to equal completed ring counts. |
| Required core | `dev/building_part_detection.osm` | One tagged building outline and nine tagged parts; grouping, heights and no accidental duplicate building shell. |
| Required core | `multipolygon/madrid_palacio_de_nuenavista_hipped_roof.osm` | Hipped-tagged multipolygon with inner members; courtyard handling and renderer-specific roof behavior. |
| Required core | `multipolygon/magdeburg_university_building_14_hipped_roof.osm` | Multipolygon, four tagged building parts and hipped/gabled roof tags. |
| Required core | `streets/berlin_frankfurter_allee.osm` | Release road/rail curve renderer, feature toggles, grouping and terrain projection. |
| Required core | `railways/kuzmice.osm` | Railway-focused import, category toggles and terrain projection. |
| Extended | `facade_visibility/bratislava_old_town.osm` | Larger urban regression; mansard, pyramidal, hipped and other roof tags, materials and multipolygons. |
| Extended | `facade_visibility/moscow_leningradski_prospekt.osm` | Broader roof mixture including saltbox, gambrel, round, skillion and dome tags; building relations and area layers. |

Karl Marx Allee must load its adjacent companion
`facade_visibility/berlin_karl_marx_allee_extra.osm`. Inspection found 150 elements
tagged `building`, 52 tagged `building:part`, 101 with `building:levels`, and five
multipolygon relations. These are XML tag counts, not expected Blender object
counts. The companion resolves the one missing way referenced by its multipolygon
or building relations. This makes it especially useful for both the complete
offline import and a controlled missing-member test.

The two extended fixtures also use their exact companions:

- `facade_visibility/bratislava_old_town_extra.osm`.
- `facade_visibility/moscow_leningradski_prospekt_extra.osm`.

Those extended snapshots remain incomplete for some multipolygon/building
relations even after their companions are read (10 and 154 absent way IDs,
respectively, across those relation types). Tests must assert documented behavior
for the supplied snapshot and must not treat it as complete geographic coverage.

The catalog enables `loadMissingMembers` for the three paired inputs and disables
it for the other fixtures. The test runner must reject unexpected network access.
For paired inputs, fail fixture preflight if the companion is missing; Blosm's
current loader would otherwise try to download it. A deliberate missing-member
test should use a working copy under artifacts and controlled local responses,
preserving the original repository inputs. Ordinary file-import tests should
consume the original files directly without modifying the extract checkout.

The three selected inputs under `dev/` and `multipolygon/` are currently untracked
in the OSM repository. Their hashes are pinned in the catalog, but another runner
must be supplied with those exact files before it can pass the full core suite.
Do not silently skip them. No changes were made to the OSM repository.

The existing roof samples do not establish exhaustive coverage of every supported
roof type or every malformed-input case. Map each roof tag to the appropriate
renderer/fallback and retain the additional edge-case requirements elsewhere in
this plan. In particular, a hipped tag on a courtyard does not imply that the
simple renderer must produce a hipped roof; it currently uses a flat courtyard
fallback. Fix random seeds and establish independent expected geometry before
turning snapshots into regression baselines.

## Execution design

An executable implementation is now available in `tests/run_release_tests.py`.
See `tests/README.md` for the implemented inventory, fixture configuration,
commands and coverage limits. The broader cases below remain the qualification
roadmap; the runner does not claim that a selected passing subset certifies them.

Use one release runner with explicit test inventories. Merely invoking unittest
discovery will run the 23 standalone tests but miss the `check_*.py` Blender
scripts and the PML functions. The runner should select Blender executables,
the candidate ZIP/commit, test families and one artifacts directory. It should
record pass/fail/skip, missing fixtures, duration and the exact environment for
every case. A skipped Pro/live/real-fixture test is not a pass.

Run the full deterministic suite on Blender 4.5 and 5.2, and a reduced package,
import, Geometry Nodes and lifecycle suite on 5.0 and 5.1. Include actual
import-for-export texture generation on 5.0 to cover the compositor transition,
and small diffuse/emission bakes across all four versions. Those four version
families are installed locally; this is a proposed test matrix rather than a
new claim of support. The source currently declares Blender 4.5 as its minimum.
Use separate GUI Blender processes for timers, context, Undo/Redo and modal
overlay behavior. Headless stepping cannot establish those UI guarantees.
Apply explicit process timeouts and inspect GUI JSON status: a zero GUI exit
code alone does not establish that timer assertions passed.

For every platform promised in the release, run at least installation, direct
import, helper startup/IPC/cancellation, Unicode paths and disable/re-enable
checks. A full Windows run does not establish macOS or Linux compatibility.
Keep user profiles and preferences isolated under the artifacts directory.

Use fixture generators in `tests/` for compact OSM, HGT, image and tile
datasets. Write generated payloads, scenes, ZIPs, renders, baselines, reports and
logs under `../artifacts/blosm_release_tests/<commit>/<run>/`.
Keep test source and this documentation in the repository. Do not stage generated
scenes, binaries or JSON reports. Existing downloaded real fixtures may be reused
from their artifacts directories with explicit hashes and provenance; a required
fixture directory must contain the expected fixtures, not silently skip an empty
glob. Production 3D Tiles tests must retain their no-dataset-cache guards.

Use explicit expected geometry/control points for synthetic tests, alongside
parity with the historical renderer. Agreement between two import modes can
preserve the same underlying error. Compare full geometry on small fixtures;
use documented scale-aware tolerances for real geographic coordinates and
floating-point joining differences. Require appropriate mesh invariants, not
universal manifoldness: open footprints, terrain and curves are intentional.

## Recommended implementation order

1. Wrap and rerun the existing standalone, direct-import, concurrent, progressive
   and GUI suites on the selected release candidate. Add explicit artifact status,
   timeouts, fixture checks and Blender-version recording.
2. Add distribution installation/registration and save/reopen checks. The actual
   candidate ZIP needs a smoke path that works without a Git checkout: the current
   historical renderer comparator reads commit `4adff002...` using `git show`.
3. Use the selected inputs in `tests/fixtures/osm_catalog.toml`, starting with
   the required Karl Marx Allee sample and its companion;
   implement parsing, simple buildings, roof cases, 2D layers and Geometry Nodes
   checks using the supplied Buildify reference.
4. Add synthetic terrain and cross-import geographic/elevation checks, then Pro
   asset checks if the Pro edition is included. Include the dedicated import-for-
   export and 3D Tiles baking suites in that Pro release gate.
5. Add imagery, GPX, invalid-input and background execution cases. Run the bounded
   live-provider checks close to publication.
6. Establish reviewed visual baselines and repeatable performance/resource tests.

The extension distribution now uses `blender_manifest.toml` for version and
minimum Blender version. Package checks verify the installed metadata and runtime
identity. Pro capability is detected from the presence of the `realistic`
directory, so tests exercise each actual packaged edition rather than infer
edition behavior from this full developer checkout.

The proposed release gate is all applicable P0 checks passing on the exact
candidate, with reviewed P1 results, no unexplained crashes/hangs or leaked
helpers, and clearly recorded external-service/fixture limitations. Any change
after qualification requires the relevant tests to run again against the new
candidate. This document proposes the work; it does not certify the release.
