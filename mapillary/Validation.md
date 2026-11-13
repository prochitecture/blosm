# Mapillary import validation

Date: 2026-11-10. Development branch: `mapillary_import_release`, based on local
`release` at `72286b3e77e8cd3c6021384e856291d7315a6a07`.
Runtime: Blender 5.2.0 LTS on Windows; standalone unit tests use Python 3.13.
All Blender processes used factory startup, disabled autoexec and the chat's
isolated resource directories under
`C:/dev/workspace/artifacts/blender_profiles/01a0edf2-527a-7e10-be5c-d4e614d29008/`.
Tests enabled online access for their own processes only. No normal profile,
installed addon, bundled Python installation or existing Blender session changed.

## Correctness and regression results

| Check | Result |
| --- | --- |
| `tests/test_mapillary.py` | 9 initial tests passed; 10 passed after the live dense-tile regression was added |
| `tests/check_mapillary_import.py` | Background Blender checks passed, including a 5,000-camera run |
| `tests/check_mapillary_ui.py` | 9 real UI scenarios passed using the production HTTP helper |
| `tests/test_progressive_transport.py` | 7 existing scheduler/transport regression tests passed |
| `tests/test_download_resilience.py` | 17 existing resilience/redaction tests passed |
| `tests/check_progressive_import.py` | 16 existing 3D Tiles integration/lifecycle checks passed |

Standalone tests cover cardinal compass directions, inverse angle-axis rotations
including zero and nonzero pitch/roll, portrait EXIF 3/6/8, missing EXIF/altitude,
zero altitude, normalized focal conversion, spherical and fisheye representation,
computed position preference, MVT extent and dictionary order, exact IDs greater
than 2^53, malformed input and cooperative consumer pause/cancel cleanup.

Background tests check camera object/data types, provenance, source IDs and
geographic world positions against Blosm's standard projection; selection,
active object/collection/cursor/render camera preservation; blocking/stepped
parity with reversed responses; object rename/delete; destination
rename/unlink/delete; partial cancellation; missing metadata fallback;
cross-tile buffered-ID deduplication and explicit deduplication limits.

Real UI tests exercise the actual import operator, timer controller and shared
transport against a separate localhost fixture. They verify ordinary blocking
Undo/Redo, grouped progressive Undo/Redo without restarting downloads,
intervening user edits with independent history, Cancel with retained cameras
and one aggregate result, destination deletion errors, scene/Edit Mode
pause/resume, HTTP metadata failure with explicit approximate fallback, and
add-on shutdown. They include the 5,000-camera measurement below and successfully
Undo/Redo that completed import. The test launcher checks the report's status:
a GUI process exit code alone is not interpreted as a passing test.

## Scale and responsiveness

Synthetic inputs are generated in memory, with deliberately delayed/out-of-order
localhost HTTP responses and IDs above 2^53. These are individual runs, not
statistical benchmarks or real-provider performance claims.

Final instrumented UI run with 5,000 entries, four HTTP slots and Global Undo on:

| Measurement | Observed |
| --- | ---: |
| Cameras successfully created | 5,000 |
| Import wall time, including helper startup and checkpoints | 59.70 s |
| Cumulative MVT decoder execution | 0.162 s |
| Cumulative metadata parse/validation execution | 0.041 s |
| Cumulative camera-commit execution | 2.916 s |
| Maximum individual decoder segment | 0.416 ms |
| Maximum metadata parse | 0.687 ms |
| Maximum individual camera commit | 2.792 ms |
| Maximum `Job.step()` turn | 10.53 ms |
| Cumulative undo-checkpoint execution | 40.26 s |
| Maximum undo checkpoint | 163.29 ms |
| Maximum heartbeat gap while the import remained active | 194.63 ms |
| Peak observed HTTP concurrency | 4 |
| Retained raw tile payload in this fixture | 348,570 bytes |
| Maximum decoded metadata batch | 32 entries |
| Maximum deduplication ID set | 5,000 IDs |
| Blender RSS before scale import | 587.85 MiB |
| Blender RSS with cameras and undo history | 938.79 MiB |
| Blender RSS after Undo, with Redo history retained | 939.81 MiB |

Undo checkpoints dominate this case even though data processing is cooperative.
The scene remained interactive between turns, but the approximately 195 ms
worst heartbeat gap is not a smooth-frame-rate guarantee. The reported process
RSS includes the existing test scene, earlier checks, allocator retention and
undo state; it is neither a standalone camera allocation measurement nor a
proof of a stable long-running RSS ceiling. Undo did not return RSS to baseline.

The whole UI scenario also recorded a roughly 26-second maximum heartbeat gap
outside the active import, including bulk test-scene cleanup. That value must
not be confused with the separately instrumented active-import heartbeat gap.
Mass deletion, history/dependency-graph work and other indivisible Blender
operations cannot be preempted by this adapter's scheduler budget.

Earlier background runs created 5,000 cameras in roughly 3-8 seconds without GUI
undo, depending on retained scene data and test instrumentation. Their timing
is not directly comparable with the real UI/HTTP run. In particular, the in-memory
test transport can spend time synthesizing responses synchronously. The real UI
run is the relevant measurement for interactive imports with Global Undo.

Practical guidance: start with a small area and a few thousand cameras per import.
Five thousand cameras were validated, but memory and checkpoint latency already
matter at that size. The 100,000-ID bookkeeping ceiling is not a validated usable
camera-count target. Larger areas require additional scale qualification and
possibly a separately designed undo/checkpoint strategy. Network/metadata work
is bounded independently of that remaining Blender-scene and undo growth.

## Evidence locations

Generated reports/logs are under `C:/dev/workspace/artifacts/mapillary_import/`:

- `background/background.json`: background correctness and scale results.
- `background/blender.log`: background check output.
- `ui/ui.json`: final UI scenarios, timings, RSS and aggregate result messages.
- `ui/blender.log`: final UI process output.
- `ui/http-2.json`: final scale fixture request count and peak concurrency.
- `3d_tiles_regression/checks.json`: existing 3D Tiles regression.

Synthetic fixture generators and test code stay in this repository under
`tests/`; no generated fixture payloads or scenes are tracked. Blender's test
session recovery files, if created, also remain in the artifact directories.

## Live-data validation, 2026-11-11

A development client token was read only from the user-created local file
`C:/dev/workspace/blosm/mapillary/mapillary_token.txt`. It was not copied into
this checkout, preferences, command lines, reports or scenes. The shared local
Git info/exclude protects the filename across local branches/worktrees; the
adapter's `.gitignore` also excludes it. Never stage or commit this credential.

The test exercised the real Mapillary coverage and Graph endpoints in an isolated
Blender 5.2 process. Final bounds were west 24.9395, south 60.1695, east 24.9405,
north 60.1705 (central Helsinki). The coverage tile contained 7,073,835 decoded
bytes. Discovery selected 550 candidate IDs; after the final metadata-position
bounds check, both runs created 546 cameras, all with computed rotation and
altitude and with no resource/decode errors.

| Final real UI run | Blocking | Progressive (real timer) |
| --- | ---: | ---: |
| Cameras with Graph metadata | 546 | 546 |
| Import time | 8.49 s | 28.59 s |
| Perspective models | 512 | 512 |
| Equirectangular / spherical models | 11 / 15 | 11 / 15 |
| Fisheye models (optical proxy) | 8 | 8 |
| Maximum metadata batch | 32 | 32 |

The IDs and poses matched between modes. The maximum camera heading difference
against live `computed_compass_angle`, after accounting for projected grid axes,
was approximately 0.0000101 degrees (546 comparisons). This is an API-convention
consistency check, not a survey/ground-truth accuracy claim or optical-image
alignment verification. Photographs were not downloaded. Fisheye distortion and
panoramic renderer limitations remain as documented in README.md.

The real-timer test recorded 839 viewport heartbeats, a maximum active heartbeat
gap of 76.54 ms and a maximum job turn of 8.43 ms. A separate background live
comparison also passed with 546 cameras in both blocking and stepped modes
(12.10 / 10.99 seconds). Smaller-area real UI tests passed with six cameras.
These are individual live runs with variable network/server latency, not a
statistical speed comparison.

Live data exposed limits absent from the original small synthetic fixtures:

- The 4 MiB response cap rejected the real coverage tile. It is now a bounded
  32 MiB per resource; the resulting payload budget is documented in README.md.
- Dense dictionaries exceeded the initial limits. The bound is now 1,048,576
  entries, and irrelevant layers are skipped before dictionary allocation.
- Sparse ROI filtering in a dense tile now checkpoints after at most 16 records
  or roughly 2 ms, as well as every decoder checkpoint, avoiding excessive timer
  turns for rejected points. Individual work still cannot be preempted.
- A new in-memory regression fixture exceeds the old 4 MiB limit and exercises
  dense dictionaries and skipped oversized irrelevant-layer dictionaries.

A broader diagnostic area intentionally stopped at the live test's 2,000-ID
safety ceiling, retaining 1,856 already committed cameras. It had no network or
metadata errors and confirmed the small-area test needed a narrower selection;
this test-only limit is separate from the adapter's 100,000-ID ceiling.

Evidence under `artifacts/mapillary_import`:

- `live/live.json`, `live/blender.log`: passed larger background live comparison.
- `live/live_cameras.blend`: its imported real cameras and metadata.
- `live/dense_area_limit.json`: the bounded broader-area diagnostic.
- `live_ui/live.json`, `live_ui/blender.log`: passed 546-camera real UI comparison.
- `live_ui/live_cameras.blend`: real UI import scene for inspection.
- `live_ui/small_area.json`: earlier six-camera real UI report.
- `live_regression/background.json`: 19 post-fix synthetic Blender checks passed.

All live JSON/log outputs were scanned for the exact token without printing it;
the scan passed. Persisted scene properties were checked before scene saving.

## Session credential fallback, 2026-11-12

The import operator now prefers the current preferences token and otherwise uses
`mapillary/mapillary_token.txt` alongside the adapter. The fallback is read lazily
at most once per Blender session, including missing/empty/unreadable outcomes.
Later preference changes override it immediately. File changes after the first
lookup require restarting Blender. No credential is copied into preferences,
scene data or reports.

All 13 standalone tests passed, including preference priority, lazy read counts,
cached missing/empty files, cached unreadable-file diagnostics and overrides.
`tests/check_mapillary_credentials.py` also passed in factory-startup isolated
Blender 5.2 with the real local token file and HTTP disabled. It exercised the
actual import operator with a capture-only adapter, checking that five imports
across preference/fallback changes accessed the file exactly once. Its sanitized
report/log are in `artifacts/mapillary_import/credentials/`; the token was not
printed, persisted or committed. Existing live HTTP results above remain valid;
this change selects the credential and does not alter the transport or poses.

## Reported upside-down cameras, 2026-11-12

Rechecked the saved real-data scene `live_ui/live_cameras.blend` in an isolated
Blender 5.2 background process using `tests/check_mapillary_upright.py`. All 546
cameras had positive world-Z components for their local +Y (up) axes: 526 had
EXIF orientation 1 and 20 had orientation 3. The smallest up-Z component was
0.877569. None was upside down in this dataset.

The saved Blender world transforms, including the matrix-to-quaternion
conversion, were compared with MapillaryJS `Transform.upVector()` and the
inverse rotation's viewing direction. Maximum vertical-component differences
were 1.20e-7 (up) and 3.33e-7 (forward). Comparing vertical components avoids
using the adapter's grid-convergence calculation as its own reference. All 13
standalone tests also passed, including portrait EXIF 6/8 and zero rotation.

Evidence is in `artifacts/mapillary_import/upright/upright.json` and
`upright/blender.log`. The pose implementation was not changed: this check
does not reproduce the reported problem and does not establish correctness
for an unidentified affected image or scene. It does not test alignment with
terrain elevation. An affected image ID or scene is needed to investigate
that specific import.

## Ground-relative camera placement, 2026-11-12

The user clarified that cameras below Z=0 were the reported problem, and requires
terrain-relative placement when terrain is present, otherwise placement above
Z=0. This supersedes the original absolute-altitude and missing-altitude policy.
The old implementation used provider altitude minus the scene/terrain offset;
live Karl Marx Allee metadata includes near-zero and negative computed altitudes.

New imports use the evaluated terrain surface plus an explicitly estimated 2 m
camera height, or Z=2 m without terrain. Original/computed altitude values remain
metadata; they are not interpreted as above-ground heights. All resulting camera
poses are marked approximate. Source rotations and geographic X/Y are unchanged.
Out-of-coverage terrain rays fall back to the scene plane with explicit provenance
and one aggregate count. Renaming the terrain is supported; deleting/unlinking it
stops work and retains completed cameras. Terrain data/transforms are not mutated.

The 14 standalone tests and 31 background Blender checks passed. The latter
include 5,000 synthetic cameras, negative provider altitudes, transformed/sloping
terrain, evaluated Solidify geometry, terrain below Z=0, missed rays, terrain
renaming/deletion, cancellation, destination deletion and blocking/stepped parity.
The final 5,000-camera no-terrain run took 9.25 s with a maximum job turn of
13.36 ms. A second 5,000-camera import over a 10,000-quad terrain took 17.76 s,
with a maximum job turn of 61.09 ms. Terrain ray misses were zero and geometry
was unchanged. These background fixtures have no GUI undo overhead; elapsed
times vary with other local work and are not controlled speed benchmarks.

All ten real-timer UI scenarios passed, including terrain-relative blocking and
progressive imports with grouped Undo/Redo, terrain preservation, user edits,
cancellation, destination deletion, pause/resume, simulated HTTP failure and
addon disable. The 1,000-camera scale phase took 4.66 s, with Global Undo enabled,
123 heartbeats, maximum active-import heartbeat gap 53.68 ms and maximum job
turn 9.27 ms. Evidence is in `ground_placement_ui/ui.json` and its logs. This is
not a qualification of importing the full dense area with tens of thousands of
cameras and undo history.

A fresh live Karl Marx Allee subarea test imported 750 cameras in each mode,
including 46 with negative provider altitudes. All camera Z values were exactly
2 m, with matching blocking/stepped IDs and poses, no import errors, and no
credential present in persisted metadata. Reports and the corrected scene are
in `artifacts/mapillary_import/live_ground_placement/`; background synthetic
evidence is in `ground_placement/`.

The earlier full-OSM-bounds metadata-only investigation reached 39,483 entries
before its explicit ten-minute diagnostic deadline. It was cancelled and is
not a complete-area count or a camera/undo scale test. Its retained bounded
samples and sanitized report are in `karl_marx_allee_full/`. The user's issue is
negative vertical placement, so no blanket camera rotation correction was made.

## Per-tile progress estimates, 2026-11-13

Reference `release` was synchronized at `e5fe4ac8b0e3a66e56ac9019dd90f99a8b434702`
before this change. Per-tile cooperative counting precedes each tile's existing
streaming import pass. The same bounded coverage/metadata request windows remain;
no whole-area prefetch, full decoded tile list or shared-core change was added.
Status uses "Estimating…" before a usable count, estimated percentages during
processing, and an exact cumulative "cameras created" counter. Finished candidates
include filtered/rejected entries; only normal job completion displays 100%.

Sixteen standalone tests passed. The background suite passed its original 43
checks, including 5,000 cameras without terrain and 5,000 over terrain, plus the
subsequently added valid-empty-tile case in a final 96-camera regression run.
Progress checks cover single-tile totals, cross-tile deduplication, clipping,
partial-tile area weighting, downward revisions, rejected metadata, malformed
tiles, early bounded fetching, and cancellation while counting. ID sets are
released on completion/cancellation. Existing parity and lifecycle checks pass.

The hidden, isolated Blender 5.2 UI suite passed ten scenarios, including grouped
undo/redo, terrain, cancellation, user edits, pause/resume, HTTP failure and addon
disable. Recorded examples were `Mapillary: Estimating… · 0 cameras created`,
`Mapillary: ~11% · 11 cameras created`, and `Mapillary: 100% · 96 cameras created`.
The 1,000-camera timer phase had 189 heartbeats, a maximum active heartbeat gap
of 92.29 ms, and a maximum job turn of 9.51 ms. Counting itself accumulated
0.038 s of decoder work in this synthetic fixture.

A fresh live Karl Marx Allee check passed with 750 candidates processed and 750
cameras created in each mode. Both reported one counted tile, an exact candidate
total of 750 and `Mapillary: 100% · 750 cameras created`. IDs and poses matched,
with no import errors or persisted credential. The dense 12,015,618-byte tile
required about 13 s of additional counting decoder work per run; maximum single
count-decoder calls were 5.76 ms (blocking) and 4.89 ms (stepped). This live test
used background stepping, not the UI timer; UI timing above uses HTTP fixtures.
Counting is real extra work and delays the first camera until its tile is counted.

Evidence is under `artifacts/mapillary_import/progress_messaging/`: `background/`,
`background_final/`, `ui/` and `live/`. Timings are individual runs, not controlled
performance comparisons. Estimated percentages describe candidate work, not
time remaining; spatially uneven density can still give poor early estimates.

## Hidden terrain evaluation, 2026-11-13

The user reported `Object 'Terrain' has no evaluated mesh data` for
`C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend`. A factory-startup
isolated Blender process reproduced the native `Object.ray_cast` error: the
object has a valid 1,298-vertex Mesh datablock, but viewport visibility is
disabled (`hide_viewport`, `hide_get` and `hide_render` all true).

The sampler now evaluates a private temporary copy and builds a per-import BVH
snapshot in the first camera turn. Source visibility, mesh and transforms are
unchanged. Temporary object/collection IDs are removed in a synchronous finally
scope before any cooperative yield or undo checkpoint. The independent snapshot
survives dependency-graph rebuilds and is released on completion/cancellation.
Geometry and world transform are captured once; terrain edits require a new
import to affect subsequent camera placement.

`tests/check_mapillary_hidden_terrain.py` reads the supplied scene without saving
it and uses synthetic metadata over 16 independently known terrain triangle
centroids. Blocking and stepped imports each created 16 cameras at terrain+2 m.
It verified mesh/transform/visibility, selection, active object/collection,
cursor and rendering camera preservation, no temporary collection remaining,
and an unchanged hash of the user's original scene file. No HTTP or real token
is needed for this regression.

The background suite additionally covers hidden terrain with a Solidify modifier,
terrain excluded from the view layer, translated/rotated/scaled terrain and
temporary-ID cleanup. A 5,000-camera run over a 10,000-quad surface completed;
the first evaluated snapshot turn reached about 200 ms, while later rays reused
the BVH. As with other native evaluation calls, that initial synchronous work
cannot be preempted by the scheduler's soft time budget. All ten isolated UI
scenarios also passed, including hidden-terrain placement in both modes and
undo/redo with source visibility preserved.

Evidence is under `artifacts/mapillary_import/terrain_evaluation/`: original
`scene_probe.json`/`probe.log`, `scene_regression.json`/`regression.log`, plus the
background and UI subdirectories. The original scene and normal Blender profile
were never written or reconfigured.

## Counting latency, 2026-11-13

Reference `release` was checked at `e5fe4ac8b0e3a66e56ac9019dd90f99a8b434702`,
with no missing committed changes before this optimization. The earlier count
pass decoded every feature's attributes, including those outside the requested
bounds. Its frequent decoder/record checkpoints also exhausted the scheduler's
32-operation limit before using its available time budget, adding idle timer
intervals to the initial "Estimating…" phase.

The reader now tests geometry before attributes and provides an ID-only count
path. Single-byte protobuf integers have a fast path, retaining uint64 overflow
and truncation checks. Counting combines decoder checkpoints into approximately
2 ms chunks under the existing job budget. No shared-core change, whole-area
prefetch or larger request/decoded-record window was introduced.

A freshly downloaded Berlin tile was reused offline for identical-input
comparisons: 12,015,652 bytes, 165,861 image features and 750 unique entries in
the selected subarea. Standalone Python counting fell from 11.65 s to 2.04 s
(5.7 times faster), with exactly the same selected IDs. Bounded full-record
decoding took 1.98 s. These comparisons exclude download and metadata requests.

`tests/benchmark_mapillary_count.py` drives the production Job through Blender
UI timers with a cached-tile transport. The baseline reader/count method is
loaded read-only from commit `9dd2509`. Time until the first counted tile fell
from 60.98 s to 14.21 s (77% less). Decoder time was 12.88 s versus 3.66 s;
maximum job turns were 10.43 ms versus 13.97 ms, and maximum timer heartbeat
gaps were 36.52 ms versus 39.24 ms. Both returned 750 candidates without errors.
This is a count-phase timer benchmark, not a complete network/import benchmark
or a measurement of the shared Controller's undo handling. Timings are individual
local runs and will vary with coverage density and machine load.
A repeat with the corrected artifact-local temporary directory took 16.71 s,
returned the same 750 candidates and left the shared recovery-file timestamp
unchanged. The observed optimized timer range is therefore 14–17 s here.

Nineteen standalone tests pass, including count/full-record parity across three
extents and three bounding regions, duplicate IDs, exact large IDs, invalid
selected tags, truncated/overflowing integers and out-of-bounds attribute
skipping. All 49 background Blender checks pass, including 5,000-camera runs,
terrain, count cancellation, parity and lifecycle behavior. All ten isolated UI
scenarios pass, including grouped undo/redo and intervening user edits. Its
1,000-camera fixture completed in 4.25 s with Global Undo enabled, a maximum
active-import heartbeat gap of 59.25 ms and a maximum job turn of 8.70 ms.

Evidence is in `artifacts/mapillary_import/count_optimization/`, including the
cached coverage, baseline/optimized profiles and reports, `timer_baseline/`,
`timer_optimized/`, `background/` and `ui/`. The initial timer benchmark omitted
Blender's temporary-directory override and wrote the shared temporary
`quit.blend` recovery file, potentially replacing a previous last-session
recovery file. The benchmark now sets its temporary directory to its artifact
folder. The corrected repeat is in `timer_optimized_isolated/`. Saved source
scenes and existing user Blender sessions were unchanged.

## Remaining limits

- The live test validates current fields and 32-image requests with this client
  token and area, not every account/provider state or an API maximum. The main
  current API documentation was login-gated during initial implementation.
- Camera frames/poses and metadata are imported; photos, image viewers, accurate
  fisheye projection and radial-distortion rendering are subsequent work.
- EXIF mirrors are approximate. Camera heights are estimated at 2 m above terrain
  or the scene plane; actual capture heights and cross-provider vertical-datum
  alignment are not recovered from provider altitude.
- Main-thread work is cooperative and size-bounded, not hard-real-time. An
  individual Blender call, small JSON parse or undo checkpoint can exceed the
  nominal eight-millisecond job budget. The first terrain ray may build Blender's
  evaluated BVH snapshot; very large evaluated terrain meshes can
  cause a longer initial camera turn.
- Decoder support intentionally targets documented single-point image features;
  malformed/oversized resources produce aggregate errors/incomplete coverage.
  Dateline-crossing or polar extents and datasets beyond stated bounds are not
  silently accepted.
- Validation used Blender 5.2 on Windows. Other supported Blender versions and
  operating systems still need release qualification.
