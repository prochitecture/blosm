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

## Remaining limits

- The live test validates current fields and 32-image requests with this client
  token and area, not every account/provider state or an API maximum. The main
  current API documentation was login-gated during initial implementation.
- Camera frames/poses and metadata are imported; photos, image viewers, accurate
  fisheye projection and radial-distortion rendering are subsequent work.
- EXIF mirrors are approximate. Missing altitude uses the explicitly marked
  scene-plane fallback. Cross-provider vertical-datum alignment is not solved.
- Main-thread work is cooperative and size-bounded, not hard-real-time. An
  individual Blender call, small JSON parse or undo checkpoint can exceed the
  nominal eight-millisecond job budget.
- Decoder support intentionally targets documented single-point image features;
  malformed/oversized resources produce aggregate errors/incomplete coverage.
  Dateline-crossing or polar extents and datasets beyond stated bounds are not
  silently accepted.
- Validation used Blender 5.2 on Windows. Other supported Blender versions and
  operating systems still need release qualification.
