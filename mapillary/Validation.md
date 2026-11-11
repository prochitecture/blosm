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
| `tests/test_mapillary.py` | 9 standalone tests passed |
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

## Remaining limits

- No Mapillary token was available, so a real-area import and current live API
  responses have not been validated. Public SDK/viewer source confirms the
  requested fields and multi-image query; the main current API page required
  login. Do not claim a verified API maximum for the internal 32-ID batch size.
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
