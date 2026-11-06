# Concurrent download validation, 2026-11-06

This validation covers the independent stalled-worker supervisor and retries of
incomplete HTTP transfers. All downloaded manifests and tile payloads remained in
memory. Generated reports and logs are under
`C:/dev/workspace/artifacts/download_live_validation`; existing Helsinki benchmark
files were also reused for controlled fault injection. No Google requests were
made; the live sources were public endpoints from the
[Blosm dataset list](https://github.com/vvoovv/blosm/wiki/Data-sets-with-OGC-3D-Tiles)
and the user's Gaussian-splatting example.

## Live imports

Each supported dataset below passed six Blender 5.2 configurations: blocking and
progressive imports with one and four workers, plus joined blocking/progressive
imports with four workers. Every run fetched its resources from the real provider.

| Dataset | Tiles per run | Vertices | Largest sampled position difference across modes |
| --- | ---: | ---: | ---: |
| Helsinki 2024 | 3 | 34,674 | 0.504 mm |
| Aarhus | 14 | 2,062 | 0.129 mm |
| AGI headquarters | 43 | 384,437 | 0.102 mm |
| San Francisco demo | 11 | 24,354 | 0.150 mm |

The small differences occur when joining changes float32 arithmetic order. Tests
require identical source hashes, vertex/face/loop counts, packed texture bytes,
tile counts and errors. World-space geometry is compared against a baseline tree,
sampling up to approximately 20,000 vertices per subsequent run, with a 2 mm
minimum tolerance scaled to coordinate magnitude. These checks establish parity
between modes; they are not an independent geodetic accuracy assessment.

AGI's 43-tile region also passed blocking and progressive four-worker imports in
Blender 4.5, 5.0 and 5.1. The Svirnas GLB dataset passed its expected-failure tests:
all 18 selected tiles were rejected for `KHR_gaussian_splatting`, with zero partial
mesh objects, in blocking and progressive modes on 5.2.

These are 32 background live-import runs: 684 successfully imported tile instances
and 36 expected unsupported-tile errors, across four Blender versions. Reports:

- `live_5_2_north/live_checks.json`: Helsinki's six passing runs. The initial Aarhus
  run here exposed a test-harness sorting error for unpacked images; the corrected
  Aarhus result is in the next report.
- `live_5_2_aarhus/live_checks.json`: Aarhus's six passing runs.
- `live_5_2_west/live_checks.json`: AGI and San Francisco, six runs each.
- `live_5_2_glb_errors/live_checks.json`: the two expected-failure GLB runs.
- `live_4_5_compat`, `live_5_0_compat`, `live_5_1_compat`: two AGI runs per version.

Live timings are recorded for diagnosis, not a controlled performance claim:
provider/CDN conditions vary and some independent tests ran concurrently.

## Real UI and lifecycle

`live_ui_4_5/live_ui_checks.json` and `live_ui_5_2/live_ui_checks.json` both passed
with eight download workers. Each ran a nine-tile Helsinki region with separate
objects, a joined import, a partial import cancelled after the first tile, and
the 18-tile unsupported GLB dataset. Viewport timer activity and intermediate tile
counts were observed. One Undo/Redo removed/restored a completed import, final
joining shared that undo step, and cancellation retained committed geometry.
The 18 glTF failures produced one aggregate error popup. All helpers exited.

The adapter checks in `adapter_4_5/checks.json` and `adapter_5_2/checks.json` also
verify that a fatal helper error preserves already committed geometry and releases
glTF patches/transport resources, for blocking/progressive and joined/separate
imports. Existing checks cover destination deletion, per-tile failure cleanup,
context preservation and add-on/scene lifecycle.

## Fault injection

The standalone resilience suite covers both ordinary and chunked responses that
end early, exhaustion of the three-attempt limit, and the original deadline
expiring during retry backoff. A real subprocess with a deliberately stuck DNS
resolver must terminate with the specific supervisor exit code despite continuing
helper heartbeats. A completed 4 MiB response held in IPC beyond its short test
deadline and grace period must still succeed when the caller resumes polling.

A live attempt exposed unclassified Windows IPC reset reporting. The regression
test models a reset becoming visible before the process exit status and verifies
that the specific stalled-worker diagnosis is preserved. Other tests cover
cooldowns, credential scopes/redaction, response limits, cancellation and broken
HTTP error bodies.

`check_concurrent_tiles.py` additionally truncates an actual Helsinki payload and
the root manifest on their first local HTTP attempt. Recovered blocking and
progressive imports must have the same geometry as complete-response imports,
and each truncated resource must be requested exactly twice. The local server
runs outside Blender and keeps its payloads in memory.

Final results: all 16 resilience tests passed with Blender 4.5's bundled Python
and standalone Python 3.13; all seven existing transport/scheduler tests passed.
The Blender 5.2 local integration matrix passed all 28 runs, including four
truncation-recovery cases with zero geometry difference. Both adapter suites
passed 17 checks each. The existing Blender 5.2 UI regression suite passed all
14 checks, including intervening user edits, Edit Mode pause/resume, Undo during
import, scene changes and add-on disable.

## Live limitations observed

The original 90-second resource deadline was retained during live testing:

- **Aurillac:** a nested manifest exceeded the deadline. An independent ordinary
  urllib probe also read slowly: 179,892 of 570,812 declared bytes in about 25 s.
- **Clermont-Ferrand:** the root manifest timed out. An independent initial probe
  also timed out.
- **Strasbourg:** a nested manifest exceeded the deadline on recheck. Independent
  urllib probing also encountered a TLS handshake timeout. An earlier run exposed
  the Windows IPC reset reporting case described above.

These datasets are recorded as incomplete live tests, not successful imports or
proof of permanent provider unavailability. Their reports are in `inventory`,
`strasbourg_retry` and the `provider_*probe` logs. Tests with network failures do
not substitute for successful full import coverage. No timeout default was raised
to hide these observations.

## Reproduction

Use factory-startup background Blender for `tests/check_live_tiles.py`, with
`--python-exit-code 1`, followed by `-- --artifacts OUTPUT --datasets helsinki aarhus
agi san_francisco`. Add `--quick` for just four-worker blocking/progressive runs,
or `--inventory` to traverse manifests without downloading selected tile content.
Public URLs, bounds and geometric-error thresholds are explicit in
`tests/live_tile_datasets.py`.

Run `tests/check_live_tiles_ui.py` in a separate factory-startup GUI Blender
process with `-- --artifacts OUTPUT`. Inspect its JSON status as well as the exit
code, since exceptions in Blender timers do not control process exit codes.
Run `tests/test_download_resilience.py` and `tests/test_progressive_transport.py`
with ordinary Python, including Blender's bundled Python for compatibility checks.
