# Performance preference during Import

## User setting

Open Blender **Preferences → Add-ons → Blosm → Advanced settings** and use
**Prefer performance during Import**. The checkbox is enabled by default and
shown on Windows. It is an addon preference, shared by Import operations and
independent of the selected street-generation settings or source file.

When enabled, Blosm temporarily requests Windows HighQoS while performing
Import work. This can help when Blender loses focus during a lengthy Import.
Turning it off leaves execution policy entirely to the existing process/OS
settings. An existing explicit EcoQoS power-saving policy is always respected.
On other platforms, or when the API is unavailable, Import proceeds normally.

The request changes only the current Blender process's execution-speed QoS
control. It does not set CPU affinity, raise scheduling priority, change a
system power plan or save an OS setting. It can influence CPU core selection
and power management; the effect depends on the machine and execution context.
See Microsoft's [QoS guidance](https://learn.microsoft.com/en-us/windows/win32/procthread/quality-of-service)
and [SetProcessInformation documentation](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-setprocessinformation).

## Lifetime and failure handling

The synchronous Import operator owns a context manager around setup and
dispatch for OSM, terrain, imagery, GPX and 3D Tiles. Returning `FINISHED` or
`CANCELLED`, or raising an exception, exits the same scope and restores the
exact previous control/state masks. Unrelated policy bits are preserved.

Progressive 3D Tiles and modal imagery imports perform work after the initial
operator returns. Each progressive processing step, cancellation/finalization
and imagery timer step has its own scope. The previous policy is restored
between callbacks, while paused and when the job ends. No persistent request,
cleanup timer or OS-policy token is retained for an idle job. Blosm applies
the API request to the Blender process, without separately tagging download
helper processes.

Nested requests reuse an existing HighQoS state without taking ownership of it.
An explicit EcoQoS state is left unchanged. Missing or denied APIs cause a
quiet fallback, with details available through debug logging. A partially
applied/unconfirmed request is rolled back. Restoration errors are logged as
warnings and do not replace the Import result or its original exception.
Restoration therefore remains best effort if Windows itself rejects the
restoring API call.

## Implementation

[`app/import_performance.py`](../app/import_performance.py) contains the shared,
Blender-independent policy class and the preference-resolving factory.
The Windows API is loaded lazily through Python's standard `ctypes` module;
no extra package is required.

- `ImportPerformance` owns the reversible request and failure handling.
- `_WindowsProcessPolicy` reads/writes only the current process policy.
- `import_performance(context, addon_name)` resolves the addon preference,
  supporting both legacy and nested extension package names.

The call sites are the Import and imagery operators in
[`__init__.py`](../__init__.py), and the progressive controller in
[`progressive/blender.py`](../progressive/blender.py). New timer-driven Import
adapters should scope each processing/finalization step with the same factory.
Synchronous adapters invoked by the normal Import operator are already covered.

## Validation — 11 November 2026

Development used the isolated `blosm_streets_dev` checkout and factory Blender
profile, synchronized with `buildings_gn` at
`3a777d98f99d82b92f7072fb9992644ff91b1092`. The main checkout, normal Blender
installation and user sessions were not modified.

[`tests/test_import_performance.py`](../tests/test_import_performance.py)
passes 12 host tests: opt-out/unsupported platforms, explicit Eco/High,
unrelated bits, nesting, exceptions and interruption, API/read/write failures,
unconfirmed policy changes, restoration failures and preference resolution.
It uses fake policy APIs and changes no real process policy.

[`tests/check_import_performance.py`](../tests/check_import_performance.py)
passes in native Blender 5.2 using the nested extension namespace
`bl_ext.import_performance_check.blosm`. It checks registered RNA defaults,
actual `bpy.ops.blosm.import_data()` success/cancellation/error paths, an actual
missing-file parser failure, opt-out and API fallback, existing Eco/High,
imagery timer scopes and progressive work/pause/cancel/failure scopes.
The original masks were `(0, 0)`, active requests used `(1, 0)`, and all
completed scopes restored `(0, 0)`. Process and thread priorities were unchanged.
The five extension-contract and three vendor-import tests also pass. The
existing [progressive import regression checks](../tests/check_progressive_import.py)
pass with the production scopes enabled, including fixture geometry parity,
incremental publication, joining, error handling and lifecycle shutdown.

The complete Karl-Marx-Allee scene was imported in a fresh background Blender
process through the production preference, without an external test QoS
override. Whole Import took **68.865 s**, including **68.370 s** in street
generation. Coverage was COMPLETE: two networks, 328 meshes, 361,422 quads.
This is one validation run, not a new controlled speedup estimate.

The independent saved-scene checker passes against the pre-integration
automatic-policy result: normalized source data, prepared and terrain-deformed
float32 positions, quad topology, UVs, material names and construction metadata
are identical. See the earlier
[controlled experiment](STREETS_IMPORT_HIGHQOS_TEST.md) for timing comparisons
and its limits.

Artifacts are under
`C:/dev/workspace/artifacts/streets/performance/production_import_highqos_20261111/`:

- `native/checks.json`: native production integration checks.
- `karl_marx_allee/run_01/timings.json`: whole Import timing and settings.
- `equivalence/report.json`: independent saved-scene comparison.
- `progressive/checks.json`: existing progressive adapter regression checks.
- Per-test `process/` directories: resource-guard receipts and native logs.

Native runs use factory settings, autoexec disabled, 16 Blender threads,
one numeric-library thread, a 3 GiB process memory cap and 1 GiB system reserve.
No working-set throttle or CPU affinity is used. Generated JSON, scenes and
logs stay in artifacts and must not be staged.

## Repeating the checks

Run the host tests from the source checkout:

```powershell
python tests/test_import_performance.py
```

Launch `tests/check_import_performance.py` in isolated factory-startup Blender,
with `--python-exit-code 1 --disable-autoexec` and `-- --output <fresh-artifact-directory>`.
Use the chat's `BLENDER_USER_*` profile overrides and the resource guard.

For the full production Import benchmark, add `--operator`, `--scene <scene.blend>`
and `--import-performance` to [`tools/benchmark_streets.py`](../tools/benchmark_streets.py).
Use `--qos system`, the default, so there is no external test-process override.
For example, adapting paths/profile to the current development chat:

```powershell
python tools/benchmark_streets.py `
  --osm C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm `
  --scene C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend `
  --operator --import-performance --runs 1 --threads 16 --numeric-threads 1 `
  --output C:/dev/workspace/artifacts/streets/performance/production_qos_new_run `
  --profile-root C:/dev/workspace/artifacts/blender_profiles/<chat-id>
```

Controlled benchmark and GUI comparison helpers deliberately disable the
production preference unless `benchmark_streets.py --import-performance` is
supplied. This preserves their system-managed control and explicit `--qos`
comparisons. Their receipts distinguish the production preference from an
external test override; those modes cannot be combined. Normal user Imports
retain the enabled-by-default addon preference.
