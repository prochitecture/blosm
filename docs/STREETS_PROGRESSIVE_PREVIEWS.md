# Early mesh previews during street import

Development record, 13 November 2026. This extends the
[cooperative Python street import](STREETS_PROGRESSIVE_IMPORT.md), using the
existing shared OSM scheduler and allocation ownership protocol. Development
started from `f1eca6bb1c876faeddfd42bda875a9cbe7a4e93f` on this chat's
`streets_development` branch, synchronized with `buildings_gn`.

## Visible behavior

Enable **Progressive import** and import streets. Complete mesh parts now appear
while a connected network is still being calculated. A first available part is
shown immediately; subsequent parts appear in batches of up to eight. The
progress status distinguishes these **preview meshes** from completed import
units. Intersections, carriageways and separators can become visible at different
times according to their dependencies.

These meshes are provisional. Complete-network overlap checks and final terrain
placement still have to succeed. Early previews use prepared surface elevations;
terrain fitting can subsequently move them into place. Once validation finishes,
the same objects and mesh datablocks become the accepted output. There is no
second duplicate set of final street meshes.

**Cancel import** removes previews from the incomplete component and retains
previously accepted components. An expected construction failure removes the
failed component's previews but retains its normal diagnostic sources/report.
Source invalidation aborts the import instead of accepting stale geometry.

This first implementation applies to **new progressive imports only**. Blocking
imports and the explicit **Update** button keep their existing behavior. An
existing network's accepted output is not replaced by partial updates.

## Scheduling and validation boundaries

The acceptance unit remains a complete connected component. Showing a part does
not commit its numerical state, increment the completed-component count or make
its provisional surface valid for downstream editing.

Preparation still reads the complete input, resolves physical ownership and
retained boundary support, and prepares dependencies before the first suitable
part exists. Internal plans, trial boundary-support calculations and intermediate
arrays are not previewed. Some junction/Bundle tasks depend on multiple streets
or all junction surfaces. Long individual builders remain indivisible, so early
previews improve visible progress without guaranteeing a fixed maximum UI pause.

The construction formulas, task order, tolerances and final mesh topology are
unchanged. The pipeline is:

1. Prepare sources and dependencies with the existing revision guards.
2. As a top-level mesh part completes, allocate its owned Blender object/mesh,
   preserve its UV charts (or apply the standard fallback), assign its material,
   and show it in a separate provisional collection at a display checkpoint.
3. Finish the numerical update, UV preparation and full overlap validation.
4. Fit the complete terrain corridor and validate final stored coordinates,
   including the existing quad-quality and overlap checks.
5. Refine each preview mesh in place to the validated final coordinates and UVs;
   attach the usual terrain attributes/modifiers. Other output types use the
   existing detached publisher.
6. Recheck source/terrain revisions and preview integrity. In one non-yielding
   transaction, link outputs into the canonical collection, update the manifest,
   commit the raw numerical state and remove provisional links/tags.

The engine still commits raw, unplaced numerical products. Terrain offsets cannot
accumulate through repeated Updates. Preview publication never mutates the
upstream numerical products.

## Source layout and APIs

| Location | Responsibility |
| --- | --- |
| `way/street_network/python_backend/updates.py` | Frozen `ProductReady(key, product)` notification after a rebuilt product and its digest are ready; opt-in `emit_products=True` on `UpdateCoordinator.prepare_steps` |
| `way/street_network/python_backend/backend.py` | Forwards opt-in product events before UV fallback and complete-update validation; default producer contract still yields ordinary checkpoints |
| `way/street_network/python_backend/previews.py` | `PendingPreviews`: explicit allocation ownership, provisional collection, batching, integrity checks, promotion and rollback |
| `way/street_network/python_backend/scene.py` | Optional `NetworkScene.update_steps(..., on_preview=callback)` for an initial publication; consumes product events and reuses previews in `SceneWriter.publish_steps` |
| `way/street_network/python_backend/renderer.py` | Enables previews only for a progressive OSM session |
| `app/osm_import.py` | `preview_changed(count)` updates visible revision/status and preview timing metrics separately from accepted-component progress |

`on_preview(count)` reports currently displayed provisional meshes. It is also
called after an in-place display refinement and with zero after promotion or
cleanup. Visible changes advance the session revision so the shared driver can
record an appropriate grouped Undo checkpoint. `YieldToUI()` requests a redraw
opportunity after display batches; plain `None` checkpoints retain their existing
budgeted scheduling behavior.

The new final report fields are `first_preview_seconds`, `preview_batches` and
`peak_preview_meshes`. A preview batch includes an in-place visual refinement.
`first_result_seconds` and `committed` continue to mean accepted results; they do
not count provisional surfaces. Measurements from a UI heartbeat represent
visibility in Blender's scene, not a measurement of pixels presented on screen.

## Ownership, cancellation and history

Preview objects carry `st_py_preview=True` and ordinary controller/key ownership
metadata. They are linked outside `st_py_outputs` in a separate pending collection.
The collection is tagged as transient; all allocations join both the private
publication ledger and, during import, the parent session ledger. No global
allocation snapshot spans a yield.

Closing or failing the producer rolls back only explicitly owned pending IDs.
Unrelated objects and Text datablocks created by a user between checkpoints are
preserved. The common Undo/load cleanup removes tagged pending/transient IDs from
restored intermediate snapshots. Successful promotion strips the tags and removes
the now-unneeded provisional collection before the next yield.

Existing source, destination, terrain and pending-ID guards still apply. Editing
the sources while generation is pending invalidates the result; previews are
display output, not an editing interface. Final publication also verifies the
reused mesh topology and other preview state so external changes cannot silently
replace the validated output.

## Qualification and reproduction

Targeted checks live in:

- `tests/test_python_street_product_events.py`: opt-in/default producer behavior,
  event ordering, cancellation, timing and numerical result parity.
- `tests/check_street_previews.py`: early visibility within one connected
  network, pointer reuse, exact final parity with blocking generation, terrain
  modes, invalidation and rollback.
- `tests/check_street_progressive.py`: existing import/session regression checks.
- `tests/check_street_progressive_ui.py`: actual Blender timer, preview display,
  cancellation, scene-switch pause and native Undo/Redo; real-extract timing.

Run Blender tests with this chat's isolated factory profile, `--factory-startup
--disable-autoexec`, and a resource guard. Use `--background --python-exit-code 1`
for native checks; the UI test needs a separate GUI process. All scripts accept
`--artifacts <output-directory>` after Blender's `--`. The numerical test uses
Blender's bundled Python/NumPy without changing that installation. Keep outputs
under `C:/dev/workspace/artifacts/streets/early_previews_20261113/`.

For ordinary real-extract GUI timing, pass `--no-navigation
--performance-request --osm <path>` to `check_street_progressive_ui.py`. Its
default continuous-navigation stress mode intentionally adds viewport redraw
load and is not comparable to an idle import benchmark. The optional
`--expected-hash <digest>` asserts parity with an earlier exact final-output
record. Do not run multiple Blender qualification jobs concurrently.

### Results

All checks below passed with Blender 5.2.0 LTS, using the isolated factory
profile. Native runs and the GUI runs were sequential, with four Blender threads,
one numerical-library thread, a 3 GiB process cap and 1 GiB system reserve.

| Check | Evidence under `artifacts/streets` | Result |
| --- | --- | --- |
| Numerical producer behavior and existing numerical regressions | `early_product_previews_20261113/` | 35 tests passed: 9 product-event, 10 cooperative numerical and 16 core tests |
| Preview scene integration | `early_previews_20261113/native_03/` | 10 test methods passed, including 14 mutation/promotion-failure subcases |
| Existing progressive street integration | `early_previews_20261113/regression_01/` | 17 tests passed, including mixed native buildings/streets, terrain, source invalidation and failed-commit recovery |
| Actual GUI timer and native history | `early_previews_20261113/gui_01/` | Early visibility, scene-switch pause, cancellation, preserved user object/Text and repeated Undo/Redo passed |
| Karl-Marx-Allee production GUI import | `early_previews_20261113/gui_berlin_01/` | Success, two accepted components, exact prior final-output hash preserved |

The preview-specific suite compares blocking and preview-driven output in the
flat, PROFILE and DIRECT placement modes, checking coordinates, faces, UVs,
material state and metadata. It records Blender object and mesh pointers during
the preview and verifies those same pointers after acceptance. It also verifies
that existing-network Update ignores the preview opt-in.

The GUI tests use the actual shared driver, not a manually drained substitute.
The small three-component test additionally holds an early preview checkpoint
for deterministic Cancel and native history checks. Real-extract timing uses the
unmodified production scheduler, without that hold or simulated navigation.

Karl-Marx-Allee measurements (`gui_berlin_01`):

| Measurement | Result |
| --- | --- |
| First preview notification | 13.969 s |
| First preview visible to the GUI heartbeat | **14.028 s** |
| First accepted component visible | 73.878 s |
| Preview lead over first accepted component | **59.849 s** |
| Complete import elapsed time | **74.248 s** |
| Peak pending preview meshes | 327 |
| Display batches, including final refinements | 85 |
| Longest import callback | 3.141 s |
| Longest GUI heartbeat gap | 3.156 s |
| Peak resident / sampled private memory | 1.109 / 1.576 GiB |

The heartbeat observed the pending mesh count growing as 1, 9, 17, and so on up
to 327 before component acceptance. The first preview appeared during
`Building junction_grass:19`, while numerical construction was still in progress.
The final exact-output SHA-256 was
`5a9b3400076a78fb1d27a204d86a579920b0ba603f0e55ab07ec60f57098c4aa`,
matching the pre-preview qualification record. The previous comparable GUI run
took 77.657 s; these are individual runs, so the difference is not evidence of a
speed improvement. The established improvement is earlier visible progress,
with unchanged final output. Occasional multi-second synchronous work remains.
