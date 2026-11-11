# Mapillary camera import

Choose **Mapillary cameras**, set the geographic bounds and your **Mapillary
client access token** in Blosm preferences, then import. Use the existing
**Progressive import** checkbox for interactive loading; disable it for a
synchronous result. Command-line mode blocks. Progressive mode requires an
interactive window and Object Mode. No photograph or thumbnail is downloaded.

Each image is represented by an actual Blender camera in a dedicated collection.
The scene's rendering camera is not changed. Selection, active object, active
layer collection and cursor are preserved. The shared controller provides pause,
grouped undo, cancellation and lifecycle cleanup. Cancellation retains completed
cameras. Deleting/renaming a camera does not stop remaining work. Removing the
destination collection stops with an error; renaming it preserves its identity.

## Adapter structure

- `tiles.py`: zoom-14 tile enumeration, cooperative protobuf/MVT point decoding,
  exact IDs and conversion of tile-local geometry to geographic coordinates.
- `pose.py`: pure-Python pose, projection-basis and camera-intrinsics conversion.
- `progressive.py`: discovery producer, metadata batches and cooperative camera
  consumer using production `Job`, `Download`, `Drain`, `Controller` and
  `DownloadTransport`. Both modes share this processing path.

The existing hidden `concurrentDownloads` / `downloadWorkers` settings apply.
Turning concurrency off uses one outstanding request. No concurrency UI or
dependency on the earlier downloader mockup is introduced.

The geographic projection comes from `BlenderApp.setProjection()`; scene origin
selection matches other Blosm imports. Source coordinates and pose provenance
are stored as `mapillary_*` object properties. Image IDs remain decimal strings,
including protobuf integers greater than 2^53; float IDs are rejected. Sequence
ID, capture time, intrinsics, original/computed metadata and creator information
when returned are retained. `mapillary_image_url` links to the source image.
Credentials belong in preferences/requests, never in scene properties or logs.

## API evidence and authentication

Computed coverage is requested at:

```text
https://tiles.mapillary.com/maps/vtp/mly1_computed_public/2/{z}/{x}/{y}
```

The `access_token` query is added at request construction. Metadata requests use
`https://graph.mapillary.com/images?image_ids=...&fields=...` with an
`Authorization: OAuth ...` header. Transport authentication scopes contain only
a token hash. URLs and tokens travel via shared IPC, not command-line arguments
or dataset files. Existing online-access checks and TLS verification apply.

Public MapillaryJS `GraphQueryCreator` / `GraphDataProvider` source, inspected on
2026-11-10, confirms the multi-image query and requested fields. A batch of 32 is
an internal memory/processing bound, not a verified provider maximum. Fields are
`id`, `geometry`, `computed_geometry`, `altitude`, `computed_altitude`,
`computed_rotation`, `compass_angle`, `computed_compass_angle`, `camera_type`,
`camera_parameters`, `width`, `height`, `exif_orientation`, `sequence`,
`captured_at`, and `creator`.

The current main API documentation required login during initial implementation.
Public implementation evidence was checked then; subsequent credentialed results
are recorded in `Validation.md`. No
undocumented Mapillary error-reason codes were added to shared authentication
suppression. HTTP 401 stops as authentication failure. HTTP 400/403 and Graph
error envelopes stop as rejected requests without classifying every rejection
as a credential error. Transport retries, deadlines, rate-limit cooldowns and
sanitized diagnostics are reused.

## Discovery and processing bounds

The MVT reader respects the layer's actual extent, top-left origin and Y-down
coordinates, including buffered points outside tile borders. It first scans
dictionaries, then streams image features, yielding checkpoints throughout.
There is no full-area decoded feature list. Coverage entries are clipped against
the selected bounds and deduplicated by image ID. The chosen metadata position
is clipped again; discovery uses computed coverage, not a union of computed and
original locations.

At most one worker-window of raw coverage tiles is retained. One tile is decoded
at a time, with bounded 32-image metadata requests. Camera consumers yield after
each entry. Network requests and unconsumed replies share the production window,
including during pauses. Deduplication state is capped and released on cleanup.

| Resource | Explicit bound |
| --- | --- |
| Coverage area | At most 4096 tiles; larger selections rejected before downloads |
| Downloaded/decoded response | 32 MiB per resource, enforced by transport |
| JSON metadata parsing | 256 KiB and 32 entities per batch |
| MVT dictionaries | 1,048,576 keys and values each; unused layers skipped before dictionary allocation |
| MVT feature or dictionary string | 64 KiB |
| Global image ID set | 100,000 IDs; stop with a smaller-area error at the limit |

The ID limit is a bookkeeping ceiling, not a recommendation to create 100,000
cameras with undo enabled. Four workers reserve at most 128 MiB of raw tile bytes
and another 128 MiB in the HTTP window; an active consumer, dictionaries, copies,
Blender data and undo add overhead. This is not a total process-RSS guarantee.
See `Validation.md` for measured practical limits.

Metadata request/absent-entity failures fall back to available coverage data and
mark the camera approximate/missing metadata. Coverage failures are counted as
incomplete coverage, never interpreted as proof that no imagery exists. Error
details are bounded to 20 sanitized strings; one aggregate result reports totals.
Oversized/malformed payloads are errors. Dateline-crossing selections and
latitudes outside Web Mercator coverage are rejected in this version.

## Documented pose conversion

### Position and altitude

Prefer valid `computed_geometry`, then original `geometry`, then coverage
position. Source selection is recorded independently for position, altitude and
orientation. Latitude/longitude use WGS84 degrees and Blosm's chosen projection.

Prefer finite `computed_altitude`, then original `altitude`; zero/negative values
are valid. Subtract a finite scene `height_offset` if present, otherwise the
selected terrain object's `height_offset`, otherwise zero. Record offset/source;
do not introduce another scene height datum or change an existing offset.

When altitude is absent, place the camera at **scene Z = 0** and mark
`altitude_source=missing_scene_plane`, incomplete and approximate. This is a
visualization fallback, not sea level or height above terrain. No terrain raycast,
fixed eye height, geoid adjustment or assumed ellipsoid/MSL conversion is used.
Original EXIF altitude is documented as above sea level, but datum compatibility
with a terrain provider is not established. Every camera records
`vertical_datum=provider_unspecified; no geoid conversion`.

### Rotation and EXIF

Mapillary `Transform.createWorldToCamera()` constructs world-to-camera `R` from
the angle-axis vector: length is radians, normalized vector is axis. Zero is the
identity, not an assumed north heading. OpenSfM georeferenced world coordinates
are east/north/up; camera coordinates are right/down/forward (+X/+Y/+Z).
Blender camera coordinates are right/up/back (+X/+Y/+Z), viewing down -Z.

`Transform._createBasicWorldToCamera()` left-multiplies `R` by image orientation
`E`: identity for EXIF 1, 180 degrees for 3, 90 degrees for 6 and 270 degrees for
8, around camera Z. The Blender camera-to-world matrix is:

```text
B = G * transpose(E * R) * diag(1, -1, -1)
```

`G` is the orthonormal local east/north basis in Blosm's projected grid, calculated
from the projection's east derivative to account for grid convergence. Projection
scale does not shear camera axes. Scene Z stays vertical; this is a projected-map
import, not a curved-Earth reconstruction import.

Tests compare the inverse direction against Mapillary `Spatial.viewingDirection()`
and the up axis against `Transform.upVector()` for portrait EXIF orientations.
The raw computed orientation receives EXIF correction once. Mirrored EXIF
2/4/5/7 are retained and flagged approximate with unsupported optical mirroring:
a reflection cannot be represented by a proper camera rotation. Missing/invalid
EXIF defaults to 1; absent/invalid values are marked assumed and approximate.

Without computed rotation, prefer computed heading, original heading, then
coverage heading. Heading is clockwise from geographic north. A horizontal,
upright camera is assumed; pitch and roll are unknown. Missing heading assumes
north and records that assumption. These fallbacks are always approximate and
already upright, so no second portrait correction is applied.

`pose_complete` means full rotation, altitude and valid EXIF orientation are
present without unsupported mirroring. It is not an accuracy or vertical-datum assertion. Optical projection
limitations are recorded separately.

### Intrinsics and projection limits

OpenSfM focal length is normalized to the **larger original image dimension**.
For the adapter's nominal 36 mm horizontal sensor:

```text
f_pixels = focal_normalized * max(original_width, original_height)
lens_mm = 36 * f_pixels / upright_width
sensor_height_mm = 36 * upright_height / upright_width
```

Upright dimensions swap for EXIF 5-8. Horizontal sensor fit is explicit. The
36 mm width is a representation convention, not a claim about the real sensor.
Missing/invalid/out-of-range intrinsics retain a marked 50 mm / 36x24 mm proxy.
The scene render aspect is preserved; a later rendering feature must use stored
upright dimensions for the intended perspective aspect.

Radial distortion parameters are preserved but not applied by Blender's native
perspective camera. Fisheye/unknown types use marked perspective proxies.
Spherical/equirectangular types create a PANO/EQUIRECTANGULAR camera where
supported; renderer/viewport behavior varies. Pose representation does not imply
an accurate optical projection or a completed photograph viewer.

## Shared-core additions

A `Download.consume(data, error)` callback may return a checkpoint iterator.
`Job.step()` advances one active consumer within the same time/step budget as its
producer. Paused `poll()` never advances it; `Drain` includes its completion.
Cleanup closes consumer and producer generators, then transport. Consumer yields
are checkpoints only; network commands still come from the producer. Existing
ordinary callbacks returning `None` retain their behavior.

This reusable extension lets Mapillary parse and create cameras cooperatively
without a separate scheduler. It cannot preempt one Blender API call or JSON
parse. Parse sizes/feature work are bounded and timings recorded. The shared
diagnostic redactor also skips empty secret suffixes from incomplete Authorization
headers, avoiding runaway redaction-marker insertion on errors.

## Validation and local development

Checkout: `C:/dev/workspace/mapillary_dev/blosm`, branch `mapillary_import_release`,
reference `release`, initially `72286b3e77e8cd3c6021384e856291d7315a6a07`.
The main checkout and normal Blender installation/profile are not modified.
Generated outputs belong under `C:/dev/workspace/artifacts/mapillary_import/`.

Run `tests/test_mapillary.py` in standalone Python for pose/decoding/core checks.
Run `tests/run_mapillary_checks.py` for background Blender; add `--ui` for real UI
timers, HTTP helper, undo and scale checks. `--count` selects synthetic image count.
The launcher uses factory startup, disabled autoexec and this chat's five isolated
resource paths. Online access is enabled only for the test process. It checks the
UI JSON report status, not just Blender's GUI exit code. No preferences are saved.
See `Validation.md` for results, measurements and remaining limits.

### Live development check

`tests/run_mapillary_checks.py --live` reads the credential from
`C:/dev/workspace/blosm/mapillary/mapillary_token.txt` without copying it into
the checkout or printing it. Add `--ui` for the real timer controller. Use
`--bounds west south east north` to select a small test area. Reports and an
imported camera scene go under `artifacts/mapillary_import/live` or `live_ui`.
The test has a 2,000-image safety ceiling and compares two live imports for
blocking/progressive ID and pose parity; changing provider data can affect parity.

The token filename is excluded by `mapillary/.gitignore`; the shared local Git
info/exclude also protects it across local branches. Never stage or commit a
token, even when collecting development fixtures or publishing a test scene.

## Primary sources inspected

- [Mapillary vector coverage](https://mapillary.github.io/mapillary-python-sdk/docs/mapillary.config.api/mapillary.config.api.vector_tiles/)
- [Mapillary ImageEnt metadata](https://mapillary.github.io/mapillary-js/api/interfaces/api.ImageEnt/)
- [GraphQueryCreator: fields and batching](https://github.com/mapillary/mapillary-js/blob/main/src/api/provider/GraphQueryCreator.ts)
- [GraphDataProvider](https://github.com/mapillary/mapillary-js/blob/main/src/api/provider/GraphDataProvider.ts)
- [Transform: rotation and EXIF](https://github.com/mapillary/mapillary-js/blob/main/src/geo/Transform.ts)
- [Spatial: viewing direction](https://github.com/mapillary/mapillary-js/blob/main/src/geo/Spatial.ts)
- [GraphCalculator: compass convention](https://github.com/mapillary/mapillary-js/blob/main/src/graph/GraphCalculator.ts)
- [GeoCoords: ENU](https://github.com/mapillary/mapillary-js/blob/main/src/geo/GeoCoords.ts)
- [OpenSfM geometry](https://opensfm.org/docs/geometry.html)
- [Current API docs (login-gated)](https://www.mapillary.com/developer/api-documentation/)
