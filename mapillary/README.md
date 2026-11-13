# Mapillary camera import

The [native viewport API](../viewport/README.md) implements immersive panorama
viewing and reusable controls. The
[framework proposal](../docs/VIEWPORT_FRAMEWORK_PROPOSAL.md) describes future map
selection, content panels and smooth transitions.

Choose **Mapillary cameras**, set the geographic bounds and your **Mapillary
client access token** in Blosm preferences, then import. Use the existing
**Progressive import** checkbox for interactive loading; disable it for a
synchronous result. Command-line mode blocks. Progressive mode requires an
interactive window and Object Mode. Camera import itself downloads no photographs
or thumbnails.

### View a panorama

Select an imported spherical Mapillary camera and click **View Mapillary panorama**.
Drag with the left mouse button to look around, use the wheel to zoom, and adjust
**Photo opacity** to reveal the Blender scene. **Esc** or **Exit panorama** removes
the overlay and keeps the final viewpoint, orientation and field of view.
The imported camera, scene rendering camera and render settings are unchanged.

The first view downloads and packs the usual 2048-pixel preview; later views reuse
the image. If an import is running, the viewer queues its request and remains
closable. Failed downloads show a retry control. Opacity is session-local; it
starts from an existing flat-photo background's opacity when available.
Controls adapt to viewport size and sidebars. Mirrored/ambiguous image layouts,
quad view, locked views and Local View are currently unsupported.

The photograph and 3D scene use the same view matrices. Geographic accuracy still
depends on the capture pose and assumed 2 m camera height. No depth or parallax is
reconstructed. Smooth transitions and navigation between captures are next-step
work. Use the existing **View Mapillary photo** button for a flat preview.

### View a camera's flat photo

Select an imported Mapillary camera and click **View Mapillary photo** in the
Blosm panel. The button appears only for the active Mapillary camera in the
current scene, regardless of the panel's selected import type. It sets the
scene's rendering camera and switches the 3D View into that camera. The photo
fills the camera frame as an opaque front background; the scene's render
resolution/aspect is set to the photo dimensions with square pixels. Camera
position, rotation and lens remain unchanged. Background photos are viewport
references, not rendered geometry or a compositor background.
Use **Photo opacity** below the button to see the scene through the photo:
1 is opaque and 0 is fully transparent. Each camera retains its opacity when
its photo is reopened, including after saving the `.blend` file.

The first click downloads a preview of up to 2048 pixels from Mapillary. Later
clicks reuse the image in Blender's image data, identified by the Mapillary image
ID even if its camera or image is renamed. Images are packed into the `.blend`
and kept as cache data, so saving/reopening works without another download or a
token. Delete an image datablock to discard its cache. Many viewed photos increase
scene memory and saved-file size; there is no automatic cache eviction.

Photos download using the shared HTTP helper and import controller, including
pause, cancellation, grouped undo and lifecycle cleanup. Cached photos can be
viewed during an import. An uncached request waits until the current Blosm job
finishes or is cancelled; only the latest queued request is retained. Changing
the selected camera cancels a stale queued request. A finished download is cached
without changing the view if the user has selected another object meanwhile.
Undo, file loading and disabling the addon clear pending work. The panel offers
**Cancel queued photo** while a request is waiting.

The image's contributor (when available), Mapillary/CC-BY-SA attribution and a
source link appear below the button. OAuth credentials are sent only to the
Graph API, never to the photo CDN. Signed photo URLs and credentials are not
stored in scene/image metadata. Temporary download files are packed and removed.

The downloaded file's embedded EXIF orientation is applied once; the original
capture's Graph EXIF tag is not reapplied to an already upright preview. Flat
panorama previews from this flat-photo action show the complete equirectangular
image. Use the separate panorama action for an interactive 360-degree view.
The flat-photo action does not correct fisheye distortion.
Downloads are capped at 32 MiB, decoded images at 8 MiPixels and each side at
8192 pixels. Native image decoding/packing and optional pixel rotation run on
Blender's main thread and can cause a short pause despite asynchronous HTTP.

If the preference token is empty, the importer falls back to
`mapillary/mapillary_token.txt` next to the installed adapter. In the main
development checkout this is
`C:/dev/workspace/blosm/mapillary/mapillary_token.txt`. The file is read lazily,
at most once per Blender session, and its trimmed contents are cached in memory.
Missing, empty and unreadable outcomes are cached too. Restart Blender to use an
edited/new file after its first lookup. A nonempty preference token takes priority
immediately and does not read the file; clearing preferences reuses the cached
fallback. File credentials are never copied into preferences or scene data.

Each image is represented by an actual Blender camera in a dedicated collection.
The scene's rendering camera is not changed. Selection, active object, active
layer collection and cursor are preserved. The shared controller provides pause,
grouped undo, cancellation and lifecycle cleanup. Cancellation retains completed
cameras. Deleting/renaming a camera does not stop remaining work. Removing the
destination collection stops with an error; renaming it preserves its identity.

### Progress messages

The status starts as `Mapillary: Estimating… · 0 cameras created`, then displays
an estimated percentage, for example `Mapillary: ~37% · 1,850 cameras created`.
The creation counter is exact and cumulative, including cameras the user later
deletes. Error counts appear when nonzero. The final successful job state shows
100%; cancellation and fatal errors never promote a partial estimate to 100%.
The aggregate completion message uses "camera objects created".

Tile fetching keeps the same bounded worker-sized windows; there is no initial
download of all tiles. Before importing a downloaded tile, a cooperative first
pass counts its unique, in-bounds image IDs. A second pass streams its records
into the existing 32-entry metadata/camera batches. The per-tile `ImageTile`
reader retains compact byte offsets for dictionary values and selected features,
so the import pass reuses the count's selection instead of rebuilding dictionaries
and scanning all features again. It retains no whole-tile decoded record list or
camera metadata. These offsets and the raw tile are released after replay, on
cancellation and on failure.

Points outside the selected bounds do not allocate property dictionaries.
Dictionary scalars are decoded only when a selected feature references them;
counting validates all referenced attributes, while retaining only the ID.
Outer protobuf framing, point geometry, selected-feature tags/values and resource
caps remain checked. Unreferenced scalar contents and out-of-bounds attributes
are intentionally skipped.

The estimated total is the unique candidate count multiplied by the ratio of
selected area to successfully counted area. Each tile contributes only its
intersection with the selected bounds, using geographic rectangle area. Empty
valid tiles contribute area; failed/malformed tiles do not become empty samples.
Before a positive usable sample, the message remains "Estimating…". Once every
tile is counted successfully, the candidate total needs no extrapolation.

Percentage tracks finished candidate processing, including rejected entries;
it is distinct from the number of cameras created and from elapsed-time progress.
Density varies across the area, so subsequent samples may revise the percentage
downward. Running percentages are capped at 99 until the job finishes. Normal
completion can still report per-resource errors or incomplete coverage.

Counting adds a tile scan before that tile's first cameras appear. Its reduced
decoder and roughly 2 ms processing chunks avoid spending UI timer intervals
on tiny checkpoint operations. During this read-only counting phase, the shared
controller uses a 5 ms timer interval while keeping its 8 ms processing budget.
Camera commits and pauses keep the normal 30 ms interval. Network latency and
dense tiles still delay the first estimate; an empty sample cannot predict image
density. Checkpoints keep
the scan cooperative, but an individual parsing/Blender call can exceed the
scheduler's soft budget. See `Validation.md` for measured count-phase timings.

## Adapter structure

- `tiles.py`: zoom-14 tile enumeration, cooperative protobuf/MVT point decoding,
  exact IDs and conversion of tile-local geometry to geographic coordinates.
- `pose.py`: pure-Python pose, projection-basis and camera-intrinsics conversion.
- `placement.py`: non-destructive evaluated-terrain sampling on Blender's main thread.
- `photo.py`: conditional photo action, packed-image cache and shared-controller download adapter.
- `photo_data.py`: bounded Graph/CDN URL, image-header and embedded EXIF validation.
- `progress.py`: pure-Python overlap-area density estimates and image-work progress.
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
dictionaries into compact offsets, then streams image features, yielding
checkpoints throughout. The count pass indexes selected features for replay.
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
| MVT dictionaries | 1,048,576 keys and values each across the active tile; unused layers skipped before dictionary indexing |
| MVT replay indices | 1,048,576 selected feature spans; value/feature offset arrays together at most 16 MiB, excluding allocation overhead |
| MVT image layers | 256 per tile |
| MVT feature or dictionary string | 64 KiB |
| Counted/imported image ID sets | Each capped at 100,000 IDs; per-tile temporary IDs share the counted-ID ceiling |

The ID limit is a bookkeeping ceiling, not a recommendation to create 100,000
cameras with undo enabled. Four workers reserve at most 128 MiB of raw tile bytes
and another 128 MiB in the HTTP window; an active consumer, dictionaries, copies,
up to 16 MiB of compact tile indices, Blender data and undo add overhead. This is
not a total process-RSS guarantee.
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

Camera elevation follows the local scene ground, using an **estimated 2 m camera
height**. Select the mesh with Blosm's existing **Terrain** object selector. Cast
a vertical ray onto its evaluated world-space surface, including modifiers and
object transforms, and place the camera 2 m above the uppermost intersection.
Terrain can be hidden or excluded from the active view layer. At the first camera,
the sampler evaluates a private temporary object copy and builds an independent
BVH snapshot. The temporary object and collection are removed synchronously
before any yield, redraw or undo checkpoint; the source is not unhidden or edited.
The snapshot fixes geometry and transforms for that import. Restart the import
after terrain edits if new cameras must use the changed surface. Completion,
cancellation and failures release the snapshot.
Without a selected terrain, place cameras at **scene Z = 2 m**. If a ray misses
the selected terrain, also use Z=2 m, record `scene_plane_outside_terrain`, and
include the number of these fallbacks in the aggregate result.

This policy does not convert provider altitude into height above terrain.
Preserve finite `computed_altitude`, otherwise original `altitude`, as
`source_altitude` with `altitude_source`; both original values also remain in the
requested metadata JSON. Zero and negative source altitudes are valid data and
are never changed to their absolute value. A scene/terrain `height_offset` is
recorded for diagnostics, together with `provider_altitude_scene_z` when possible,
but is **not subtracted from terrain surface coordinates** or the assumed height.
`height_offset_applied=false` makes this explicit.

Every camera records `vertical_position_source` (`terrain_surface`, `scene_plane`
or `scene_plane_outside_terrain`), `ground_elevation`, `camera_height=2` and
`camera_height_source=assumed`. This estimated height makes the resulting pose
approximate, even when full computed orientation is available. A terrain below
world Z=0 can legitimately produce a negative camera Z: the camera remains 2 m
above that terrain. The scene's datum, terrain geometry/transforms and rendering
camera are not changed. The selected terrain is captured by object identity;
renaming it is supported, while deleting/unlinking it stops remaining import work.

Original EXIF altitude is documented as above sea level; the primary API's
computed-altitude description does not establish a compatible terrain datum.
Every camera retains `vertical_datum=provider_unspecified; no geoid conversion`.
This ground placement is a visualization estimate, unsuitable for recovering
actual capture heights on bridges, elevated platforms or aerial images.

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

`orientation_complete` means full rotation and valid, non-mirrored EXIF orientation
are available. `provider_pose_complete` additionally requires provider altitude;
it is not an accuracy or vertical-datum assertion. Ground-anchored cameras have
`pose_complete=false` and `pose_approximate=true` because their height is estimated.
Optical projection limitations are recorded separately.

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
See [worldwide altitude samples](AltitudeSampling.md) for measured distributions
of provider `computed_altitude` across ten areas and the limits of interpreting it.

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

### Orientation diagnostics

`tests/check_mapillary_upright.py --scene <blend> --artifacts <directory>` inspects
actual saved camera world transforms against Mapillary's up/viewing conventions.
It reports downward local +Y axes; it does not assume every source pose is upright.

For a dense area, `tests/check_mapillary_area_pose.py --osm <fixture.osm>
--token-file <local-token-file> --artifacts <directory>` checks the OSM file's
bounds through the production discovery/metadata/pose pipeline without allocating
thousands of camera objects. It uses the production request and ID limits,
a ten-minute test deadline, and retains at most 32 downward metadata samples plus
four upright references. This measures metadata processing rather than viewport,
camera-object or undo overhead. Its `area_pose.json` can be replayed offline with
`check_mapillary_upright.py --pose-report <area_pose.json> --artifacts <directory>`
to create `sample_cameras.blend` through the real camera-creation path.

Run these Blender scripts only in a separate factory-startup process, with
autoexec disabled and all five isolated resource paths described above. Keep
all generated reports/scenes under `artifacts/mapillary_import/`.

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
