# Native viewport API, version 1

The shared session, drawing and controls are implemented. The first production
consumer is [Mapillary panoramas](../mapillary/panorama.py). Map selection, feeds,
transition animation and cross-service download admission remain separate work.
See the [framework proposal](../docs/VIEWPORT_FRAMEWORK_PROPOSAL.md) for that plan.

## Entry points and ownership

```python
from blosm import viewport
from blosm.viewport.example import Example

session = viewport.open(context, Example())
session.invalidate()  # schedule a redraw of the owner only
session.request_close()  # idempotent
```

In an installed extension, use relative package imports rather than assuming
`blosm` is its top-level package. `viewport.current()` returns the active session
or `None`. There is one interactive session across the addon; opening another
replaces the previous one. The registry is the only place enforcing that limit.

Call from the target 3D View, including its sidebar. The host resolves that
area's WINDOW region and never chooses another editor. Background Blender,
quad view and non-3D contexts are rejected. No network or scene transaction is
required by a surface.

Subclass `viewport.session.Surface`. Its default hooks are no-ops:

| Hook | Contract |
| --- | --- |
| `attach(session)` | Main-thread setup, capabilities validation and bounded preparation; call `super().attach(session)` |
| `layout(frame)` | Allocate logical rectangles; called when available size or UI scale changes |
| `handle(event)` | Return `ui.HANDLED` or `ui.PASS_THROUGH`; never block on IO |
| `tick(budget_seconds)` | Bounded preparation; return next polling interval, clamped to 50–250 ms; default 200 ms |
| `draw(frame, painter)` | Draw prepared resources; no decoding, network, IDs or undo writes |
| `cancel_gesture()` | Release canvas/widget capture on resizing or focus loss |
| `detach(reason)` | Release owned resources; called once, including failed partial attachment |

Set `requires_object_mode=True` for an adapter that depends on Object Mode.
Expose `capturing` when a canvas/widget drag owns input outside its rectangle.
Temporary view state and feature-specific target validation belong to the
adapter. Public `session.context()` provides a fresh validated owner override;
do not retain a `bpy.context` object between callbacks.

The host binds window, screen, area, space, region, scene and view layer. It closes
on owner replacement, editor/scene/layer change, undo/redo, file loading and addon
disable. Esc closes the session; it does not stop an unrelated import. An adapter
exception is reported once and closes its surface. Provider exceptions must be
sanitized before they are raised into this layer.

Closed modal instances finish on the next input event; graphics, assets, loading
subscriptions and the registry are cleaned immediately. The lightweight global
lifecycle timer stops at its next tick. No perpetual redraw is performed at idle.

## Layout, controls and composition

`Frame.region` is the full WINDOW rectangle. `Frame.available` excludes overlapping
sidebars, tool shelves and headers. All UI coordinates are logical pixels with a
bottom-left origin; `frame.scale` converts them to physical pixels. The host
normalizes mouse coordinates and increments `frame.revision` on layout changes.

`ui.Container` owns ordered children, clips drawing/hit testing and disposes its
children. Later children draw above and receive hits before earlier children.
`ui.Router` supplies pointer capture, hover, keyboard focus and Tab navigation.
Buttons support Enter/Space; sliders support arrows. A resize cancels an active
gesture at its last value. Use a surface's own geometry for domain interactions;
the router does not know about maps or cameras.

`Button`, `Slider`, `Rect`, `InputEvent` and `shell_layout` are Blender-independent
and tested with ordinary Python. Layout is deliberately small: assign rectangles
in `layout()`, show/hide optional children at meaningful widths, and use nested
containers where useful. There is no fixed list of feature modes, HTML renderer,
text-input implementation or general constraint solver.

[example.py](example.py) demonstrates a custom canvas, toolbar, opacity slider
and optional side card that collapses in a narrow viewport. Pass an existing
Blender image to display it instead of the plain canvas. It creates no scene IDs
and makes no requests. Run it from a 3D View context; it is not registered or
displayed automatically by the addon.

## Drawing and borrowed images

The host creates `draw.Painter` while a graphics context is available. It provides:

* `rect(rect, color)` and `text(text, x, y, max_width=..., size=...)`.
* `image(image_ref, rect, fit='contain', opacity=1)`; `cover` and explicit
  `stretch` are also supported.
* Nested `clip(rect)` and `pixel_clip()` for custom shaders.

`assets.ImageRef(image)` borrows a Blender image/texture and records its dimensions
and alpha convention. Construct it during preparation, not during drawing.
`release()` drops references and never deletes the image, packed data or fake user.
Transient image allocation/eviction is not provided by this first release.

The painter restores blend/depth/write-mask states. Geometry, fragment-shader
clipping and bounded text are used instead of altering Blender's scissor-enable
bit, whose getter is unavailable in the supported API. A custom shader must use
the current `pixel_clip()` (framebuffer pixel coordinates) to respect containers.
Text is single-line with ellipsis; future multiline/scrolling controls can extend
the UI layer. Color conversion explicitly displays scene-linear image texture
values as sRGB, respecting premultiplied alpha. It does not apply a second scene
view transform to the reference photograph.

## Loading and future services

`session.ticket(key)` returns an `AssetTicket(session_id, generation, key)`.
`session.accepts(ticket)` rejects publication after replacement/close or owner
invalidation. A feature advancing to a new content identity increments its session
generation. Layout revision is separate and does not invalidate unchanged images.

The panorama adapter borrows existing packed photos or runs a finite
`PhotoImport(activate=False, report_result=False)` through the shared import
controller. That task publishes image data without selecting a camera or changing
render settings. The surface polls for completion after the controller releases
its slot/undo checkpoint. It never starts a new task from a finish callback.
Closing cancels only the exact controller owned by the viewer. A queued panorama
does not allocate a second helper while another Blosm task is active.

The viewport host has no HTTP implementation. A future map/feed service must
implement the coordinated admission policy in the proposal before running
independent transports. `Job` owns and closes its transport; raw transports
cannot be shared among independently cancellable jobs. The current HTTP limits
are per helper, not global. Do not add one helper per card.

## Panorama-specific behavior

The Mapillary adapter owns image interpretation, pose validation, spherical shader
and the fixed-eye perspective view. The host knows none of those details. Dragging
uses local quaternion rotations (no forced horizon leveling); wheel zoom controls
vertical FOV from 20 to 110 degrees. Resizing recalibrates the actual viewport
projection to preserve that FOV. Opacity is session-local, initially inherited
from an existing photo background if one exists.

Esc or Exit panorama retains the final eye position, rotation and FOV in normal
perspective navigation. Camera objects, scene camera, render dimensions, selection
and cursor are unchanged. Failed startup can restore the pre-entry view. Lifecycle
teardown never writes a snapshot into the incoming scene or undo state.

The imported camera already includes the basic EXIF rotation. The adapter first
converts its local ray back to raw SfM axes, projects it spherically, and then
rotates UVs into the upright downloaded image. Rotating a 3D bearing is not the
same operation as rotating a 2D equirectangular image. See
[panorama_math.py](../mapillary/panorama_math.py), Mapillary's
[SphericalCamera](https://github.com/mapillary/mapillary-js/blob/main/src/geometry/camera/SphericalCamera.ts)
and [Transform](https://github.com/mapillary/mapillary-js/blob/main/src/geo/Transform.ts).
The shader explicitly blends wrapped longitude samples and clamps the poles.

Supported: spherical/equirectangular metadata with a full 2:1 source (or an
upright portrait equivalent), EXIF 1/3/6/8, and heading-only EXIF-1 poses marked
approximate by the importer. Mirrored EXIF, oriented images without computed
rotation, mismatched image dimensions, local view, locked views and quad view
are rejected clearly. Higher-resolution images, fisheye correction, navigation
between captures and all smooth transitions are subsequent work.

## Validation and handoff

Tests live under [tests](../tests): `test_viewport.py`, `check_panorama_ui.py`,
`check_panorama_loading.py` and `check_panorama_real.py`. Use the isolated Blender
profile/environment required by workspace `AGENTS.md`; UI tests need
`--enable-event-simulate`. Only an opt-in real HTTP test uses `--online-mode` and
`--token-file`. All outputs go to the supplied `--artifacts` directory.

Use the committed Mapillary implementation as the dependency for the area-selector
chat, preserving that chat's own checkout/reference/profile. Start with the
provider-free example and public APIs above. Do not copy internal modal or GPU
ownership code into another feature. Request additional public hooks when a
concrete consumer needs them. A transition service is intentionally not shipped
yet; neither photo nor panorama transitions are part of this implementation.
