# Upright panorama controls

Reviewed and implemented 14 November 2026. The previous input code multiplied
incremental camera-local quaternions. Moving around several axes accumulated
roll, and vertical movement could pass a pole and invert the view. The new
navigation state is heading, bounded pitch and bounded field of view. The
quaternion used by Blender is reconstructed from these values on every change.

## Primary implementations reviewed

These are pinned source revisions, rather than claims that every browser viewer
has the same settings. No third-party implementation was copied into Blosm.

| Viewer | Rotation and dragging | Zoom and motion | Sources |
| --- | --- | --- | --- |
| Pannellum | Separate yaw/pitch state; pointer displacement is scaled by FOV and viewport dimensions; render clamps pitch. Configurable horizon correction is separate. | Bounded horizontal FOV, wheel adjustments and damped motion; explicit friction setting. | [Mouse drag](https://github.com/mpetroff/pannellum/blob/818b698d2b60d9cd534a9395ea838d4d33ccb06d/src/js/pannellum.js#L871), [wheel](https://github.com/mpetroff/pannellum/blob/818b698d2b60d9cd534a9395ea838d4d33ccb06d/src/js/pannellum.js#L1137), [pitch bounds](https://github.com/mpetroff/pannellum/blob/818b698d2b60d9cd534a9395ea838d4d33ccb06d/src/js/pannellum.js#L1638), [reference](https://pannellum.org/documentation/reference/) |
| Photo Sphere Viewer | Yaw/pitch dynamics; pitch limited to ±90°. Drag deltas use horizontal/vertical FOV over the corresponding viewport dimension. | Default vertical FOV 30–90°, wheel zoom steps, normalized zoom-to-FOV mapping and configurable movement damping. | [Dynamics](https://github.com/mistic100/Photo-Sphere-Viewer/blob/eba89dedea7358f0278a5fa5e8e25bcc36def813/packages/core/src/services/ViewerDynamics.ts), [drag/wheel](https://github.com/mistic100/Photo-Sphere-Viewer/blob/eba89dedea7358f0278a5fa5e8e25bcc36def813/packages/core/src/services/EventsHandler.ts#L577), [defaults](https://github.com/mistic100/Photo-Sphere-Viewer/blob/eba89dedea7358f0278a5fa5e8e25bcc36def813/packages/core/src/data/config.ts), [FOV conversion](https://github.com/mistic100/Photo-Sphere-Viewer/blob/eba89dedea7358f0278a5fa5e8e25bcc36def813/packages/core/src/services/DataHelper.ts#L53) |
| Marzipano | Separate yaw/pitch/roll; default roll zero. Drag controls scale input by viewport size and view FOV. Its traditional limiter restricts pitch to ±90°. | Combined horizontal/vertical FOV caps and an image-resolution limiter; drag-release velocity with friction is supported. | [View/control conversion](https://github.com/google/marzipano/blob/7936a4e0791fcfe55d199c03080be4f89ccf71bd/src/views/Rectilinear.js#L568), [traditional limiter](https://github.com/google/marzipano/blob/7936a4e0791fcfe55d199c03080be4f89ccf71bd/src/views/Rectilinear.js#L978), [drag controls](https://github.com/google/marzipano/blob/7936a4e0791fcfe55d199c03080be4f89ccf71bd/src/controls/Drag.js), [limiter reference](https://www.marzipano.net/reference/RectilinearView.limit.html) |

The common pattern is constrained heading/pitch, FOV-aware dragging and bounded
zoom. Roll/horizon metadata and interactive look controls are distinct. This
research concerns the public source of these viewers; it does not establish how
Google's proprietary Street View implementation is written.

## Blosm behavior

- World +Z stays up. Opening preserves the capture's forward direction and
  removes its roll from the **view controls**. The imported camera's complete
  rotation remains the shader's world-to-capture reference. The photograph and
  scene continue to use the same rays; leveling is not a pixel rotation or a
  change to camera metadata. Incorrect source pose metadata can still cause
  incorrect photo/scene alignment.
- Pitch stops at ±89°, leaving a margin from the pole singularity. Reversing
  direction at a limit responds immediately. Heading wraps through full turns.
  The camera's right axis stays horizontal and its up axis always has a positive
  world-Z component, regardless of gesture history. For an exactly vertical
  source direction, its right axis supplies a deterministic opening heading.
- Dragging **right moves the image right** (look left); dragging up looks down.
  This matches the grabbing interaction in the reviewed browser viewers. Drag
  sensitivity follows current visible FOV, so zooming in slows angular movement.
  Blender's click/drag threshold and Shift-click object selection remain active.
- Opening vertical FOV is 60°. Its normal zoom range is 30–90°, with a 120°
  horizontal cap to avoid extreme wide-screen stretching. At unusually wide,
  short viewport sizes the horizontal cap reduces the allowed vertical range;
  the range can collapse when even 30° vertical would exceed that cap.
- Each wheel notch moves the zoom target by 3° (5% of the normal range).
  Numpad Plus/Minus use the same controls. FOV approaches that target with
  time-based exponential easing (90 ms time constant), with no overshoot and a
  final 0.01° snap tolerance. Frame rate does not determine the zoom distance.
- Dragging follows the pointer directly and stops at release. There is no
  rotational coasting after release: precise object picking and the retained
  exit orientation take precedence over optional browser-style release inertia.
  Starting a selection press stops a zoom tail at the displayed projection.
- Resize/focus loss stops the active gesture/zoom at its current value. Resize
  recalibrates the native projection and only changes FOV when a bound requires
  it. Exit retains the actual current view, including a partly completed zoom.
- Native orbit/frame/view/walk shortcuts exit before proceeding. A viewport
  gizmo or custom operator can bypass key handling: a read-only draw guard then
  suppresses photo sampling until the next event/tick closes the viewer. The
  viewer never draws the photo through an externally inverted projection or
  overwrites the user's resulting view. Object-mode exit and lifecycle cleanup
  remain handled by the shared controller.

For logical viewport width `W`, height `H`, current vertical FOV `v` and drag
delta `(dx,dy)`:

```
h = 2 atan(tan(v/2) W/H)
heading -= dx h/W
pitch    = clamp(pitch - dy v/H, -89°, +89°)
```

Heading is clockwise from +Y; pitch is positive up. If `s`/`c` denote sine/cosine
of heading and pitch respectively, the camera basis is constructed from a
horizontal right vector, a forward vector at the bounded elevation, and their
cross product. It is a right-handed orthonormal basis with no accumulated roll.

## Shared framework additions

`viewport/navigation.py` contains Blender-independent `PanoramaControls`; its
initial direction and up convention have no Mapillary dependency. It supplies
angles, effective projection bounds, dragging, zoom targets/easing and the final
camera-to-world rotation matrix.

The existing shared session timer now supports an opt-in animation hook:
`Surface.animating`, `Surface.animate(elapsed_seconds)` and
`session.request_animation()`. Animation aims for 60 Hz only while requested.
Ordinary bounded `tick()` work retains its requested 50–250 ms polling cadence.
There is one timer, no decoder/download worker per frame and no perpetual idle
redraw. Wall-clock steps are capped at 100 ms after stalls. Starting animation
from a tick does not register another timer; closing/replacing a session removes
its drawing and stops its animation through the established lifecycle.

Other surfaces keep their previous behavior unless they opt in. This is a frame
hook, not an implementation of photo transitions, transition loading or a new
scene transaction controller. The latter remain separate work.

## Verification

Use `tests/test_panorama_navigation.py` for numerical invariants, including
10,000 mixed gestures, complete turns, pole reversal, FOV/size scaling, zoom
limits, frame-rate-independent easing and cancellation.

`tests/check_panorama_controls.py` checks actual native camera axes/projection,
rolled source poses, both pitch limits, smooth shared-timer zoom, native wheel
bursts, resize, focus loss, selection, external navigation guards, exiting during
zoom and mode/lifecycle teardown. See [validation](Validation.md) for measured
results and artifacts. Launch Blender only with the isolated profile and startup
flags required by workspace instructions.
