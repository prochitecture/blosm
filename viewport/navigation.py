"""Upright perspective panorama controls, independent of Blender and providers.

Heading is clockwise from world +Y, pitch is positive upward, and world +Z
is up. Keep view controls separate from the capture pose used to sample an image.
See mapillary/PanoramaControls.md for the browser-viewer research and choices.
"""
import math

PITCH_LIMIT = math.radians(89.)
MIN_FOV = math.radians(30.)
MAX_FOV = math.radians(90.)
MAX_HORIZONTAL_FOV = math.radians(120.)
ZOOM_TIME_CONSTANT = .09
ZOOM_EPSILON = math.radians(.01)
WHEEL_STEP = (MAX_FOV-MIN_FOV)*.05


def clamp(value, low, high):
    return max(low, min(high, value))


class PanoramaControls:
    def __init__(self, forward=(0., 1., 0.), right=(1., 0., 0.)):
        x, y, z = forward
        length = math.hypot(x, y, z)
        if not math.isfinite(length) or length == 0:
            raise ValueError('A finite viewing direction is required')
        if math.hypot(x, y) > length*1e-8:
            self.yaw = math.atan2(x, y)
        else:
            # At a pole the forward vector supplies no heading. The capture's
            # right axis gives a deterministic opening heading instead.
            self.yaw = math.atan2(-right[1], right[0])
        self.pitch = clamp(math.asin(clamp(z/length, -1., 1.)), -PITCH_LIMIT, PITCH_LIMIT)
        self.width = self.height = 1.
        self.fov = self.target_fov = math.radians(60.)

    @property
    def maximum_fov(self):
        horizontal_bound = 2*math.atan(math.tan(MAX_HORIZONTAL_FOV/2)*self.height/self.width)
        return min(MAX_FOV, horizontal_bound)

    @property
    def minimum_fov(self):
        # In an extremely wide, short region the horizontal bound wins. A
        # collapsed zoom range is preferable to violating that projection cap.
        return min(MIN_FOV, self.maximum_fov)

    @property
    def animating(self):
        return self.fov != self.target_fov

    def resize(self, width, height):
        self.width, self.height = max(1., float(width)), max(1., float(height))
        self.fov = clamp(self.fov, self.minimum_fov, self.maximum_fov)
        self.target_fov = clamp(self.target_fov, self.minimum_fov, self.maximum_fov)

    def drag(self, dx, dy):
        previous = (self.yaw, self.pitch)
        hfov = 2*math.atan(math.tan(self.fov/2)*self.width/self.height)
        # Grab the image: dragging right looks left; dragging up looks down.
        # Angular speed follows the visible FOV, rather than fixed radians/pixel.
        self.yaw = (self.yaw-dx/self.width*hfov+math.pi) % (2*math.pi)-math.pi
        self.pitch = clamp(self.pitch-dy/self.height*self.fov, -PITCH_LIMIT, PITCH_LIMIT)
        return previous != (self.yaw, self.pitch)

    def wheel(self, steps):
        """Positive steps zoom in. Change the target without a lens jump."""
        self.target_fov = clamp(self.target_fov-steps*WHEEL_STEP, self.minimum_fov, self.maximum_fov)

    def advance(self, seconds):
        """Time-based, monotone easing; independent of mouse/event/frame rate."""
        if not self.animating or seconds <= 0:
            return False
        self.fov += (self.target_fov-self.fov)*(-math.expm1(-seconds/ZOOM_TIME_CONSTANT))
        if abs(self.target_fov-self.fov) <= ZOOM_EPSILON:
            self.fov = self.target_fov
        return True

    def stop_zoom(self):
        self.target_fov = self.fov

    def matrix(self):
        """Camera-to-world rotation: columns are right, up, and backward.

        Up's world-Z component is cos(pitch), always strictly positive. No
        sequence of gestures can accumulate roll or cross a pitch pole.
        """
        sy, cy = math.sin(self.yaw), math.cos(self.yaw)
        sp, cp = math.sin(self.pitch), math.cos(self.pitch)
        return ((cy, -sy*sp, -sy*cp), (-sy, -cy*sp, -cy*cp), (0., cp, -sp))
