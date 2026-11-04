"""Shared profile catalogue and footprint roof overrides, independent of Blender.

The compact FACE modes are an importer binding, not additional PADL syntax.
Mode zero preserves the authored roof settings. Both geometry backends consume
the same normalized profile points and the same directed-axis convention.
"""
from dataclasses import dataclass
import math

from .errors import PadlError


DIRECTION_MODES = {0: 'inherit', 1: 'longest_along', 2: 'longest_perpendicular', 3: 'vector2d'}


def normalized_profile(points):
    try:
        points = tuple(tuple(p) for p in points)
        if len(points) < 3 or any(len(p) != 2 or any(type(v) not in (int, float) or not math.isfinite(v)
                                                   or not 0 <= v <= 1 for v in p) for p in points):
            raise ValueError('Expected at least three finite normalized [x, height] points')
        if points[0][0] != 0 or points[-1][0] != 1 or any(a[0] >= b[0] for a, b in zip(points, points[1:])):
            raise ValueError('Profile x must increase strictly from 0 to 1')
        low, high = min(points[0][1], points[-1][1]), max(p[1] for p in points)
        if high - low <= 1e-6 or any(p[1] < low for p in points):
            raise ValueError('Profile needs a peak above its lower eave and cannot extend below that eave')
        return tuple((float(x), (y-low)/(high-low)) for x, y in points)
    except (TypeError, ValueError, OverflowError) as exc:
        raise PadlError(f'Invalid roof profile: {exc}', code='invalid_profile') from exc


@dataclass(frozen=True)
class RoofProfile:
    mode: int
    name: str
    points: tuple
    osm_shape: str


class RoofProfileRegistry:
    """Stable preset codes; extensions can supply additional profile records."""
    def __init__(self, profiles=None):
        if profiles is None:
            arc = tuple(((1-math.cos(i*math.pi/16))/2, math.sin(i*math.pi/16)) for i in range(17))
            profiles = (
                RoofProfile(1, 'gabled', ((0, 0), (.5, 1), (1, 0)), 'gabled'),
                RoofProfile(2, 'gambrel', ((0, 0), (.2, .6), (.5, 1), (.8, .6), (1, 0)), 'gambrel'),
                RoofProfile(3, 'round', arc, 'round'),
                # Equal pitch, offset ridge and unequal eaves: the OSM catslide
                # shape, not the trapezoid historically called saltbox in Blosm.
                RoofProfile(4, 'saltbox', ((0, 0), (.625, 1), (1, .4)), 'saltbox'),
            )
        self.profiles = {p.mode: RoofProfile(p.mode, p.name, normalized_profile(p.points), p.osm_shape) for p in profiles}

    def get(self, mode):
        if type(mode) is not int or mode not in self.profiles:
            raise PadlError(f'Unknown roof_profile_mode {mode!r}', code='invalid_category')
        return self.profiles[mode]


PROFILES = RoofProfileRegistry()


@dataclass(frozen=True)
class RoundProfilePolicy:
    """Execution quality, in metres; never a random architectural parameter."""
    chord_error: float = .025
    min_segments: int = 12
    max_segments: int = 512

    def __post_init__(self):
        if (not math.isfinite(self.chord_error) or self.chord_error <= 0 or
                type(self.min_segments) is not int or type(self.max_segments) is not int or
                self.min_segments < 4 or self.min_segments % 2 or self.max_segments % 2 or
                self.max_segments < self.min_segments):
            raise ValueError('Round profile quality requires positive error and ordered even segment limits')

    def segments(self, span, height):
        radius = max(span/2, height)
        if not math.isfinite(radius) or radius <= 0 or span <= 0 or height <= 0:
            raise PadlError('Round profile needs a positive finite span and height', code='invalid_geometry')
        # 1-cos(x) <= x*x/2: conservative sagitta bound without acos near 1.
        count = max(self.min_segments, 2*math.ceil(math.pi/2*math.sqrt(radius/(8*self.chord_error))))
        if count > self.max_segments:
            raise PadlError('Round profile exceeds the configured segment limit at the requested chord error', code='geometry_limit')
        return count

    def points(self, span, height):
        count = self.segments(span, height)
        return tuple(((1-math.cos(i*math.pi/count))/2, 0. if i in (0,count) else math.sin(i*math.pi/count))
                     for i in range(count+1))


def ridge_vector(value):
    if not isinstance(value, (list, tuple)) or len(value) not in (2, 3):
        raise PadlError('roof.ridge_vector requires finite nonzero XY', code='invalid_input')
    try:
        valid = all(type(v) in (int, float) and math.isfinite(v) for v in value)
    except OverflowError:
        valid = False
    if not valid or (len(value) == 3 and value[2] != 0):
        raise PadlError('roof.ridge_vector requires finite nonzero XY and Z=0', code='invalid_input')
    length = math.hypot(*value[:2])
    if not math.isfinite(length) or length == 0:
        raise PadlError('roof.ridge_vector requires finite nonzero XY', code='invalid_input')
    return value[0]/length, value[1]/length


class RoofInputResolver:
    def resolve(self, kind, authored, inputs):
        options = dict(authored)
        if kind not in ('skillion', 'gabled'):
            return options
        mode = inputs.get('roof_direction_mode', 0)
        if type(mode) is not int or mode not in DIRECTION_MODES:
            raise PadlError(f'Unknown roof_direction_mode {mode!r}', code='invalid_category')
        if mode in (1, 2):
            options.update(reference_source='longest_edge', alignment='along' if mode == 1 else 'perpendicular')
        elif mode == 3:
            options.update(reference_source='vector2d', alignment='along',
                           ridge_vector=ridge_vector(inputs.get('roof_ridge_vector')))
        if kind == 'gabled':
            profile_mode = inputs.get('roof_profile_mode', 0)
            if type(profile_mode) is not int or profile_mode < 0:
                raise PadlError('Invalid roof_profile_mode', code='invalid_category')
            points = PROFILES.get(profile_mode).points if profile_mode else options.get('roof_profile', PROFILES.get(1).points)
            options['roof_profile'] = normalized_profile(points)
            options['roof_profile_mode'] = profile_mode
        return options


def direction(footprint, options):
    """Return the section/uphill axis, ridge cross +Z, in footprint coordinates."""
    reference = options.get('reference_source', 'longest_edge')
    alignment = options.get('alignment', 'along')
    orientation = options.get('orientation', 'ccw')
    if alignment not in ('along', 'perpendicular') or reference not in ('longest_edge', 'edge_index', 'vector2d'):
        raise PadlError('Unsupported roof direction configuration', code='unsupported_feature')
    if reference == 'vector2d':
        dx, dy = ridge_vector(options.get('ridge_vector'))
    else:
        if orientation not in ('ccw', 'cw'):
            raise PadlError('Unsupported roof orientation', code='invalid_input')
        ring = footprint.outer
        lengths = [math.dist(a[:2], b[:2]) for a, b in zip(ring, ring[1:] + ring[:1])]
        if reference == 'longest_edge':
            maximum = max(lengths)
            index = next(i for i, length in enumerate(lengths) if length >= maximum-.0001)
        else:
            index = options.get('edge_index', 0)
            if type(index) not in (int, float) or not math.isfinite(index) or int(index) != index:
                raise PadlError('roof.edge_index must be integral', code='invalid_input')
            index = int(index) % len(ring)
        a, b = ring[index], ring[(index+1) % len(ring)]
        dx, dy = ridge_vector((b[0]-a[0], b[1]-a[1]))
        if orientation == 'cw':
            dx, dy = -dx, -dy
    if alignment == 'perpendicular':
        dx, dy = -dy, dx
    return dy, -dx
