"""Extensible OSM -> logical PADL inputs, independent of Blender and backends.

Rules are ordered and replaceable. Values remain sparse: FootprintBindings is
responsible for writing the corresponding per-face presence flags.
"""
from dataclasses import dataclass, field
import math
import re


@dataclass(frozen=True)
class MappingIssue:
    code: str
    message: str
    fatal: bool = False


@dataclass
class MappingResult:
    attributes: dict = field(default_factory=lambda: {'building_use': 'unknown', 'min_elevation': 0.0})
    issues: list = field(default_factory=list)

    @property
    def accepted(self):
        return not any(issue.fatal for issue in self.issues)


def length_metres(value):
    """OSM length syntax: metres, explicit m, feet, or feet and inches.

    Reject ranges, lists, approximate values and nonfinite values, rather than
    setting an authoritative height to a guessed or partially parsed value.
    """
    text = str(value).strip().lower()
    number = r'(?:\d+(?:\.\d*)?|\.\d+)'
    metric = re.fullmatch(r'(' + number + r')\s*(m|ft|feet|foot)?', text)
    imperial = re.fullmatch(r'(' + number + r")\s*'\s*(?:([0-9]+(?:\.[0-9]*)?)\s*\")?", text)
    if metric:
        result = float(metric[1]) * (0.3048 if metric[2] in ('ft', 'feet', 'foot') else 1.0)
    elif imperial:
        inches = float(imperial[2] or 0)
        if inches >= 12:
            raise ValueError('Inches must be less than 12')
        result = float(imperial[1]) * 0.3048 + inches * 0.0254
    else:
        raise ValueError('Expected a length in metres or feet/inches')
    if not math.isfinite(result) or result > 100000:
        raise ValueError('Length is outside the supported range')
    return result


class TagRule:
    def apply_to(self, tags, result, footprint=None):
        self.apply(tags, result)

    def apply(self, tags, result):
        raise NotImplementedError


class BuildingUseRule(TagRule):
    """Current explicit use wins over form; ambiguous residential stays unknown."""
    uses = {
        'apartments': 'apartments', 'house': 'house', 'detached': 'house',
        'semidetached_house': 'house', 'terrace': 'house', 'bungalow': 'house',
        'office': 'office',
    }

    def __init__(self, uses=None):
        self.uses = dict(self.uses if uses is None else uses)

    def apply(self, tags, result):
        if tags.get('building:use'):
            value = self.uses.get(tags['building:use'], 'unknown')
        else:
            form = tags.get('building', tags.get('building:part', 'yes'))
            value = self.uses.get(form, 'unknown')
            # An office POI inside an apartment block doesn't reclassify it.
            if form == 'yes' and tags.get('office') not in (None, '', 'no', 'vacant'):
                value = 'office'
            # These tags identify functions outside our three initial uses.
            if tags.get('amenity') not in (None, '', 'no') and form == 'yes':
                value = 'unknown'
        result.attributes['building_use'] = value


class RoofShapeRule(TagRule):
    def __init__(self, shapes=None):
        self.shapes = (dict(shapes) if shapes is not None else
                       {**{name: name for name in ('flat', 'skillion', 'gabled', 'hipped')},
                        **{name: 'gabled' for name in ('gambrel', 'round', 'saltbox')}})
        self.profiles = {'gabled': 1, 'gambrel': 2, 'round': 3, 'saltbox': 4}

    def apply(self, tags, result):
        shape = tags.get('roof:shape')
        if shape:
            if shape in self.shapes:
                result.attributes['roof_type'] = self.shapes[shape]
                if self.shapes[shape] == 'gabled' and shape in self.profiles:
                    result.attributes['roof_profile_mode'] = self.profiles[shape]
            else:
                result.issues.append(MappingIssue('unsupported_roof', f'Unsupported roof:shape={shape}', True))


def bearing_degrees(raw):
    text = str(raw).strip().upper()
    compass = ('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
               'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW')
    if text in compass:
        return compass.index(text)*22.5
    if not re.fullmatch(r'(?:\d+(?:\.\d*)?|\.\d+)', text):
        raise ValueError('Expected a compass point or bearing from 0 to 360 degrees')
    value = float(text)
    if not math.isfinite(value) or not 0 <= value <= 360:
        raise ValueError('Bearing must be from 0 to 360 degrees')
    return value % 360


class ProjectedRoofBearing:
    """Convert geographical downhill bearings to the import projection's XY."""
    def __init__(self, projection):
        self.projection = projection

    def __call__(self, degrees, footprint=None):
        theta = math.radians(degrees)
        downhill = math.sin(theta), math.cos(theta)
        if self.projection is not None and footprint is not None:
            # Local finite difference handles meridian convergence away from
            # the projection origin. Input owners are created in this frame.
            ring = footprint.outer
            x, y = (sum(p[i] for p in ring)/len(ring) for i in (0,1))
            lat, lon = self.projection.toGeographic(x, y)
            delta = math.degrees(10/6378137.)
            a = self.projection.fromGeographic(lat, lon)
            b = self.projection.fromGeographic(lat+delta*downhill[1], lon+delta*downhill[0]/math.cos(math.radians(lat)))
            downhill = b[0]-a[0], b[1]-a[1]
        length = math.hypot(*downhill)
        if not math.isfinite(length) or length == 0:
            raise ValueError('Cannot project the roof bearing')
        # Ridge cross +Z is uphill, opposite the OSM downhill direction.
        return downhill[1]/length, -downhill[0]/length, 0.


class RoofDirectionRule(TagRule):
    def __init__(self, bearing_transform=None):
        self.bearing_transform = bearing_transform or ProjectedRoofBearing(None)

    def apply(self, tags, result):
        self.apply_to(tags, result)

    def apply_to(self, tags, result, footprint=None):
        if result.attributes.get('roof_type') in ('flat', 'hipped'):
            return
        orientation = tags.get('roof:orientation')
        if orientation is not None:
            modes = {'along': 1, 'across': 2}
            if orientation in modes:
                result.attributes['roof_direction_mode'] = modes[orientation]
            else:
                result.issues.append(MappingIssue('invalid_tag', f'roof:orientation={orientation!r}: expected along or across; input left absent'))
        if 'roof:direction' in tags:
            try:
                bearing = bearing_degrees(tags['roof:direction'])
                vector = self.bearing_transform(bearing, footprint)
                result.attributes.update(roof_direction_mode=3, roof_ridge_vector=vector)
                if orientation in ('along', 'across'):
                    result.issues.append(MappingIssue('roof_direction_priority', 'roof:direction takes precedence over roof:orientation'))
            except (ValueError, OverflowError, ZeroDivisionError) as exc:
                result.issues.append(MappingIssue('invalid_tag', f'roof:direction={tags["roof:direction"]!r}: {exc}; using valid roof:orientation or PADL settings'))
        if result.attributes.get('roof_type') == 'skillion' and result.attributes.get('roof_direction_mode') in (1,2):
            result.issues.append(MappingIssue('incomplete_roof_direction', 'roof:orientation supplies a skillion axis but no downhill side; using PADL orientation'))


class LengthTagRule(TagRule):
    def __init__(self, tag, attribute, *, allow_zero=True):
        self.tag, self.attribute, self.allow_zero = tag, attribute, allow_zero

    def apply(self, tags, result):
        if self.tag not in tags:
            return
        if self.tag == 'roof:height' and result.attributes.get('roof_type') == 'flat':
            return
        try:
            value = length_metres(tags[self.tag])
            if not self.allow_zero and value == 0:
                raise ValueError('Height must be positive')
            result.attributes[self.attribute] = value
        except ValueError as exc:
            result.issues.append(MappingIssue('invalid_tag', f'{self.tag}={tags[self.tag]!r}: {exc}; input left absent'))


class WallLevelsRule(TagRule):
    def apply(self, tags, result):
        if 'building:levels' in tags:
            raw = tags['building:levels']
            if re.fullmatch(r'[0-9]+', raw.strip()) and 0 < int(raw) < 2**31:
                result.attributes['levels_wall_count'] = int(raw)
            else:
                result.issues.append(MappingIssue('invalid_tag', f'building:levels={raw!r} is not a positive integer; input left absent'))
        # Skipped storeys cannot be represented by merely clipping at min_height:
        # the base and exposed floor ranges need a dedicated future adapter.
        if tags.get('building:min_level', '0').strip() not in ('', '0'):
            result.issues.append(MappingIssue('unsupported_min_level', 'Nonzero building:min_level requires a building-part/raised-floor adapter', True))


class OsmFaceAttributeMapper:
    def __init__(self, rules=None, *, projection=None):
        self.rules = tuple(rules) if rules is not None else (
            BuildingUseRule(), RoofShapeRule(), RoofDirectionRule(ProjectedRoofBearing(projection)),
            LengthTagRule('height', 'max_elevation', allow_zero=False),
            LengthTagRule('min_height', 'min_elevation'),
            LengthTagRule('roof:height', 'roof_height'), WallLevelsRule(),
        )

    def map(self, tags, *, footprint=None):
        result = MappingResult()
        for rule in self.rules:
            rule.apply_to(tags, result, footprint)
        if result.attributes.get('roof_type') == 'flat':
            result.attributes.pop('roof_height', None)  # OSM roof height is immaterial for flat roofs.
        return result
