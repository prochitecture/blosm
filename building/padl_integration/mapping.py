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
                       {name: name for name in ('flat', 'skillion', 'gabled', 'hipped')})

    def apply(self, tags, result):
        shape = tags.get('roof:shape')
        if shape:
            if shape in self.shapes:
                result.attributes['roof_type'] = self.shapes[shape]
            else:
                result.issues.append(MappingIssue('unsupported_roof', f'Unsupported roof:shape={shape}', True))


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
    def __init__(self, rules=None):
        self.rules = tuple(rules) if rules is not None else (
            BuildingUseRule(), RoofShapeRule(),
            LengthTagRule('height', 'max_elevation', allow_zero=False),
            LengthTagRule('min_height', 'min_elevation'),
            LengthTagRule('roof:height', 'roof_height'), WallLevelsRule(),
        )

    def map(self, tags):
        result = MappingResult()
        for rule in self.rules:
            rule.apply(tags, result)
        if result.attributes.get('roof_type') == 'flat':
            result.attributes.pop('roof_height', None)  # OSM roof height is immaterial for flat roofs.
        return result
