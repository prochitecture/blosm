"""Declared versus allocated height resolution for the initial wall-level profile."""
from dataclasses import dataclass, field
from typing import Mapping

from .evaluation import is_number
from .model import Value


def close(a, b):
    return abs(a - b) <= max(0.0001, 0.000001 * max(abs(a), abs(b)))


@dataclass(frozen=True)
class ResolvedVolume:
    operation_path: tuple
    operation_name: str
    footprint_identity: str
    roof_type: str
    base_elevation: float
    min_elevation: float
    max_elevation: float
    wall_min_height: float
    roof_height: float
    levels_wall_count: int | None = None
    level_height: float | None = None
    ground_floor_height: float | None = None
    top_height: float | None = None
    level_overrides: Mapping = field(default_factory=dict)
    roof_options: Mapping = field(default_factory=dict)

    @property
    def roof_min_elevation(self):
        return self.base_elevation + self.wall_min_height


class HeightResolver:
    def __init__(self, source):
        self.source = source

    def resolve(self, operation, resolver, context):
        spec = operation.params
        p = operation.path + ('params',)

        def numeric(value, minimum=None, positive=False, integer=False):
            result = resolver.numeric(resolver.resolve(value, context), value)
            if (minimum is not None and result < minimum) or (positive and result <= 0):
                raise self.source.error(value.path, 'Height/count is outside its permitted range', 'invalid_height')
            if integer and result != int(result):
                raise self.source.error(value.path, 'Count/index must be integral; use an explicit floor/ceil expression', 'invalid_count')
            return int(result) if integer else result

        def scalar(name, default=0):
            return numeric(spec.values.get(name, Value('literal', default, p + (name,))))

        def fit(target, total, policy, path):
            if policy == 'fixed':
                if not close(target, total):
                    raise self.source.error(path, f'Fixed height {target} does not match declared sum {total}', 'height_conflict')
                return 1
            if total == 0 and target == 0:
                return 1
            if total == 0 or not is_number(target / total) or target / total <= 0:
                raise self.source.error(path, 'Proportional fitting requires a finite positive factor', 'height_conflict')
            return target / total

        roof = {k: resolver.resolve(v, context) for k, v in spec.roof.items()}
        kind = roof.pop('type')
        if kind not in ('flat', 'skillion', 'gabled', 'hipped'):
            raise self.source.error(p + ('roof', 'type'), f'Unsupported roof {kind!r}', 'unsupported_roof')
        b, m = scalar('base_elevation'), scalar('min_elevation')
        if 'roof_height' not in spec.values and kind != 'flat':
            raise self.source.error(p, 'A non-flat roof requires roof_height', 'unresolved_height')
        r = scalar('roof_height')
        if r < 0 or (kind == 'flat' and r != 0):
            raise self.source.error(p + ('roof_height',), 'Roof height must be nonnegative and zero for a flat roof', 'invalid_height')
        n = h = g = t = None
        overrides = {}
        total = None
        if spec.wall is not None:
            wall = spec.wall
            n = numeric(wall.count, minimum=0, integer=True)
            h = numeric(wall.height, positive=True)
            t = numeric(wall.top_height, minimum=0)
            for index, value in wall.overrides:
                i = numeric(index, minimum=0, integer=True)
                if i >= n or i in overrides:
                    raise self.source.error(index.path, 'Duplicate or out-of-range wall level index', 'invalid_count')
                overrides[i] = numeric(value, positive=True)
            total = n * h + sum(v - h for v in overrides.values()) + t
        w = scalar('wall_min_height') if 'wall_min_height' in spec.values else total
        if not is_number(w) or w < 0:
            raise self.source.error(p + ('wall_min_height',), 'Wall range must be a finite nonnegative height', 'invalid_height')
        if total is not None:
            s = fit(w, total, spec.wall.height_type, spec.wall.path)
            h, t = h * s, t * s
            overrides = {i: v * s for i, v in overrides.items()}
            g = overrides.get(0, h) if n else None
        total_elevation = b + w + r
        if not is_number(total_elevation):
            raise self.source.error(p, 'Elevation sum is not finite', 'invalid_height')
        scale = fit(scalar('max_elevation'), total_elevation, spec.height_type, p) if 'max_elevation' in spec.values else 1
        b, m, w, r = (x * scale for x in (b, m, w, r))
        if total is not None:
            h, t = h * scale, t * scale
            g = g * scale if g is not None else None
            overrides = {i: v * scale for i, v in overrides.items()}
        maximum = b + w + r
        if not all(is_number(x) for x in (b, m, w, r, maximum)):
            raise self.source.error(p, 'Allocated dimensions must be finite', 'invalid_height')
        if m > b + w:
            raise self.source.error(p + ('min_elevation',), 'min_elevation above roof_min_elevation requires unsupported clipping', 'unsupported_clipping')
        if kind != 'flat' and r <= 0:
            raise self.source.error(p + ('roof_height',), 'The initial pitched-roof backend requires a positive rise', 'invalid_height')
        return ResolvedVolume(operation.path, operation.name, context.identity, kind, b, m, maximum, w, r,
                              n, h, g, t, overrides, roof)
