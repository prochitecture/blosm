"""Backend-independent, persistent identities and authored physical inputs."""
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math


def fingerprint(value):
    """A versioned source digest, independent of dict and object iteration order."""
    if hasattr(value, '__dataclass_fields__'):
        value = asdict(value)
    raw = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(('street-python-1:' + raw).encode('utf8')).hexdigest()


@dataclass(frozen=True, order=True)
class Key:
    kind: str
    id: int

    def __post_init__(self):
        if not self.kind or ':' in self.kind or type(self.id) is not int or self.id < 0:
            raise ValueError('An entity needs a kind and nonnegative persistent ID')

    def text(self):
        return f'{self.kind}:{self.id}'

    @classmethod
    def parse(cls, text):
        kind, identity = text.split(':')
        return cls(kind, int(identity))


@dataclass(frozen=True)
class Section:
    id: int
    street_id: int
    points: tuple
    width: float = 6.
    offset: float = 0.
    rise: float = .08
    interpolation: int = 3
    datum: float = 0.
    category: str = 'residential'
    selection_ends: tuple = ()

    def __post_init__(self):
        Key('section', self.id)
        Key('street', self.street_id)
        points = tuple(tuple(float(v) for v in p) for p in self.points)
        if len(points) < 2 or any(len(p) != 3 or not all(map(math.isfinite, p)) for p in points):
            raise ValueError('A Section needs at least two finite XYZ controls')
        if any(math.dist(a, b) <= 1e-7 for a, b in zip(points, points[1:])):
            raise ValueError('A Section has a zero-length segment')
        if (not all(math.isfinite(v) for v in (self.width, self.offset, self.rise, self.datum))
                or self.width <= 0 or self.rise < 0 or self.interpolation not in (0, 3)):
            raise ValueError('Invalid Section dimensions, profile or interpolation')
        object.__setattr__(self, 'points', points)
        if (len(set(self.selection_ends)) != len(self.selection_ends) or
                any(type(end) is not int or end not in (0, 1) for end in self.selection_ends)):
            raise ValueError('Selection cuts must name distinct directed Section ends')
        object.__setattr__(self, 'selection_ends', tuple(self.selection_ends))

    @property
    def key(self):
        return Key('section', self.id)


@dataclass(frozen=True, order=True)
class End:
    section_id: int
    end: int

    def __post_init__(self):
        Key('section', self.section_id)
        if type(self.end) is not int or self.end not in (0, 1):
            raise ValueError('An endpoint is the directed Section start or end')


@dataclass(frozen=True)
class Intersection:
    id: int
    ends: tuple
    radius: float = 2.
    radii: tuple = ()
    radius_modes: tuple = ()

    def __post_init__(self):
        Key('intersection', self.id)
        if len(self.ends) < 2 or len(set(self.ends)) != len(self.ends):
            raise ValueError('An intersection needs distinct endpoint incidences')
        if not math.isfinite(self.radius) or self.radius <= 0:
            raise ValueError('Corner radius must be positive')
        if self.radii and (len(self.radii) != len(self.ends) or
                           any(not math.isfinite(r) or r <= 0 for r in self.radii)):
            raise ValueError('Per-approach corner radii must match endpoint ownership')
        if self.radius_modes and (len(self.radius_modes) != len(self.ends) or
                                  any(mode not in ('STRICT', 'STYLE') for mode in self.radius_modes)):
            raise ValueError('Corner radius modes must be STRICT or STYLE for each endpoint')

    def corner_radius(self, end):
        """The outgoing counterclockwise corner belongs to this incidence."""
        index = self.ends.index(end)
        return (self.radii[index] if self.radii else self.radius,
                self.radius_modes[index] if self.radius_modes else 'STRICT')

    def with_corner_radius(self, end, radius, *, mode='STRICT'):
        """An explicit override is strict unless fitting is expressly selected."""
        values = [self.corner_radius(e) for e in self.ends]
        values[self.ends.index(end)] = (radius, mode)
        return replace(self, radii=tuple(v[0] for v in values), radius_modes=tuple(v[1] for v in values))

    @property
    def key(self):
        return Key('intersection', self.id)

    @property
    def layout_key(self):
        return Key('mouths', self.id)

    @property
    def mesh_key(self):
        return Key('junction_mesh', self.id)


@dataclass(frozen=True)
class Connection:
    """Explicit Section incidence; differing dimensions or frames need a join."""
    id: int
    ends: tuple

    def __post_init__(self):
        Key('connection', self.id)
        if len(self.ends) != 2 or len(set(self.ends)) != 2:
            raise ValueError('A Section seam needs two distinct endpoint incidences')

    @property
    def key(self):
        return Key('connection', self.id)

    @property
    def layout_key(self):
        return Key('seam_layout', self.id)

    @property
    def mesh_key(self):
        return None


@dataclass(frozen=True)
class Traversal:
    section_id: int
    forward: bool = True

    def __post_init__(self):
        Key('section', self.section_id)
        if type(self.forward) is not bool:
            raise ValueError('Section traversal orientation must be explicit')


@dataclass(frozen=True)
class BundleMember:
    path: tuple
    complete: bool = True

    def __post_init__(self):
        if not self.path or len({t.section_id for t in self.path}) != len(self.path):
            raise ValueError('Bundle member requires a distinct ordered Section traversal')


@dataclass(frozen=True)
class Bundle:
    id: int
    members: tuple
    interior_section_ids: tuple = ()
    profile_rise: float | None = None

    def __post_init__(self):
        Key('bundle', self.id)
        if len(self.members) < 2 or len(set(self.sections)) != len(self.sections):
            raise ValueError('A Bundle needs distinct declared member and interior Sections')
        if self.profile_rise is not None and (not math.isfinite(self.profile_rise) or self.profile_rise < 0):
            raise ValueError('Bundle profile rise must be finite and nonnegative')

    @property
    def sections(self):
        return tuple(t.section_id for m in self.members for t in m.path) + self.interior_section_ids

    @property
    def key(self):
        return Key('bundle', self.id)


@dataclass(frozen=True)
class Placement:
    """Source-space anchor, never a generated vertex or face reference.

    offset is measured from the selected physical carriageway edge; positive
    means farther outside. Station is metres along the directed Section.
    """
    id: int
    section_id: int
    station: float
    side: str = 'RIGHT'
    offset: float = .6
    height: float = 0.
    asset: str = 'placeholder'

    def __post_init__(self):
        Key('furniture', self.id)
        Key('section', self.section_id)
        if (self.side not in ('LEFT', 'RIGHT') or self.station < 0
                or not all(math.isfinite(v) for v in (self.station, self.offset, self.height))):
            raise ValueError('Invalid furniture anchor')

    @property
    def key(self):
        return Key('furniture', self.id)


@dataclass(frozen=True)
class Settings:
    cell_size: float = .7
    resolution: int = 1
    exaggeration: float = 10.
    curve_tolerance: float = .005
    pattern_cache_bytes: int = 8 * 1024 * 1024
    continuation_angle_degrees: float = 4.5

    def __post_init__(self):
        if (not all(math.isfinite(v) for v in (self.cell_size, self.exaggeration, self.curve_tolerance, self.continuation_angle_degrees))
                or self.cell_size <= 0 or self.curve_tolerance <= 0 or self.exaggeration < 0
                or not 0. <= self.continuation_angle_degrees <= 15.
                or type(self.resolution) is not int or not 1 <= self.resolution <= 4
                or type(self.pattern_cache_bytes) is not int or self.pattern_cache_bytes < 0):
            raise ValueError('Invalid generation settings')

    @property
    def spacing(self):
        return self.cell_size / self.resolution


@dataclass(frozen=True)
class Network:
    uuid: str
    sections: tuple
    intersections: tuple = ()
    bundles: tuple = ()
    furniture: tuple = ()
    settings: Settings = Settings()
    connections: tuple = ()

    def __post_init__(self):
        if not self.uuid:
            raise ValueError('A network needs a persistent UUID')
        index = self.section_index
        attached = set()
        for collection in (self.sections, self.intersections, self.connections, self.bundles, self.furniture):
            if len({v.key for v in collection}) != len(collection):
                raise ValueError('Duplicate persistent entity IDs')
        for junction in (*self.intersections, *self.connections):
            for end in junction.ends:
                if end.section_id not in index or end in attached:
                    raise ValueError('Missing or multiply attached Section endpoint')
                if end.end in index[end.section_id].selection_ends:
                    raise ValueError('A selection cut cannot also own an explicit junction incidence')
                attached.add(end)
        for bundle in self.bundles:
            if any(sid not in index for sid in bundle.sections):
                raise ValueError('Bundle references an absent Section')
        if any(item.section_id not in index for item in self.furniture):
            raise ValueError('Furniture references an absent Section')

    @property
    def section_index(self):
        return {s.id: s for s in self.sections}

    def records(self):
        return {item.key: item for group in (self.sections, self.intersections, self.connections, self.bundles, self.furniture)
                for item in group}

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        return cls(data['uuid'], tuple(Section(**s) for s in data['sections']),
            tuple(Intersection(j['id'], tuple(End(**e) for e in j['ends']), j['radius'],
                               tuple(j.get('radii', ())), tuple(j.get('radius_modes', ())))
                  for j in data.get('intersections', ())),
            tuple(Bundle(b['id'], tuple(BundleMember(tuple(Traversal(**t) for t in m['path']),
                  m.get('complete', True)) for m in b['members']), tuple(b.get('interior_section_ids', ())), b.get('profile_rise'))
                  for b in data.get('bundles', ())),
            tuple(Placement(**p) for p in data.get('furniture', ())), Settings(**data.get('settings', {})),
            tuple(Connection(c['id'], tuple(End(**e) for e in c['ends'])) for c in data.get('connections', ())))
