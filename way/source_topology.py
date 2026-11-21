"""Keep OSM endpoint identity separate from geometric coincidence.

Most preparation code uses frozen XY vectors as graph keys. Distinct OSM nodes
at one XY need distinct keys without moving the geometry. Only those collisions
use SourcePoint; all arithmetic deliberately returns ordinary geometry vectors.
Shared OSM node IDs remain connected, including at bridge/tunnel transitions.
"""
from collections import defaultdict

from mathutils import Vector


class SourcePoint(Vector):
    """An immutable-position graph key with an explicit source identity."""
    def __new__(cls, position, identity):
        point = super().__new__(cls, position)
        point._source_identity = identity
        return point.freeze()

    def __eq__(self, other):
        identity = getattr(self, '_source_identity', None)
        if identity is None:
            if isinstance(other, SourcePoint) and getattr(other, '_source_identity', None) is not None:
                return False
            return Vector.__eq__(self, other)
        return (isinstance(other, SourcePoint)
                and identity == getattr(other, '_source_identity', None)
                and Vector.__eq__(self, other))

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        identity = getattr(self, '_source_identity', None)
        return (Vector.__hash__(self) if identity is None
                else hash((tuple(self), identity)))

    def __add__(self, other):
        return Vector(self) + Vector(other)

    def __radd__(self, other):
        return Vector(other) + Vector(self)

    def __sub__(self, other):
        return Vector(self) - Vector(other)

    def __rsub__(self, other):
        return Vector(other) - Vector(self)

    def __mul__(self, value):
        return Vector(self) * value

    def __rmul__(self, value):
        return value * Vector(self)

    def __truediv__(self, value):
        return Vector(self) / value

    def __neg__(self):
        return -Vector(self)

    def copy(self):
        identity = getattr(self, '_source_identity', None)
        return Vector(self) if identity is None else SourcePoint(self, identity)


def graph_point(value):
    """Do not discard identity when a graph API receives an existing key."""
    return (value.freeze() if isinstance(value, SourcePoint)
            and getattr(value, '_source_identity', None) is not None else Vector(value).freeze())


class SourceEndpoints:
    """Resolve collected endpoint incidences without geometric welding."""
    def __init__(self):
        self.identities = defaultdict(set)
        self.positions = {}
        self.points = {}

    def add(self, position, identity):
        position = tuple(position)
        if identity[0] == 'osm_node':
            previous = self.positions.setdefault(identity, position)
            if previous != position:
                raise ValueError('One OSM node has conflicting projected coordinates')
        self.identities[position].add(identity)

    def point(self, position, identity):
        position = tuple(position)
        key = position, identity
        if key not in self.points:
            self.points[key] = (SourcePoint(position, identity)
                if len(self.identities[position]) > 1 else Vector(position).freeze())
        return self.points[key]

    @staticmethod
    def identity(segment, position, way_id, segment_index, end):
        """Clipped endpoints identify a source edge cut, never a new junction."""
        position = tuple(position)
        for node, original in ((segment.id1, segment.v1), (segment.id2, segment.v2)):
            if position == tuple(original) and node is not None:
                return 'osm_node', str(node)
        nodes = segment.id1, segment.id2
        if all(node is not None for node in nodes):
            return 'source_cut', tuple(sorted(str(node) for node in nodes)), position
        return 'unidentified_end', str(way_id), segment_index, end
