"""Repeatable source edits for native real-network qualification/measurement."""
import json
import math
import bpy
from mathutils import Vector

from blosm.way.street_network.python_backend.model import Network, End


class SectionProbeEdit:
    """Extend one unattached terminal end and change width in physical metres."""
    def __init__(self, controller, section_id=None, *, distance=.05, width=.05):
        self.network = Network.from_dict(json.loads(controller['st_py_source_data']))
        bundled = {sid for b in self.network.bundles for sid in b.sections}
        attached = {e for j in (*self.network.intersections, *self.network.connections) for e in j.ends}
        candidates = [s for s in self.network.sections if s.id not in bundled and
                      sum(End(s.id, end) in attached for end in (0, 1)) == 1]
        if section_id is not None:
            candidates = [s for s in candidates if s.id == section_id]
        if not candidates:
            raise ValueError('Probe needs an unbundled Section with one unattached terminal end')
        self.section = max(candidates, key=lambda s: sum(math.dist(a, b) for a, b in zip(s.points, s.points[1:])))
        section = self.section
        self.free = 0 if End(section.id, 0) not in attached else 1
        self.terminal_kind = 'SELECTION_CUT' if self.free in section.selection_ends else 'SOURCE_ENDPOINT'
        self.source = next(o for o in controller['st_py_sources'].objects if
            section.id in {v.value for v in o.data.attributes['st_section_id'].data})
        assert self.source.mode == 'OBJECT'
        ids = self.source.data.attributes['st_section_id']
        degrees = {}
        for edge, identity in zip(self.source.data.edges, ids.data):
            if identity.value == section.id:
                for vertex in edge.vertices:
                    degrees[vertex] = degrees.get(vertex, 0) + 1
        scale = bpy.context.scene.unit_settings.scale_length
        inverse = self.source.matrix_local.inverted()
        expected = inverse @ (Vector(section.points[0 if self.free == 0 else -1]) / scale)
        self.vertex = min((self.source.data.vertices[i] for i, degree in degrees.items() if degree == 1),
                          key=lambda v: (v.co - expected).length)
        assert (self.vertex.co - expected).length * scale < 1e-4
        self.position = self.vertex.co.copy()
        a, b = (section.points[0], section.points[1]) if self.free == 0 else (section.points[-1], section.points[-2])
        tangent = Vector(a) - Vector(b)
        tangent.normalize()
        self.delta = inverse.to_3x3() @ (distance * tangent / scale)
        self.width = width
        self.widths = self.source.data.attributes['st_width']
        self.original_widths = [(i, v.value) for i, v in enumerate(self.widths.data) if ids.data[i].value == section.id]

    def apply(self):
        self.vertex.co = self.position + self.delta
        for i, value in self.original_widths:
            self.widths.data[i].value = value + self.width
        self.source.data.update()

    def restore(self):
        self.vertex.co = self.position
        for i, value in self.original_widths:
            self.widths.data[i].value = value
        self.source.data.update()


class SectionSeamProbeEdit:
    """Move all authored copies of a real Section seam and widen both sides.

    Identity comes from Connection incidences and edge Section IDs, rather than
    source-object grouping or generated mesh indices. No topology is changed.
    """
    def __init__(self, controller, connection_id=None, *, distance=.05, width=.05):
        self.network = Network.from_dict(json.loads(controller['st_py_source_data']))
        sections = {s.id: s for s in self.network.sections}
        bundled = {sid for b in self.network.bundles for sid in b.sections}
        candidates = [c for c in self.network.connections if
                      all(e.section_id not in bundled for e in c.ends)]
        if connection_id is not None:
            candidates = [c for c in candidates if c.id == connection_id]
        if not candidates:
            raise ValueError('Seam probe needs an explicit unbundled Section connection')

        def reach(connection):
            return min(math.dist(sections[e.section_id].points[0 if e.end == 0 else -1],
                                 sections[e.section_id].points[1 if e.end == 0 else -2])
                       for e in connection.ends)

        self.connection = max(candidates, key=lambda c: (reach(c), -c.id))
        self.sections = tuple(sections[e.section_id] for e in self.connection.ends)
        self.section, self.free, self.terminal_kind = self.sections[0], None, None
        first = self.connection.ends[0]
        points = self.section.points
        tangent = Vector(points[1 if first.end == 0 else -2]) - Vector(points[0 if first.end == 0 else -1])
        tangent.normalize()
        displacement = distance * Vector((-tangent.y, tangent.x, 0.))
        scale = bpy.context.scene.unit_settings.scale_length
        self.positions, self.widths, self.sources = {}, {}, {}
        self.width = width
        for endpoint in self.connection.ends:
            section = sections[endpoint.section_id]
            source = next(o for o in controller['st_py_sources'].objects if
                          section.id in {v.value for v in o.data.attributes['st_section_id'].data})
            assert source.mode == 'OBJECT'
            identities = source.data.attributes['st_section_id']
            degrees = {}
            for edge, identity in zip(source.data.edges, identities.data):
                if identity.value == section.id:
                    for vertex in edge.vertices:
                        degrees[vertex] = degrees.get(vertex, 0) + 1
            inverse = source.matrix_local.inverted()
            expected = inverse @ (Vector(section.points[0 if endpoint.end == 0 else -1]) / scale)
            vertex = min((source.data.vertices[i] for i, degree in degrees.items() if degree == 1),
                         key=lambda v: (v.co - expected).length)
            assert (vertex.co - expected).length * scale < 1e-4
            address = (source.as_pointer(), vertex.index)
            self.positions[address] = (vertex, vertex.co.copy(), inverse.to_3x3() @ (displacement / scale))
            values = source.data.attributes['st_width']
            for i, identity in enumerate(identities.data):
                if identity.value == section.id:
                    self.widths[source.as_pointer(), i] = (values.data[i], values.data[i].value)
            self.sources[source.as_pointer()] = source

    def apply(self):
        for vertex, position, delta in self.positions.values():
            vertex.co = position + delta
        for value, original in self.widths.values():
            value.value = original + self.width
        for source in self.sources.values():
            source.data.update()

    def restore(self):
        for vertex, position, delta in self.positions.values():
            vertex.co = position
        for value, original in self.widths.values():
            value.value = original
        for source in self.sources.values():
            source.data.update()
