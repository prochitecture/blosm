"""Address canonical Section polylines through whole-Street solver views.

Editable Sections keep their separate endpoint roles. The private solver view
shares consecutive Section ends by explicit identity, never by proximity.
Object grouping is deliberately absent from these disposable lookup indices.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PortControlReference:
    vertex: int
    alias: int = -1


@dataclass(frozen=True)
class PortEdgeReference:
    edge: int
    first: int
    last: int
    section_id: int


@dataclass(frozen=True)
class PortSourceLayout:
    section_ids: tuple
    source_vertex_count: int
    source_edge_count: int
    controls: tuple[PortControlReference, ...]
    edges: tuple[PortEdgeReference, ...]

    @classmethod
    def prepare(cls, snapshot, plan):
        """Bind already validated original spans to canonical Section order.

        Stable live indices follow plan.source_streets, then original Section
        and authored control order. Managed regrouping can keep those indices;
        changing topology requires a new plan and layout. No coordinates are
        fitted or generated here, and aliases must still agree during GN edits.
        """
        controls, edges, sections = [], [], []
        vertex = edge = 0
        if not plan.source_streets or len(set(plan.source_streets)) != len(plan.source_streets):
            raise ValueError('Editable port layout needs unique original Streets')
        for sid in plan.source_streets:
            source = snapshot.streets.get(sid)
            if source is None or not plan.source_supported(source):
                raise ValueError('Editable port layout needs supported original Streets')
            items = [s for s in source['items'] if s['kind'] == 'Section']
            spans = [s for s in plan.spans if s.street_id == sid]
            if not items or tuple(s['id'] for s in items) != tuple(sec for span in spans for sec in span.section_ids):
                raise ValueError('Original Section order changed; rebuild the port plan')
            expected = sum(len(s['points']) - 1 for s in items) + 1
            if not spans or any(s.source_control_count != expected for s in spans):
                raise ValueError('Original control topology changed; rebuild the port plan')
            ranges = {};at = 0
            for item in items:
                end = at + len(item['points']) - 1
                ranges[item['id']] = (at, end);at = end
            if any((span.first_control, span.last_control) !=
                   (ranges[span.section_ids[0]][0], ranges[span.section_ids[-1]][1]) for span in spans):
                raise ValueError('Original Section seam changed; rebuild the port plan')
            for order, section in enumerate(items):
                count = len(section['points']);section_id = section['id']
                if count < 2 or section_id in sections:
                    raise ValueError('Editable Sections need distinct IDs and at least two controls')
                sections.append(section_id)
                if order:
                    previous = controls[-1]
                    controls[-1] = PortControlReference(previous.vertex, vertex)
                controls.extend(PortControlReference(vertex + i) for i in range(1 if order else 0, count))
                edges.extend(PortEdgeReference(edge + i, vertex + i, vertex + i + 1, section_id)
                             for i in range(count - 1))
                vertex += count;edge += count - 1
        return cls(tuple(sections), vertex, edge, tuple(controls), tuple(edges))
