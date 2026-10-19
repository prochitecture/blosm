"""Declare mouth-free carriageway spans without relabeling category endpoints.

This is topology and dimensional admission only. Geometry, interpolation and
live endpoint checks belong to GN. Excluded Sections remain in the complete
editable source Street; no replacement Street or Section identity is allocated.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class BoundaryRun:
    span_identity: tuple
    roles: tuple
    owners: tuple
    next_sections: tuple


class OriginalBoundaryRuns:
    def prepare(self, snapshot, ports, network, *, dimensions=False):
        runs, issues = [], []
        for route in network.routes:
            if route.kind != 'BOUNDARY_ONLY':
                continue
            span = route.span
            try:
                if not span.in_scope or route.start_port is not None or route.end_port is not None:
                    raise ValueError('A boundary run cannot own a junction mouth')
                source = snapshot.streets[span.street_id]
                if source.get('bundle_id') is not None:
                    raise ValueError('Boundary-bounded Bundle members need the Bundle surface branch')
                sections = [s for s in source['items'] if s['kind'] == 'Section' and s['id'] in span.section_ids]
                if tuple(s['id'] for s in sections) != span.section_ids:
                    raise ValueError('Boundary run lost its ordered original Sections')
                sizes = [(float(s['width']), float(s.get('offset') or 0.)) for s in sections]
                if not sizes or any(not math.isfinite(w) or w <= 0 or not math.isfinite(o) for w, o in sizes):
                    raise ValueError('Boundary run has invalid physical dimensions')
                if dimensions:
                    from .physical_dimensions import PreparedPhysicalDimensions
                    PreparedPhysicalDimensions().read_spans(snapshot,(span,),source_plan=ports)
                elif any(d != sizes[0] for d in sizes):
                    raise ValueError('Boundary run dimensional transitions need the tapered branch')
                roles, owners, peers = [], [], []
                for end in (span.start, span.end):
                    if end.junction_id is not None:
                        raise ValueError('Boundary run still has an original junction owner')
                    if end.owner_kind == 'FREE':
                        roles.append(1); owners.append(-1); peers.append(-1)
                    elif end.owner_kind == 'SELECTION_BOUNDARY':
                        ports.selection_boundary(end.endpoint)
                        roles.append(6); owners.append(end.endpoint.section_id); peers.append(-1)
                    elif end.owner_kind == 'CATEGORY_BOUNDARY':
                        boundary = ports.category_boundary(end.endpoint)
                        peer, = [e for e in boundary.endpoints if e != end.endpoint]
                        cross = boundary.source_kind == 'OriginalStreetSeam'
                        roles.append(5 if cross else 4)
                        owners.append(boundary.source_id if cross else span.street_id)
                        peers.append(peer.section_id)
                    else:
                        raise ValueError('Boundary run has an unresolved original endpoint')
                if all(role == 1 for role in roles):
                    raise ValueError('A genuine standalone Street uses its own declaration')
                runs.append(BoundaryRun(span.identity, tuple(roles), tuple(owners), tuple(peers)))
            except ValueError as exc:
                issues.append(dict(kind='BOUNDARY_RUN', id=span.street_id,
                                   section_ids=list(span.section_ids), reason=str(exc)))
        return tuple(runs), tuple(issues)
