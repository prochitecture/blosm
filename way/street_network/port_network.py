"""Once-owned surface routes over original Section spans.

Topology preparation does not fit curves or solve geometry. A short span may
belong to a declared compound; every other span retains its original owners.
Unselected junctions remain explicit pending connections, never free ends.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PortSpanRoute:
    span: object
    kind: str
    start_port: object | None
    end_port: object | None


@dataclass(frozen=True)
class PortNetworkPlan:
    ordinary_site_ids: tuple
    compounds: tuple
    routes: tuple[PortSpanRoute, ...]
    deferred: tuple
    logical_clusters: tuple = ()
    boundary_runs: tuple = ()

    @property
    def unresolved(self):
        declared = {run.span_identity for run in self.boundary_runs}
        return tuple(r for r in self.routes if r.kind in ('REMOTE_CONNECTION', 'UNSELECTED')
                     or (r.kind == 'BOUNDARY_ONLY' and r.span.identity not in declared))


class PreparedPortNetwork:
    def __init__(self, compound_reader):
        self.compound_reader = compound_reader

    def prepare(self, plan, *, compound_spans=(), logical_clusters=()):
        logical_clusters = tuple(logical_clusters)
        sites = {s.junction_id: s for s in plan.sites}
        spans = {s.identity: s for s in plan.spans}
        if len(sites) != len(plan.sites) or len(spans) != len(plan.spans):
            raise ValueError('Network sites and original spans must have unique identities')
        endpoint_spans = {}
        for span in plan.spans:
            for end in (span.start, span.end):
                if end.endpoint in endpoint_spans:
                    raise ValueError('An original Section end belongs to multiple spans')
                endpoint_spans[end.endpoint] = span
        ports = {}
        for site in plan.sites:
            for port in site.ports:
                if port.junction_id != site.junction_id or port.endpoint in ports:
                    raise ValueError('An original port has conflicting site ownership')
                if port.endpoint not in endpoint_spans:
                    raise ValueError('An original port has no source span')
                ports[port.endpoint] = port
        compounds = []; absorbed = {}; used_sites = set()
        for identity in compound_spans:
            compound = self.compound_reader.prepare(plan, identity)
            if compound.span_identity in absorbed or used_sites.intersection(compound.site_ids):
                raise ValueError('Compound declarations overlap an existing surface owner')
            absorbed[compound.span_identity] = ('COMPOUND', set(compound.internal_ports))
            used_sites.update(compound.site_ids); compounds.append(compound)
        logical_ids = set()
        for cluster in logical_clusters:
            ids = set(cluster.physical_junction_ids)
            owners = tuple(cluster.owner_ids)
            if (not owners or len(set(owners)) != len(owners) or
                    any(type(jid) is not int or jid < 0 for jid in owners) or
                    cluster.logical_junction_id != min(owners)):
                raise ValueError('Logical surface needs its complete distinct original owner IDs')
            group_owners = {getattr(group, 'logical_junction_id', cluster.logical_junction_id)
                            for group in (*cluster.groups, *cluster.internal_groups)}
            if group_owners != set(owners):
                raise ValueError('Logical owner IDs disagree with original connector memberships')
            if (logical_ids.intersection(owners) or used_sites.intersection(ids) or
                    len(ids) != len(cluster.physical_junction_ids)):
                raise ValueError('Logical cluster declarations overlap an existing surface owner')
            if len(ids) < 2 or not ids.issubset(sites):
                raise ValueError('A logical surface must own its complete selected physical cluster')
            internal = {s.identity: s for s in plan.spans
                        if s.start.junction_id in ids and s.end.junction_id in ids}
            if (len(cluster.internal_spans) != len(internal) or
                    {s.identity: s for s in cluster.internal_spans} != internal):
                raise ValueError('Logical surface must retain every exact original internal span')
            outer = {p for sid in ids for p in sites[sid].ports
                     if endpoint_spans[p.endpoint].identity not in internal}
            if len(cluster.external_ring) != len(outer) or set(cluster.external_ring) != outer:
                raise ValueError('Logical surface must retain every exact original external port')
            for identity, span in internal.items():
                owned = {ports.get(span.start.endpoint), ports.get(span.end.endpoint)}
                if None in owned or len(owned) != 2:
                    raise ValueError('Logical surface needs both original ends of every internal span')
                absorbed[identity] = ('LOGICAL_CLUSTER', owned)
            used_sites.update(ids); logical_ids.update(owners)
        routes = []
        for span in plan.spans:
            ends = (span.start, span.end)
            owned = tuple(ports.get(end.endpoint) for end in ends)
            for end, port in zip(ends, owned):
                if end.owner_kind in ('FREE', 'CATEGORY_BOUNDARY', 'SELECTION_BOUNDARY'):
                    if end.junction_id is not None or port is not None:
                        raise ValueError('A free or category end cannot have a junction owner')
                    if end.owner_kind == 'CATEGORY_BOUNDARY':plan.category_boundary(end.endpoint)
                    if end.owner_kind == 'SELECTION_BOUNDARY':plan.selection_boundary(end.endpoint)
                elif end.junction_id is None:
                    raise ValueError('A non-free end needs an explicit original owner')
                if port is not None and port.junction_id != end.junction_id:
                    raise ValueError('Span endpoint and port name different junction owners')
                if end.junction_id in sites and port is None and span.in_scope:
                    raise ValueError('A selected junction is missing an original span incidence')
            count = sum(p is not None for p in owned)
            if not span.in_scope:
                if count or span.identity in absorbed:
                    raise ValueError('An excluded span cannot have a carriageway surface owner')
                kind = 'DEFERRED_CATEGORY'
            elif span.identity in absorbed:
                kind, expected = absorbed[span.identity]
                if count != 2 or set(owned) != expected:
                    raise ValueError('Surface does not own both original ends of its span')
            elif count == 2:
                kind = 'CORRIDOR'
            elif count == 1:
                remote = ends[0 if owned[0] is None else 1]
                kind = ('APPROACH' if remote.owner_kind == 'FREE' else 'CATEGORY_APPROACH'
                        if remote.owner_kind == 'CATEGORY_BOUNDARY' else 'SELECTION_APPROACH'
                        if remote.owner_kind == 'SELECTION_BOUNDARY' else 'REMOTE_CONNECTION')
            elif all(end.owner_kind == 'FREE' for end in ends):
                kind = 'STANDALONE'
            elif all(end.owner_kind in ('FREE', 'CATEGORY_BOUNDARY', 'SELECTION_BOUNDARY') for end in ends):
                kind = 'BOUNDARY_ONLY'
            else:
                kind = 'UNSELECTED'
            routes.append(PortSpanRoute(span, kind, *owned))
        return PortNetworkPlan(tuple(s.junction_id for s in plan.sites if s.junction_id not in used_sites),
                               tuple(compounds), tuple(routes), tuple(plan.deferred), tuple(logical_clusters))
