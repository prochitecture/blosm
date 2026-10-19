"""Once-owned routes for explicit aggregate connectors on whole Streets.

No geometry is inferred here. Complete terminal and Bundle declarations are
required; the shared native provider qualifies their actual surface joins.
"""
from dataclasses import dataclass

from .aggregate_ports import PreparedAggregatePorts
from .port_network import PreparedPortNetwork
from .physical_dimensions import PreparedPhysicalDimensions
from .original_span_separators import (
    ClosedBundleSeparator, FreeEndBundleSeparator, SelectionEndBundleSeparator, SeparatorPreparation,
)


@dataclass(frozen=True)
class AggregateBundleMembers:
    bundle_id: int
    street_ids: tuple
    forward: tuple


@dataclass(frozen=True)
class AggregateNetworkPlan:
    ports: object
    network: object
    separators: SeparatorPreparation
    dimensions: tuple


class PreparedAggregateNetwork:
    def prepare(self, snapshot, junctions, *, free_endpoints,
                selection_boundaries=(), bundle_members=()):
        ports = PreparedAggregatePorts().complete(snapshot, junctions,
            free_endpoints=free_endpoints, selection_boundaries=selection_boundaries)
        network = PreparedPortNetwork(None).prepare(ports)
        if network.unresolved:
            raise ValueError('Aggregate network retains an unresolved original connection')
        dimensions = PreparedPhysicalDimensions().read(snapshot, ports.source_streets)
        members = tuple(bundle_members)
        if len({b.bundle_id for b in members}) != len(members):
            raise ValueError('Aggregate Bundle declaration has duplicate ownership')
        expected = {s['bundle_id'] for s in snapshot.streets.values() if s.get('bundle_id') is not None}
        if {b.bundle_id for b in members} != expected:
            raise ValueError('Every original Bundle needs its complete member declaration')
        routes = {r.span.street_id: r for r in network.routes}
        design = {p.endpoint: d for p, d in zip(
            (p for site in ports.sites for p in site.ports), ports.designs)}
        closed, free, clipped, pending = [], [], [], []
        for bundle in members:
            ids = bundle.street_ids
            if (type(bundle.bundle_id) is not int or bundle.bundle_id < 0 or
                    len(ids) != 2 or len(set(ids)) != 2 or
                    any(type(sid) is not int or sid < 0 for sid in ids) or
                    len(bundle.forward) != 2 or any(type(v) is not bool for v in bundle.forward)):
                raise ValueError('Aggregate Bundle needs two distinct, explicitly directed lateral members')
            original = {sid for sid, s in snapshot.streets.items() if s.get('bundle_id') == bundle.bundle_id}
            if set(ids) != original:
                raise ValueError('Aggregate Bundle members disagree with original Street ownership')
            oriented = []
            for rank, (sid, forward) in enumerate(zip(ids, bundle.forward)):
                route = routes[sid]
                ends = (route.span.start, route.span.end)
                if not forward:
                    ends = ends[::-1]
                for end_index, end in enumerate(ends):
                    if end.owner_kind == 'AGGREGATE':
                        d = design[end.endpoint]
                        expected_half = (1 if rank == 0 else -1)*(1 if end_index == 0 else -1)
                        if d.bundle_id != bundle.bundle_id or d.half != expected_half:
                            raise ValueError('Aggregate Bundle direction/order disagrees with its mouth profile')
                oriented.append(ends)
            owners = [tuple(e.junction_id for e in ends) for ends in oriented]
            if owners[0] != owners[1]:
                raise ValueError('Aggregate Bundle members name different terminal owners')
            attached = [i for i, owner in enumerate(owners[0]) if owner is not None]
            if len(attached) == 2 and owners[0][0] != owners[0][1]:
                closed.append(ClosedBundleSeparator(bundle.bundle_id, ids,
                    (owners[0][0],), (owners[0][1],)))
            elif len(attached) == 1 and all(ends[1-attached[0]].owner_kind == 'FREE' for ends in oriented):
                at = attached[0]
                free.append(FreeEndBundleSeparator(bundle.bundle_id, ids, (owners[0][at],), at == 0,
                    tuple(ends[1-at].endpoint for ends in oriented)))
            elif len(attached) == 1 and all(ends[1-attached[0]].owner_kind in
                                            ('FREE','SELECTION_BOUNDARY') for ends in oriented):
                at=attached[0]
                terminals=tuple(ends[1-at] for ends in oriented)
                for end in terminals:
                    if end.owner_kind=='SELECTION_BOUNDARY':ports.selection_boundary(end.endpoint)
                clipped.append(SelectionEndBundleSeparator(bundle.bundle_id,ids,(owners[0][at],),at==0,terminals))
            else:
                pending.append(dict(bundle_id=bundle.bundle_id, status='PENDING',
                    reason='Aggregate separator needs two distinct junction owners or one owner and two explicit free/selection ends'))
        return AggregateNetworkPlan(ports, network,
            SeparatorPreparation(tuple(closed), tuple(pending), tuple(free),selection_ended=tuple(clipped)), tuple(dimensions))
