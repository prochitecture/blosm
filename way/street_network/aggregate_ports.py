"""Original path bindings for an explicitly declared aggregate junction mouth.

An aggregate reference point is not a coincident physical OSM node. Its paired
members keep independent original endpoints and one shared profile/connector.
This plan binds paths only; it is not a surface or coverage certificate.
"""
from dataclasses import dataclass, replace

from .physical_ports import SectionEnd, SpanEnd, PhysicalSpan


@dataclass(frozen=True)
class AggregatePort:
    junction_id: int
    endpoint: SectionEnd
    connector_id: int
    source_kind: str = 'AGGREGATE'


@dataclass(frozen=True)
class AggregateSite:
    junction_id: int
    location: tuple
    ports: tuple


@dataclass(frozen=True)
class AggregateDesign:
    radius: float
    rise: float
    datum: float
    half: int
    bundle_id: int
    partner: int  # Local ring index, -1 for an ordinary mouth.
    bridge: bool  # The following side crosses this connector's separator.


@dataclass(frozen=True)
class AggregatePortPlan:
    sites: tuple
    spans: tuple
    source_streets: tuple
    designs: tuple
    section_scopes: tuple = ()
    category_connections: tuple = ()
    selection_boundaries: tuple = ()
    deferred: tuple = ()

    def selection_boundary(self, endpoint):
        rows = [b for b in self.selection_boundaries if b.endpoint == endpoint]
        if len(rows) != 1:
            raise ValueError('Selection boundary needs one exact original endpoint declaration')
        return rows[0]

    def selection_peer(self, endpoint):
        self.selection_boundary(endpoint)
        # This bridge owns whole Streets, so accepted cuts are terminal only.
        # Interior Section cuts require the physical span partitioner.
        return None

    @staticmethod
    def source_supported(source):
        return source.get('status') == 'SUPPORTED'

    def span_for(self, port):
        found=[s for s in self.spans if port.endpoint in (s.start.endpoint,s.end.endpoint)]
        if len(found)!=1:raise ValueError('Aggregate incidence needs one original complete Street')
        return found[0]


class PreparedAggregatePorts:
    def complete(self, snapshot, junctions, *, free_endpoints, selection_boundaries=()):
        """Resolve complete network extents from explicit terminal declarations.

        A junction-local UNBOUND end never becomes free by omission. The
        caller must supply the free-end inventory and every selected junction.
        This method only binds original topology; native surfaces still have
        to qualify the actual cuts, tapers and separator openings.
        """
        plan = self.prepare(snapshot, junctions)
        if set(plan.source_streets) != set(snapshot.streets):
            raise ValueError('Aggregate network needs a route for every declared original Street')
        free_endpoints = tuple(free_endpoints)
        selection_boundaries = tuple(selection_boundaries)
        free = set(free_endpoints)
        cuts = {b.endpoint: b for b in selection_boundaries}
        if len(free) != len(free_endpoints) or len(cuts) != len(selection_boundaries):
            raise ValueError('Duplicate original terminal declaration')
        if free.intersection(cuts):
            raise ValueError('Selection boundary cannot also be a free end')
        terminal = {e.endpoint: e for s in plan.spans for e in (s.start, s.end)}
        if not free.union(cuts).issubset(terminal):
            raise ValueError('Terminal declaration is not an original complete Street end')
        spans = []
        for span in plan.spans:
            resolved = []
            source = snapshot.streets[span.street_id]
            owners = source.get('endpoints')
            if owners is None or len(owners) != 2:
                raise ValueError('Complete Street requires its original endpoint ownership inventory')
            for end in (span.start, span.end):
                owner = owners[end.endpoint.end]
                if end.junction_id is not None:
                    if end.endpoint in free or end.endpoint in cuts or owner != end.junction_id:
                        raise ValueError('Junction port conflicts with original terminal ownership')
                    resolved.append(end)
                else:
                    if owner is not None:
                        raise ValueError('An unselected original junction cannot become a free end')
                    kind = ('FREE' if end.endpoint in free else 'SELECTION_BOUNDARY'
                            if end.endpoint in cuts else None)
                    if kind is None:
                        raise ValueError('Unbound original end needs an explicit terminal declaration')
                    resolved.append(replace(end, owner_kind=kind))
            spans.append(replace(span, start=resolved[0], end=resolved[1]))
        return replace(plan, spans=tuple(spans), selection_boundaries=selection_boundaries)

    def prepare(self, snapshot, junctions):
        junctions=tuple(junctions)
        if not junctions or len({j.source_id for j in junctions})!=len(junctions):
            raise ValueError('Aggregate junction IDs must be nonempty and unique')
        sites=[];designs=[];used={};sources=set()
        for junction in junctions:
            ports=[];local=[]
            for connector in junction.connectors:
                profile=junction.profiles[connector.profile_index]
                paired=len(connector.members)==2
                if len(connector.members) not in (1,2) or profile.separator!=paired:
                    raise ValueError('Aggregate mouths need an ordinary profile or an explicit paired separator')
                if profile.half_profile or profile.preset!=1:
                    raise ValueError('Aggregate mouth construction requires full smooth incoming profiles')
                if paired and (connector.item_kind!=1 or connector.bundle_source_id<0):
                    raise ValueError('Paired aggregate mouth has no original Bundle owner')
                entries=[]
                for member in connector.members:
                    street=junction.streets[member.street_index];sid=street.source_id
                    source=snapshot.streets.get(sid)
                    if source is None or not AggregatePortPlan.source_supported(source):
                        raise ValueError('Aggregate member needs a supported complete original Street')
                    sections=[s for s in source['items'] if s['kind']=='Section']
                    expected=tuple((s['id'],tuple(s['points']),s['width'],s['offset']) for s in sections)
                    actual=tuple((s.source_id,tuple(s.centerline),s.width,s.offset) for s in street.sections)
                    if actual!=expected:raise ValueError('Aggregate incidence changed its original Street chain')
                    end=int(not member.leaving);section=sections[0 if end==0 else -1]
                    endpoint=SectionEnd(sid,section['id'],end)
                    if endpoint in used:raise ValueError('Original endpoint has multiple aggregate owners')
                    used[endpoint]=junction.source_id;sources.add(sid)
                    half=0
                    if paired:
                        if member.head_order not in (0,1) or member.head_order!=member.tail_order:
                            raise ValueError('Aggregate Bundle must retain its declared lateral order')
                        half=(1 if member.head_order==0 else -1)*(1 if connector.leaving else -1)
                        if source.get('bundle_id')!=connector.bundle_source_id:
                            raise ValueError('Aggregate connector and complete Street Bundle owners disagree')
                    entries.append((half,AggregatePort(junction.source_id,endpoint,connector.source_id)))
                # Walk across the first member, then the separator, then the
                # second member. Half signs describe the outward mouth frame.
                entries.sort(key=lambda row:row[0])
                if paired and [h for h,p in entries]!=[-1,1]:
                    raise ValueError('Paired aggregate mouth needs two complementary half profiles')
                start=len(ports)
                for i,(half,port) in enumerate(entries):
                    ports.append(port)
                    local.append(AggregateDesign(connector.radius,profile.rise,junction.datum,half,
                        connector.bundle_source_id if paired else -1,start+1-i if paired else -1,paired and i==0))
            sites.append(AggregateSite(junction.source_id,junction.location,tuple(ports)));designs.extend(local)
        spans=[]
        for sid in sorted(sources):
            items=[s for s in snapshot.streets[sid]['items'] if s['kind']=='Section']
            count=1+sum(len(s['points'])-1 for s in items)
            ends=[]
            for end,section in ((0,items[0]),(1,items[-1])):
                endpoint=SectionEnd(sid,section['id'],end);jid=used.get(endpoint)
                # An unselected far end remains unbound. It is not declared
                # free by this junction-local input bridge.
                ends.append(SpanEnd(endpoint,'AGGREGATE' if jid is not None else 'UNBOUND',jid))
            spans.append(PhysicalSpan(sid,tuple(s['id'] for s in items),*ends,0,count-1,count))
        return AggregatePortPlan(tuple(sites),tuple(spans),tuple(sorted(sources)),tuple(designs))
