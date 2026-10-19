"""Original-ID declarations for a Bundle region ending at a crossing road.

No coordinates or generated mesh are accepted here. The first branch has one
complete crossing Street and a common attached cap at the opposite end. A
declaration is not separator coverage: leftover member extents are reported
explicitly, and native boundary/surface checks remain mandatory.
"""
from dataclasses import dataclass

from .original_span_separator_chains import ChainedSeparatorPlanner
from .physical_ports import SectionEnd


@dataclass(frozen=True)
class CrossingSeparatorRegion:
    bundle_id: int
    member_street_ids: tuple[int, int]
    head_sites: tuple[int, ...]
    tail_sites: tuple[tuple[int, ...], tuple[int, ...]]
    member_span_ids: tuple[tuple[tuple[int, int, int], ...], ...]
    head_ends: tuple[SectionEnd, SectionEnd]
    tail_ends: tuple[SectionEnd, SectionEnd]
    intermediate_sites: tuple[tuple[tuple[int, ...], ...], ...]
    crossing_span_id: tuple[int, int, int]
    crossing_ends: tuple[SectionEnd, SectionEnd]
    remaining_member_span_ids: tuple[tuple[tuple[int, int, int], ...], ...]
    deferred_interiors: tuple[dict, ...]


@dataclass(frozen=True)
class _MemberChain:
    street_id: int
    routes: tuple
    ends: tuple

    def reversed(self):
        return _MemberChain(self.street_id, self.routes[::-1],
                            tuple(pair[::-1] for pair in self.ends[::-1]))


@dataclass(frozen=True)
class CrossingRegionPreparation:
    regions: tuple[CrossingSeparatorRegion, ...]
    pending: tuple[dict, ...]


class CrossingSeparatorRegionPlanner(ChainedSeparatorPlanner):
    """Prepare an attached region without mistaking it for the whole Bundle."""

    def prepare(self, snapshot, component, *, bundle_ids):
        """Explicit region inventory; never return these as complete chains."""
        ids = tuple(bundle_ids)
        if len(ids)!=len(set(ids)) or any(type(b) is not int or b not in snapshot.bundles for b in ids):
            raise ValueError('Region inventory requires unique original Bundle identities')
        regions=[];pending=[]
        for bid in ids:
            try: regions.append(self.candidate(snapshot.bundles[bid],snapshot,component))
            except ValueError as exc: pending.append(dict(bundle_id=bid,reason=str(exc)))
        return CrossingRegionPreparation(tuple(regions),tuple(pending))

    def interiors(self, bundle, snapshot, component):
        ids = tuple(bundle.get('topology', {}).get('interior_streets', ()))
        if self.categories is None:
            raise ValueError('Crossing regions require an explicit feature category scope')
        if len(ids) != len(set(ids)) or any(type(sid) is not int or sid < 0 for sid in ids):
            raise ValueError('Interior Street identities must be unique nonnegative integers')
        active = []; deferred = []
        generated = {r.span.street_id for r in component.network.routes if r.span.in_scope}
        for sid in ids:
            street = snapshot.streets.get(sid)
            if street is None: raise ValueError('Missing original interior Street')
            sections = [s for s in street.get('items', ()) if s.get('kind') == 'Section']
            categories = {s.get('category') for s in sections}
            if not categories or any(not isinstance(c, str) or not c for c in categories):
                raise ValueError('Interior Street lacks complete category evidence')
            if categories.intersection(self.categories) or sid in generated:
                active.append(sid)
            else:
                deferred.append(dict(street_id=sid, section_ids=tuple(s['id'] for s in sections),
                                     categories=tuple(sorted(categories))))
        return tuple(active), tuple(deferred)

    def deferred_interiors(self, bundle, snapshot, component):
        # Only this region planner permits active interiors. candidate() below
        # must prove their exact original connectivity before returning a region.
        return self.interiors(bundle, snapshot, component)[1]

    def _member_chain(self, member, snapshot, component, declared, used, owners):
        sid = member['street_id']
        routes = sorted((r for r in component.network.routes if r.span.street_id == sid),
                        key=lambda r:r.span.first_control)
        if not routes or any(r.kind not in ('CORRIDOR','APPROACH') or not r.span.in_scope for r in routes):
            raise ValueError('Region member needs accepted in-scope original routes throughout')
        sections = tuple(s['id'] for s in snapshot.streets[sid]['items'] if s['kind']=='Section')
        actual = tuple(s for r in routes for s in r.span.section_ids)
        count = routes[0].span.source_control_count
        if (actual != sections or len(set(actual)) != len(actual) or
                routes[0].span.first_control != 0 or routes[-1].span.last_control != count-1 or
                any(r.span.source_control_count != count or
                    not 0 <= r.span.first_control < r.span.last_control < count for r in routes) or
                any(a.span.last_control != b.span.first_control for a,b in zip(routes,routes[1:]))):
            raise ValueError('Region member must retain its complete original Street coverage')
        ends = tuple((self.route_end(r,0,allow_free=True),self.route_end(r,1,allow_free=True)) for r in routes)
        for i,(route,pair) in enumerate(zip(routes,ends)):
            free = [at for at,end in enumerate(pair) if end.owner_kind=='FREE']
            if route.kind != ('APPROACH' if free else 'CORRIDOR') or len(free)>1:
                raise ValueError('Only a genuine free terminal may use an approach')
            if free and (i,free[0]) not in ((0,0),(len(routes)-1,1)):
                raise ValueError('A free member end cannot occur inside a region chain')
        for previous,following in zip(ends,ends[1:]):
            a,b = previous[1],following[0]
            key = (sid,a.junction_id); row = declared.get(key)
            if a.junction_id != b.junction_id or row is None or key in used or a.junction_id not in owners:
                raise ValueError('Region chain needs every original intermediate junction exactly once')
            through = [tuple(v) for v in row.get('through',())]
            if len(through)!=2 or set(through)!={(a.endpoint.section_id,1),(b.endpoint.section_id,0)}:
                raise ValueError('Region through endpoints differ from the original chain')
            used.add(key)
        chain = _MemberChain(sid,tuple(routes),ends)
        if not member['forward']: chain=chain.reversed()
        if len(member['spans'])!=len(chain.routes):
            raise ValueError('Bundle extents disagree with the region member chain')
        for i,(extent,route,pair) in enumerate(zip(member['spans'],chain.routes,chain.ends)):
            expected = route.span.section_ids if member['forward'] else route.span.section_ids[::-1]
            rows = extent.get('sections',())
            if (extent.get('ordinal')!=i or tuple(v['start']['section_id'] for v in rows)!=expected or
                    any(v['start']['section_id']!=v['end']['section_id'] for v in rows)):
                raise ValueError('Region extent lost its directed original Section inventory')
            for at,label in enumerate(('start','end')):
                edge = (i==0 and at==0) or (i==len(chain.routes)-1 and at==1)
                endpoint = extent.get(label,{})
                if (endpoint.get('bundle_end') != (('HEAD' if at==0 else 'TAIL') if edge else None) or
                        endpoint.get('junction_id') != (None if edge else pair[at].junction_id)):
                    raise ValueError('Region extent endpoints disagree with original member joins')
        return chain

    def candidate(self, bundle, snapshot, component, owners=None):
        owners = self.owners(component) if owners is None else owners
        members,ids = self.members(bundle,component,snapshot,split_extents=True)
        active,deferred = self.interiors(bundle,snapshot,component)
        if len(active)!=1:
            raise ValueError('This region branch requires one complete active crossing Street')
        crossing_id, = active
        if crossing_id in ids: raise ValueError('A crossing cannot be its own Bundle member')
        crossing_routes = [r for r in component.network.routes if r.span.street_id==crossing_id]
        if len(crossing_routes)!=1:
            raise ValueError('Crossing must be one complete original corridor')
        crossing, = crossing_routes
        sections = tuple(s['id'] for s in snapshot.streets[crossing_id]['items'] if s['kind']=='Section')
        if (crossing.kind!='CORRIDOR' or not crossing.span.in_scope or
                crossing.span.section_ids!=sections or len(set(sections))!=len(sections) or
                crossing.span.first_control!=0 or crossing.span.last_control!=crossing.span.source_control_count-1):
            raise ValueError('Crossing must retain its complete in-scope original Street')
        cross_ends = (self.route_end(crossing,0),self.route_end(crossing,1))
        if cross_ends[0].junction_id==cross_ends[1].junction_id:
            raise ValueError('Crossing endpoints must be distinct original junctions')
        junctions = bundle['topology'].get('junctions',())
        declared = {(j['street_id'],j['junction_id']):j for j in junctions}
        if len(declared)!=len(junctions): raise ValueError('Duplicate original Bundle junction')
        used=set()
        chains=tuple(self._member_chain(m,snapshot,component,declared,used,owners) for m in members)
        if used!=set(declared): raise ValueError('Region inventory omitted an original Bundle junction')
        def common_head(items):
            sites = [owners.get(c.ends[0][0].junction_id) for c in items]
            return sites[0] if sites[0] and sites[0]==sites[1] else None
        head = common_head(chains)
        if head is None:
            chains=tuple(c.reversed() for c in chains); head=common_head(chains)
        if head is None: raise ValueError('Crossing-capped region needs a common attached opposite end')
        identities=[]; remainders=[]; first=[]; last=[]; middle=[]; caps=[]; cross_order=[]
        for chain in chains:
            sites=[chain.ends[0][0].junction_id]+[pair[1].junction_id for pair in chain.ends]
            matches=[(i,end) for i,site in enumerate(sites) for end in cross_ends if site==end.junction_id]
            if len(matches)!=1 or matches[0][0]==0:
                raise ValueError('Crossing must meet each member exactly once beyond its attached cap')
            at,end=matches[0]
            prefix=chain.routes[:at]
            if any(r.kind!='CORRIDOR' for r in prefix): raise ValueError('Bounded region needs attached member corridors')
            tail=owners.get(end.junction_id)
            if not tail or tail==head: raise ValueError('Crossing terminal needs a separate accepted original owner')
            mids=tuple(owners[pair[1].junction_id] for pair in chain.ends[:at-1])
            identities.append(tuple((r.span.street_id,r.span.section_ids[0],r.span.section_ids[-1]) for r in prefix))
            remainders.append(tuple((r.span.street_id,r.span.section_ids[0],r.span.section_ids[-1]) for r in chain.routes[at:]))
            first.append(chain.ends[0][0].endpoint);last.append(chain.ends[at-1][1].endpoint)
            middle.append(mids);caps.append(tail);cross_order.append(end.endpoint)
        if len(set(cross_order))!=2 or caps[0]==caps[1]:
            raise ValueError('Crossing must connect two distinct member-side owners')
        all_middle=[site for row in middle for site in row]
        if len(all_middle)!=len(set(all_middle)) or set(all_middle).intersection((head,*caps)):
            raise ValueError('Shared intermediate owners require a different region partition')
        return CrossingSeparatorRegion(bundle['id'],ids,head,tuple(caps),tuple(identities),
            tuple(first),tuple(last),tuple(middle),(crossing_id,sections[0],sections[-1]),
            tuple(cross_order),tuple(remainders),deferred)
