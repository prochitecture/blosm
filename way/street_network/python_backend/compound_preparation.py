"""Prepare one physical exterior when close source junctions cannot stand alone.

A short connecting road is not an entrance to the final compound. Requiring a
complete isolated mouth on it can fail before ordinary footprint-union planning
is possible. This module defers that failure without inventing a provisional
polygon or extrapolating the road, then certifies one directly constructed rim.
"""
from dataclasses import dataclass, replace
import math

from .compounds import CompoundGroup, SharedClaim, sample_footprint
from .model import End, Key
from .routes import RouteTopology
from .regions import RegionUnion
from .updates import drain_steps


class CompoundExpansionRequired(ValueError):
    """More original source intervals overlap before the rim is final."""
    def __init__(self,neighbors):
        self.neighbors=tuple(sorted(neighbors))
        super().__init__('Compound exterior claims an opposite connected junction mouth')


@dataclass(frozen=True)
class DeferredCompoundFootprint:
    site: object
    approaches: tuple
    neighbors: tuple
    failure: str
    absorbed_transitions: tuple = ()

    @classmethod
    def from_failure(cls, context, junction, approaches, failure):
        """Only explicit short routes to another physical junction may defer."""
        from .intersections import ApproachExtentError
        topology = RouteTopology(context.network)
        physical = {j.key for j in context.network.intersections if len(j.ends) > 2}
        failed = {End(a.source.id, a.end) for a in failure.approaches}
        neighbors = set()
        for approach in approaches:
            end = End(approach.source.id, approach.end)
            if end not in failed:
                continue
            candidate = topology.trace(end)
            if (candidate.terminal_site not in physical or candidate.terminal_site == junction.key
                    or candidate.terminal_selection or candidate.terminal_scope
                    or getattr(approach.path, 'absorbed_sites', ())
                    or tuple(p.section.id for p in approach.path.pieces)
                       != tuple(v.section_id for v in candidate.traversals)):
                continue
            # A failed numerical contact on a long road must not pull a remote
            # junction into this one. A width/radius-scaled shortness bound only
            # admits a candidate; final containment and surface tests certify it.
            radius, _ = junction.corner_radius(end)
            reach = max(a.width + 2. * abs(a.offset) for a in approaches) + 2. * radius
            overlaps=False
            if not isinstance(failure,ApproachExtentError) and approach.length>reach:
                from .overlapping_approaches import persistent_overlap
                from .corner_radii import CornerRadiusPolicy
                first,second=failure.approaches
                requested,mode=junction.corner_radius(End(first.source.id,first.end))
                minimum=min(CornerRadiusPolicy.contact_candidates(first,second,requested,mode))
                # A long neighboring pair can remain too close for even its
                # smallest permitted fillet all the way to an actual junction.
                # Prove that bounded construction corridor before considering
                # the remote incidence; final region/ownership guards still run.
                overlaps=any(other is not approach and persistent_overlap(approach,other,clearance=2.*minimum)
                    for other in failure.approaches)
            if isinstance(failure, ApproachExtentError) or approach.length <= reach or overlaps:
                neighbors.add(candidate.terminal_site)
        if not neighbors:
            return None
        # Cached plans must not retain an exception traceback and its entire
        # build context/mesh inventory during large progressive imports.
        return cls(junction.key, tuple(approaches), tuple(sorted(neighbors)), str(failure))


class DeferredCompoundPlanner:
    maximum_members = 8

    def __init__(self, context, junctions, layouts):
        self.context, self.layouts = context, layouts
        self.sites = {j.key:j for j in junctions}
        self.topology = RouteTopology(context.network)

    def build(self, ordinary_groups):
        return drain_steps(self.build_steps(ordinary_groups))

    def build_steps(self, ordinary_groups):
        ordinary_groups=tuple(ordinary_groups)
        pending = {key for key, layout in self.layouts.items() if isinstance(layout,DeferredCompoundFootprint)}
        from .compound_ports import CompoundPortPlanner
        # Individually valid junctions can still clip one another's mouths.
        # Pure source-Section compounds need the same complete exterior policy
        # as deferred contacts; a partial profile is not a valid chart input.
        # Absorbed dimensional seams retain their established dedicated port
        # builder below OwnershipPlan's deferred pass.
        pending.update(key for group in ordinary_groups
            if CompoundPortPlanner.incomplete_ports(group)
            and not any(getattr(self.layouts[member],'absorbed_transitions',()) for member in group.members)
            for key in group.members)
        width_planner=CompoundPortPlanner(self.context,tuple(self.sites.values()),self.layouts)
        for group in ordinary_groups:
            yield
            if (width_planner.incomplete_source_widths(group) and not any(
                    getattr(self.layouts[member], 'absorbed_transitions', ()) for member in group.members)):
                pending.update(group.members)
        if not pending:
            return tuple(ordinary_groups)
        parent = {key:key for key in self.layouts}
        join_budget=max(0,len(parent)-1)
        joins=0
        def root(key):
            while parent[key] != key:
                parent[key] = parent[parent[key]]
                key = parent[key]
            return key
        def join(first,last):
            nonlocal joins
            if first not in parent or last not in parent:
                raise ValueError('Deferred expansion requires an available original physical neighbor')
            a,b = sorted((root(first),root(last)))
            if a==b:return False
            if joins>=join_budget:
                raise ValueError('Deferred source closure exceeded its finite original graph join budget')
            parent[b] = a
            joins+=1
            return True
        for group in ordinary_groups:
            for key in group.members[1:]:join(group.members[0],key)
        for key in pending:
            for neighbor in getattr(self.layouts[key],'neighbors',()):
                if neighbor not in self.layouts:
                    raise ValueError('Deferred junction requires an available physical neighbor')
                join(key,neighbor)
        # A paired Bundle whose two original terminals are common source
        # junctions is a closed incidence cycle. Its overlap cannot be owned
        # by only one side's short-road compound while the other is rendered
        # independently. Retain every original member junction in the claim.
        from .compound_cycles import shared_bundle_sites
        for bundle in self.context.network.bundles:
            complete=shared_bundle_sites(self.context.network,bundle)
            if complete and pending.intersection(complete):
                if any(key not in parent for key in complete):
                    raise ValueError('Closed source Bundle needs all original physical junctions')
                for key in complete[1:]:join(complete[0],key)
        # Context/network/layouts stay fixed throughout this closure. A
        # certified member tuple remains valid while unrelated owners grow.
        # Failed attempts are never cached. Every restart must reduce the
        # finite original forest; there is no network-wide eight-pass cap.
        exteriors={}
        while True:
            members = {}
            for key in sorted(parent):members.setdefault(root(key),[]).append(key)
            active={tuple(keys) for keys in members.values()}
            exteriors={token:group for token,group in exteriors.items() if token in active}
            rebuilt = {}
            expanded=False
            for owner,keys in members.items():
                if pending.intersection(keys):
                    if len(keys)>self.maximum_members:
                        from .source_pool_budget import MAXIMUM_ORIGINAL_SITES
                        if len(keys)>MAXIMUM_ORIGINAL_SITES:
                            raise ValueError('Deferred junction compound exceeds its local member budget')
                    token=tuple(keys)
                    try:
                        if token not in exteriors:
                            yield
                            candidate=yield from self.exterior_steps(token)
                            if len(keys)>self.maximum_members:
                                from .source_pool_budget import certify_steps as certify_source_pool
                                try:yield from certify_source_pool(self.context,candidate)
                                except ValueError as failure:
                                    raise ValueError('Deferred junction compound exceeds its local member budget: '+str(failure)) from failure
                            exteriors[token]=candidate
                        rebuilt[owner] = exteriors[token]
                    except CompoundExpansionRequired as required:
                        if any(neighbor not in self.sites or neighbor not in parent for neighbor in required.neighbors):
                            raise ValueError('Deferred expansion requires an available original physical neighbor') from required
                        if not any(root(owner)!=root(neighbor) for neighbor in required.neighbors):
                            raise ValueError('Deferred source expansion must add a new original neighbor') from required
                        for neighbor in required.neighbors:join(owner,neighbor)
                        expanded=True
                        break
            if expanded:continue
            # A newly solved rim may consume another *connected* source mouth.
            # Expand ownership before constructing any faces, never by distance.
            mouth_index = {}
            for key,layout in self.layouts.items():
                owner=root(key)
                if owner in rebuilt:
                    if key != owner:continue
                    layout=rebuilt[owner].physical_layout
                for (sid,end),mouth in layout.mouths.items():
                    mouth_index.setdefault(sid,{0:[],1:[]})[end].append((owner,mouth.station))
            changed=False
            for ends in mouth_index.values():
                for first,lo in ends[0]:
                    for last,hi in ends[1]:
                        if root(first)!=root(last) and lo>hi+1e-7 and (root(first) in rebuilt or root(last) in rebuilt):
                            join(first,last);changed=True
            if not changed:
                untouched=[g for g in ordinary_groups if root(g.members[0]) not in rebuilt]
                return (*untouched,*tuple(rebuilt[key] for key in sorted(rebuilt)))

    def exterior(self, members, *, retain_source_curbs=False):
        return drain_steps(self.exterior_steps(members, retain_source_curbs=retain_source_curbs))

    def exterior_steps(self, members, *, retain_source_curbs=False):
        from .compound_ports import CompoundPortPlanner, ExteriorSpecification, SourceContainmentError
        from .intersections import IntersectionProvider, Mouth
        from .transitions import required
        sites = set(members)
        internal_ends=set();consumed=set();owners={};seams=set()
        for key in members:
            for end in self.sites[key].ends:
                route=self.topology.trace(end)
                if route.terminal_site not in sites or route.terminal_site==key:
                    continue
                if route.terminal_scope or route.terminal_selection:
                    raise ValueError('A compound cannot cross a declared source-scope boundary')
                internal_ends.add(end)
                seams.update(route.seams)
                for visit in route.traversals:
                    consumed.add(visit.section_id)
                    owners.setdefault(visit.section_id,{})[0 if visit.forward else 1]=key
        seam_sources={s.key:s for s in self.context.network.connections}
        if any(required(self.context.network,seam_sources[key]) for key in seams):
            raise ValueError('Deferred junction compound needs a separate dimensional-seam policy')
        if not consumed or any(set(ends)!={0,1} for ends in owners.values()):
            raise ValueError('Deferred compound needs both original incidences of every consumed route')
        approaches=[];radii={}
        for key in members:
            for approach in self.layouts[key].approaches:
                end=End(approach.source.id,approach.end)
                if end in internal_ends:continue
                if end in radii:raise ValueError('Deferred compound repeats an exterior source incidence')
                approaches.append(approach);radii[end]=self.sites[key].corner_radius(end)
        prefix_claims=tuple(a.path.claim for key in members for a in self.layouts[key].approaches
                            if hasattr(a.path,'claim'))
        seams.update(claim.seam for claim in prefix_claims)
        if len(approaches)<3:
            raise ValueError('Deferred junction compound needs at least three external approaches')
        profiles=[self.context[Key('profile',sid)] for sid in consumed]+[a.profile for a in approaches]
        bundle_ports,partial_claims=self.check_bundle_terminals(consumed,approaches,members)
        if max(p.datum for p in profiles)-min(p.datum for p in profiles)>1e-8:
            raise ValueError('Deferred junction compound needs a common minimum height')
        order,chains,cycles=self.exterior_order(members,internal_ends,include_chains=True,include_cycles=True)
        if any(claim.bundle_end=='WHOLE' for claim in partial_claims) and not cycles:
            raise ValueError('Whole Bundle ownership needs the physical source-cycle boundary policy')
        def source_union(failure=None):
            from .compound_cycles import source_cycle_group_steps
            group=yield from source_cycle_group_steps(self,members,consumed,internal_ends,owners,seams)
            if failure is not None:
                group.physical_layout.source_union_fallback=dict(
                    policy='FINITE_ORIGINAL_SOURCE_UNION', failed_curb_policy=type(failure).__name__,
                    failure=str(failure), original_members=[key.text() for key in members],
                    complete_source_sections=sorted(consumed))
            if not getattr(group.physical_layout,'source_cycle_full_width_proof',None):
                from .cycle_dimensional_retry import certify_layout
                group.physical_layout.source_cycle_full_width_proof=certify_layout(self.context,group,consumed)
            from .owned_prefixes import certify
            certify(self.context,self.sites[members[0]],group.physical_layout,prefix_claims,group=group)
            from .partial_bundle_claims import certify as certify_partial,certify_whole_separation
            certify_whole_separation(self.context,group,partial_claims)
            group=certify_partial(self.context,group,partial_claims)
            from .wholly_owned_bundle_preparation import certify as certify_free_counterpart
            return certify_free_counterpart(self.context,group)
        if cycles:return (yield from source_union())
        by_end={End(a.source.id,a.end):a for a in approaches}
        median_pairs=tuple((first,last) for first,last in zip(order,order[1:]+order[:1])
            if first in bundle_ports and last in bundle_ports
            and bundle_ports[first][::2]==bundle_ports[last][::2]
            and bundle_ports[first][1]!=bundle_ports[last][1]
            and by_end[first].profile.half==-1 and by_end[last].profile.half==1)
        specification=ExteriorSpecification(members[0].id,radii,order,median_pairs)
        from .compound_curbs import (retained_layout,RetainedCurbOverlapError,RetainedSourceRouteError,
                                    RetainedCornerContactError)
        try:
            if retain_source_curbs:
                physical=retained_layout(self,specification,approaches,chains,retain_all=True)
            else:
                physical=IntersectionProvider.curved_boundary(self.context,specification,approaches)
        except (RetainedCurbOverlapError,RetainedSourceRouteError,RetainedCornerContactError) as failure:
            return (yield from source_union(failure))
        except Exception as failure:
            from .intersections import CornerConstructionError
            if not isinstance(failure,CornerConstructionError):raise
            try:
                physical=retained_layout(self,specification,approaches,chains)
            except (RetainedCurbOverlapError,RetainedSourceRouteError,RetainedCornerContactError) as retained_failure:
                return (yield from source_union(retained_failure))
        mouths=dict(physical.mouths)
        claims=[]
        for sid in sorted(consumed):
            length=self.context[Key('path',sid)].length
            mouths[sid,0]=Mouth(length,2*physical.radial)
            mouths[sid,1]=Mouth(0.,2*physical.radial)
            claims.append(SharedClaim(sid,owners[sid][0],owners[sid][1],length,0.))
        physical=replace(physical,mouths=mouths,absorbed_transitions=tuple(sorted(seams)))
        # A provisional external profile can overlap another source site's
        # opposite claim. Resolve that connected ownership before asking the
        # unfinished rim to be simple; no provisional faces are generated.
        neighbors=set()
        for key,other in self.layouts.items():
            if key in sites or isinstance(other,DeferredCompoundFootprint):continue
            for (sid,end),mouth in physical.mouths.items():
                opposite=other.mouths.get((sid,1-end))
                if opposite is not None and (mouth.station>opposite.station+1e-7 if end==0
                                             else mouth.station<opposite.station-1e-7):
                    neighbors.add(key)
        if neighbors:raise CompoundExpansionRequired(neighbors)
        footprint=sample_footprint(members[0],physical)
        try:
            loops=RegionUnion().build((footprint.region,))
        except ValueError as failure:
            if str(failure)!='A physical footprint boundary intersects or touches itself':raise
            from .compound_port_overlap import shared_port_neighbors
            neighbors=shared_port_neighbors(self,set(members),physical)
            if neighbors:raise CompoundExpansionRequired(neighbors) from failure
            neighbors=self.next_source_hops(members,approaches)
            if not neighbors:raise
            raise CompoundExpansionRequired(neighbors) from failure
        if len(loops)!=1 or loops[0].hole:
            raise ValueError('Deferred junction compound needs one simple physical exterior')
        try:
            CompoundPortPlanner(self.context,tuple(self.sites.values()),{}).check_support(consumed,loops[0])
        except SourceContainmentError:
            if retain_source_curbs:raise
            return (yield from self.exterior_steps(members,retain_source_curbs=True))
        from .owned_prefixes import certify
        certify(self.context,self.sites[members[0]],physical,prefix_claims)
        from .partial_bundle_claims import certify as certify_partial
        group=CompoundGroup(members,tuple(claims),(footprint,),loops,physical)
        from .cycle_dimensional_retry import certify_layout
        from .source_strip_containment import SourceStripContainmentError
        try:
            physical.source_cycle_full_width_proof=certify_layout(self.context,group,consumed)
        except SourceStripContainmentError as failure:
            # Centerline clearance alone does not retain complete original
            # curbs, particularly beside a displaced returning source road.
            # Rebuild from the same finite strips and original local corners;
            # its unchanged source/root/profile guards must qualify again.
            from .compound_curbs import RetainedSourceWidthError
            return (yield from source_union(RetainedSourceWidthError(consumed,failure)))
        group=certify_partial(self.context,group,partial_claims)
        from .wholly_owned_bundle_preparation import certify as certify_free_counterpart
        return certify_free_counterpart(self.context,group)

    def next_source_hops(self,members,approaches):
        """Consider one more real junction beside a persistently close chain.

        A consumed internal road may continue through an already-owned source
        junction. Another external approach can remain too close to that whole
        original chain for a fillet, even though it extends beyond the first
        consumed road. The chain supplies a planning certificate only. No
        source site is removed, and the next pass repeats the incidence-cycle,
        simple-boundary, hole and complete-containment checks.
        """
        import numpy as np
        from .intersections import Approach
        from .routes import Route,Piece
        from .overlapping_approaches import persistent_overlap
        from .corner_radii import CornerRadiusPolicy
        sites=set(members)
        incidence={end:key for key in members for end in self.sites[key].ends}
        graph={key:[] for key in members}
        for end,key in incidence.items():
            route=self.topology.trace(end)
            if route.terminal_site in sites and route.terminal_site!=key:
                graph[key].append((route.terminal_site,route.traversals))

        def chain(first,last):
            queue=[(first,())];seen={first}
            for key,visits in queue:
                if key==last:return visits
                for neighbor,following in graph[key]:
                    if neighbor not in seen:
                        seen.add(neighbor);queue.append((neighbor,visits+following))
            return None

        neighbors=set()
        for first in approaches:
            end=End(first.source.id,first.end)
            route=self.topology.trace(end)
            if (first.profile.half or abs(first.offset)>1e-8
                    or getattr(first.path,'absorbed_sites',())
                    or route.terminal_site not in self.sites or route.terminal_site in sites
                    or route.terminal_selection or route.terminal_scope
                    or tuple(p.section.id for p in first.path.pieces)
                       !=tuple(v.section_id for v in route.traversals)):
                continue
            for other in approaches:
                other_end=End(other.source.id,other.end)
                if other is first or other.profile.half or getattr(other.path,'absorbed_sites',()):continue
                prefix=chain(incidence[end],incidence[other_end])
                if not prefix:continue
                pieces=tuple(Piece(self.context[Key('path',v.section_id)],v.forward) for v in prefix)+other.path.pieces
                if any(abs(p.section.offset)>1e-8 or self.context[Key('profile',p.section.id)].half for p in pieces):continue
                if any(np.linalg.norm(a.frame(a.path.length)[0]-b.frame(0.)[0])>1e-8
                       for a,b in zip(pieces,pieces[1:])):continue
                support=Approach(Route(pieces),0,self.context[Key('profile',pieces[0].section.id)])
                support.width=min(p.section.width for p in pieces)
                radius,mode=self.sites[incidence[other_end]].corner_radius(other_end)
                minimum=min(CornerRadiusPolicy.contact_candidates(other,first,radius,mode,displaced=True))
                if persistent_overlap(first,support,clearance=2.*minimum):
                    neighbors.add(route.terminal_site)
                    break
        return neighbors

    def exterior_order(self,members,internal,*,include_chains=False,include_cycles=False):
        """Walk the source junction rotations across wholly internal routes.

        Sorting only the final tangent directions is not a planar ordering
        when the mouths originate at different junctions. Nearly parallel
        Bundle arms can swap, selecting their outer instead of facing curbs.
        The incidence walk preserves adjacency without inventing one origin.
        A caller explicitly requesting closed cycles must retain their physical
        coverage and remaining holes before this order can become a surface.
        """
        following={}
        for key in members:
            ordered=sorted(self.layouts[key].approaches,
                key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
            ends=[End(a.source.id,a.end) for a in ordered]
            following.update(zip(ends,ends[1:]+ends[:1]))
        exterior=set(following)-internal
        if not exterior:raise ValueError('Compound has no external source incidence')
        first=min(exterior);current=first;order=[];visited=set();chains={}
        while current not in order:
            start=current;order.append(start);chain=[]
            while True:
                last=following[current]
                chain.append((current,last))
                if last not in internal:
                    current=last
                    break
                if last in visited:
                    raise ValueError('Compound incidence walk contains an unexposed internal cycle')
                visited.add(last)
                route=self.topology.trace(last)
                visit=route.traversals[-1]
                opposite=End(visit.section_id,1 if visit.forward else 0)
                if opposite not in internal:
                    raise ValueError('Compound route lacks its opposite original incidence')
                current=opposite
            chains[start,current]=tuple(chain)
        if current!=first or set(order)!=exterior or (visited!=internal and not include_cycles):
            raise ValueError('Compound exterior needs one complete source-incidence cycle')
        cycles=[]
        while visited!=internal:
            route=self.topology.trace(min(internal-visited));visit=route.traversals[-1]
            start=End(visit.section_id,1 if visit.forward else 0)
            current=start;chain=[];seen=set()
            while current not in seen:
                seen.add(current)
                last=following[current]
                if last not in internal or last in visited:
                    raise ValueError('Compound internal cycle has conflicting source incidences')
                visited.add(last);chain.append((current,last))
                route=self.topology.trace(last);visit=route.traversals[-1]
                current=End(visit.section_id,1 if visit.forward else 0)
                if current not in internal:
                    raise ValueError('Compound internal cycle lacks its opposite source incidence')
            if current!=start:
                raise ValueError('Compound internal incidence cycle does not close')
            cycles.append(tuple(chain))
        if include_cycles:return tuple(order),chains,tuple(cycles)
        return (tuple(order),chains) if include_chains else tuple(order)

    def check_bundle_terminals(self, consumed, approaches, members):
        """Retain complete exterior half-profiles and explicit Bundle claims.

        The established corner construction already supports full and half
        crowns. Global member terminals and a separately proved contiguous
        internal interval can use it here. A retained adjacent through pair at
        one owned original side junction keeps both complete half ports without
        consuming a Bundle body interval. Physical source coverage and exposed
        profiles qualify the latter before any mesh; accepted junction caps,
        visible regions and overlap checks still qualify the separator.
        """
        from .profiles import memberships
        bundled=memberships(self.context.network)
        ports={};partial_claims=[];checked=set();partial_cuts={};through_ports={};through_checked=set()
        from .compound_cycles import shared_bundle_sites
        for sid in consumed:
            record=bundled.get(sid)
            if record is None:continue
            bundle=record[0]
            if bundle.id in checked:continue
            checked.add(bundle.id)
            from .wholly_owned_bundle_preparation import candidate as free_counterpart
            if free_counterpart(self.context,bundle,consumed,set(members)) is not None:
                # The opposite original member has no physical terminal owner
                # to acquire. It stays independent; the completed candidate
                # must prove its entire ribbon outside the physical exterior.
                continue
            complete=shared_bundle_sites(self.context.network,bundle)
            if not complete or not set(bundle.sections)<=consumed:
                from .bundle_terminal_expansion import missing_terminals
                missing=missing_terminals(self.context,bundle,consumed,set(members),self.layouts)
                if missing:
                    raise CompoundExpansionRequired(missing)
                from .bundle_partial_terminal_expansion import missing_terminal
                missing=missing_terminal(self.context,bundle,consumed,set(members),self.layouts)
                if missing:
                    raise CompoundExpansionRequired(missing)
                from .partial_bundle_claims import prepare as prepare_partial
                claims=prepare_partial(self.context,bundle,consumed,approaches,set(members))
                partial_claims.extend(claims)
                for claim in claims:
                    partial_cuts.update({end:(bundle.id,claim.role,
                        ('partial',claim.bundle_end) if claim.bundle_end else ('internal',claim.role,end))
                        for end in claim.cuts})
        for approach in approaches:
            if not approach.profile.half:continue
            piece=getattr(approach.path,'mouth_piece',approach.path.pieces[0])
            endpoint=End(piece.section.id,0 if piece.forward else 1)
            record=bundled.get(endpoint.section_id)
            if record is None:
                raise ValueError('Compound half-profile needs its original Bundle membership')
            bundle,role,_=record
            member=bundle.members[role]
            first,last=member.path[0],member.path[-1]
            terminals={End(first.section_id,0 if first.forward else 1),
                       End(last.section_id,1 if last.forward else 0)}
            if not member.complete or endpoint not in terminals:
                if endpoint not in partial_cuts:
                    if bundle.id not in through_checked:
                        from .retained_bundle_through import prepare as prepare_through
                        through_checked.add(bundle.id)
                        through=prepare_through(self.context,bundle,consumed,approaches,set(members))
                        partial_claims.extend(through)
                        through_ports.update({end:claim.token for claim in through for end in claim.cuts})
                    if endpoint not in through_ports:
                        raise ValueError('Deferred compound requires a complete Bundle terminal profile')
                    ports[End(approach.source.id,approach.end)]=through_ports[endpoint]
                    continue
                ports[End(approach.source.id,approach.end)]=partial_cuts[endpoint]
                continue
            terminal=0 if endpoint==End(first.section_id,0 if first.forward else 1) else 1
            ports[End(approach.source.id,approach.end)]=(bundle.id,role,terminal)
        return ports,tuple(partial_claims)
