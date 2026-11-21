"""Actual separated side-junction gaps in original Bundle source spans.

The field may cross an already owned junction for evaluation. That connector is
never a pavement strip. Mesh consumers must retain the real near-side profiles
and curb, and leave the opposite continuous source interval body-owned.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np

from .bundle_axis import AxialChartUnavailable
from .bundle_angular_curves import AngularCrown
from .bundle_angular_internal import shared_frame
from .bundle_angular_open_islands import gap_intervals
from .bundle_ends import PathRange
from .bundle_visible_runs import BundleRun,MemberRun
from .merged_run_inventory import (physical_source_spans,hydrate_physical_spans,sample_sources,
    _SourceBoundary,_inner_bridge,InnerBridge)
from .model import BundleMember,End,Key
from .ownership import PLAN,CompoundProduct,surface_key
from .paired_source_sampling import verify_shared_source_phase
from .source_strip_containment import certify
from .updates import BuildContext


class GappedSourceUnavailable(AxialChartUnavailable):
    pass


@dataclass(frozen=True)
class OwnedSourceSpan:
    role: int
    visits: tuple
    owner: object
    intervals: tuple
    width_proof: dict


@dataclass(frozen=True)
class SeparatedSourceGap:
    role: int
    before: object
    after: object
    opposite: object
    bridge: object
    minimum_gap: float


def complete_run(bundle):
    """A source descriptor only; its endpoints are actual authored Ends."""
    boundaries=(('bundle_end',0),('bundle_end',1))
    return BundleRun(tuple(MemberRun(tuple(member.path),(
        End(member.path[0].section_id,int(not member.path[0].forward)),
        End(member.path[-1].section_id,int(member.path[-1].forward))),boundaries)
        for member in bundle.members),boundaries)


def _owned(context,span):
    plan=context.products.get(PLAN)
    if plan is None:return None
    owners=tuple(plan.owner(end.site) for end in (span.head,span.tail))
    if owners[0] is None or owners[0]!=owners[1]:return None
    owner=owners[0];product=context.products.get(surface_key(owner))
    if not isinstance(product,CompoundProduct) or product.pavement.mesh is None:
        return None
    claimed={claim.section_id for claim in product.group.claims}
    ids={visit.section_id for visit in span.visits}
    if not ids<=claimed:return None
    proof=certify(context,product.group,ids)
    intervals=tuple(dict(section=visit.section_id,forward=visit.forward,first=0.,
        last=context[Key('path',visit.section_id)].length) for visit in span.visits)
    return OwnedSourceSpan(span.role,tuple(span.visits),owner,intervals,proof)


def _bridge(context,role,before,after,frame,owned_between):
    if not owned_between:return _inner_bridge(context,role,before,after,frame)
    a,b=before.tail,after.head
    if (a.kind!='PHYSICAL' or b.kind!='PHYSICAL' or a.owner!=b.owner
            or any(span.owner!=a.owner for span in owned_between)):
        raise GappedSourceUnavailable('Skipped original source intervals need one actual physical curb owner')
    product=context[surface_key(a.owner)]
    if not isinstance(product,CompoundProduct) or product.physical_rims is None:
        raise GappedSourceUnavailable('An owned source gap needs its actual completed compound curb')
    from dataclasses import replace
    physical=replace(product,rims=product.physical_rims,exaggeration=1.)
    xyz=physical.cap((a.endpoint.section_id,a.endpoint.end),(b.endpoint.section_id,b.endpoint.end))
    if not np.array_equal(xyz[0],a.physical_xyz[-1]) or not np.array_equal(xyz[-1],b.physical_xyz[-1]):
        raise ValueError('Owned source curb changed its actual complete profile peaks')
    stations=frame.local(xyz[:,:2])[:,0]
    if np.any(np.diff(stations)<=1e-9):
        raise GappedSourceUnavailable('An owned source curb is not ordered in the certified body frame')
    xyz=xyz.copy();stations=stations.copy();xyz.flags.writeable=stations.flags.writeable=False
    return InnerBridge(role,a.owner,a.owner,a.endpoint,b.endpoint,xyz,stations)


def _limits(span,frame):
    values=[]
    for end in (span.head,span.tail):
        points=end.prepared_mouth.xyz if end.prepared_mouth is not None else end.physical_xyz
        values.append(frame.local(points[:,:2])[:,0])
    return float(values[0].max()),float(values[1].min())


def source_ledger(context,run,spans_by_role,owned_spans):
    """Account ordered original visits, including certified moved-mouth prefixes.

    A real junction's original Route can consume a complete first/last Section
    before exposing its mouth on the next Section. That interval remains the
    junction's pavement; it is neither a new body visit nor a SharedClaim.
    """
    if len(spans_by_role)!=len(run.members):
        raise ValueError('Gapped source ledger lost an original member role')
    all_ids=[v.section_id for member in run.members for v in member.visits]
    if len(all_ids)!=len(set(all_ids)):
        raise ValueError('Gapped source ledger repeats an original source visit')
    prefixes=[];roles=[]
    sites={site.key:site for site in context.network.intersections}
    plan=context.products.get(PLAN)
    for role,(member,spans) in enumerate(zip(run.members,spans_by_role)):
        original=tuple(member.visits);positions={v.section_id:i for i,v in enumerate(original)}
        blocks=[];last_position=-1
        for span in spans:
            lineage=tuple(span.original_visits or span.visits);visible=tuple(span.visits)
            if not lineage or not visible or span.role!=role:
                raise ValueError('Gapped source ledger lost its ordered original visits')
            start=positions.get(lineage[0].section_id,-1)
            if start<0 or original[start:start+len(lineage)]!=lineage or start<=last_position:
                raise ValueError('Gapped source ledger reorders or repeats an original source visit')
            last_position=start+len(lineage)-1;blocks.append((start,lineage))
            at=next((i for i,v in enumerate(lineage) if v==visible[0]),-1)
            if at<0 or lineage[at:at+len(visible)]!=visible:
                raise ValueError('Gapped visible visits changed their original source lineage')
            expected=[]
            for end_index,node in enumerate((span.head,span.tail)):
                if node.original_terminal is None:continue
                first_visit=lineage[0 if end_index==0 else -1]
                terminal=End(first_visit.section_id,int(not first_visit.forward) if end_index==0 else int(first_visit.forward))
                site=sites.get(node.site)
                owner=(plan.owner(node.site) if plan is not None else None) or node.site
                cut_visit=visible[0 if end_index==0 else -1]
                cut=End(cut_visit.section_id,int(not cut_visit.forward) if end_index==0 else int(cut_visit.forward))
                if (node.kind!='PHYSICAL' or node.registry_key is not None
                        or node.original_terminal!=terminal or node.endpoint!=cut
                        or site is None or terminal not in site.ends or node.owner!=owner):
                    raise ValueError('Owned source prefix changed its actual original junction owner or endpoint')
                stop=at if end_index==0 else at+len(visible)-1
                indices=range(stop+1) if end_index==0 else range(stop,len(lineage))
                for index in indices:
                    visit=lineage[index];path=context[Key('path',visit.section_id)]
                    first,last=0.,path.length
                    if index==stop:
                        station=node.source_station
                        if station is None or not 0.<=station<=path.length:
                            raise ValueError('Owned source prefix leaves its original finite Path')
                        if (end_index==0)==visit.forward:last=station
                        else:first=station
                    if last-first>1e-9:expected.append((visit,owner,float(first),float(last)))
            actual=tuple(span.owned_source_intervals)
            if len(actual)!=len(expected):
                raise ValueError('Gapped source ledger has missing or duplicate owned route prefixes')
            remaining=list(expected)
            for claim in actual:
                sid=claim['section']
                if sid not in positions:
                    raise ValueError('Owned source prefix invents an original member Section')
                path=context[Key('path',sid)]
                matches=[i for i,(v,owner,first,last) in enumerate(remaining)
                         if (v.section_id,v.forward,owner)==(sid,claim['forward'],claim['owner'])
                         and abs(first-claim['first'])<=1e-9 and abs(last-claim['last'])<=1e-9]
                proof=claim.get('containment',{});sections=proof.get('sections',())
                certified=[row for row in sections if row['section']==sid
                           and np.max(abs(np.asarray(row['interval'])-(claim['first'],claim['last'])))<=1e-9
                           and row.get('maximum_combined_error',float('inf'))<=.002]
                if (len(matches)!=1 or claim.get('source_path') is not path
                        or claim.get('proof')!='complete_original_approach_route_prefix'
                        or proof.get('maximum_error',float('inf'))>.002 or len(certified)!=1):
                    raise ValueError('Owned source prefix lost its original Path, owner, interval or width proof')
                remaining.pop(matches[0]);prefixes.append(dict(claim,role=role))
            dropped={v.section_id for v in lineage}-{v.section_id for v in visible}
            for sid in dropped:
                full=[row for row in actual if row['section']==sid and row['first']==0.
                      and row['last']==context[Key('path',sid)].length]
                if len(full)!=1:
                    raise ValueError('A missing visible visit needs one complete original owner-prefix proof')
        for span in owned_spans:
            if span.role!=role:continue
            visits=tuple(span.visits);start=positions.get(visits[0].section_id,-1) if visits else -1
            if start<0 or original[start:start+len(visits)]!=visits:
                raise ValueError('Wholly owned source spans changed their original ordered visits')
            blocks.append((start,visits))
        represented=tuple(v for _,visits in sorted(blocks,key=lambda row:row[0]) for v in visits)
        if represented!=original:
            raise ValueError('Gapped source inventory lost or repeated an original source visit')
        roles.append(tuple(dict(section=v.section_id,forward=v.forward) for v in represented))
    return tuple(prefixes),dict(original_ordered_visits=tuple(roles),
        owner_certified_route_prefixes=tuple(dict(section=row['section'],forward=row['forward'],role=row['role'],
            first=row['first'],last=row['last'],owner=row['owner'].text(),
            proof=row['proof'],containment=row['containment']) for row in prefixes),
        original_source_visits_preserved=True)


@dataclass
class GappedSourceInventory:
    bundle: object
    run: object
    frame: object
    field: object
    spans_by_role: tuple
    gaps: tuple
    owned_spans: tuple
    computational_outer_bridges: tuple
    source_sections: tuple
    source_visits: tuple
    source_samples: dict
    proofs: dict
    owned_prefixes: tuple = ()

    @classmethod
    def build(cls,context,bundle,*,run=None,frame=None,max_samples=200000,maximum_gaps=8):
        if (len(bundle.members)!=2 or not all(member.complete for member in bundle.members)
                or type(maximum_gaps) is not int or not 1<=maximum_gaps<=8):
            raise GappedSourceUnavailable('Gapped source inventory needs two complete original members and bounded gaps')
        run=complete_run(bundle) if run is None else run
        frame=shared_frame(context,bundle) if frame is None else frame
        original=physical_source_spans(context,bundle,run,frame=frame)
        active=[];owned=[];positions=[]
        for role,spans in enumerate(original):
            visible=[];indices=[]
            for index,span in enumerate(spans):
                claim=_owned(context,span)
                if claim is None:visible.append(span);indices.append(index)
                else:owned.append(claim)
            if not visible:
                raise GappedSourceUnavailable('A gapped body needs visible original intervals on both roles')
            active.append(tuple(visible));positions.append(indices)
        spans=hydrate_physical_spans(context,bundle,tuple(active),frame=frame,max_samples=max_samples)
        from .gapped_source_validity import certify_visible_centerlines
        centerline_proof=certify_visible_centerlines(context,spans)
        entrances=tuple({id(end.merged_entrance):end.merged_entrance
            for own in spans for span in own for end in (span.head,span.tail)
            if end.merged_entrance is not None}.values())
        samples=sample_sources(context,spans,None,max_samples=max_samples,additional_entrances=entrances)
        bridges=[]
        for role,own in enumerate(spans):
            row=[]
            for i,(before,after) in enumerate(zip(own,own[1:])):
                skipped=original[role][positions[role][i]+1:positions[role][i+1]]
                covered=[claim for claim in owned if claim.role==role
                         and any(claim.visits==tuple(span.visits) for span in skipped)]
                row.append(_bridge(context,role,before,after,frame,covered))
            bridges.append(tuple(row))
        curves=[];outer=[];source_ids=[]
        for role,own in enumerate(spans):
            parts=[[],[]]
            for index,span in enumerate(own):
                products=dict(context.products)
                for at,visit in enumerate(span.visits):
                    path=context[Key('path',visit.section_id)];low,high=0.,path.length
                    # Outer endpoint reserve supports the complete real mouth
                    # plane. Only actual intermediate owner gaps trim the field.
                    for end,trim in ((span.head,index>0 and at==0),
                                     (span.tail,index<len(own)-1 and at==len(span.visits)-1)):
                        if not trim:continue
                        if end.kind!='PHYSICAL':
                            raise GappedSourceUnavailable('An intermediate paired owner needs its prepared paired-source inventory')
                        if end.endpoint.end==0:low=max(low,end.source_station)
                        else:high=min(high,end.source_station)
                    if high<=low+1e-7:raise GappedSourceUnavailable('An actual physical owner exhausts a field-support interval')
                    if low>0. or high<path.length:products[Key('path',visit.section_id)]=PathRange(path,low,high)
                    source_ids.append(visit.section_id)
                boundary=_SourceBoundary(BuildContext(context.network,products,context.patterns),
                    BundleMember(span.visits),frame,samples,join_tolerance=span.route.join_tolerance)
                values=tuple(boundary.points(side) for side in (-1.,1.))
                for end,at,trim in ((span.head,0,index>0),(span.tail,-1,index<len(own)-1)):
                    if not trim:continue
                    for side,node in (((0,0),(1,-1)) if role==0 else ((1,0),(0,-1))):
                        if np.linalg.norm(values[side][at]-end.physical_xyz[node,:2])>1e-8:
                            raise ValueError('Actual source curb and its published profile endpoint disagree')
                        values[side][at]=end.physical_xyz[node,:2]
                inner=1 if role==0 else 0
                for side in (0,1):
                    if parts[side]:
                        bridge=bridges[role][index-1]
                        joined=(bridge.physical_xyz[:,:2] if side==inner else
                                np.array((parts[side][-1][-1],values[side][0])))
                        if side!=inner:
                            outer.append(dict(role=role,owner=bridge.owner.text(),
                                policy='COMPUTATIONAL_ONLY_UNMESHED_JUNCTION_FIELD',xy=joined.tolist()))
                        if not np.array_equal(parts[side][-1][-1],joined[0]) or not np.array_equal(joined[-1],values[side][0]):
                            raise ValueError('An actual source gap lost its complete curb/profile seam')
                        parts[side].append(joined[1:])
                    parts[side].append(values[side] if not parts[side] else values[side][1:])
            curves.extend(np.concatenate(part) for part in parts)
        if sum(map(len,curves))>max_samples:raise GappedSourceUnavailable('Gapped field exceeds its finite original sample budget')
        profiles=[end.port.profile if end.port is not None else end.prepared_mouth.approach.profile
                  for own in spans for span in own for end in (span.head,span.tail)]
        datum,rise=profiles[0].datum,profiles[0].rise
        if any(abs(p.datum-datum)+abs(p.rise-rise)>1e-9 for p in profiles):
            raise GappedSourceUnavailable('Gapped source fields need the original common datum and crown rise')
        field=AngularCrown(frame,curves,datum=datum,rise=rise)
        gaps=[]
        for role,own in enumerate(spans):
            for before,after,bridge in zip(own,own[1:],bridges[role]):
                first,last=map(float,bridge.angular_stations[[0,-1]])
                opposite=[span for span in spans[1-role] if _limits(span,frame)[0]<first
                          and _limits(span,frame)[1]>last]
                if len(opposite)!=1:
                    raise GappedSourceUnavailable('A separated junction gap needs one continuous actual opposite source span')
                regions,roots=gap_intervals(SimpleNamespace(frame=frame,field=field,rays=field.curves),first,last)
                if roots or not regions or any(not positive for _,_,positive in regions):
                    raise GappedSourceUnavailable('An overlapping side junction needs actual paired physical preparation')
                stations=np.unique(np.concatenate((bridge.angular_stations,
                    *(ray.stations[(ray.stations>first)&(ray.stations<last)] for ray in field.curves))))
                _,a,b,_=field.edges(stations);minimum=float(np.min(a-b))
                gaps.append(SeparatedSourceGap(role,before,after,opposite[0],bridge,minimum))
        if not 1<=len(gaps)<=maximum_gaps:
            raise GappedSourceUnavailable('Gapped source construction needs its bounded actual separated junctions')
        gaps.sort(key=lambda gap:float(gap.bridge.angular_stations[0]))
        if any(a.bridge.angular_stations[-1]>=b.bridge.angular_stations[0] for a,b in zip(gaps,gaps[1:])):
            raise GappedSourceUnavailable('Overlapping original junction gaps need a simultaneous source preparation')
        visits=tuple(tuple(v for span in own for v in span.visits) for own in spans)
        prefixes,ledger=source_ledger(context,run,spans,owned)
        proofs=dict(original_source_visits_preserved=True,actual_physical_profiles=True,
            visible_centerline_separation=centerline_proof,
            actual_inner_curbs=True,opposite_source_body_owned=True,new_source_aliases=0,
            original_source_paths=True,phase=verify_shared_source_phase(entrances),
            complete_owned_source_width=[span.width_proof for span in owned],
            original_source_ledger=ledger,
            minimum_separated_gap=min(gap.minimum_gap for gap in gaps))
        return cls(bundle,run,frame,field,spans,tuple(gaps),tuple(owned),tuple(outer),
                   tuple(source_ids),visits,samples,proofs,prefixes)
