"""Finite source-curb envelopes across one wholly owned dimensional link.

This failed-only family keeps explicit original incidence walks.  Contact
curves never represent a changed-width link as a uniform authored road.
Original paths, profiles, endpoints and radius acceptance remain unchanged.
"""
from dataclasses import replace
from types import SimpleNamespace
import math
import numpy as np

from .model import End,Key
from .mesh import cross
from .paths import right
from .routes import RouteTopology
from .continuations import WidthContinuation,ContinuationReachError,_evaluate
from .curb_envelope import CurbEnvelope
from .compound_root_envelopes import RootEnvelope
from .compound_dimensional_corners import DimensionalCompositeCurb,_CertifiedContinuation
from .compound_displaced_contacts import SourceContactWalk,_ordinary,_centerline
from .corner_patch import _PreparedCurve
from .corner_radii import CornerRadiusPolicy
from .free_pads import EnvelopeBoundary
from .intersections import IntersectionProvider
from .regions import Region,RegionUnion,signed_area,_classify
from .compounds import SampledFootprint,RimSegment,rim_label
from .sampling import CertifiedCurbSampler
from .compound_contact_routes import regular_contact


POLICY='FINITE_ORIGINAL_SOURCE_CONTACT_WALK_ENVELOPE'


def source_curvature_vector(approach,station,distance):
    """Analytical geometric offset curvature; reversing a visit keeps its vector."""
    path=approach.path
    if not 0.<=station<=approach.length:raise ValueError('Source contact jet exceeds its finite original approach')
    if getattr(path,'outward',False):
        index=int(np.clip(np.searchsorted(path.stations,station,side='right')-1,0,len(path.pieces)-1))
        piece=path.pieces[index];station-=float(path.stations[index]);path=piece.path
        if not piece.forward:station=path.length-station;distance=-distance
    elif approach.end==1:
        station=path.length-station;distance=-distance
    if path.exact:return np.zeros(2)
    _,velocity,acceleration=path.jet(path.parameter(station))
    speed=float(np.linalg.norm(velocity));tangent=velocity/speed
    curvature=float(cross(velocity,acceleration))/speed**3
    denominator=1.+distance*curvature
    if denominator<=0.:raise ValueError('Source contact jet needs an unfolded original offset')
    return curvature/denominator*np.array((-tangent[1],tangent[0]))


class SourceJetContinuation(WidthContinuation):
    """Quintic Hermite curb with actual geometric endpoint curvature jets."""
    def __init__(self,first,last,start_tangent,end_tangent,minimum_radius,source_curvatures):
        values=tuple(np.asarray(v,float) for v in (first,last,start_tangent,end_tangent,*source_curvatures))
        if (len(values)!=6 or any(v.shape!=(2,) or not np.isfinite(v).all() for v in values)
                or not math.isfinite(minimum_radius) or minimum_radius<=0.):
            raise ValueError('Invalid finite source-curvature Hermite request')
        first,last,a,b,ka,kb=values
        if min(np.linalg.norm(a),np.linalg.norm(b))<=1e-12:raise ValueError('Collapsed source tangent')
        a=a/np.linalg.norm(a);b=b/np.linalg.norm(b);axis=a+b
        if np.linalg.norm(axis)<=1e-8:raise ContinuationReachError('Source contact jets turn back')
        axis/=np.linalg.norm(axis);length=float((last-first)@axis)
        if length<=1e-8:raise ContinuationReachError('Source contact jets need positive finite reach')
        p1=first+length*a/5.;p4=last-length*b/5.
        self.controls=np.array((first,p1,2.*p1-first+length**2*ka/20.,
            2.*p4-last+length**2*kb/20.,p4,last))
        self.first_derivative=5.*np.diff(self.controls,axis=0)
        self.second_derivative=20.*np.diff(self.controls,n=2,axis=0)
        if np.min(self.first_derivative@axis)<=1e-8:
            raise ContinuationReachError('Source contact jets need forward finite reach')
        numerator=np.zeros(8)
        for i,velocity in enumerate(self.first_derivative):
            for j,acceleration in enumerate(self.second_derivative):
                numerator[i+j]+=math.comb(4,i)*math.comb(3,j)/math.comb(7,i+j)*cross(velocity,acceleration)
        self.curvature_bound,self.curvature_intervals=self._refined_curvature(numerator,minimum_radius)
        self.maximum_speed=float(np.linalg.norm(self.first_derivative,axis=1).max())
        self.maximum_second=float(np.linalg.norm(self.second_derivative,axis=1).max())
        nodes,weights=np.polynomial.legendre.leggauss(8)
        self.length=float(.5*np.sum(weights*np.linalg.norm(_evaluate(self.first_derivative,.5*(nodes+1.)),axis=1)))
        for array in (self.controls,self.first_derivative,self.second_derivative):array.flags.writeable=False


def contact_patch(owner,corner,curve,stations,settings,walks,maximum=8192):
    """Close a certified curb only through the two actual source walks."""
    chains=[_centerline(w,s,min(.3,settings.spacing),settings.curve_tolerance,maximum)
            for w,s in zip(walks,stations)]
    if len(chains)!=2 or np.linalg.norm(chains[0][0]-chains[1][0])>1e-8:
        raise ValueError('Source contact envelope lacks its original common incidence')
    parameters=CertifiedCurbSampler(settings.curve_tolerance,settings.spacing,maximum).sample(curve,0.,1.)
    curb=curve(parameters);curve.contact_reference_points=curb.copy()
    curves=[];points=[];labels=[]
    for chain,kind,actual,values in ((np.vstack((chains[0],curb[0])),'internal',None,None),
            (curb,'curb',curve,parameters),(np.vstack((curb[-1],chains[1][::-1])),'internal',None,None)):
        keep=np.r_[True,np.linalg.norm(np.diff(chain,axis=0),axis=1)>1e-8];chain=chain[keep]
        if len(chain)<2:continue
        actual=actual if actual is not None else _PreparedCurve(chain)
        values=values[keep] if values is not None else actual.stations/actual.length
        index=len(curves);curves.append(actual);points.extend(chain[:-1])
        labels.extend(RimSegment(kind,index,float(a),float(b)) for a,b in zip(values,values[1:]))
    points=np.asarray(points)
    if len(points)>maximum or signed_area(points)<=0.:
        raise ValueError('Source contact envelope needs a simple bounded oriented patch')
    region=Region(owner,points);RegionUnion().build((region,))
    layout=SimpleNamespace(approaches=(),outward_stations=(),boundaries=tuple(curves),rows=len(points),mouths={},corners=(corner,))
    return SampledFootprint(region,tuple(labels),layout)


def phase_parameters(source,corner,stations,*,maximum_edges=16):
    """Add evaluation cuts on original paths; retain every existing source cut.

    Matching C2 curves separate cubically near a contact.  Independent chord
    phases can otherwise hide a short new curb in a coarser source chord.
    Project only the first/last sixteen prepared curb edges onto their own
    original source approach; these cuts make the source approximation finer.
    """
    if type(maximum_edges) is not int or not 1<=maximum_edges<=16:
        raise ValueError('Source contact phase exceeds its finite edge budget')
    edges=[i for i,label in enumerate(source.segments) if label.kind=='curb'];result={}
    for side,approach in enumerate((corner.first,corner.second)):
        selected=edges[:maximum_edges] if side==0 else edges[-maximum_edges:]
        vertices=sorted({i for edge in selected for i in (edge,(edge+1)%len(source.region.xy))})
        station=float(stations[side]);offset=(-.5 if side==0 else .5)*approach.width
        for vertex in vertices:
            point=source.region.xy[vertex];current=station
            for _ in range(16):
                position,velocity=approach.jet(current,offset)
                denominator=float(velocity@velocity)
                if denominator<=1e-16:raise ValueError('Source phase has a collapsed original derivative')
                current=float(np.clip(current+(point-position)@velocity/denominator,0.,station))
            for (sid,end),cut in approach.path.trims(current).items():result.setdefault(sid,set()).add(float(cut))
        for (sid,end),cut in approach.path.trims(station).items():result.setdefault(sid,set()).add(float(cut))
    return {sid:tuple(sorted(values)) for sid,values in result.items()}


def _outer_half(context,approach,original,side,members):
    from .profiles import memberships
    if (not approach.profile.half or approach.profile.rise!=original.profile.rise
            or approach.profile.datum!=original.profile.datum
            or abs(float(approach.profile.value(side)))>1e-10):return False
    membership=memberships(context.network).get(approach.source.id)
    if membership is None:return False
    bundle,role,*_=membership
    if len(bundle.members)!=2 or not all(member.complete for member in bundle.members):return False
    opposite={v.section_id for v in bundle.members[1-role].path}
    from .ownership import footprint_key
    for key in members:
        layout=context.products.get(footprint_key(key))
        if layout is not None and any(a.source.id in opposite and a.profile.half for a in layout.approaches):return True
    return False


def _walks(context,original,side,members,consumed):
    """All real incidences define order; only compatible selected curbs are used."""
    from .ownership import footprint_key
    topology=RouteTopology(context.network);sites={j.key:j for j in context.network.intersections}
    result=[SourceContactWalk(original,original)];prefix=[];seen=set();current=original
    for _ in range(8):
        if not _ordinary(current) or any(p.section.id not in consumed or p.section.id in seen for p in current.path.pieces):break
        route=topology.trace(End(current.source.id,current.end))
        if (route.terminal_site not in members or route.terminal_site not in sites
                or route.terminal_scope or route.terminal_selection
                or tuple((p.section.id,p.forward) for p in current.path.pieces)
                   !=tuple((v.section_id,v.forward) for v in route.traversals)):break
        last=route.traversals[-1];reverse=End(last.section_id,1 if last.forward else 0)
        layout=context.products.get(footprint_key(route.terminal_site))
        if layout is None:break
        ordered=sorted(layout.approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
        ends=[End(a.source.id,a.end) for a in ordered]
        if len(ends)!=len(set(ends)) or set(ends)!=set(sites[route.terminal_site].ends) or reverse not in ends:break
        following=ordered[(ends.index(reverse)+(1 if side>0 else -1))%len(ordered)]
        local_first=reverse if side>0 else End(following.source.id,following.end)
        if sites[route.terminal_site].corner_radius(local_first)[1]!='STYLE':break
        half_ok=_outer_half(context,following,original,side,members)
        if (not (_ordinary(following) or half_ok) or following.source.id in seen
                or following.source.id in {p.section.id for p in current.path.pieces}
                or not (following.profile==original.profile or half_ok)
                or np.linalg.norm(current.frame(current.length)[0]-following.origin)>1e-8):break
        prefix.extend(current.path.pieces);seen.update(p.section.id for p in current.path.pieces);current=following
        result.append(SourceContactWalk(original,current,tuple(prefix)))
    return tuple(result)


def _pair(corner):return tuple(End(a.source.id,a.end) for a in (corner.first,corner.second))


def _caps(loops,sources):
    inventory={s.region.owner:s for s in sources};result=[]
    for loop in loops:
        for span in loop.spans:
            if rim_label(span,inventory).kind!='internal':continue
            source=inventory[span.owner];end=getattr(source.layout,'original_root_caps',{}).get(span.segment)
            if end is None:return None
            xy=source.region.xy;a=xy[span.segment];d=xy[(span.segment+1)%len(xy)]-a
            result.append((end,a+np.array((span.first,span.last))[:,None]*d))
    return result


def _strip(owner,path,at):
    from .compound_cycles import _footprint
    section=path.section
    points=np.vstack((path.point(at,section.offset+.5*section.width),path.point(at,section.offset-.5*section.width)[::-1]))
    source=_footprint(owner,points,internal_edges=(len(at)-1,2*len(at)-1))
    if source is None:raise ValueError('Source contact proof has no original strip area')
    caps={tuple(sorted((tuple(points[i]),tuple(points[(i+1)%len(points)])))):End(section.id,end)
          for i,end in ((len(at)-1,1),(2*len(at)-1,0))}
    source.layout.original_root_caps={i:caps[token] for i,(a,b) in enumerate(zip(source.region.xy,np.roll(source.region.xy,-1,axis=0)))
        if (token:=tuple(sorted((tuple(a),tuple(b))))) in caps}
    return source


def _operands(planner,sources,corners,consumed,replaced,candidate,settings):
    from .compound_cycles import _canonical_samples,_port_footprint,_sampling_intervals
    sources=tuple(s for s in sources if not any(_pair(c)==replaced for c in s.layout.corners))
    sources=(*sources,candidate.source);paths={sid:planner.context[Key('path',sid)] for sid in consumed};cuts={}
    contacts={}
    for source in sources:
        for a in source.layout.approaches:
            if not _ordinary(a) and not a.profile.half:raise ValueError('Source contact proof needs original finite ports')
            for piece in a.path.pieces:paths[piece.section.id]=piece.path
        for c in source.layout.corners:
            for a,distance in zip((c.first,c.second),c.stations):
                if distance<0.:continue
                end=End(a.source.id,a.end);contacts[end]=max(contacts.get(end,0.),distance)
                for (sid,side),at in a.path.trims(distance).items():cuts.setdefault(sid,set()).add(float(at))
            record=getattr(c,'source_contact_envelope',None)
            if record is not None:
                for sid,values in record.report['source_phase_parameters'].items():cuts.setdefault(sid,set()).update(values)
    port_stations={}
    for source in sources:
        if not source.layout.approaches:continue
        if len(source.layout.approaches)!=1:raise ValueError('Source contact proof cannot replace a merged physical port')
        a=source.layout.approaches[0];station=source.layout.outward_stations[0]
        station=max(station,contacts.get(End(a.source.id,a.end),0.)+min(.25,.05*a.width))
        if not IntersectionProvider.mouth_fits(a,station):raise ValueError('Source contact proof exhausts an original finite port')
        port_stations[source.region.owner]=station
        for (sid,side),at in a.path.trims(station).items():cuts.setdefault(sid,set()).add(float(at))
    ports=[source for source in sources if source.layout.approaches]
    partitions=_canonical_samples(paths,cuts,settings.spacing,settings.curve_tolerance,
        intervals=_sampling_intervals(paths,consumed,[s.layout.approaches[0] for s in ports],
            [port_stations[s.region.owner] for s in ports]))
    rebuilt=[]
    for source in sources:
        roots=getattr(source.layout,'original_root_caps',{})
        if source.region.owner.kind=='cycle_strip':
            sid=source.region.owner.id
            if sid not in consumed:raise ValueError('Source contact proof cannot add an unclaimed whole strip')
            source=_strip(source.region.owner,paths[sid],partitions[sid])
        elif source.layout.approaches:
            a=source.layout.approaches[0]
            source=_port_footprint(source.region.owner,a,port_stations[source.region.owner],source.layout.radial,
                settings.spacing,settings.curve_tolerance,partitions)
        rebuilt.append(source)
    if sum(len(s.region.xy) for s in rebuilt)>20000:raise ValueError('Source contact envelope exceeds the source sampling budget')
    return tuple(rebuilt)


def _audit(planner,sources,corners,consumed,replaced,candidate,targets,settings,holes):
    operands=_operands(planner,sources,corners,consumed,replaced,candidate,settings)
    xy=candidate.source.region.xy;edge_ids=[i for i,label in enumerate(candidate.source.segments) if label.kind=='curb']
    probe_ids=sorted(set(edge_ids[:16]+edge_ids[-16:]))
    probes=np.array([.5*(xy[i]+xy[(i+1)%len(xy)]) for i in probe_ids])
    for source in operands:
        if source.region.owner==candidate.source.region.owner:continue
        inside,boundary=_classify(probes,source.region.xy,1e-8)
        if np.any(inside&(boundary<0)):return False
    loops=RegionUnion(max_segments=20000).build(s.region for s in operands)
    if sum(not loop.hole for loop in loops)!=1 or sum(loop.hole for loop in loops)!=holes:return False
    amount=0.
    for loop in loops:
        for cap in targets:
            axis=cap[1]-cap[0];length=float(np.linalg.norm(axis));delta=loop.xy-cap[0]
            close=abs(cross(delta,axis))/length<=1e-8;t=delta@axis/(axis@axis)
            low=np.maximum(0.,np.minimum(t,np.roll(t,-1)));high=np.minimum(1.,np.maximum(t,np.roll(t,-1)))
            if np.any(close&np.roll(close,-1)&((high-low)*length>1e-7)):return False
        if loop.hole:continue
        for span in loop.spans:
            if span.owner!=candidate.source.region.owner:continue
            label=candidate.source.segments[span.segment]
            if label.kind=='curb':amount+=abs((label.last-label.first)*(span.last-span.first))
    return abs(amount-1.)<=1e-8


def _fit(planner,sources,corners,consumed,original_pair,walks,template,owner,targets,root_points,settings,holes):
    site=next(j for j in planner.sites.values() if original_pair[0] in j.ends)
    preferred,mode=site.corner_radius(original_pair[0])
    if mode!='STYLE':return None
    attempts=[];a,b=(w.approach for w in walks);width=min(a.width,b.width)
    if (any(abs(x.width-x.source.width)>1e-8 or abs(x.offset)>1e-8 for x in (a,b))
            or abs(a.profile.datum-b.profile.datum)>1e-10):return None
    gap=math.atan2(float(cross(a.direction,b.direction)),float(a.direction@b.direction))%math.tau
    for radius in CornerRadiusPolicy.candidates(preferred,mode):
        attempts.append(float(radius));increment=max(.125,max(radius,.5*max(x.width for x in (a,b)))/8.)
        for actual_jets in (False,True):
            for step in range(17):
                stations=(.25+step*increment,)*2
                if not all(IntersectionProvider.mouth_fits(x,s) for x,s in zip((a,b),stations)):break
                if not all(regular_contact(x,s) for x,s in zip((a,b),stations)):continue
                first,va=a.jet(stations[0],-.5*a.width);last,vb=b.jet(stations[1],.5*b.width)
                try:
                    curvatures=(source_curvature_vector(a,stations[0],-.5*a.width),source_curvature_vector(b,stations[1],.5*b.width))
                    if not actual_jets and (float(curvatures[0]@(-a.frame(stations[0])[2]))>1e-10
                            or float(curvatures[1]@b.frame(stations[1])[2])>1e-10):
                        # A zero-curvature curb would immediately enter this
                        # outward-bending original ribbon at its contact.
                        continue
                    if actual_jets:
                        base=SourceJetContinuation(first,last,-va,vb,radius,curvatures)
                    else:base=WidthContinuation(first,last,-va,vb,radius,refine_curvature=True)
                except ContinuationReachError:continue
                outward=right(last-first);outward/=np.linalg.norm(outward);axis=(last-first)/np.linalg.norm(last-first)
                required=0.
                for point in root_points:
                    target=float((point-first)@axis)
                    if not 0.<target<float((last-first)@axis):continue
                    low,high=0.,1.
                    for _ in range(60):
                        middle=.5*(low+high)
                        if float((_evaluate(base.controls,middle)-first)@axis)<target:low=middle
                        else:high=middle
                    t=.5*(low+high);required=max(required,float((point-_evaluate(base.controls,t))@outward)/(64.*t**3*(1.-t)**3))
                margin=settings.curve_tolerance
                distances=sorted({0.,required+margin,required+2.*margin,required+4.*margin,
                    required+8.*margin,required*1.01+margin,required*1.05+margin,required*1.1+margin})
                for distance in distances:
                    if not 0.<=distance<=.25*width:continue
                    try:
                        envelope=CurbEnvelope(base,outward,distance,radius)
                        corner=replace(template,first=a,second=b,gap=gap,center=None,radius=radius,
                            requested_radius=preferred,radius_mode=mode,radius_attempts=tuple(attempts),
                            continuation=True,stations=stations,fit_policy=POLICY)
                        curve=DimensionalCompositeCurb(corner,(_CertifiedContinuation(EnvelopeBoundary(corner,envelope)),),settings.curve_tolerance)
                        source=contact_patch(owner,corner,curve,stations,settings,walks)
                        report=dict(policy=POLICY,original_pair=[dict(section=e.section_id,end=e.end) for e in original_pair],
                            original_junction=site.key.text(),requested_radius=preferred,accepted_radius=radius,
                            radius_attempts=list(attempts),stations=list(stations),source_widths=[a.width,b.width],
                            source_profile_halves=[a.profile.half,b.profile.half],
                            original_facing_profile_derivatives=[float(a.profile.derivative(-1.)),float(b.profile.derivative(1.))],
                            source_contact_walks=[w.report(s) for w,s in zip(walks,stations)],
                            actual_source_curvature_jets=actual_jets,curvature_bound=envelope.curvature_bound,
                            outward_envelope_m=distance,maximum_outward_envelope_m=.25*width,
                            reach_steps=step,maximum_reach_steps=17,maximum_candidates_per_reach=8,
                            original_sources_unchanged=True,endpoint_extrapolation=False,
                            source_phase_parameters=phase_parameters(source,corner,stations),
                            source_phase_maximum_edges_per_contact=16)
                        record=RootEnvelope(source,curve,stations,report,corner)
                        corner.source_contact_envelope=record;corner.source_contact_walks=walks
                        corner.source_contact_support=report
                        if _audit(planner,sources,corners,consumed,original_pair,record,targets,settings,holes):
                            report.update(whole_curb_exposed=True,actual_root_caps_internal=True,original_hole_count=holes)
                            return record
                    except ContinuationReachError:continue
                    except ValueError:continue  # A candidate must still pass the exact arrangement proof.
    return None


def prepare_required(planner,loops,sources,corners,unresolved,consumed,members,settings):
    """Two real curb walks may support a typed exhausted narrow owned link."""
    exposed=_caps(loops,sources)
    if not exposed:return None
    if any(hasattr(a,'original_approaches') for s in sources for a in s.layout.approaches):return None
    by_end={End(a.source.id,a.end):a for key in members for a in planner.layouts[key].approaches}
    holes=sum(loop.hole for loop in loops)
    for width_corner in corners:
        record=getattr(width_corner,'root_envelope',None)
        if record is None or record.report.get('admission')!='TYPED_FINITE_NARROW_WIDTH_LEAD_EXHAUSTION':continue
        link=record.report['exhausted_original_section']
        if link not in consumed:continue
        link_end=next((e for e in _pair(width_corner) if e.section_id==link),None)
        if link_end is None:continue
        original=by_end.get(link_end);opposite=End(link,1-link_end.end)
        if original is None or len(original.path.pieces)!=1 or opposite not in by_end:continue
        for first,last,failure in unresolved:
            upper_pair=(first,last)
            if opposite not in upper_pair or any(by_end[e].profile.half for e in upper_pair):continue
            upper_index=upper_pair.index(opposite)
            upper_walks=[SourceContactWalk(by_end[e],by_end[e]) for e in upper_pair]
            wide=next(e for e in _pair(width_corner) if e!=link_end)
            options=_walks(planner.context,by_end[opposite],-1 if upper_index==0 else 1,set(members),consumed)
            extension=next((w for w in options if tuple(p.section.id for p in w.prefix)==(link,)
                and End(w.approach.source.id,w.approach.end)==wide),None)
            if extension is None:continue
            upper_walks[upper_index]=extension
            if min(w.approach.width for w in upper_walks)<=original.width+1e-6:continue
            upper_targets=[];lower_targets=[]
            for end,cap in exposed:
                selected=False
                for side,w in enumerate(upper_walks):
                    if end!=End(w.approach.source.id,w.approach.end):continue
                    a=w.approach;low=a.point(0.,-.5*a.width);axis=a.point(0.,.5*a.width)-low
                    u=(cap-low)@axis/(axis@axis)
                    if (side==0 and u.max()<=.5+1e-8) or (side==1 and u.min()>=.5-1e-8):selected=True
                (upper_targets if selected else lower_targets).append(cap)
            if not upper_targets or not lower_targets:continue
            donor=next((c for c in corners if c is not width_corner and link_end in _pair(c)
                and c.center is not None and not c.first.profile.half and not c.second.profile.half),None)
            if donor is None:continue
            donor_pair=_pair(donor);donor_index=donor_pair.index(link_end)
            lower_walks=[SourceContactWalk(a,a) for a in (donor.first,donor.second)]
            lower_options=_walks(planner.context,by_end[link_end],-1 if donor_index==0 else 1,set(members),consumed)
            lower_extension=next((w for w in lower_options if tuple(p.section.id for p in w.prefix)==(link,)
                and End(w.approach.source.id,w.approach.end) not in upper_pair),None)
            if lower_extension is None:continue
            lower_walks[donor_index]=lower_extension
            donor_source=next((s for s in sources if any(c is donor for c in s.layout.corners)),None)
            if donor_source is None:raise ValueError('Source contact envelope lost its actual original donor')
            lower=_fit(planner,sources,corners,consumed,donor_pair,tuple(lower_walks),donor,donor_source.region.owner,
                lower_targets,tuple(p for cap in lower_targets for p in cap),settings,holes)
            if lower is None:continue
            remaining=tuple(s for s in sources if s is not donor_source)+(lower.source,)
            points=[p for cap in upper_targets for p in cap]
            points.extend(record.source.region.xy[i] for i,label in enumerate(record.source.segments) if label.kind=='curb')
            serial=1+max((s.region.owner.id for s in sources if s.region.owner.kind=='cycle_corner'),default=0)
            upper=_fit(planner,remaining,corners,consumed,upper_pair,tuple(upper_walks),width_corner,Key('cycle_corner',serial),
                tuple(cap for end,cap in exposed),tuple(points),settings,holes)
            if upper is None:continue
            marker=wide
            return (((donor_pair,marker),lower),((upper_pair,marker),upper))
    return None
