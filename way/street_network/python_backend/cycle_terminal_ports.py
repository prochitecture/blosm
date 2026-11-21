"""Finite original-strip unions at a complete paired Bundle terminal.

The two source incidences remain distinct. Only their physical incoming
profile is shared on one transverse fibre of a certified finite source chart.
No region between disjoint source strips is added.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np

from .model import End, Key
from .regions import Region, RegionUnion, signed_area
from .compounds import RimSegment, SampledFootprint, rim_label
from .bundle_visible_runs import BundleRun, MemberRun


@dataclass(frozen=True)
class TerminalPair:
    bundle_id: int
    run: object
    end: int
    originals: tuple

    @property
    def endpoints(self):
        return tuple(member.ends[self.end] for member in self.run.members)

    def build(self, context, cuts):
        settings=context.network.settings
        tolerance=min(.00025,.0125*settings.curve_tolerance)
        parameters={}
        for endpoint in self.endpoints:
            parameters[endpoint.section_id]=tuple(sorted(cuts.get(endpoint.section_id,())))
        from .finite_terminal_entrances import FinitePairedEntrance
        result=FinitePairedEntrance(context,self.run,self.end,
            source_tolerance=tolerance,source_parameters=parameters)
        result.bundle_id=self.bundle_id
        result.original_approaches=self.originals
        result.source_curve_tolerance=tolerance
        return result


@dataclass(frozen=True)
class OriginalSectionGap:
    """An actual separated fibre on both finite original terminal Sections."""
    bundle_id: int
    end: int
    station: float
    gap: float
    allowance: float
    margin: float
    sources: tuple

    def report(self):
        return dict(policy='REQUIRED_FINITE_FIBRE_HAS_SEPARATED_ORIGINAL_SOURCE_CURBS',
            bundle=self.bundle_id,bundle_end=self.end,station=self.station,gap=self.gap,
            combined_source_curve_allowance=self.allowance,numerical_margin=self.margin,
            sources=self.sources,original_separate_profiles_retained=True)


def original_section_gap(approach,station):
    """Invert the four ORIGINAL curbs; sampled fields only propose this cut.

    A merged minimum cut does not establish that a farther corner-required
    cut still overlaps. Declining that physical merge requires a clear gap
    beyond both source sampling allowances, on complete finite source fibres.
    """
    from .bundle_axis import AxialChartUnavailable
    if not minimum_station(approach)<=station<approach.length-1e-6:
        raise AxialChartUnavailable('Required terminal cut leaves its finite original source support')
    frame=approach.frame_chart;target=float(approach.chart_station(station))
    scale=max(1.,abs(target),float(abs(frame.origin).max()))
    margin=max(1e-7,512.*np.finfo(float).eps*scale)
    records=[];edges=[]
    for path,visit,endpoint,interval in zip(approach.paths,approach.visits,
            approach.original_ends,approach.source_ranges):
        first,last=map(float,interval);direction=1. if visit.forward else -1.
        if not 0.<=first<last<=path.length:
            raise ValueError('Terminal separation lost its finite original source interval')
        lower=None;regularity=None
        if not path.exact:
            from .bundles import CenterlineRange
            from .bundle_boundaries import _certify_increasing
            regularity=1.-(abs(path.section.offset)+.5*path.section.width)*path.curvature_bound
            if regularity<=0.:
                raise AxialChartUnavailable('Original terminal offset curb lacks a regularity certificate')
            lower=min(_certify_increasing(3.*np.diff(curve,axis=0)@frame.axis*direction,
                coordinate='original terminal separation axial')
                for curve in CenterlineRange.build(path,first,last).curves)
        curbs=[]
        for side in (-1.,1.):
            offset=path.section.offset+side*direction*.5*path.section.width
            if path.exact:
                knots=np.r_[first,path.stations[(path.stations>first)&(path.stations<last)],last]
                coordinates=frame.local(path.point(knots,offset))[:,0]
                if np.any(direction*np.diff(coordinates)<=1e-9):
                    raise AxialChartUnavailable('Original terminal offset curb is not axially monotone')
            endpoints=frame.local(path.point(np.array((first,last)),offset))[:,0]
            if not min(endpoints)+margin<target<max(endpoints)-margin:
                raise AxialChartUnavailable('Required terminal fibre leaves a complete original transverse profile')
            low,high=first,last
            for _ in range(64):
                middle=.5*(low+high)
                if middle==low or middle==high:break
                value=float(frame.local(path.point(middle,offset))[0])
                if direction*(value-target)<0.:low=middle
                else:high=middle
            at=.5*(low+high);point=path.point(at,offset);local=frame.local(point)
            residual=float(abs(local[0]-target))
            if residual>margin:
                raise ValueError('Original terminal offset-curb inversion did not converge')
            edges.append(float(local[1]))
            curbs.append(dict(original_station=at,offset=offset,point=point.tolist(),
                transverse=float(local[1]),inversion_residual=residual))
        records.append(dict(section=path.section.id,end=endpoint.end,forward=visit.forward,
            original_interval=[first,last],axial_derivative_lower=lower,
            offset_tangent_scale_lower=regularity,curbs=curbs))
    ao,ai,bi,bo=edges
    if not ao>ai or not bi>bo:
        raise ValueError('Original terminal transverse profiles lost their positive widths')
    gap=max(ai,bo)-min(ao,bi)
    allowance=2.*approach.source_curve_tolerance
    if gap<=allowance+margin:return None
    return OriginalSectionGap(approach.bundle_id,approach.end_index,float(station),
        float(gap),float(allowance),margin,tuple(records))


def prepare_requested(pair,context,cuts,stations):
    """Keep old merged success first; otherwise require analytic separation."""
    from .bundle_axis import AxialChartUnavailable
    approach=pair.build(context,cuts);station=requested_station(approach,stations)
    try:approach.cross_section(station)
    except AxialChartUnavailable:
        certificate=original_section_gap(approach,station)
        if certificate is None:raise
        return None,station,certificate
    return approach,station,None


def discover(planner,members,consumed,approaches):
    """Use only both declared, complete, original external terminal owners."""
    from .bundle_axis import AxialChartUnavailable
    context=planner.context
    by_end={End(a.source.id,a.end):a for a in approaches}
    owners={end:site for site in context.network.intersections for end in site.ends}
    pairs=[];used=set()
    for bundle in context.network.bundles:
        if len(bundle.members)!=2 or not all(member.complete for member in bundle.members):continue
        runs=[]
        for member in bundle.members:
            first,last=member.path[0],member.path[-1]
            endpoints=(End(first.section_id,int(not first.forward)),End(last.section_id,int(last.forward)))
            tokens=(('bundle_end',0),('bundle_end',1))
            runs.append(MemberRun(member.path,endpoints,tokens))
        run=BundleRun(tuple(runs),tokens)
        for end in (0,1):
            endpoints=tuple(member.ends[end] for member in runs)
            if any(e not in by_end or e.section_id in consumed or e in used for e in endpoints):continue
            if any(e not in owners or owners[e].key not in members for e in endpoints):continue
            originals=tuple(by_end[e] for e in endpoints)
            if any(e.end in context.network.section_index[e.section_id].scope_ends
                    or e.end in context.network.section_index[e.section_id].selection_ends for e in endpoints):continue
            if any(getattr(a.path,'absorbed_sites',()) or hasattr(a.path,'mouth_piece')
                    or hasattr(a.path,'deferred_profile_reserve') or abs(a.offset)>1e-8
                    or abs(a.width-context.network.section_index[a.source.id].width)>1e-8 for a in originals):continue
            # A terminal may have a finite source Scope at its remote end, but
            # the selected incidence itself must be an actual junction owner.
            try:
                pair=TerminalPair(bundle.id,run,end,originals)
                candidate=pair.build(context,{})
            except AxialChartUnavailable:continue
            minimum=minimum_station(candidate)
            if minimum>=candidate.length-1e-6:continue
            try:candidate.cross_section(minimum)
            except AxialChartUnavailable:continue
            pairs.append(pair);used.update(endpoints)
    return tuple(pairs)


def minimum_station(approach):
    """The aliases must intersect both complete finite source centrelines.

    A staggered original root cap can cross the common fibre. Its retained
    portion stays internal; the remainder belongs to the finite Bundle body.
    The footprint separately checks the *actual* union profile endpoints.
    """
    if hasattr(approach,'minimum_station'):return approach.minimum_station
    origins=[path.frame(0. if end.end==0 else path.length)[0]
             for path,end in zip(approach.paths,approach.original_ends)]
    coordinates=approach.frame_chart.local(origins)[:,0]
    value=np.max(approach.sign*(coordinates-approach.origin_station))
    return max(0.,float(value))+1e-5


def requested_station(approach,stations):
    value=minimum_station(approach)
    for original in approach.original_approaches:
        station=stations[End(original.source.id,original.end)]
        xy=original.point(np.asarray(station),original.offset+.5*original.width*np.array((-1.,1.)))
        local=approach.frame_chart.local(xy)[:,0]
        value=max(value,float(np.max(approach.sign*(local-approach.origin_station))))
    return value


def _clip(points,center,normal):
    """Exact finite half-plane clipping of an authored curb-chord polygon."""
    result=[]
    for a,b in zip(points,np.roll(points,-1,axis=0)):
        da,db=float((a-center)@normal),float((b-center)@normal)
        if da<=0.:result.append(a)
        if (da<0.<db) or (db<0.<da):result.append(a+(b-a)*(da/(da-db)))
    values=np.asarray(result,float)
    if len(values)<3:raise ValueError('A terminal fibre removes an entire finite source strip')
    keep=np.r_[True,np.linalg.norm(np.diff(values,axis=0),axis=1)>1e-8]
    values=values[keep]
    if np.linalg.norm(values[0]-values[-1])<=1e-8:values=values[:-1]
    if signed_area(values)<0.:values=values[::-1]
    return values


def footprint(owner,approach,station,radial):
    """Clip and union the two original strips; expose one complete profile."""
    from .compound_cycles import PolylineCurb
    if not minimum_station(approach)<=station<approach.length-1e-6:
        raise ValueError('Merged terminal profile is outside its finite original curb support')
    approach.profile=approach.cross_section(station)
    aliases=approach.source_mouths(station)
    if hasattr(approach,'cut_plane'):
        center,normal=approach.cut_plane(station)
    else:
        center=approach.frame_chart.center
        ray=approach.point(station,0.)-center;ray/=np.linalg.norm(ray)
        normal=approach.sign*approach.frame_chart.direction*np.array((-ray[1],ray[0]))
    from .terminal_curb_ancestry import original_strip,clip as clip_ancestry,validate as validate_ancestry
    operands=[];source_polygons=[];original_banks={}
    for role in (0,1):
        first,last=approach.field.curves[2*role:2*role+2]
        points=np.vstack((first.points,last.points[::-1]))
        clipped=_clip(points,center,normal)
        original,tags=original_strip(approach,role)
        if not np.array_equal(original,points):
            raise ValueError('Terminal curb ancestry changed its actual original bank points')
        traced,tags=clip_ancestry(original,tags,center,normal)
        if not np.array_equal(traced,clipped):
            raise ValueError('Terminal curb ancestry changed the existing clipped source geometry')
        # Retain the exact finite portion of the original root cap.
        cap=np.array((first.points[0 if approach.sign>0. else -1],
                      last.points[0 if approach.sign>0. else -1]))
        labels=[]
        for a,b in zip(clipped,np.roll(clipped,-1,axis=0)):
            if max(abs((a-center)@normal),abs((b-center)@normal))<1e-7:kind='mouth'
            elif _on_segment(a,cap) and _on_segment(b,cap):kind='internal'
            else:kind='curb'
            labels.append(RimSegment(kind,0,0.,1.))
        source=SampledFootprint(Region(Key('terminal_strip',2*approach.bundle_id+role),clipped),tuple(labels))
        original_banks[source.region.owner]=tags
        operands.append(source);source_polygons.append(clipped)
    loops=RegionUnion().build(source.region for source in operands)
    if len(loops)!=1 or loops[0].hole:
        raise ValueError('Merged terminal prefixes need one actual overlapping strip union')
    inventory={source.region.owner:source for source in operands}
    loop=loops[0];kinds=[rim_label(span,inventory).kind for span in loop.spans]
    mouth=[i for i,kind in enumerate(kinds) if kind=='mouth']
    if not mouth:raise ValueError('Merged terminal union has no physical incoming profile')
    starts=[i for i in mouth if kinds[i-1]!='mouth']
    if len(starts)!=1:raise ValueError('Merged terminal source strips do not share one complete profile')
    start=starts[0];order=np.roll(np.arange(len(loop.xy)),-start)
    points=loop.xy[order];kinds=[kinds[i] for i in order]
    count=next((i for i,kind in enumerate(kinds) if kind!='mouth'),len(kinds))
    u=approach.mouth_parameters(radial,approach.width/(2*radial))
    profile=approach.point(station,approach.width*(u-.5))
    if max(np.linalg.norm(points[0]-profile[0]),np.linalg.norm(points[count]-profile[-1]))>1e-7:
        raise ValueError('Merged terminal union disagrees with its exact physical profile endpoints')
    rim=np.vstack((profile[:-1],points[count:]))
    labels=[RimSegment('mouth',0,float(a),float(b)) for a,b in zip(u,u[1:])]
    curves=[];caps={};curb_ancestry={}
    for i in range(count,len(points)):
        curve=PolylineCurb(np.array((points[i],points[(i+1)%len(points)])))
        index=len(curves);curves.append(curve)
        labels.append(RimSegment(kinds[i],index,0.,1.))
        if kinds[i]=='curb':
            span=loop.spans[int(order[i])]
            tag=original_banks[span.owner][span.segment]
            if tag.kind!='curb' or tag.bank is None:
                raise ValueError('Terminal union curb is missing its explicit original bank ancestry')
            curb_ancestry[index]=tag.bank.subspan(span.first,span.last)
            validate_ancestry(approach,curb_ancestry[index],curve.points)
        if kinds[i]=='internal':
            # Both original source identities are retained separately below;
            # this is an original root cap, never an invented node bridge.
            caps[len(labels)-1]=tuple(a.endpoint for a in aliases)
    layout=SimpleNamespace(approaches=(approach,),outward_stations=(station,),boundaries=tuple(curves),
        rows=len(rim),radial=radial,mouths={},corners=(),original_terminal_root_caps=caps)
    layout.terminal_merged_port=dict(mode='FINITE_ORIGINAL_TERMINAL_STRIP_UNION',bundle=approach.bundle_id,
        bundle_end=approach.end_index,station=float(station),source_curve_tolerance=approach.source_curve_tolerance,
        aliases=[dict(section=a.endpoint.section_id,end=a.endpoint.end,station=a.station) for a in aliases],
        source_sections=list(approach.source_run_sections),source_prefix_proof=approach.source_prefix_proof)
    layout.terminal_source_polygons=tuple(source_polygons)
    layout.terminal_curb_ancestry=curb_ancestry
    return SampledFootprint(Region(owner,rim),tuple(labels),layout)


def _on_segment(point,ends):
    delta=ends[1]-ends[0]
    t=float((point-ends[0])@delta/(delta@delta))
    return -1e-9<=t<=1.+1e-9 and np.linalg.norm(point-(ends[0]+t*delta))<=1e-7


def certify(group):
    """Continuously retain both clipped finite strips, excluding every hole."""
    from .finite_containment import contains_polygon,segment_on_polygon_side
    from .source_strip_containment import _distance
    from .regions import _classify
    outer=[loop.xy for loop in group.loops if not loop.hole]
    holes=[loop.xy for loop in group.loops if loop.hole]
    if len(outer)!=1:raise ValueError('Merged terminal strips require one physical exterior')
    reports=[]
    for source in group.footprints:
        record=getattr(source.layout,'terminal_merged_port',None)
        if record is None:continue
        band=.002-record['source_curve_tolerance']
        for polygon in source.layout.terminal_source_polygons:
            if not contains_polygon(polygon,outer[0],tolerance=band):
                raise ValueError('Compound exterior clips a merged terminal original source strip')
            for hole in holes:
                if not all(segment_on_polygon_side(a,b,hole,inside=False,tolerance=band)
                        for a,b in zip(polygon,np.roll(polygon,-1,axis=0))):
                    raise ValueError('Compound hole crosses a merged terminal original source strip')
            sides=np.stack((polygon,np.roll(polygon,-1,axis=0)),axis=1)
            for boundary in (outer[0],*holes):
                occupied,near=_classify(boundary,polygon,1e-9)
                if np.any((occupied | (near>=0)) & (_distance(boundary,sides)>.002)):
                    raise ValueError('Compound boundary enters a merged terminal original source strip')
        reports.append(dict(bundle=record['bundle'],bundle_end=record['bundle_end'],
            aliases=record['aliases'],source_chord_error=record['source_curve_tolerance'],
            boundary_band=band,maximum_combined_error=.002,
            proof='continuous_clipped_original_strip_boundaries_with_explicit_hole_exclusion'))
    return reports
