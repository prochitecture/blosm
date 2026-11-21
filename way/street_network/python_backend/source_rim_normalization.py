"""Bounded resampling of a derived source-union curb, with provenance."""
from copy import copy
from dataclasses import replace
import numpy as np
from .finite_free_caps import PreparedCapCurb
from .compounds import rim_label
from .regions import RegionUnion
from .model import Key


class ResampledSourceCurb(PreparedCapCurb):
    """Project onto a source-union polyline within one micrometre.

    This records additional *polyline* error. Original analytical primitives
    remain authoritative in the layout and source references; no analytical
    circle/offset equivalence is claimed for the projected chord chain.
    """
    def __init__(self,points,source_spans,source_curves,tolerance=1e-6):
        original=np.asarray(points,dtype=float).copy()
        if (original.ndim!=2 or original.shape[1]!=2 or len(original)<2
                or not np.isfinite(original).all() or not np.isfinite(tolerance)
                or not 0.<tolerance<=1e-6):
            raise ValueError('Source curb normalization needs finite points and at most one micrometre projection')
        keep={0,len(original)-1};pending=[(0,len(original)-1)]
        while pending:
            first,last=pending.pop()
            if last-first<2:continue
            edge=original[last]-original[first];squared=float(edge@edge)
            if squared<=1e-20:
                index=(first+last)//2;keep.add(index);pending.extend(((first,index),(index,last)));continue
            fraction=np.clip((original[first+1:last]-original[first])@edge/squared,0.,1.)
            error=np.linalg.norm(original[first+1:last]-original[first]-fraction[:,None]*edge,axis=1)
            at=int(np.argmax(error))
            if error[at]>tolerance:
                index=first+1+at;keep.add(index);pending.extend(((first,index),(index,last)))
        indices=sorted(keep);maximum=0.
        for first,last in zip(indices,indices[1:]):
            edge=original[last]-original[first]
            if float(edge@edge)<=1e-20:raise ValueError('Source curb normalization collapses an original finite interval')
            fraction=np.clip((original[first:last+1]-original[first])@edge/(edge@edge),0.,1.)
            maximum=max(maximum,float(np.linalg.norm(original[first:last+1]-original[first]-fraction[:,None]*edge,axis=1).max()))
        if maximum>tolerance+1e-12:raise ValueError('Source curb normalization exceeds its polyline projection allowance')
        super().__init__(original[indices])
        self.original_polyline=original;self.original_polyline.flags.writeable=False
        self.source_spans,self.source_curves=tuple(source_spans),tuple(source_curves)
        self.maximum_projection_error,self.projection_allowance=maximum,tolerance
        parameters=self.stations/self.length
        if (len(parameters)-1)%2:
            at=int(np.argmax(np.diff(self.stations)))
            parameters=np.sort(np.r_[parameters,.5*(parameters[at]+parameters[at+1])])
        self.initial_parameters=parameters

    def report(self):
        return dict(reference='sampled_original_ground_union_polyline',original_points=len(self.original_polyline),
            retained_points=len(self.points),maximum_projection_error=self.maximum_projection_error,
            additional_error_allowance=self.projection_allowance,profile_and_scope_endpoints_exact=True,
            original_analytical_primitives_retained=True,
            source_curves=[dict(owner=owner.text(),index=index,curve_type=type(curve).__name__)
                for owner,index,curve in self.source_curves])


def normalize(group,*,maximum_projection=1e-6):
    """Consolidate one proven source rim; preserve complete profile samples.

    Logical claims, all original proof operands, dynamic layout metadata and
    original analytical primitives survive. The caller must re-run complete
    source-width containment against the resulting physical rim.
    """
    if not np.isfinite(maximum_projection) or not 0.<maximum_projection<=1e-6:
        raise ValueError('Source rim normalization allows at most one micrometre projection')
    if group.physical_layout is None or len(group.loops)!=1 or group.loops[0].hole:
        raise ValueError('Source rim normalization requires one prepared physical exterior')
    loop=group.loops[0];layout=group.physical_layout
    previous=getattr(layout,'source_rim_normalization',None)
    if previous:
        if any(row['additional_error_allowance']>maximum_projection for row in previous):
            raise ValueError('An already normalized source rim needs its original projection budget')
        return group
    inventory={source.region.owner:source for source in group.footprints}
    labels=[rim_label(span,inventory) for span in loop.spans]
    if any(label.kind not in ('curb','mouth') for label in labels):
        raise ValueError('Source rim normalization cannot hide an internal source seam')
    by_identity={id(a):i for i,a in enumerate(layout.approaches)}
    tokens=[]
    for span,label in zip(loop.spans,labels):
        a=inventory[span.owner].layout.approaches[label.index] if label.kind=='mouth' else None
        tokens.append(by_identity.get(id(a)) if a is not None else None)
        if label.kind=='mouth' and tokens[-1] is None:
            raise ValueError('Source rim normalization has an unmapped original physical profile')
    starts={}
    for i,token in enumerate(tokens):
        if token is not None and token!=tokens[i-1]:
            if token in starts:raise ValueError('Source rim normalization splits a complete profile')
            starts[token]=i
    if set(starts)!=set(range(len(layout.approaches))):
        raise ValueError('Source rim normalization loses a complete incoming profile')
    start=starts[0];xy=np.roll(loop.xy,-start,axis=0)
    spans=list(loop.spans[start:])+list(loop.spans[:start]);labels=labels[start:]+labels[:start]
    tokens=tokens[start:]+tokens[:start]
    starts=sorted((i,t) for i,t in enumerate(tokens) if t is not None and t!=tokens[i-1])
    order=[t for _,t in starts];curves=[];profile_rows=[]
    for side,(begin,token) in enumerate(starts):
        finish=starts[side+1][0] if side+1<len(starts) else len(tokens)
        index=begin;value=1.;samples=[]
        while index<finish and tokens[index]==token:
            label=labels[index]
            if abs(label.first-value)>1e-8 or label.last>=label.first:
                raise ValueError('Source rim normalization breaks an original profile partition')
            samples.append(label.first);value=label.last;index+=1
        if abs(value)>1e-8 or index>=finish:
            raise ValueError('Source rim normalization clips a profile or loses its outer curb')
        samples.append(0.);profile_rows.append(np.asarray(samples))
        points=np.vstack((xy[index:finish],xy[finish%len(xy)]));provenance={}
        for span,label in zip(spans[index:finish],labels[index:finish]):
            provenance[span.owner,label.index]=inventory[span.owner].layout.boundaries[label.index]
        curves.append(ResampledSourceCurb(points,spans[index:finish],
            [(owner,index,curve) for (owner,index),curve in sorted(provenance.items())],maximum_projection))
    physical=copy(layout)
    physical.approaches=tuple(layout.approaches[i] for i in order)
    physical.outward_stations=tuple(layout.outward_stations[i] for i in order)
    physical.guide_weights=[0.]*len(order);physical.boundaries=tuple(curves)
    physical.rows=max(len(curve.points)-1 for curve in curves)
    physical.source_rim_normalization=tuple(curve.report() for curve in curves)
    # sample_footprint queries the exact same profile partition contract.
    # Preserve an explicit row if a source provider supplied a nonstandard one.
    from .compounds import SampledFootprint,RimSegment
    from .regions import Region
    rim=[];segments=[]
    for side,(approach,station,parameters,curve) in enumerate(zip(physical.approaches,
            physical.outward_stations,profile_rows,curves)):
        points=approach.point(station,approach.offset+approach.width*(parameters-.5))
        old=xy[starts[side][0]:starts[side][0]+len(parameters)-1]
        if not np.allclose(points[:-1],old,atol=1e-7,rtol=0.):
            raise ValueError('Source rim normalization changed original incoming profile samples')
        rim.extend(old)
        segments.extend(RimSegment('mouth',side,float(a),float(b)) for a,b in zip(parameters,parameters[1:]))
        t=curve.initial_parameters;rim.extend(curve(t)[:-1])
        segments.extend(RimSegment('curb',side,float(a),float(b)) for a,b in zip(t,t[1:]))
    owner=Key('normalized_source_rim',group.members[0].id)
    if owner in inventory:raise ValueError('Source rim normalization repeats a derived physical owner')
    normalized=SampledFootprint(Region(owner,np.asarray(rim)),tuple(segments),physical)
    loops=RegionUnion().build((normalized.region,))
    if len(loops)!=1 or loops[0].hole:raise ValueError('Source curb normalization changed physical topology')
    # Original operands supply immutable proof polygons and analytical source
    # provenance, even when only the consolidated rim owns boundary spans.
    return replace(group,physical_layout=physical,footprints=(normalized,*group.footprints),loops=loops)


def _profiles(group):
    from .conditions import BoundaryConditions
    layouts={source.region.owner:source.layout for source in group.footprints}
    profiles={}
    for loop in group.loops:
        boundary=BoundaryConditions().build(group,loop,layouts)
        for port in boundary.ports:
            key=(port.section_id,port.end,port.station)
            if key in profiles:raise ValueError('Source exterior has a duplicate physical profile')
            profiles[key]=(boundary.xy[list(port.vertices)],boundary.heights[list(port.vertices)])
    return profiles


def normalize_exterior(group,*,maximum_projection=1e-6):
    """Normalize only the derived exterior; retain every physical hole exactly.

    The new footprint supplies span/curve lookup only. Its Region is explicitly
    boundary-only and RegionUnion rejects it as an area operand. Original area
    operands remain available through group.area_footprints; the authoritative
    physical domain is group.loops (outer minus holes), never the carrier disk.
    """
    from .finite_containment import contains_segment
    if (not group.loops or group.loops[0].hole
            or any(not loop.hole for loop in group.loops[1:])
            or len(group.loops)>9 or sum(len(loop.xy) for loop in group.loops)>20000):
        raise ValueError('Source exterior normalization needs one bounded outer rim and its physical holes')
    old=_profiles(group)
    candidate=normalize(replace(group,loops=(group.loops[0],)),maximum_projection=maximum_projection)
    if candidate is not group and candidate.physical_layout is group.physical_layout:
        # normalize() returned its already-normalized exterior view.
        if not getattr(group.footprints[0].region,'boundary_only',False):
            raise ValueError('An already normalized exterior needs explicit boundary-only area provenance')
        return group
    holes=group.loops[1:]
    exterior=candidate.loops[0].xy
    for hole in holes:
        if not all(contains_segment(a,b,exterior,tolerance=1e-8)
                for a,b in zip(hole.xy,np.roll(hole.xy,-1,axis=0))):
            raise ValueError('Source exterior normalization crosses an unchanged physical hole')
    derived=candidate.footprints[0]
    carrier=replace(derived,region=replace(derived.region,boundary_only=True))
    candidate=replace(candidate,footprints=(carrier,*candidate.footprints[1:]),
                      loops=(candidate.loops[0],*holes))
    current=_profiles(candidate)
    if set(old)!=set(current) or any(not np.array_equal(a,b)
            for key in old for a,b in zip(old[key],current[key])):
        raise ValueError('Source exterior normalization changed an exact physical profile XY or height')
    if len(candidate.area_footprints)!=len(group.area_footprints) or any(
            a is not b for a,b in zip(candidate.area_footprints,group.area_footprints)):
        raise ValueError('Source exterior normalization changed original physical area operands')
    candidate.physical_layout.source_exterior_sampling=dict(
        boundary_only_owner=carrier.region.owner.text(),original_holes_exact=True,
        original_profiles_exact=True,original_area_operands_unchanged=True,
        physical_area_reference='group.loops: exterior minus original holes')
    return candidate
