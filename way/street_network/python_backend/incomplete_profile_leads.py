"""Typed evidence for an existing finite complete-profile lead search.

This does not move a source curve or relax a full-profile test. It records the
ordinary source-union builder's existing outward step when the previous full
profile was actually clipped by the other original source operands.
"""
from dataclasses import dataclass
import hashlib
import math
import numpy as np
from .compounds import rim_label
from .regions import RegionUnion


def _exposed(port,sources,loops):
    inventory={source.region.owner:source for source in sources}
    if inventory.get(port.region.owner) is not port:
        raise ValueError('Incomplete-profile evidence lost its actual previous source port')
    return sum(label.first-label.last for loop in loops for span in loop.spans
               if span.owner==port.region.owner and (label:=rim_label(span,inventory)).kind=='mouth')


def _source_fingerprints(sources):
    return tuple((source.region.owner,np.asarray(source.region.xy).dtype.str,source.region.xy.shape,
                  hashlib.sha256(np.ascontiguousarray(source.region.xy).tobytes()).digest())
                 for source in sources)


@dataclass(frozen=True)
class IncompleteProfileLead:
    approach: object
    previous: float
    candidate: float
    previous_port: object
    previous_sources: tuple
    source_fingerprints: tuple
    exposed_measure: float
    iteration: int

    def validate(self,context):
        from .intersections import Approach,IntersectionProvider
        from .routes import Route
        a=self.approach
        if (type(a) is not Approach or type(a.path) is not Route
                or len(self.previous_port.layout.approaches)!=1
                or self.previous_port.layout.approaches[0] is not a
                or self.previous_port.layout.outward_stations!=(self.previous,)
                or type(self.iteration) is not int or not 1<=self.iteration<=31
                or not all(math.isfinite(v) for v in (self.previous,self.candidate,self.exposed_measure))
                or not 0.<=self.previous<self.candidate<a.length
                or not 0.<=self.exposed_measure<1.-1e-8
                or not IntersectionProvider.mouth_fits(a,self.candidate)):
            raise ValueError('Incomplete-profile evidence needs an actual finite ordinary source lead')
        if self.candidate!=IntersectionProvider.mouth_station(a,
                self.previous+max(.25,context.network.settings.spacing),context.network.settings):
            raise ValueError('Incomplete-profile evidence changed the existing bounded lead step')
        from .model import Key
        if any(context[Key('path',piece.section.id)] is not piece.path for piece in a.path.pieces):
            raise ValueError('Incomplete-profile evidence changed its original Path identity')
        if (getattr(a.path,'mouth_piece',None) is not None
                or a.profile!=context[Key('profile',a.source.id)].oriented(a.end)
                or a.width!=a.source.width):
            raise ValueError('Incomplete-profile evidence changed its original complete Profile')
        old,new=a.path.trims(self.previous),a.path.trims(self.candidate)
        if set(old)!=set(new) or any(new[key]<value if key[1]==0 else new[key]>value
                                     for key,value in old.items()):
            raise ValueError('Incomplete-profile evidence acquired another original source or moved inward')
        if sum(len(source.region.xy) for source in self.previous_sources)>20000:
            raise ValueError('Incomplete-profile evidence exceeds the existing source-union sample bound')
        if (len({source.region.owner for source in self.previous_sources})!=len(self.previous_sources)
                or _source_fingerprints(self.previous_sources)!=self.source_fingerprints):
            raise ValueError('Incomplete-profile evidence changed its original source operands')
        for curve in self.previous_port.layout.boundaries:
            coordinates=getattr(curve,'source_coordinates',None)
            if coordinates is not None:coordinates.validate(curve.points,context)
        loops=RegionUnion(max_segments=20000).build(source.region for source in self.previous_sources
            if not getattr(source.region,'boundary_only',False))
        actual=_exposed(self.previous_port,self.previous_sources,loops)
        if actual!=self.exposed_measure:
            raise ValueError('Incomplete-profile evidence no longer proves the actual missing profile')
        return True


def capture(approach,previous,candidate,port,sources,loops,iteration):
    from .intersections import Approach
    from .routes import Route
    if type(approach) is not Approach or type(approach.path) is not Route:return None
    measure=_exposed(port,sources,loops)
    if not 0.<=measure<1.-1e-8:return None
    return IncompleteProfileLead(approach,float(previous),float(candidate),port,tuple(sources),_source_fingerprints(sources),
                                 float(measure),iteration)


def admission(context,group,a,s,b,t):
    """Validate an actual typed chain; a dictionary or moved final cut is insufficient."""
    if a is not b or not 0.<=s<t<a.length:return None
    # A separate root-support construction retains its own donor-plane proof.
    # This policy cannot supersede that construction's explicit fixed bounds.
    if any(getattr(group.physical_layout,name,()) for name in
           ('source_cycle_root_envelopes','source_cycle_root_support_envelopes')):
        return None
    records=getattr(group.physical_layout,'source_cycle_profile_lead_advances',())
    selected=[record for record in records if type(record) is IncompleteProfileLead and record.approach is a]
    if not selected:return None
    current=s;proof=[]
    for record in selected:
        if record.previous!=current:return None
        record.validate(context)
        proof.append(dict(previous_station=record.previous,candidate_station=record.candidate,
                          actual_exposed_profile_measure=record.exposed_measure,
                          normal_union_iteration=record.iteration))
        current=record.candidate
    if current!=t:return None
    return dict(policy='TYPED_EXISTING_INCOMPLETE_PROFILE_LEAD_SEARCH',section=a.source.id,end=a.end,
        original_route_station=s,candidate_route_station=t,advance=t-s,steps=proof,
        old_original_cuts=[dict(section=sid,end=end,station=value) for (sid,end),value in a.path.trims(s).items()],
        new_original_cuts=[dict(section=sid,end=end,station=value) for (sid,end),value in b.path.trims(t).items()],
        original_path_profile_and_incidence_unchanged=True,actual_old_profile_incomplete=True)
