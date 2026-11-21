"""An independently necessary width-transition plane remains authoritative.

A root repair may retain a previously certified recipient plane. It cannot
create a new plane advance or borrow the admission from another source graph.
"""
from dataclasses import dataclass
import hashlib
import numpy as np

from .model import End,Key


class StagedRecipientUnavailable(ValueError):
    pass


def _dimension(group,approach):
    values=[c.dimensional_corner for c in group.physical_layout.corners
        if hasattr(c,'dimensional_corner') and any(a is approach for a in (c.first,c.second))]
    if len(values)!=1:raise StagedRecipientUnavailable('Staged recipient needs one actual certified dimensional corner')
    return values[0]


def _dimension_token(record):
    digest=hashlib.sha256(record.source.region.xy.tobytes())
    for piece in record.curve.pieces:
        continuation=getattr(getattr(piece,'curve',None),'continuation',None)
        if continuation is not None:
            digest.update(np.asarray(continuation.controls).tobytes())
    corner=record.curve.corner
    return (record.stations,record.report.get('policy'),record.report.get('minimum_radius'),
        record.report.get('curvature_bound'),corner.radius,corner.requested_radius,corner.radius_mode,
        tuple((a.source.id,a.end,a.width,a.offset,a.profile.rise,a.profile.datum,a.profile.half)
              for a in (corner.first,corner.second)),digest.hexdigest())


def _group_token(group):
    digest=hashlib.sha256()
    for source in group.area_footprints:
        digest.update(source.region.owner.text().encode());digest.update(source.region.xy.tobytes())
        digest.update(repr(source.segments).encode())
    for loop in group.loops:digest.update(loop.xy.tobytes());digest.update(repr(loop.spans).encode())
    layout=group.physical_layout
    return (group.members,tuple(layout.outward_stations),tuple(sorted(layout.mouths.items())),digest.hexdigest())


@dataclass(frozen=True,eq=False)
class StagedRecipientWitness:
    network_uuid: str
    site: object
    approach: object
    previous: float
    plane: float
    contact: float
    original_layout: object
    pre_group: object
    pre_group_token: tuple
    dimensional_token: tuple
    source_objects: tuple
    proof: dict

    @classmethod
    def capture(cls,context,group,bounds):
        """Call the existing recipient proof before any root support exists."""
        from .ordinary_singleton_sources import OriginalSingletonRootBounds,_ends,_same_mouth_representation,_dimensional_recipient_lead
        from .source_strip_containment import SourceStripContainmentError
        from .cycle_dimensional_retry import certify_layout
        if (type(bounds) is not OriginalSingletonRootBounds or not isinstance(bounds.failure,SourceStripContainmentError)
                or group.members!=(bounds.site,) or group.claims
                or any(hasattr(c,'root_support_envelope') for c in group.physical_layout.corners)):
            raise StagedRecipientUnavailable('Staged recipient needs its actual pre-root original single-junction group')
        owned={piece.section.id for a in group.physical_layout.approaches for piece in a.path.pieces}
        if bounds.failure.section_id not in owned:
            raise StagedRecipientUnavailable('Staged recipient lost its actual failing original Section')
        original=_ends(bounds.layout);current=_ends(group.physical_layout)
        if set(original)!=set(current):raise StagedRecipientUnavailable('Staged recipient changed original End ownership')
        result=[]
        for key,(a,s) in original.items():
            b,t=current[key]
            if _same_mouth_representation(a,s,b,t):continue
            proof=_dimensional_recipient_lead(context,group,a,s,b,t)
            if proof is None:
                raise StagedRecipientUnavailable('A pre-root plane advance has no independent width-transition necessity proof')
            certify_layout(context,group,())
            dependencies=tuple((Key(kind,piece.section.id),context[Key(kind,piece.section.id)])
                for piece in a.path.pieces for kind in ('path','profile'))
            result.append(cls(context.network.uuid,bounds.site,a,float(s),float(t),
                float(proof['prepared_narrow_contact']),bounds.layout,group,_group_token(group),
                _dimension_token(_dimension(group,a)),dependencies,dict(proof)))
        return tuple(result)

    @property
    def endpoint(self):return End(self.approach.source.id,self.approach.end)

    def validate_pre(self,context):
        from .ordinary_singleton_sources import _dimensional_recipient_lead,_ends
        from .source_identity import matches_current_section
        if (context.network.uuid!=self.network_uuid or self.pre_group.members!=(self.site,)
                or any(context[key] is not value for key,value in self.source_objects)
                or any(not matches_current_section(context,piece.section) for piece in self.approach.path.pieces)
                or _group_token(self.pre_group)!=self.pre_group_token
                or _dimension_token(_dimension(self.pre_group,self.approach))!=self.dimensional_token):
            raise StagedRecipientUnavailable('Staged recipient witness is stale or belongs to another source graph')
        original=_ends(self.original_layout).get((self.approach.source.id,self.approach.end))
        if original is None or original[0] is not self.approach or original[1]!=self.previous:
            raise StagedRecipientUnavailable('Staged recipient changed its actual original incoming plane')
        repeated=_dimensional_recipient_lead(context,self.pre_group,self.approach,self.previous,self.approach,self.plane)
        if repeated is None or repeated!=self.proof:
            raise StagedRecipientUnavailable('Staged recipient lost its independently necessary original profile proof')

    def validate_final(self,context,group,a,s,b,t):
        from .ordinary_singleton_sources import _ends
        from .cycle_dimensional_retry import certify_layout
        from .compound_ports import CompoundPortPlanner
        from .source_root_templates import NativeRootEnvelope
        self.validate_pre(context)
        if (group.members!=(self.site,) or group.claims or a is not self.approach or b is not a
                or s!=self.previous or t!=self.plane
                or _dimension_token(_dimension(group,a))!=self.dimensional_token):
            raise StagedRecipientUnavailable('Root construction changed its independently certified recipient')
        before,after=_ends(self.pre_group.physical_layout),_ends(group.physical_layout)
        if (set(before)!=set(after) or any(before[key][0] is not after[key][0] or before[key][1]!=after[key][1] for key in before)):
            raise StagedRecipientUnavailable('Root construction advanced an existing complete original plane')
        roots=[c.root_support_envelope for c in group.physical_layout.corners if hasattr(c,'root_support_envelope')]
        matching=[r for r in roots if type(r) is NativeRootEnvelope
            and any(value is self for value in r.staged_witnesses)]
        if not matching:raise StagedRecipientUnavailable('Final root construction does not own this staged witness')
        for record in matching:
            record.reference.validate(context)
            if record.contact_bounds!=record.reference.contact_bounds:
                raise StagedRecipientUnavailable('Final root construction has forged contact bounds')
            if (len(record.staged_witnesses)!=len(record.reference.staged_witnesses)
                    or any(a is not b for a,b in zip(record.staged_witnesses,record.reference.staged_witnesses))):
                raise StagedRecipientUnavailable('Final root construction borrowed another staged witness')
            for arm,station,limit in zip((record.corner.first,record.corner.second),record.stations,record.contact_bounds):
                if station>limit:raise StagedRecipientUnavailable('Root construction exceeds its original certified contact plane')
                if End(arm.source.id,arm.end)==self.endpoint and station>self.contact:
                    raise StagedRecipientUnavailable('Root construction advances the staged recipient contact')
        if CompoundPortPlanner.incomplete_ports(group):
            raise StagedRecipientUnavailable('Final staged group loses an original incoming profile')
        full=certify_layout(context,group,())
        return dict(self.proof,admission_stage='INDEPENDENT_PRE_ROOT_WIDTH_TRANSITION',
            root_support_did_not_advance_original_recipient=True,
            final_dimensional_source_geometry_exact=True,fresh_final_source_width=full)

    def signature(self,context):
        self.validate_pre(context)
        return (self.site.text(),self.endpoint.section_id,self.endpoint.end,self.previous,self.contact,self.plane,
            self.dimensional_token,tuple((key.text(),context.digest(key)) for key,_ in self.source_objects),
            self.proof['minimum_radius'],self.proof['curvature_bound'],
            self.proof['actual_old_exposed_profile_measure'],self.proof['actual_new_exposed_profile_measure'])


def admission(context,group,a,s,b,t):
    """No witness can authorize a different recipient or any further advance."""
    seen=set()
    for corner in group.physical_layout.corners:
        for witness in getattr(getattr(corner,'root_support_envelope',None),'staged_witnesses',()):
            if type(witness) is not StagedRecipientWitness:
                raise StagedRecipientUnavailable('Final root construction has a forged staged witness')
            if id(witness) in seen:continue
            seen.add(id(witness))
            if witness.approach is a:
                return witness.validate_final(context,group,a,s,b,t)
    return None
