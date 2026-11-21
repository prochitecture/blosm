"""Reserve finite curb contacts through an explicitly consumed original road.

Logical corner policy and actual contact profiles are separate witnesses. A
changed-width connector is an internal centreline closure, never a uniform
contact Route or an extrapolated successor bank.
"""
from dataclasses import dataclass,replace
import numpy as np
from .model import End,Key
from .routes import Route,RouteTopology
from .intersections import Approach
from .corner_radii import CornerRadiusPolicy
from .generated_junction_curbs import GeneratedJunctionArea,is_generated_junction_area
from .compound_displaced_contacts import SourceContactWalk,patch
from .source_bank_reconciliation import (SourceBankDiscovery,NativeBankWitness,
    CircularCornerWitness,BankOrigin,BankSwitch,owned_ranges)
from .compounds import rim_label

SOURCE_OWNED_CONNECTOR_POLICY='2-visible-successor-bank-contacts'


@dataclass(frozen=True)
class NativeCurbContact:
    path: object
    side: int
    station: float
    owned_intervals: tuple=()


@dataclass(frozen=True)
class OwnedConnectorWalk:
    original_origin: BankOrigin
    through_site: object
    contact_origin: BankOrigin
    pieces: tuple
    internal_banks: tuple
    sector_policies: tuple


@dataclass(frozen=True)
class PreparedOwnedConnectorFillet(GeneratedJunctionArea):
    proposal: BankSwitch
    original_corner: CircularCornerWitness
    connector_walk: OwnedConnectorWalk
    corner: object
    source_contact_walks: tuple
    contacts: tuple

    @property
    def first_origin(self):
        return self.original_corner.first.origins[0]

    @property
    def token(self):
        p=self.proposal;w=self.connector_walk
        return (p.site,p.first.token,p.last.token,w.original_origin.endpoint,
            w.through_site,w.contact_origin.endpoint)


def _source_approach(planner,origin,path):
    layout=planner.layouts.get(origin.site)
    matches=[arm for arm in getattr(layout,'approaches',())
        if type(arm) is Approach and End(arm.source.id,arm.end)==origin.endpoint]
    if len(matches)!=1:return None
    arm=matches[0]
    if (type(arm.path) is not Route or len(arm.path.pieces)!=1
            or arm.path.pieces[0].path is not path
            or arm.path.pieces[0].forward!=(origin.endpoint.end==0)):
        return None
    return arm


def _policy(records):
    strict={radius for site,end,radius,mode in records if mode=='STRICT'}
    style=[radius for site,end,radius,mode in records if mode=='STYLE']
    if len(strict)>1:return None
    if strict:
        radius=next(iter(strict))
        return (radius,'STRICT') if not style or radius<=min(style) else None
    return (min(style),'STYLE') if style else None


def discover(planner,group):
    """Adjacent boundary ancestry, one real consumed connector, no XY pairing."""
    context=planner.context;d=SourceBankDiscovery(context,group)
    scoped={End(section.id,end) for section in context.network.sections
        for end in (*section.scope_ends,*section.selection_ends)}
    internal_ids={end.section_id for end in d.internal};found=[]
    for li,loop in enumerate(group.loops):
        for index,span in enumerate(loop.spans):
            prior=loop.spans[index-1]
            if prior.owner==span.owner or is_generated_junction_area(d.sources[prior.owner]):continue
            first,last=d.describe(prior,True),d.describe(span,False)
            if not isinstance(first,CircularCornerWitness) or not isinstance(last,NativeBankWitness):continue
            if (d.same_bank(first.last,last) or not d.visible_anchor(loop,index-1,first,entry=True)
                    or d.policy(first.first,last,'direct') is not None):continue
            origins=[]
            for origin in first.last.origins:
                start=origin.endpoint;far=End(start.section_id,1-start.end)
                if (start not in d.internal or far not in d.internal or start in scoped or far in scoped):continue
                through=d.by_end.get(far)
                observed=[value for value in last.origins if value.site==through]
                if len(observed)==1:origins.append((origin,far,through))
            if len(origins)!=1:continue
            origin,far,through=origins[0]
            connector=_source_approach(planner,origin,first.last.path)
            if connector is None:continue
            route=planner.topology.trace(origin.endpoint)
            if (route.terminal_site!=through or route.terminal_scope or route.terminal_selection
                    or tuple((p.section.id,p.forward) for p in connector.path.pieces)
                       !=tuple((v.section_id,v.forward) for v in route.traversals)):continue
            internal=[];target=None;previous=None;owners=0
            for step in range(len(loop.spans)):
                current=loop.spans[(index+step)%len(loop.spans)]
                label=rim_label(current,d.sources);identity=(current.owner,label.kind,label.index)
                if identity==previous:continue
                previous=identity;owners+=1
                if owners>32 or label.kind!='curb':break
                bank=d.describe(current,False)
                if not isinstance(bank,NativeBankWitness):break
                if bank.path.section.id in internal_ids:
                    if not any(value.site==through for value in bank.origins):break
                    internal.append(bank.token);continue
                choices=[value for value in bank.origins if value.site==through
                    and value.endpoint not in scoped and value.outward_side==1]
                if len(choices)!=1:break
                target=bank;contact_origin=choices[0];break
            if target is None or not internal:continue
            first_origins=[value for value in first.first.origins
                if value.site==first.site and value.outward_side==-1]
            if len(first_origins)!=1:continue
            a=_source_approach(planner,first_origins[0],first.first.path)
            b=_source_approach(planner,contact_origin,target.path)
            if a is None or b is None:continue
            original_policy=d.policy(first.first,first.last,'logical_original_sector')
            far_bank=replace(first.last,station=0. if far.end==0 else first.last.path.length,
                origins=(BankOrigin(through,far,-origin.outward_side),))
            through_policy=d.policy(far_bank,target,'owned_connector_exit_sector')
            if original_policy is None or through_policy is None:continue
            sector_policies=tuple((site,end,radius,mode) for site,_,_,rows in (original_policy,through_policy)
                for end,radius,mode in rows)
            fitted_policy=_policy(sector_policies)
            if fitted_policy is None:continue
            radius,mode=fitted_policy
            proposal=BankSwitch(li,index,'CLIPPED_CORNER_OWNED_CONNECTOR_SUCCESSOR',first.first,
                target,first.site,radius,mode,original_policy[3],tuple(loop.xy[index]))
            witness=OwnedConnectorWalk(origin,through,contact_origin,tuple(connector.path.pieces),
                tuple(internal),sector_policies)
            found.append((proposal,first,witness,a,b,connector))
            if len(found)>16:raise ValueError('Owned connector discovery exceeds its finite original-sector budget')
    return tuple(found)


def required(planner,group,previous):
    """Return one reservation; existing owner extents are not changed here."""
    used={record.token for record in previous}
    for proposal,original,witness,a,b,connector in discover(planner,group):
        arms=[]
        for bank,end in zip((proposal.first,proposal.last),(a.end,b.end)):
            arms.append(Approach(bank.path,end,planner.context[Key('profile',bank.path.section.id)]))
        def accept(corner):
            return corner.center is not None and all(0.<=at<=arm.length for arm,at in zip(arms,corner.stations))
        corner=CornerRadiusPolicy.solve(*arms,proposal.radius,proposal.mode,displaced=True,accept=accept)
        walks=(SourceContactWalk(a,arms[0]),SourceContactWalk(connector,arms[1],witness.pieces))
        corner.source_contact_walks=walks
        corner.source_contact_support=dict(policy='EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS',
            original_junction=proposal.site.text(),through_site=witness.through_site.text(),
            original_connector_sections=tuple(piece.section.id for piece in witness.pieces),
            actual_contact_sources=(proposal.first.token,proposal.last.token))
        contacts=tuple(NativeCurbContact(bank.path,bank.side,float(at if arm.end==0 else arm.length-at))
            for bank,arm,at in zip((proposal.first,proposal.last),arms,corner.stations))
        record=PreparedOwnedConnectorFillet(proposal,original,witness,corner,walks,contacts)
        if record.token in used:
            raise ValueError('Owned connector reservation failed to remove its original clipped-corner switch')
        return record
    return None


def insert(planner,members,consumed,sources,stations,records,spacing,tolerance):
    """Reserve actual native contacts before the normal source port search."""
    if len(records)>16:raise ValueError('Owned connector reservations exceed their finite retry budget')
    for index,record in enumerate(records):
        w=record.connector_walk
        if (record.proposal.site not in members or w.through_site not in members
                or any(piece.section.id not in consumed for piece in w.pieces)):
            raise ValueError('Owned connector reservation lost its original fully consumed incidences')
        source=patch(Key('source_owned_connector_fillet',index),record.corner,spacing,tolerance)
        source.layout.original_junction=record.proposal.site
        source.layout.generated_junction_area=record;source.layout.source_owned_connector_fillet=record
        sources.append(source)
        origins=(record.first_origin,w.contact_origin)
        for contact,arm,origin in zip(record.contacts,(record.corner.first,record.corner.second),origins):
            actual=_source_approach(planner,origin,contact.path)
            if actual is None or actual.end!=arm.end:
                raise ValueError('Owned connector contact lost its exact original native Approach')
            endpoint=origin.endpoint;at=contact.station if actual.end==0 else actual.length-contact.station
            stations[endpoint]=max(stations.get(endpoint,0.),at)


def contact_cuts(cuts,records):
    for record in records:
        for contact in record.contacts:cuts.setdefault(contact.path.section.id,set()).add(contact.station)


def _endpoint(end):
    return end.section_id, end.end


def _sector_rows(record):
    return tuple((site.text(), _endpoint(end), radius, mode)
                 for site,end,radius,mode in record.connector_walk.sector_policies)


def bind_contacts(context,group):
    """Seal reservations only against the final, normally prepared mouths."""
    from copy import copy
    ranges=owned_ranges(context,group);sources=[];records=[]
    for source in group.footprints:
        record=getattr(source.layout,'source_owned_connector_fillet',None)
        if not isinstance(record,PreparedOwnedConnectorFillet):sources.append(source);continue
        contacts=[]
        for contact in record.contacts:
            sid=contact.path.section.id;owned=ranges.get(sid,())
            if (contact.path is not context[Key('path',sid)] or contact.side not in (-1,1)
                    or not any(a<=contact.station<=b for a,b in owned)):
                raise ValueError('Owned connector contact lacks final finite native ownership')
            contacts.append(replace(contact,owned_intervals=owned))
        record=replace(record,contacts=tuple(contacts));layout=copy(source.layout)
        layout.generated_junction_area=record;layout.source_owned_connector_fillet=record
        sources.append(replace(source,layout=layout));records.append(record)
    if not records:return group
    layout=copy(group.physical_layout);layout.source_owned_connector_fillets=tuple(records)
    layout.source_owned_connector_reconciliation=dict(policy=SOURCE_OWNED_CONNECTOR_POLICY,
        records=[dict(original_site=r.proposal.site.text(),original_corner_ends=(
            _endpoint(r.first_origin.endpoint),_endpoint(r.connector_walk.original_origin.endpoint)),
            through_site=r.connector_walk.through_site.text(),actual_contacts=[dict(
                section=v.path.section.id,side=v.side,station=v.station,owned_intervals=v.owned_intervals)
                for v in r.contacts],radius=r.corner.radius,radius_mode=r.proposal.mode,
            sector_policies=_sector_rows(r),
            original_connector_sections=tuple(piece.section.id for piece in r.connector_walk.pieces),
            original_corners_unchanged=True,profile_leads_recomputed=True) for r in records])
    return replace(group,footprints=tuple(sources),physical_layout=layout)


def signature(group):
    records=getattr(group.physical_layout,'source_owned_connector_fillets',())
    if not records:return None
    return (SOURCE_OWNED_CONNECTOR_POLICY,tuple((record.proposal.kind,record.proposal.site.text(),
        record.proposal.first.token,record.proposal.last.token,_endpoint(record.connector_walk.original_origin.endpoint),
        record.connector_walk.through_site.text(),_endpoint(record.connector_walk.contact_origin.endpoint),
        record.proposal.radius,record.proposal.mode,_sector_rows(record),
        tuple((piece.section.id,piece.forward) for piece in record.connector_walk.pieces),
        tuple((contact.path.section.id,contact.side,contact.station,contact.owned_intervals)
              for contact in record.contacts),record.corner.radius) for record in records))
