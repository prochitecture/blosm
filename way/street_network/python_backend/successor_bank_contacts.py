"""Fit an acute exterior corner beyond the end of a fully owned Section."""
from dataclasses import dataclass, replace
import math
import numpy as np

from .compound_displaced_contacts import SourceContactWalk
from .compounds import rim_label
from .corner_radii import CornerRadiusPolicy
from .intersections import Approach, Corner
from .local_bank_contacts import DirectedBank
from .model import End, Key
from .source_bank_reconciliation import (
    CircularCornerWitness, NativeBankWitness, SourceBankDiscovery, visible_bank_stations)
from .source_owned_connector_fillets import (
    NativeCurbContact, OwnedConnectorWalk, PreparedOwnedConnectorFillet, _source_approach)
from .source_parent_chain_fillets import ParentChainUnavailable, offset_solution


@dataclass(frozen=True)
class SuccessorBankFillet(PreparedOwnedConnectorFillet):
    start: object

    @property
    def first_origin(self):
        return self.start


def validate_exposure(group):
    records = getattr(group.physical_layout,'source_owned_connector_fillets',())
    if not any(isinstance(record,SuccessorBankFillet) for record in records):
        return
    sources = {area.region.owner:area for area in group.footprints}
    for owner, source in sources.items():
        if not isinstance(getattr(source.layout,'source_owned_connector_fillet',None),SuccessorBankFillet):
            continue
        labels = sorted((rim_label(span,sources) for loop in group.loops for span in loop.spans
                         if span.owner == owner),key=lambda label:label.first)
        if not labels or len({(label.kind,label.index) for label in labels}) != 1 or labels[0].kind != 'curb':
            raise ValueError('Successor bank must expose its circular curb and enclose its centreline closure')
        previous = 0.
        for label in labels:
            if abs(label.first-previous) > 1e-9 or label.last <= label.first:
                raise ValueError('Successor-bank arc has an interrupted contact interval')
            previous = label.last
        if abs(previous-1.) > 1e-9:
            raise ValueError('Successor-bank arc ends before its tangent contact')


def _next_source(discovery, spans, index, matches):
    for _ in spans:
        witness = discovery.describe(spans[index], False)
        if not matches(witness):
            return index, witness
        index = (index + 1) % len(spans)
    return index, None


def _regular(bank, domain, radius):
    low, high = sorted(domain)
    path = bank.path
    if path.exact:
        return not np.any((path.knots > low) & (path.knots < high))
    offset = path.section.offset + .5 * path.section.width * bank.side
    return (abs(offset) + radius) * path.visible_curvature_bound(low, high) < 1. - 1e-10


def reserve(planner, group, switch):
    loop = group.loops[switch.loop]
    if loop.hole or switch.kind != 'DIRECT_ORIGINAL_BANK_SWITCH' or switch.mode != 'STYLE':
        return None
    context = planner.context
    discovery = SourceBankDiscovery(context, group)
    a, connector = switch.first, switch.last
    origins = [[o for o in bank.origins if o.site == switch.site and o.outward_side == side]
               for bank, side in ((a,-1), (connector,1))]
    if any(len(items) != 1 for items in origins):
        return None
    start, departure = (items[0] for items in origins)
    far = End(connector.path.section.id, 1-departure.endpoint.end)
    if not {departure.endpoint, far} <= discovery.internal:
        return None
    spans = loop.spans
    i, bend = _next_source(discovery, spans, switch.vertex, lambda w:
        isinstance(w, NativeBankWitness) and discovery.same_bank(w,connector))
    if (not isinstance(bend,CircularCornerWitness) or bend.site != discovery.by_end[far]
            or bend.site == switch.site or bend.corner.radius_mode != 'STYLE'
            or bend.corner.sweep <= 0. or not discovery.same_bank(bend.first,connector)
            or not discovery.visible_anchor(loop,i,bend,entry=True)
            or not discovery.visible_anchor(loop,i,bend,entry=False)):
        return None
    j, b = _next_source(discovery, spans, i, lambda w:
        isinstance(w,CircularCornerWitness) and w.corner is bend.corner)
    if not isinstance(b,NativeBankWitness) or not discovery.same_bank(b,bend.last):
        return None
    arrival, = bend.last.origins
    if (arrival.outward_side != 1 or arrival.endpoint in discovery.internal
            or any(s.path.section.scope_ends or s.path.section.selection_ends for s in (a,connector,b))):
        return None
    sector = discovery.policy(bend.first,bend.last,'direct')
    if sector is None or sector[0] != bend.site or sector[2] != 'STYLE':
        return None
    policies = tuple((site,end,r,mode) for site,rows in ((switch.site,switch.sectors),(bend.site,sector[3]))
                     for end,r,mode in rows)
    first = _source_approach(planner,start,a.path)
    through = _source_approach(planner,departure,connector.path)
    if first is None or through is None or _source_approach(planner,arrival,b.path) is None:
        return None
    witnessed = replace(switch, vertex=j, last=b)
    da = visible_bank_stations(discovery,switch,incoming=True)
    db = visible_bank_stations(discovery,witnessed,incoming=False)
    if da is None or db is None:
        return None
    preferred = min(r for _,_,r,_ in policies)
    arms = tuple(Approach(bank.path,origin.endpoint.end,context[Key('profile',bank.path.section.id)])
                 for bank,origin in ((a,start),(b,arrival)))
    locus = DirectedBank(a,da)
    outward = db if arms[1].end == 0 else tuple(arms[1].length-s for s in db)
    attempts = []
    for radius in CornerRadiusPolicy.candidates(preferred,'STYLE'):
        attempts.append(radius)
        if not (_regular(a,da,radius) and _regular(b,db,radius)):
            continue
        try:
            hit = offset_solution(locus,(0.,1.),arms[1],outward,radius)
        except ParentChainUnavailable:
            continue
        at = float(locus.station(hit['poly_contact']))
        second = hit['bank_outward_contact']
        stations = (at if arms[0].end == 0 else arms[0].length-at, second)
        corner = Corner(*arms,stations,np.asarray(hit['center']),radius,hit['start_angle'],hit['sweep'],
                        math.pi+hit['sweep'],requested_radius=preferred,radius_mode='STYLE',
                        radius_attempts=tuple(attempts),fit_policy='VISIBLE_SUCCESSOR_BANK')
        pieces = through.path.pieces
        walks = (SourceContactWalk(first,arms[0]),SourceContactWalk(through,arms[1],pieces))
        corner.source_contact_walks = walks
        corner.source_contact_support = dict(policy='EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS',
            original_junction=switch.site.text(), through_site=bend.site.text(),
            original_connector_sections=tuple(piece.section.id for piece in pieces),
            actual_contact_sources=(a.token,b.token))
        contacts = (NativeCurbContact(a.path,a.side,at),
                    NativeCurbContact(b.path,b.side,second if arms[1].end == 0 else arms[1].length-second))
        walk = OwnedConnectorWalk(departure,bend.site,arrival,tuple(pieces),(connector.token,),policies)
        return SuccessorBankFillet(replace(switch,last=b,radius=preferred),bend,walk,corner,walks,contacts,start)
    return None
