"""Once-owned compound sources in a complete prepared Bundle component.

This is an opt-in planning layer. It neither changes the normal snapshot's
support decisions nor declares a covered source generated. A native provider
must bind the retained internal inputs and validate them with the final patch.
The existing component planner continues to reject any unowned external end.
"""
from copy import copy
from dataclasses import asdict, dataclass, replace
import math


@dataclass(frozen=True)
class CompoundComponentPlan:
    network: object
    junction_interiors: tuple

    def report(self):
        result=self.network.report()
        result['status']='PLANNED_COMPOUND_COMPONENT'
        result['junction_interiors']=[dict(junction_id=p.junction.source_id,
            original_connector_ids=p.original_connector_ids,
            external_connector_ids=tuple(c.source_id for c in p.junction.connectors),
            links=[asdict(link) for link in p.links]) for p in self.junction_interiors]
        result['owned_internal_sections']=sorted(s.source_id for p in self.junction_interiors
            for street in p.internal_streets for s in street.sections)
        return result


def candidate_snapshot(snapshot, packets):
    """Make a private contract view; all original source inventories stay intact."""
    packets=tuple(packets)
    ids=[p.junction.source_id for p in packets]
    if len(ids)!=len(set(ids)):
        raise ValueError('A compound junction has more than one candidate contract')
    view=copy(snapshot)
    view.junctions=dict(snapshot.junctions)
    for packet in packets:
        jid=packet.junction.source_id
        internal=[cid for link in packet.links for cid in link.connector_ids]
        external=[c.source_id for c in packet.junction.connectors]
        original=packet.original_connector_ids
        if (not packet.links or len(original)!=len(set(original)) or
            len(internal)!=len(set(internal)) or len(external)!=len(set(external)) or
            set(internal)&set(external) or set(internal)|set(external)!=set(original)):
            raise ValueError('Original connector incidences need exactly one internal or external role')
        if tuple(cid for cid in original if cid not in set(internal))!=tuple(external):
            raise ValueError('The external contract must retain the original cyclic connector order')
        sources={s.source_id:s for s in packet.internal_streets}
        if len(sources)!=len(packet.internal_streets) or set(sources)!={link.street_id for link in packet.links}:
            raise ValueError('Every internal Street must have exactly one original input')
        if len({link.street_id for link in packet.links})!=len(packet.links):
            raise ValueError('An internal Street has more than one patch owner')
        for link in packet.links:
            if link.junction_id!=jid or link.bundle_side not in (0,1):
                raise ValueError('An internal link must belong to this declared junction and Bundle end')
            paired=[c for c in packet.junction.connectors if c.bundle_source_id==link.bundle_id]
            if len(paired)!=1 or len(paired[0].members)!=2:
                raise ValueError('An internal link requires its original paired Bundle connector')
            if (int(not paired[0].leaving)!=link.bundle_side or
                tuple(sorted(link.endpoint_members))!=(0,1)):
                raise ValueError('Internal endpoints must name the two members at the declared Bundle end')
            source=snapshot.streets.get(link.street_id)
            if source is None or source.get('status')!='SUPPORTED' or tuple(source.get('endpoints',()))!=(jid,jid):
                raise ValueError('Both original ends of the internal Street must belong to this junction')
            original_sections=source['items']
            supplied=sources[link.street_id].sections
            if (any(s['kind']!='Section' for s in original_sections) or
                tuple(s['id'] for s in original_sections)!=link.section_ids or
                tuple(s.source_id for s in supplied)!=link.section_ids):
                raise ValueError('An internal Street must retain its complete original Section chain')
            for old,new in zip(original_sections,supplied):
                if (tuple(tuple(p) for p in old['points'])!=tuple(tuple(p) for p in new.centerline) or
                    old['width']!=new.width or (old['offset'] or 0.)!=new.offset):
                    raise ValueError('Internal source geometry and dimensions differ from the snapshot')
            endpoints={}
            for member in paired[0].members:
                order=member.head_order if link.bundle_side==0 else member.tail_order
                street=packet.junction.streets[member.street_index]
                section=street.sections[0 if member.leaving else -1]
                endpoints[order]=section.centerline[0 if member.leaving else -1]
            if (set(endpoints)!={0,1} or not supplied or
                any(math.dist(point[:2],endpoints[order][:2])>1e-5 for point,order in zip(
                    (supplied[0].centerline[0],supplied[-1].centerline[-1]),link.endpoint_members))):
                raise ValueError('Original internal endpoints must match their declared Bundle members')
        view.junctions[jid]=packet.junction
    return view


class PreparedCompoundComponentPlanner:
    """Compose the existing complete-network planner with explicit internal ownership."""
    def __init__(self, component_planner):
        self.component_planner=component_planner

    def plan(self, snapshot, seed, packets):
        packets=tuple(packets)
        view=candidate_snapshot(snapshot,packets)
        network=self.component_planner.plan(view,seed)
        junctions=set(network.terminal_junctions)
        selected=tuple(sorted((p for p in packets if p.junction.source_id in junctions),
                              key=lambda p:p.junction.source_id))
        ordinary={s.section_id for s in network.sections}
        internal=set()
        for packet in selected:
            for link in packet.links:
                if link.bundle_id not in network.bundle_ids:
                    raise ValueError('An internal Street references a Bundle outside its component')
                for sid in link.section_ids:
                    if sid in ordinary or sid in internal:
                        raise ValueError('An original Section cannot be swept and absorbed, or owned twice')
                    internal.add(sid)
        return CompoundComponentPlan(network,selected)


@dataclass(frozen=True)
class CompoundComponentInputs:
    network: object
    junction_interiors: tuple  # Complete source packets in the common local frame.


class PreparedCompoundComponentInputs:
    """Localize complete internal sources without turning them into sweep paths."""
    def __init__(self, bounded_adapter):
        self.bounded_adapter=bounded_adapter

    def prepare(self, snapshot, plan):
        view=candidate_snapshot(snapshot,plan.junction_interiors)
        network=self.bounded_adapter.prepare(view,plan.network)
        junctions={j.source_id:j for j in network.junctions}
        origin=network.origin
        localized=[]
        for packet in plan.junction_interiors:
            streets=tuple(replace(street,sections=tuple(replace(section,centerline=tuple(
                (p[0]-origin[0],p[1]-origin[1],p[2] if len(p)>2 else 0.)
                for p in section.centerline)) for section in street.sections))
                for street in packet.internal_streets)
            localized.append(replace(packet,junction=junctions[packet.junction.source_id],
                                     internal_streets=streets))
        return CompoundComponentInputs(network,tuple(localized))
