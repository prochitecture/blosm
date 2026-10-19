"""Prepare once-owned physical sites behind aggregated Bundle junctions.

This opt-in adapter retains the complete original physical ring. It does not
promote a site to supported, absorb connecting Streets, flatten half profiles,
sort directions, or construct surface coordinates. The injected native
contract performs its own capability checks after this semantic transport.
"""
from dataclasses import dataclass
import math


class UnresolvedBundleProfile(ValueError):
    """Original membership is retained, but no supported lateral role exists."""


@dataclass(frozen=True)
class PhysicalJunctionInputs:
    logical_junction_ids: tuple
    junctions: tuple
    endpoint_bundles: tuple  # Physical junction, connector, longitudinal Bundle role (-1 otherwise).


@dataclass(frozen=True)
class PhysicalPortDesign:
    junction_id: int
    street_id: int
    section_id: int
    end: int
    radius: float
    rise: float
    half_profile: int | None
    datum: float = 0.
    bundle_id: int = -1
    profile_error: str | None = None


class PreparedPhysicalJunctions:
    def __init__(self, contract, section_reader, category_properties, *, rise=.08, cell_size=.7):
        if not math.isfinite(rise) or rise < 0 or not math.isfinite(cell_size) or cell_size <= 0:
            raise ValueError('Profile rise and cell size must be finite with nonnegative rise and positive size')
        self.contract=contract
        self.section_reader=section_reader
        self.category_properties=category_properties
        self.rise=rise
        self.cell_size=cell_size

    def prepare_port_designs(self, snapshot, plan, *, retain_unresolved=False):
        """Apply existing profile/category choices to original Section-end ports.

        Internal through ports retain the complete Street's Bundle side role.
        Each Section supplies its own category radius; no shorter-link radius
        or replacement Street is synthesized to make the geometry fit.
        """
        result=[]
        for site in plan.sites:
            for port in site.ports:
                endpoint=port.endpoint;source=snapshot.streets[endpoint.street_id]
                sections=[s for s in source['items'] if s['kind']=='Section' and s['id']==endpoint.section_id]
                supported=(plan.source_supported(source) if hasattr(plan,'source_supported') else source['status']=='SUPPORTED')
                if not supported or len(sections)!=1:
                    raise ValueError('Port design needs one supported original Section')
                section=sections[0]
                try:radius=self.category_properties[section['category']]['radius']
                except KeyError:raise ValueError('Physical port category has no declared radius') from None
                if not math.isfinite(radius) or radius<=0:
                    raise ValueError('Physical port radius must be positive and finite')
                try:
                    half,bid=self._half(snapshot,source,endpoint.end)
                except ValueError as exc:
                    message = (
                        f'Junction {site.junction_id}, Street {endpoint.street_id}, '
                        f'Section {endpoint.section_id} end {endpoint.end}, '
                        f'Bundle {source.get("bundle_id")}: {exc}')
                    if retain_unresolved and isinstance(exc, UnresolvedBundleProfile):
                        result.append(PhysicalPortDesign(site.junction_id, endpoint.street_id,
                            endpoint.section_id, endpoint.end, float(radius), self.rise, None,
                            bundle_id=source['bundle_id'], profile_error=message))
                        continue
                    raise ValueError(message) from exc
                result.append(PhysicalPortDesign(site.junction_id,endpoint.street_id,endpoint.section_id,
                    endpoint.end,float(radius),self.rise,half,bundle_id=bid))
        return tuple(result)

    def _half(self, snapshot, street, end):
        bid=street.get('bundle_id')
        if bid is None:
            return 0,-1
        bundle=snapshot.bundles.get(bid)
        if bundle is None:
            raise ValueError('Physical endpoint references an absent Bundle')
        topology=bundle.get('topology',{})
        members=topology.get('members',())
        roles=[m for m in members if m['street_id']==street['id']]
        interiors=tuple(topology.get('interior_streets',()))
        chains=topology.get('chains',())
        if chains:
            ids=[sid for chain in chains for sid in chain]
            if (len(chains)!=2 or any(not chain for chain in chains) or
                    len(ids)!=len(set(ids)) or len(members)!=len(ids) or
                    set(ids)!={m['street_id'] for m in members} or
                    any(m['head_order']!=i or m['tail_order']!=i
                        for i,chain in enumerate(chains) for m in members if m['street_id'] in chain)):
                raise ValueError('Physical Bundle chains need distinct complete lateral roles')
        if topology.get('status')!='MAPPED' or (not chains and len(members)!=2):
            reason = topology.get('reason') or bundle.get('reason')
            raise UnresolvedBundleProfile('Physical Bundle endpoint needs a mapped two-member profile role'
                + (': '+reason if reason else ''))
        if street['id'] in interiors:
            if roles or interiors.count(street['id'])!=1:
                raise ValueError('An interior Street cannot also own a longitudinal Bundle profile')
            # Bundle area ownership is not a lateral carriageway role. An
            # explicitly recorded interior connecting Street keeps its full
            # crown and must not be paired as a longitudinal half-profile.
            # The original source Bundle ID remains in the input snapshot.
            return 0,-1
        if len(roles)!=1 or not roles[0]['complete']:
            raise ValueError('Physical Bundle endpoint needs an explicit complete two-member profile role')
        member=roles[0]
        if member['head_order'] not in (0,1) or member['head_order']!=member['tail_order']:
            raise ValueError('Physical Bundle endpoint has inconsistent lateral ownership')
        inner_u=float((member['head_order']==0)==member['forward'])
        return int((2.*inner_u-1.)*(1. if end==0 else -1.)),bid

    def prepare(self, snapshot, logical_junction_id):
        origins=snapshot.junction_origins
        components=[c for c in origins.get('components',()) if logical_junction_id in c['logical_junction_ids']]
        if len(components)!=1:
            raise ValueError('Logical junction needs one recorded physical component')
        component=components[0]
        if component['status']!='MAPPED':
            raise ValueError('Physical component has unresolved original endpoint incidences')
        return self.prepare_sites(snapshot, component['physical_junction_ids'],
                                  logical_junction_ids=component['logical_junction_ids'])

    def prepare_sites(self, snapshot, physical_junction_ids, *, logical_junction_ids=()):
        """Read explicitly owned physical sites, including ordinary junctions.

        This transports original rings, not an assertion that a selected group
        is a complete logical component. Each endpoint and source is checked by
        the same rules used by aggregate-junction preparation.
        """
        ids=tuple(physical_junction_ids)
        if not ids or len(set(ids))!=len(ids):
            raise ValueError('Physical site selection must be nonempty and unique')
        origins=snapshot.junction_origins
        physical={j['junction_id']:j for j in origins['physical']}
        if len(physical)!=len(origins['physical']) or any(jid not in physical for jid in ids):
            raise ValueError('Physical site selection has missing or ambiguous original records')
        k=self.contract; junctions=[];bundles=[];used=set()
        for jid in ids:
            site=physical[jid];streets=[];connectors=[];profiles=[]
            for connector in site['ring']:
                if connector['item_kind']!='Street' or len(connector['endpoints'])!=1:
                    raise ValueError('Physical ring must retain individual original Street endpoint incidences')
                endpoint=connector['endpoints'][0];sid=endpoint['street_id'];end=endpoint['end']
                source=snapshot.streets.get(sid)
                if source is None or source['status']!='SUPPORTED':
                    raise ValueError('Physical ring contains a missing or unsupported original Street '+str(sid))
                sections=tuple(self.section_reader(snapshot,sid))
                if not sections or end not in (0,1) or endpoint['terminal_kind']!='Section':
                    raise ValueError('Physical endpoint must identify an original terminal Section')
                section=sections[0 if end==0 else -1]
                if section.source_id!=endpoint['terminal_id']:
                    raise ValueError('Original physical endpoint terminal changed during aggregation')
                key=(section.source_id,end)
                if key in used:
                    raise ValueError('One original Section endpoint has multiple physical surface owners')
                used.add(key)
                point=section.centerline[0 if end==0 else -1]
                if math.dist(point[:2],site['location'][:2])>1e-5:
                    raise ValueError('Original Section endpoint does not meet its recorded physical site')
                half,bid=self._half(snapshot,source,end)
                try:radius=self.category_properties[section.category]['radius']
                except KeyError:raise ValueError('Physical approach category has no declared radius') from None
                if not math.isfinite(radius) or radius<=0:
                    raise ValueError('Physical approach radius must be positive and finite')
                index=len(streets);streets.append(k.StreetInput(sid,sections))
                profiles.append(k.ProfileInput(self.rise,half_profile=half))
                connectors.append(k.ConnectorInput(connector['connector_id'],
                    (k.MemberInput(index,end==0),),index,radius=radius,leaving=end==0))
                bundles.append((jid,connector['connector_id'],bid))
            junctions.append(k.PreparedJunction(jid,site['location'],tuple(streets),
                tuple(connectors),tuple(profiles),cell_size=self.cell_size))
        return PhysicalJunctionInputs(tuple(logical_junction_ids),tuple(junctions),tuple(bundles))
