"""Translate a semantic bounded Bundle plan into native input contracts.

No bpy, curve samples, mouth stations or surface coordinates are constructed.
The contract module is injected so the candidate kernel can be qualified before
it replaces the released asset. Coordinates are only translated to a local
metric frame; original Street, Section, connector IDs and radii are retained.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class BoundedSeam:
    before_section: int
    after_section: int
    transition_id: int = -1
    minimum_half_length: float = 1.5
    maximum_edge_slope: float = .3
    source_kind: str = ''
    source_id: int = -1


@dataclass(frozen=True)
class BoundedInputs:
    junctions: tuple
    sections: dict
    section_profiles: dict
    seams: tuple[BoundedSeam, ...]
    origin: tuple
    section_bundles: dict
    free_ends: tuple = ()


class PreparedBoundedInputs:
    def __init__(self, contract, *, profile_rise=.08, datum=0.):
        if not math.isfinite(profile_rise) or profile_rise<0. or not math.isfinite(datum):
            raise ValueError('Free member profile rise and datum must be finite; rise must be nonnegative')
        self.contract = contract
        self.profile_rise=profile_rise
        self.datum=datum

    def prepare(self, snapshot, plan):
        k = self.contract
        ids = {j for j in plan.terminal_junctions if j is not None} | {p.junction_id for p in plan.junctions}
        ids.update(getattr(plan,'ordinary_junctions',()))
        if not ids <= snapshot.junctions.keys():
            raise ValueError('Every planned junction needs its prepared input contract')
        source = {s['id']: (street['id'], s) for street in snapshot.streets.values()
            for s in street['items'] if s['kind'] == 'Section'}
        first=next((j for j in plan.terminal_junctions if j is not None),min(ids) if ids else None)
        origin = tuple(snapshot.junctions[first].location[:2]) if ids else tuple(
            source[min(r.section_id for r in plan.sections)][1]['points'][0][:2])
        datums = {snapshot.junctions[jid].datum for jid in ids}
        if not datums:datums={self.datum}
        if len(datums) != 1:
            raise ValueError('Bounded native input currently needs one common flat datum')
        datum = datums.pop()
        roles = {r.section_id:r for r in plan.sections}
        sections = {}
        for sid, role in roles.items():
            owner, s = source[sid]
            if owner != role.street_id:
                raise ValueError('Prepared Section ownership differs from the plan')
            points = tuple((p[0]-origin[0], p[1]-origin[1], p[2] if len(p)>2 else 0.) for p in s['points'])
            sections[sid] = k.SectionInput(sid, points, s['width'], s['offset'] or 0.,
                s['category'], tags=tuple(sorted(s.get('tags', {}).items())))
        junctions, rises = [], {}
        for jid in sorted(ids):
            j = snapshot.junctions[jid]
            streets = [None]*len(j.streets)
            connectors, profiles = [], []
            for c in j.connectors:
                original_profile = j.profiles[c.profile_index]
                half = 0
                members = []
                for m in c.members:
                    street = j.streets[m.street_index]
                    sid = street.sections[0 if m.leaving else -1].source_id
                    role = roles[sid]
                    if role.street_id != street.source_id:
                        raise ValueError('A junction member references a different source Street')
                    owner = role.start_junction if m.leaving else role.end_junction
                    if owner != jid:
                        raise ValueError('Junction member and planned Section end owner disagree')
                    s = sections[sid]
                    points = s.centerline[:2] if m.leaving else s.centerline[-2:]
                    local = k.SectionInput(s.source_id, points, s.width, s.offset, s.category, tags=s.tags)
                    if streets[m.street_index] is not None:
                        raise ValueError('Each local Street entry must have one explicit member reference')
                    streets[m.street_index] = k.StreetInput(street.source_id, (local,))
                    members.append(k.MemberInput(m.street_index,m.leaving,m.lateral_order,
                        m.bundle_index,m.head_order,m.tail_order))
                    if len(c.members)==1 and role.half_profile:
                        half = int((2.*role.inner_u-1.) * (1. if m.leaving else -1.))
                    rises.setdefault(role.street_id,set()).add(original_profile.rise)
                profiles.append(k.ProfileInput(original_profile.rise,original_profile.preset,
                    original_profile.separator,half_profile=half))
                connectors.append(k.ConnectorInput(c.source_id,tuple(members),len(profiles)-1,
                    c.radius,c.item_kind,c.leaving,c.bundle_source_id))
            if any(s is None for s in streets):
                raise ValueError('Prepared junction contains an unreferenced Street entry')
            junctions.append(k.PreparedJunction(jid,(j.location[0]-origin[0],j.location[1]-origin[1]),
                tuple(streets),tuple(connectors),tuple(profiles),datum,j.cell_size))
        if any(len(values)!=1 for values in rises.values()):
            raise ValueError('A continuous Street currently needs one profile rise')
        bundle_owners=(dict(plan.section_bundles) if hasattr(plan,'section_bundles') else
            {r.section_id:plan.bundle_id for r in plan.sections if r.half_profile})
        bundle_rises={}
        for sid,bid in bundle_owners.items():
            bundle_rises.setdefault(bid,set()).update(rises.get(roles[sid].street_id,()))
        for role in roles.values():
            if role.street_id not in rises:
                values=bundle_rises.get(bundle_owners.get(role.section_id),set())
                if len(values)>1:raise ValueError('A free member needs one consistent Bundle profile rise')
                rises[role.street_id]=values or {self.profile_rise}
        profiles = {sid:(role.half_profile,role.inner_u,next(iter(rises[role.street_id])),datum)
            for sid,role in roles.items()}
        transitions = {(i['transition']['before_section'],i['transition']['after_section']):i
            for street in snapshot.streets.values() for i in street['items']
            if i.get('transition')}
        seams = []; next_transition = 0
        for seam in plan.seams:
            item = transitions.get((seam.before_section,seam.after_section))
            if item:
                t = item['transition']
                # Native transition IDs occupy one private table namespace.
                # Original SymLane / SideLane / Intersection IDs may coincide.
                seams.append(BoundedSeam(seam.before_section,seam.after_section,
                    next_transition,
                    t['minimum_half_length'],t['maximum_edge_slope'],t['source_kind'],t['source_id']))
                next_transition += 1
            else:
                seams.append(BoundedSeam(seam.before_section,seam.after_section))
        return BoundedInputs(tuple(junctions),sections,profiles,tuple(seams),origin,bundle_owners,plan.free_ends)
