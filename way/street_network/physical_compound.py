"""Semantic contraction of one original Street between two physical sites.

No geometric decision is made here. A native provider decides whether the
original corners can form a valid compound boundary and cover the retained
internal source. Planning never declares that source generated or discards it.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class PhysicalCompoundPlan:
    sites: tuple
    internal_street: object
    internal_indices: tuple
    external_ring: tuple  # (physical site ID, original connector ID)
    external_streets: tuple
    join_radius: float

    @property
    def internal_sections(self):
        return tuple(s.source_id for s in self.internal_street.sections)


class PreparedPhysicalCompound:
    """Keep exact endpoint roles and circular ordering when removing one link."""
    def prepare(self, junctions, street_id, *, join_radius=1.):
        junctions=tuple(junctions)
        if not math.isfinite(join_radius) or join_radius<=0:
            raise ValueError('The additional joining radius must be finite and positive')
        if len({j.source_id for j in junctions})!=len(junctions):
            raise ValueError('Physical sites must have unique identities')
        selected=[]
        for site in junctions:
            hits=[]
            for at,connector in enumerate(site.connectors):
                if len(connector.members)!=1:
                    raise ValueError('Physical rings require individual original Street incidences')
                member,=connector.members
                if not 0<=member.street_index<len(site.streets):
                    raise ValueError('Physical incidence references a missing Street')
                street=site.streets[member.street_index]
                if street.source_id==street_id:hits.append((at,member,street))
            if hits:
                if len(hits)!=1:raise ValueError('A returning Street cannot use the two-site contraction')
                selected.append((site,*hits[0]))
        if len(selected)!=2:raise ValueError('The internal Street needs exactly two distinct physical site owners')
        selected.sort(key=lambda row:row[0].source_id)
        (a,ia,ma,sa),(b,ib,mb,sb)=selected
        if ma.leaving==mb.leaving:raise ValueError('The internal Street requires opposite authored endpoint roles')
        if sa!=sb or not sa.sections:raise ValueError('Physical owners disagree on the complete original Street')
        if len({s.source_id for s in sa.sections})!=len(sa.sections):
            raise ValueError('The internal Street repeats an original Section')
        rings=[];external=[];section_owners={s.source_id:sa.source_id for s in sa.sections}
        for site,internal in ((a,ia),(b,ib)):
            if len(site.connectors)<3:raise ValueError('A physical site requires at least three incidences')
            if len({c.source_id for c in site.connectors})!=len(site.connectors):
                raise ValueError('Physical connector IDs must be unique within their site')
            for step in range(1,len(site.connectors)):
                connector=site.connectors[(internal+step)%len(site.connectors)]
                member,=connector.members;street=site.streets[member.street_index]
                if not street.sections or len({s.source_id for s in street.sections})!=len(street.sections):
                    raise ValueError('An external Street needs a nonempty unique original Section chain')
                rings.append((site.source_id,connector.source_id));external.append(street.source_id)
                for section in street.sections:
                    owner=section_owners.setdefault(section.source_id,street.source_id)
                    if owner!=street.source_id:raise ValueError('An original Section has multiple Street owners')
        if len(set(external))!=len(external):
            raise ValueError('Another shared or returning Street retains a hole; it cannot be absorbed by a single-link contraction')
        return PhysicalCompoundPlan((a.source_id,b.source_id),sa,(ia,ib),tuple(rings),tuple(external),join_radius)
