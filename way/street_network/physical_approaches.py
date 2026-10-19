"""Classify original one-ended approaches without solving their geometry.

Absence from a selected physical group never implies a free street end. The
prepared inventory is authoritative; a dimension change is an explicit handoff
before the far Street end. These records do not certify complete road coverage.
"""
from dataclasses import dataclass
from enum import IntEnum


class ApproachEndRole(IntEnum):
    FREE = 1
    REMOTE_JUNCTION = 2
    INSPECTION_BOUNDARY = 3
    DIMENSION_CHANGE = 4


@dataclass(frozen=True)
class PhysicalApproach:
    street_id: int
    leaving: bool
    section_ids: tuple
    end_role: ApproachEndRole
    owner_id: int = -1
    next_section_id: int = -1


class PreparedPhysicalApproaches:
    """Transport external-chain ownership from the complete prepared inventory."""

    def prepare(self, snapshot, junctions):
        incidences = {}
        for junction in junctions:
            for connector in junction.connectors:
                if len(connector.members) != 1:
                    raise ValueError('A physical approach requires an individual Street incidence')
                member, = connector.members
                street = junction.streets[member.street_index]
                ends = incidences.setdefault(street.source_id, [])
                if any(old[0] == member.leaving for old in ends):
                    raise ValueError('An authored Street end has multiple physical owners')
                if ends and ends[0][1] != street:
                    raise ValueError('Original Street differs between its physical owners')
                ends.append((member.leaving, street))
        result = []
        for sid, ends in sorted(incidences.items()):
            if len(ends) == 2:
                continue
            leaving, street = ends[0]
            record = snapshot.streets.get(sid)
            if record is None or record.get('status') != 'SUPPORTED':
                raise ValueError('Approach needs a supported original Street inventory record')
            endpoints = record.get('endpoints')
            if not isinstance(endpoints, (list, tuple)) or len(endpoints) != 2:
                raise ValueError('Approach needs explicit authored head and tail endpoint roles')
            sequence = street.sections if leaving else street.sections[::-1]
            if not sequence:
                raise ValueError('Approach needs original Sections')
            first = sequence[0]
            selected = []
            following = None
            for section in sequence:
                if (section.width, section.offset) != (first.width, first.offset):
                    following = section
                    break
                selected.append(section.source_id)
            if following is not None:
                role = ApproachEndRole.DIMENSION_CHANGE
                owner, next_section = sid, following.source_id
            else:
                far = endpoints[1 if leaving else 0]
                if far is not None and (type(far) is not int or far < 0):
                    raise ValueError('Far junction identity must be a nonnegative integer or explicit None')
                role = ApproachEndRole.FREE if far is None else ApproachEndRole.REMOTE_JUNCTION
                owner, next_section = (-1 if far is None else far), -1
            # The chain mapping is authored, even when the surface runs backwards.
            if not leaving:
                selected.reverse()
            result.append(PhysicalApproach(sid, leaving, tuple(selected), role, owner, next_section))
        return tuple(result)
