"""Roles of retained original interior roads in a ground-induced Bundle.

Bridge removal can leave two unrelated approaches outside the finite median.
They remain source roads, with their own scope handoffs and full profiles.
Only a complete original opposite-member crossing belongs to an extent.
"""
import numpy as np
from .model import End


def scoped_interior_roles(extents, interior_runs, sections, sites):
    assignments = [set() for _ in extents]
    records = []
    owners = {e: site for site in sites for e in site.ends}
    for path in interior_runs:
        ids = {t.section_id for t in path}
        matches = []
        for ordinal, members in enumerate(extents):
            roles = {t.section_id: role for role, member in enumerate(members) for t in member.path}
            ports = []
            for sid in sorted(ids):
                for end in (0, 1):
                    site = owners.get(End(sid, end))
                    if site is not None:
                        attached = {roles[e.section_id] for e in site.ends if e.section_id in roles}
                        ports.extend((role, sid, end, site.key.text()) for role in attached)
            if ports:
                if len(ports) != 2 or sorted(p[0] for p in ports) != [0, 1]:
                    raise ValueError('Ground interior run needs two complete opposite-member ports or an exterior proof')
                matches.append((ordinal, ports))
        if len(matches) > 1:
            raise ValueError('Ground interior run has competing finite Bundle extents')
        if matches:
            ordinal, ports = matches[0]
            assignments[ordinal].update(ids)
            records.append(dict(section_ids=sorted(ids), role='GROUND_INTERIOR_CROSSING',
                                extent_index=ordinal, ports=ports))
            continue
        proofs = [_exterior_controls(path, members, sections) for members in extents]
        if any(proof is None for proof in proofs):
            raise ValueError('Retained ground interior run intrudes a finite Bundle without complete source crossing ports')
        records.append(dict(section_ids=sorted(ids), role='RETAINED_EXTERIOR_GROUND_RUN',
            scope_ends=[dict(section=sid, end=end) for sid in sorted(ids) for end in sections[sid].scope_ends],
            source_incidence='No retained original incidence reaches a lateral member',
            finite_source_selection=proofs, numerical_geometry_qualification_required=True))
    return tuple(tuple(sorted(ids)) for ids in assignments), records


def _controls(path, sections):
    points = []
    for traversal in path:
        piece = np.asarray(sections[traversal.section_id].points)[:, :2]
        if not traversal.forward:
            piece = piece[::-1]
        if points:
            if np.linalg.norm(points[-1][-1]-piece[0]) > 1e-7:
                raise ValueError('Ground interior role proof has discontinuous original controls')
            piece = piece[1:]
        points.append(piece)
    return np.concatenate(points)


def _exterior_controls(path, members, sections):
    """Finite source-control corridor disjointness, never new incidence.

    This selects the original role only. Curved road/median surfaces retain
    their mandatory normal numerical clearance check before publication.
    """
    from .regions import _classify
    a, b = (_controls(member.path, sections) for member in members)
    rim = np.concatenate((a, b[::-1]))
    source = _controls(path, sections)
    # Bound the paired source footprint and the retained road's full width;
    # no centreline or pavement is clipped to obtain this classification.
    pad = max(sections[t.section_id].width/2 + abs(sections[t.section_id].offset)
              for member in members for t in member.path)
    pad += max(sections[t.section_id].width/2 + abs(sections[t.section_id].offset) for t in path)
    inside, boundary = _classify(source[:1], rim, 1e-7)
    if inside[0] or boundary[0] >= 0:
        return None
    if (len(source)-1)*len(rim) > 200_000:
        raise ValueError('Ground interior exterior proof exceeds its source segment-pair budget')
    minimum = float('inf')
    for p, q in zip(source, source[1:]):
        for a, b in zip(rim, np.roll(rim, -1, axis=0)):
            minimum = min(minimum, _segment_distance(p, q, a, b))
            if minimum <= pad + 1e-7:
                return None
    return dict(method='finite original-control corridor separation',
                minimum_source_distance=minimum, original_width_allowance=pad,
                source_footprint_clearance=minimum-pad)


def _segment_distance(a, b, c, d):
    def cross(v, w):
        return float(v[0]*w[1]-v[1]*w[0])
    ab, cd = b-a, d-c
    denominator = cross(ab, cd)
    if abs(denominator) > 1e-14:
        t, u = cross(c-a, cd)/denominator, cross(c-a, ab)/denominator
        if 0. <= t <= 1. and 0. <= u <= 1.:
            return 0.
    def distance(p, q, r):
        delta = r-q
        t = np.clip(float((p-q)@delta)/max(float(delta@delta), 1e-30), 0., 1.)
        return float(np.linalg.norm(p-q-t*delta))
    return min(distance(a,c,d), distance(b,c,d), distance(c,a,b), distance(d,a,b))
