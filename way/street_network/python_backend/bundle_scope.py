"""Explicit ground-only extents of original paired lateral Bundle chains."""
from .model import BundleMember


def ground_bundle_scope(identity, members, sections, excluded_grade):
    """Restrict matching outer ends, never fill gaps or infer new incidence.

    The original complete paths remain in the diagnostic. A ground member is
    complete for this provider's finite extent only; its cut must already be
    recorded as a grade-scope endpoint by the source-incidence adapter.
    """
    if len(members) != 2 or not all(member.complete for member in members):
        raise ValueError('Ground Bundle scope requires two complete original lateral chains')
    ground, records, masks, boundaries = [], [], [], []
    for role, member in enumerate(members):
        path = member.path
        missing = [t.section_id for t in path if t.section_id not in sections]
        if any(sid not in excluded_grade for sid in missing):
            raise ValueError('Bundle member includes an invalid or unexplained absent Section')
        retained = [i for i, traversal in enumerate(path) if traversal.section_id in sections]
        if not retained:
            raise ValueError('Each ground Bundle member needs a nonempty original ground run')
        start, last = retained[0], retained[-1]
        if retained != list(range(start, last + 1)):
            return _split_ground_scope(identity, members, sections, excluded_grade)
        mask = (start > 0, last < len(path) - 1)
        masks.append(mask)
        selected = path[start:last + 1]
        for is_cut, label, current, deferred, end in (
                (mask[0], 'HEAD', selected[0], path[start - 1] if mask[0] else None,
                 0 if selected[0].forward else 1),
                (mask[1], 'TAIL', selected[-1], path[last + 1] if mask[1] else None,
                 1 if selected[-1].forward else 0)):
            if not is_cut:
                continue
            if end not in sections[current.section_id].scope_ends:
                raise ValueError('Ground Bundle cut needs an existing original grade-scope endpoint')
            deferred_end = ((1 if deferred.forward else 0) if label == 'HEAD'
                            else (0 if deferred.forward else 1))
            boundaries.append(dict(member=role, bundle_end=label, section=current.section_id,
                                   end=end, deferred_section=deferred.section_id, deferred_end=deferred_end))
        ground.append(BundleMember(selected, complete=True))
        records.append(dict(member=role, original=[dict(section=t.section_id, forward=t.forward) for t in path],
                            ground=[dict(section=t.section_id, forward=t.forward) for t in selected]))
    if masks[0] != masks[1] or not any(masks[0]):
        raise ValueError('Ground Bundle scope requires corresponding HEAD/TAIL cuts on both members')
    original_ids = sorted(t.section_id for member in members for t in member.path)
    ground_ids = sorted(t.section_id for member in ground for t in member.path)
    deferred_ids = sorted(set(original_ids) - set(ground_ids))
    return tuple(ground), dict(kind='ScopeBundle', id=identity,
        reason='Paired original Bundle ends are explicitly deferred outside ground scope',
        policy='PAIRED_CONTIGUOUS_GROUND_EXTENT', complete_in_ground_scope=True,
        original_section_ids=original_ids, ground_section_ids=ground_ids,
        deferred_section_ids=deferred_ids, members=records, boundaries=boundaries,
        deferred=[dict(section=sid, reason=excluded_grade[sid]) for sid in deferred_ids])


def _split_ground_scope(identity, members, sections, excluded_grade):
    """Select one unique finite alongside pair; never bridge a removed span.

    A source member split by a bridge may have one ground run alongside the
    other member and another beyond it. The latter remains an ordinary road.
    The source controls prove this selection only; normal path, finite-cap,
    surface-quality and global pavement-clearance checks still qualify output.
    """
    scopes = _split_ground_scopes(identity, members, sections, excluded_grade)
    if len(scopes) != 1:
        raise ValueError('Split ground Bundle needs one unique finite alongside run pair')
    return scopes[0]


def ground_bundle_scopes(identity, members, sections, excluded_grade):
    """All unambiguous disconnected ground extents in original source order.

    Each run may participate once, and the paired run orders must agree. A
    one-to-many overlap remains ambiguous rather than arbitrarily partitioned.
    """
    split = False
    for member in members:
        indices = [i for i, t in enumerate(member.path) if t.section_id in sections]
        split |= bool(indices and indices != list(range(indices[0], indices[-1]+1)))
    if not split:
        return (ground_bundle_scope(identity, members, sections, excluded_grade),)
    return _split_ground_scopes(identity, members, sections, excluded_grade)


def _split_ground_scopes(identity, members, sections, excluded_grade):
    if len(members) != 2 or not all(member.complete for member in members):
        raise ValueError('Ground Bundle scope requires two complete original lateral chains')
    runs, cuts = [], []
    for role, member in enumerate(members):
        missing = [t.section_id for t in member.path if t.section_id not in sections]
        if any(sid not in excluded_grade for sid in missing):
            raise ValueError('Bundle member includes an invalid or unexplained absent Section')
        intervals, start = [], None
        for index, traversal in enumerate((*member.path, None)):
            present = traversal is not None and traversal.section_id in sections
            if present and start is None:
                start = index
            elif not present and start is not None:
                intervals.append((start, index))
                start = None
        if not intervals or len(intervals) > 32:
            raise ValueError('Split ground Bundle needs a bounded nonempty run inventory')
        role_runs, role_cuts = [], []
        for first, stop in intervals:
            path = member.path[first:stop]
            boundaries = []
            for is_cut, label, current, deferred, end in (
                (first > 0, 'HEAD', path[0], member.path[first-1] if first else None,
                 0 if path[0].forward else 1),
                (stop < len(member.path), 'TAIL', path[-1], member.path[stop] if stop < len(member.path) else None,
                 1 if path[-1].forward else 0)):
                if not is_cut:
                    continue
                if end not in sections[current.section_id].scope_ends:
                    raise ValueError('Ground Bundle cut needs an existing original grade-scope endpoint')
                boundaries.append(dict(member=role, bundle_end=label, section=current.section_id,
                    end=end, deferred_section=deferred.section_id,
                    deferred_end=(1 if deferred.forward else 0) if label == 'HEAD' else (0 if deferred.forward else 1)))
            role_runs.append(path)
            role_cuts.append(boundaries)
        runs.append(role_runs)
        cuts.append(role_cuts)
    candidates = []
    for a, first in enumerate(runs[0]):
        for b, second in enumerate(runs[1]):
            proof = _common_fibres(first, second, sections)
            if proof is not None:
                candidates.append((a, b, proof))
    if (not candidates or len({a for a, b, _ in candidates}) != len(candidates) or
            len({b for a, b, _ in candidates}) != len(candidates) or
            [b for _, b, _ in candidates] != sorted(b for _, b, _ in candidates)):
        raise ValueError('Split ground Bundle needs unique finite alongside run pairs in source order')
    retained_ids = {t.section_id for a, b, _ in candidates for path in (runs[0][a], runs[1][b]) for t in path}
    original_ids = {t.section_id for member in members for t in member.path}
    unpaired = (original_ids & sections.keys()) - retained_ids
    missing = original_ids - sections.keys()
    result = []
    for a, b, proof in candidates:
      selected = (runs[0][a], runs[1][b])
      current_ids = {t.section_id for path in selected for t in path}
      record = dict(
        kind='ScopeBundle', id=identity,
        reason='Original ground runs have unique finite alongside support; no removed span is bridged',
        policy='UNIQUE_FINITE_GROUND_RUN_PAIR' if len(candidates) == 1 else 'DISJOINT_FINITE_GROUND_RUN_PAIRS',
        complete_in_ground_scope=True, paired_extent_count=len(candidates),
        original_section_ids=sorted(original_ids), ground_section_ids=sorted(current_ids),
        deferred_section_ids=sorted(missing), unpaired_section_ids=sorted(unpaired),
        members=[dict(member=role,
            original=[dict(section=t.section_id, forward=t.forward) for t in member.path],
            ground=[dict(section=t.section_id, forward=t.forward) for t in selected[role]])
            for role, member in enumerate(members)],
        boundaries=cuts[0][a]+cuts[1][b],
        unpaired_ground_runs=[dict(member=role, path=[dict(section=t.section_id, forward=t.forward) for t in path],
            boundaries=cuts[role][index]) for role in (0, 1) for index, path in enumerate(runs[role])
            if any(t.section_id in unpaired for t in path)],
        finite_source_selection=proof,
        deferred=[dict(section=sid, reason=excluded_grade[sid]) for sid in sorted(missing)])
      result.append((tuple(BundleMember(path, complete=True) for path in selected), record))
    return tuple(result)


def _common_fibres(first, second, sections):
    """Positive finite axial overlap and strict authored lateral order.

    Original controls, not distances between centroids or normalized paths,
    provide a monotone piecewise-linear selection proof. No geometry is moved.
    Curved interpolation and its width envelopes are independently certified by
    the numerical Bundle provider before publication.
    """
    import numpy as np
    curves = []
    for path in (first, second):
        pieces = []
        for traversal in path:
            points = np.asarray(sections[traversal.section_id].points)[:, :2]
            if not traversal.forward:
                points = points[::-1]
            if pieces:
                if np.linalg.norm(pieces[-1][-1]-points[0]) > 1e-7:
                    raise ValueError('Ground Bundle run has discontinuous original controls')
                points = points[1:]
            pieces.append(points)
        curves.append(np.concatenate(pieces))
    direction = [p[-1]-p[0] for p in curves]
    if min(np.linalg.norm(d) for d in direction) <= 1e-7 or direction[0]@direction[1] <= 0.:
        return None
    axis = sum(d/np.linalg.norm(d) for d in direction)
    axis /= np.linalg.norm(axis)
    left = np.array((-axis[1], axis[0]))
    origin = curves[0][0]
    local = [np.column_stack(((p-origin)@axis, (p-origin)@left)) for p in curves]
    if any(np.min(np.diff(p[:, 0])) <= 1e-7 for p in local):
        return None
    lo, hi = max(p[0, 0] for p in local), min(p[-1, 0] for p in local)
    if hi <= lo+1e-7:
        return None
    nodes = np.unique(np.r_[lo, hi, *(p[(p[:, 0] > lo) & (p[:, 0] < hi), 0] for p in local)])
    gap = np.interp(nodes, *local[0].T)-np.interp(nodes, *local[1].T)
    if np.min(gap) <= 1e-7:
        return None
    return dict(method='finite original-control common fibres', axis=axis.tolist(),
                overlap_length=float(hi-lo), minimum_centerline_gap=float(np.min(gap)),
                laterality_preserved=True, numerical_geometry_qualification_required=True)
