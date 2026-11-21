"""Explicit, source-incidence preparation of one-sided merged junctions.

The merged-entrance provider invokes this preparation for multi-Section runs.
It returns fresh source objects before ownership or any surface is built.
"""
from dataclasses import replace
import numpy as np

from .bundle_axis import AxialChartUnavailable
from .bundle_angular_curves import RayCurb
from .bundle_merged_entrances import MergedEntranceSpecification
from .merged_run_inventory import inventory_frame
from .model import End, Key
from .ownership import footprint_key
from .paired_source_sampling import (OneSidedPreparationRequest, canonical_reprepare,
    shared_source_parameters, verify_shared_source_phase)
from .sampling import CertifiedCurbSampler


class PairedJunctionPreparationUnavailable(AxialChartUnavailable):
    """Original incidence or certified opposite support needs another strategy."""


def _runs(context, bundle, terminal):
    if terminal.bundle_id != bundle.id or len(bundle.members) != 2 or not all(m.complete for m in bundle.members):
        raise PairedJunctionPreparationUnavailable('Paired preparation needs the same complete original Bundle')
    declared = [item for item in context.network.bundles if item.id == bundle.id]
    if len(declared) != 1 or declared[0] != bundle:
        raise PairedJunctionPreparationUnavailable('Paired preparation changed the original network Bundle')
    result, seen = [], set()
    originals = tuple(tuple(m.path) for m in bundle.members)
    for approach in terminal.layout.approaches:
        run = getattr(approach, 'run', None)
        if run is None:
            continue
        signature = tuple(tuple(m.visits) for m in run.members)
        if signature in seen:
            continue
        if len(signature) != 2:
            raise PairedJunctionPreparationUnavailable('A paired terminal lost its two original member inventories')
        for role, visits in enumerate(signature):
            if not visits or not any(originals[role][i:i+len(visits)] == visits
                                     for i in range(len(originals[role])-len(visits)+1)):
                raise PairedJunctionPreparationUnavailable('A paired run is not an ordered original Bundle subsequence')
            for visit in visits:
                path = context[Key('path', visit.section_id)]
                if path.section.id != visit.section_id:
                    raise ValueError('Paired preparation lost its original Path identity')
        seen.add(signature)
        result.append((run, approach))
    if not result:
        raise PairedJunctionPreparationUnavailable('Paired preparation has no original terminal run descriptors')
    return tuple(result)


def _inner_curve(layout, near_ends, tolerance):
    indexes = {i for i, a in enumerate(layout.approaches) if End(a.source.id, a.end) in near_ends}
    if len(indexes) != 2 or len(layout.approaches) != 3:
        raise PairedJunctionPreparationUnavailable('One-sided discovery needs both original near mouths and one stem')
    candidates = [i for i in indexes if (i+1) % len(layout.approaches) in indexes]
    if len(candidates) != 1:
        raise PairedJunctionPreparationUnavailable('An original near-side junction lacks a unique direct curb')
    index = candidates[0]
    curve = layout.boundaries[index]
    endpoints = []
    for i in (index, (index+1) % len(layout.approaches)):
        a = layout.approaches[i]
        if abs(a.profile.peak) != 1.:
            raise PairedJunctionPreparationUnavailable('One-sided discovery needs the original Bundle half profiles')
        endpoints.append(a.point(layout.outward_stations[i], a.offset+.5*a.width*a.profile.peak))
    if np.max(np.linalg.norm(curve(np.array((0., 1.)))-endpoints, axis=1)) > 1e-7:
        raise PairedJunctionPreparationUnavailable('The original stem interrupts the facing near-side curb')
    parameters = CertifiedCurbSampler(tolerance, .25, 20000).sample(curve, 0., 1.)
    return curve(parameters)


def _opposite(context, visits, frame, points, tolerance, role):
    stations, transverse = frame.local(points).T
    candidates, records = [], []
    for visit in visits:
        path = context[Key('path', visit.section_id)]
        frame.certify(path, forward=visit.forward)
        peak = context[Key('profile', visit.section_id)].peak
        if abs(peak) != 1.:
            raise PairedJunctionPreparationUnavailable('An opposite original member lost its authored half profile')
        samples = path.samples(0., path.length, min(.25, context.network.settings.spacing), tolerance=tolerance)
        xy = path.point(samples, path.section.offset+.5*path.section.width*peak)
        if not visit.forward:
            xy = xy[::-1]
        ray = RayCurb(frame, xy)
        hits = (stations >= ray.stations[0]) & (stations <= ray.stations[-1])
        records.append((visit, ray, hits))
        if np.all(hits):
            candidates.append((visit, ray))
    if len(candidates) != 1:
        raise PairedJunctionPreparationUnavailable('Facing original curb needs one explicit opposite source interval')
    visit, ray = candidates[0]
    if any(np.any(hits) for other, other_ray, hits in records if other != visit):
        raise PairedJunctionPreparationUnavailable('Facing original curb overlaps ambiguous opposite source intervals')
    opposite = frame.local(ray(stations))[:, 1]
    gaps = transverse-opposite if role == 0 else opposite-transverse
    if float(gaps.min()) < -2*tolerance:
        return visit
    if float(gaps.min()) <= 2*tolerance:
        raise PairedJunctionPreparationUnavailable('Original source contact is uncertain within its certified chord budget')
    return None


def discover(context, bundle, terminal):
    """Find original intermediate incidences, never connect nearby coordinates.

    Opposite selection uses the certified common chart and complete source
    intervals. A proposed junction must exhibit actual overlap beyond the
    source approximation band. Unsupported/ambiguous cases remain errors.
    """
    by_end = {}
    for site in (*context.network.intersections, *context.network.connections):
        for endpoint in site.ends:
            if endpoint in by_end:
                raise ValueError('An original endpoint has multiple authored incidence owners')
            by_end[endpoint] = site
    tolerance = min(.00025, .0125*context.network.settings.curve_tolerance)
    requests, seen = [], set()
    for run, approach in _runs(context, bundle, terminal):
        frame = inventory_frame(context, bundle, approach)
        for role, member in enumerate(run.members):
            for before, after in zip(member.visits, member.visits[1:]):
                exit = End(before.section_id, int(before.forward))
                entry = End(after.section_id, int(not after.forward))
                site = by_end.get(exit)
                if site is None or entry not in site.ends:
                    raise PairedJunctionPreparationUnavailable('Consecutive original visits lack their authored incidence')
                if len(site.ends) == 2:
                    continue
                if len(site.ends) != 3:
                    raise PairedJunctionPreparationUnavailable('One-sided discovery needs an original three-arm junction')
                other = run.members[1-role].visits
                if any(end.section_id in {v.section_id for v in other} for end in site.ends):
                    raise PairedJunctionPreparationUnavailable('An intermediate paired authored junction needs its own preparation')
                if site.key in seen:
                    raise ValueError('One-sided discovery repeats an original junction across source runs')
                seen.add(site.key)
                layout = context[footprint_key(site.key)]
                if getattr(layout, 'one_sided_source_preparation', False):
                    raise ValueError('Discovery requires the unmodified original local footprint')
                points = _inner_curve(layout, (exit, entry), tolerance)
                opposite = _opposite(context, other, frame, points, tolerance, role)
                if opposite is not None:
                    requests.append(OneSidedPreparationRequest(site, layout, (before, after), opposite, role))
    if len(requests) > 32:
        raise PairedJunctionPreparationUnavailable('One-sided discovery exceeds the bounded original-site budget')
    return tuple(sorted(requests, key=lambda request: request.site.key))


def prepare(context, bundle, terminal, *, requests=None):
    """Return fresh terminal, side specs, pooled global knots and joint proof.

    Rebuild every corner, source union and complete physical mouth. No accepted
    Product is edited, snapped, or reused as a newly prepared physical port.
    The normal provider calls this before physical ownership and hydration.
    """
    runs = tuple(run for run, approach in _runs(context, bundle, terminal))
    fresh = MergedEntranceSpecification.build(context, bundle, terminal.crossing, runs)
    if fresh is None:
        raise PairedJunctionPreparationUnavailable('Original paired terminal no longer qualifies its source frame')
    if requests is None:
        requests = discover(context, bundle, fresh)
    else:
        requests = tuple(requests)
        for request in requests:
            if request.near_role not in (0, 1) or len(request.near_visits) != 2:
                raise PairedJunctionPreparationUnavailable('Explicit preparation changed the original near-source role')
            near = tuple(request.near_visits)
            matches = [(run, a) for run, a in _runs(context, bundle, fresh)
                       if any(tuple(run.members[request.near_role].visits[i:i+2]) == near
                              for i in range(len(run.members[request.near_role].visits)-1))
                       and request.opposite_visit in run.members[1-request.near_role].visits]
            sites = [site for site in context.network.intersections if site.key == request.site.key]
            if (len(matches) != 1 or len(sites) != 1 or sites[0] != request.site or
                    request.local_layout is not context[footprint_key(request.site.key)]):
                raise PairedJunctionPreparationUnavailable('Explicit preparation lost its original source incidence or footprint')
            run, a = matches[0]
            ends = (End(near[0].section_id, int(near[0].forward)),
                    End(near[1].section_id, int(not near[1].forward)))
            if not all(end in request.site.ends for end in ends):
                raise PairedJunctionPreparationUnavailable('Explicit preparation lost its authored near-side incidence')
            tolerance = min(.00025, .0125*context.network.settings.curve_tolerance)
            points = _inner_curve(request.local_layout, ends, tolerance)
            opposite = _opposite(context, run.members[1-request.near_role].visits,
                                 inventory_frame(context, bundle, a), points, tolerance, request.near_role)
            if opposite != request.opposite_visit:
                raise PairedJunctionPreparationUnavailable('Explicit opposite visit has no certified original-source overlap')
    active = tuple(a for a in fresh.layout.approaches if hasattr(a, 'source_curve_stations'))
    if requests:
        sides, parameters, _ = canonical_reprepare(context, requests, additional_entrances=active)
        fresh = MergedEntranceSpecification.build(context, bundle, terminal.crossing, runs,
                                                  source_parameters=parameters)
        if fresh is None:
            raise PairedJunctionPreparationUnavailable('Canonical original paired terminal lost its source frame')
    else:
        sides = ()
    active = tuple(a for a in fresh.layout.approaches if hasattr(a, 'source_curve_stations'))
    all_entrances = (*active, *(a for side in sides for a in side.entrances))
    proof = verify_shared_source_phase(all_entrances)
    parameters = shared_source_parameters(all_entrances)
    fresh = replace(fresh, paired_sources=sides, source_phase_proof=proof)
    return fresh, sides, parameters, proof
