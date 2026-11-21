"""Original side incidences prepared from an already resolved source-cycle run.

This is source ownership preparation, before any junction Product exists. It
does not move an accepted mouth or create an opposite source End. Complete
original member inventories remain the later body's mandatory coverage scope.
"""
from dataclasses import dataclass
from types import SimpleNamespace

from .bundle_merged_entrances import PairedAngularEntrance
from .paired_junction_preparation import discover, PairedJunctionPreparationUnavailable
from .paired_source_sampling import canonical_reprepare, verify_shared_source_phase


@dataclass(frozen=True)
class PreparedSourceRun:
    bundle_id: int
    group: object
    layout: object
    paired_sources: tuple
    source_phase_proof: dict
    root_owner: object
    source_sections: tuple
    expected_run_count: int = 1
    consumed_sections: tuple = ()


def prepare(context, groups, *, prepared_bundle_ids=()):
    """Add only genuinely overlapping original side junctions to ownership.

    The resolved root is never rebuilt here. Canonical side preparation must
    agree exactly with its original source phases. Shared source changes that
    would require a different root remain unavailable rather than being snapped.
    """
    groups = tuple(groups)
    excluded = set(prepared_bundle_ids)
    candidates = {}
    for group in groups:
        layout = group.physical_layout
        if layout is None:
            continue
        for approach in layout.approaches:
            if not isinstance(approach, PairedAngularEntrance) or not hasattr(approach, 'run'):
                continue
            run = tuple(tuple(member.visits) for member in approach.run.members)
            for bundle in context.network.bundles:
                if (bundle.id in excluded or len(bundle.members) != 2
                        or not all(member.complete for member in bundle.members)
                        or run != tuple(tuple(member.path) for member in bundle.members)):
                    continue
                candidates.setdefault(bundle.id, []).append((bundle, group, approach))
    added, resources = [], {}
    occupied = {site for group in groups for site in group.members}
    for bundle_id, matches in sorted(candidates.items()):
        # Multiple whole-run roots need an explicit source partition. Do not
        # quietly choose one by distance or discard an original junction.
        if len(matches) != 1:
            continue
        bundle, group, approach = matches[0]
        terminal = SimpleNamespace(bundle_id=bundle_id,
                                   layout=SimpleNamespace(approaches=(approach,)))
        try:
            requests = discover(context, bundle, terminal)
        except PairedJunctionPreparationUnavailable:
            continue
        if not requests:
            continue
        sides, parameters, _ = canonical_reprepare(context, requests, additional_entrances=(approach,))
        entrances = (approach, *(a for side in sides for a in side.entrances))
        proof = verify_shared_source_phase(entrances)
        new_sites = [site for side in sides for site in side.group.members]
        if len(new_sites) != len(set(new_sites)) or occupied.intersection(new_sites):
            raise ValueError('Paired side preparation conflicts with an existing original physical owner')
        ids = tuple(sorted({visit.section_id for member in bundle.members for visit in member.path}
                           | {end.section_id for request in requests for end in request.site.ends}))
        # These dependencies include opposite sources that are absent from a
        # side junction's original local footprint. They invalidate that local
        # layout and surface when the common source phase or width changes.
        for side in sides:
            side.group.physical_layout.paired_source_dependencies = ids
        occupied.update(new_sites)
        added.extend(side.group for side in sides)
        resources[bundle_id] = PreparedSourceRun(bundle_id, group,
            terminal.layout, tuple(sides), proof, group.members[0], ids)
    return (*groups, *added), resources
