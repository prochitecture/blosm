"""Complete ordinary source widths before publishing physical mouth trims.

Mouth advancement retains original circular contacts, radii, Paths and profiles.
A separate failed-only preferred-radius fit records its changed corner. The
complete canonical grid includes future finite source stations for search.
"""
from dataclasses import replace
from types import SimpleNamespace
import math
import numpy as np

from .compounds import CompoundGroup, sample_footprint
from .finite_containment import segment_on_polygon_side
from .model import Key
from .updates import drain_steps
from .regions import RegionUnion
from .routes import Route, RouteTopology
from .source_strip_containment import SourceStripContainmentError

ORDINARY_MOUTH_POLICY = '15-finite-width-contacts-after-mouth-clearance'


class PreparedOwnerClearanceRequired(ValueError):
    """A real completed source owner needs an explicitly supported composer."""
    def __init__(self, sites, group, family):
        self.sites, self.owner_members, self.family = tuple(sites), tuple(group.members), family
        super().__init__('Ordinary complete-mouth closure requires prepared original owner '
            +group.members[0].text()+' ('+family+'): '+', '.join(key.text() for key in group.members))


def _family(group):
    from .finite_scope_spurs import FiniteScopeSpurPort
    from .finite_scope_compounds import ExpandedScopePort
    layout = group.physical_layout
    ports = getattr(layout, 'finite_scope_ports', ())
    if ports:
        if any(isinstance(p, ExpandedScopePort) for p in ports):
            return 'composed_finite_scope'
        return ('finite_scope_spur' if all(isinstance(p, FiniteScopeSpurPort) for p in ports)
                else 'finite_scope_port')
    if getattr(layout, 'finite_free_pads', ()):
        return 'finite_free_pad'
    if getattr(layout, 'one_sided_source_preparation', None):
        return 'one_sided_paired_source'
    if getattr(layout, 'terminal_merged_port', None):
        return 'merged_terminal_source'
    if getattr(layout, 'source_region_union', False) or hasattr(layout, 'source_cycle_proofs'):
        return 'certified_source_pool'
    if any(hasattr(a, 'source_mouths') for a in getattr(layout, 'approaches', ())):
        return 'prepared_source_mouths'
    return 'original_footprint_owner'


def _owner_mouths(context, group):
    """Only actual external original aliases plus whole consumed source bounds."""
    from .compounds import rim_label
    from .intersections import Mouth
    from .source_mouths import source_mouths
    sources = {source.region.owner: source for source in group.footprints}
    if group.physical_layout is not None:
        layout = group.physical_layout
        ports = tuple(zip(layout.approaches, layout.outward_stations))
    else:
        tokens = {(span.owner, label.index) for loop in group.loops for span in loop.spans
                  for label in (rim_label(span, sources),) if label.kind == 'mouth'}
        ports = tuple((sources[owner].layout.approaches[index],
                       sources[owner].layout.outward_stations[index]) for owner, index in sorted(tokens))
    mouths = {}
    for approach, station in ports:
        for alias in source_mouths(approach, station):
            if not alias.original_incidence:
                continue
            token = alias.endpoint.section_id, alias.endpoint.end
            if token in mouths and mouths[token].station != alias.station:
                raise ValueError('Prepared original owner has ambiguous external source mouths')
            mouths[token] = Mouth(alias.station, 0)
    for claim in group.claims:
        length = context[Key('path', claim.section_id)].length
        mouths[claim.section_id, 0] = Mouth(length, 0)
        mouths[claim.section_id, 1] = Mouth(0., 0)
    return SimpleNamespace(mouths=mouths, original_owner=group.members[0], original_members=group.members)


def _authenticate_native_half_profiles(context, layout):
    """Match each native approach to its Section, Route and Bundle-oriented crown."""
    from .model import End, Section
    from .intersections import Approach
    from .paths import Path
    from .profiles import Crown, ProfileProvider, memberships
    from .routes import route_key
    from .source_identity import matches_current_section
    ownership = memberships(context.network)
    sections = context.network.section_index
    for arm in layout.approaches:
        # Prepared/derived entrances are handled by their established composers.
        if type(arm) is not Approach or type(arm.path) is not Route or hasattr(arm, 'field'):
            return False
        if not matches_current_section(context, arm.source):
            raise ValueError('Native half clearance changed its original Section ownership')
        section = sections[arm.source.id]
        if (arm.width != section.width
                or arm.offset != section.offset*(1. if arm.end == 0 else -1.)):
            raise ValueError('Native half clearance changed its original source width/offset')
        if arm.path is not context.products.get(route_key(End(section.id, arm.end))):
            raise ValueError('Native half clearance changed its original Route ownership')
        pieces = arm.path.pieces
        if (not pieces or pieces[0].section.id != section.id
                or pieces[0].forward != (arm.end == 0)):
            raise ValueError('Native half clearance changed its original Route orientation')
        for piece in pieces:
            if (type(piece.section) is not Section or type(piece.path) is not Path
                    or piece.path is not context.products.get(Key('path', piece.section.id))
                    or not matches_current_section(context, piece.section)):
                raise ValueError('Native half clearance changed its original Path ownership')
            original = sections[piece.section.id]
            expected = ProfileProvider.build(original, ownership.get(original.id), sections)
            if context.products.get(Key('profile', original.id)) != expected:
                raise ValueError('Native half clearance needs the actual membership-resolved source profile')
        expected = ProfileProvider.build(section, ownership.get(section.id), sections).oriented(arm.end)
        if type(arm.profile) is not Crown or arm.profile != expected:
            raise ValueError('Native half clearance needs its actual oriented source profile')
    return True


def _eligible(context, layout):
    if (not getattr(layout, 'corners', ()) or getattr(layout, 'absorbed_transitions', ())
            or getattr(layout, 'finite_scope_ports', ()) or getattr(layout, 'finite_free_pads', ())
            or getattr(layout, 'transition', None) is not None
            or any(not isinstance(a.path, Route)
                   or hasattr(a, 'source_mouths') for a in layout.approaches)):
        return False
    for a in layout.approaches:
        if (getattr(a.path, 'absorbed_sites', ()) or getattr(a.path, 'mouth_piece', None) is not None
                or getattr(a.path, 'deferred_profile_reserve', None) is not None
                or getattr(a.path, 'contact_extension_through', ())):
            return False
        if a.profile != context[Key('profile', a.source.id)].oriented(a.end):
            raise ValueError('Ordinary mouth clearance changed its original source profile')
        for piece in a.path.pieces:
            if piece.path is not context[Key('path', piece.section.id)]:
                raise ValueError('Ordinary mouth clearance changed its original Path ownership')
    return (not any(arm.profile.half for arm in layout.approaches)
            or _authenticate_native_half_profiles(context, layout))


def _group(site, layout):
    count = sum(2*layout.radial + (len(boundary.initial_parameters)-1
                if getattr(boundary, 'initial_parameters', None) is not None else layout.rows)
                for boundary in layout.boundaries)
    if count > 20000:
        raise ValueError('Ordinary mouth clearance exceeds its source boundary budget')
    footprint = sample_footprint(site, layout)
    if len(footprint.region.xy) > 20000:
        raise ValueError('Ordinary mouth clearance exceeds its source boundary budget')
    loops = RegionUnion(max_segments=20000).build((footprint.region,))
    return CompoundGroup((site,), (), (footprint,), loops, physical_layout=layout)


def _canonical_candidate(context,site,layout):
    """Retain certified parameters on the same curves after a raw chord fails."""
    from copy import copy
    from .sampling import CertifiedCurbSampler,AnalyticalBoundarySampler
    from .cycle_dimensional_retry import certify_layout
    settings=context.network.settings
    sampler=CertifiedCurbSampler(settings.curve_tolerance,settings.spacing,20000)
    boundaries=[];count=2*layout.radial*len(layout.approaches)
    for boundary in layout.boundaries:
        parameters=sampler.sample(boundary,0.,1.)
        count+=len(parameters)-1
        if count>20000:
            raise ValueError('Ordinary mouth clearance exceeds its source boundary budget')
        parameters.flags.writeable=False
        prepared=copy(boundary)
        prepared.initial_parameters=parameters
        boundaries.append(prepared)
    prepared=replace(layout,boundaries=tuple(boundaries))
    group=_group(site,prepared)
    certify_layout(context,group,())
    # Chart preparation samples these same analytical curves again. Its exact
    # output must retain complete ownership too, before this layout is accepted.
    sampled=AnalyticalBoundarySampler(tolerance=settings.curve_tolerance,
        spacing=settings.spacing,max_samples=20000).build(group)
    proof=certify_layout(context,sampled,())
    return prepared,proof


def _record(layout,site,proof,records,tested,canonical=False,style_fit=None):
    layout.ordinary_complete_mouth_clearance=dict(policy=ORDINARY_MOUTH_POLICY,
        original_site=site.text(),records=records,candidates_tested=tested,
        complete_source_width=proof,original_corners_contacts_radii_preserved=style_fit is None,
        original_paths_profiles_incidence_preserved=True,complete_source_grid=True,
        certified_same_curve_sampling=canonical,style_source_width_fit=style_fit)
    return layout


def correct(planner, site, layout, *, maximum_passes=16, maximum_candidates=4096):
    return drain_steps(correct_steps(planner, site, layout, maximum_passes=maximum_passes,
                                     maximum_candidates=maximum_candidates))


def record_width_transition(layout, site, proof, mouth_records, tested, transition):
    accepted = _record(layout, site, proof, mouth_records, tested + transition['candidates_tested'], True)
    accepted.ordinary_complete_mouth_clearance.update(
        local_width_transition=transition, original_corners_contacts_radii_preserved=False,
        strict_authored_radii_preserved=True,
        original_mouth_planes_preserved=not mouth_records and transition['original_mouth_planes_preserved'])
    return accepted


def correct_steps(planner, site, layout, *, maximum_passes=16, maximum_candidates=4096):
    """Return the identical valid layout, or a freshly certified correction.

    Unrelated containment, curvature and allocation failures propagate. A true
    opposite-owner conflict retains the existing CompoundExpansionRequired type.
    """
    if (type(maximum_passes) is not int or not 1 <= maximum_passes <= 16
            or type(maximum_candidates) is not int or not 1 <= maximum_candidates <= 4096):
        raise ValueError('Ordinary mouth clearance needs bounded passes and candidates')
    context = planner.context
    if not _eligible(context, layout):
        return layout
    from .cycle_dimensional_retry import certify_layout
    try:
        certify_layout(context, _group(site, layout), ())
        return layout
    except SourceStripContainmentError as failure:
        original_failure = failure
    # Uniform layout rows can put a chord across an analytical contact. Keep
    # the existing certified sampling phase before considering a moved mouth.
    try:
        prepared,proof=_canonical_candidate(context,site,layout)
    except SourceStripContainmentError as failure:
        original_failure=failure
    else:
        return _record(prepared,site,proof,[],0,True)
    from .compound_cycles import _canonical_samples, _port_footprint
    from .finite_profile_search import advance
    from .intersections import BoundaryCurve, Mouth
    paths = {piece.section.id: piece.path for a in layout.approaches for piece in a.path.pieces}
    current, records, tested = layout, [], 0
    for _ in range(maximum_passes):
        yield
        cuts = {sid: [] for sid in paths}
        for a, station in zip(current.approaches, current.outward_stations):
            for (sid, end), value in a.path.trims(station).items():
                cuts[sid].append(value)
        for corner in current.corners:
            for a, station in zip((corner.first, corner.second), corner.stations):
                if station < 0.:
                    continue
                for (sid, end), value in a.path.trims(station).items():
                    cuts[sid].append(value)
        tolerance = min(.001, .05 * context.network.settings.curve_tolerance)
        partitions = _canonical_samples(paths, cuts, context.network.settings.spacing,
                                         tolerance, max_samples=20000)
        # Both directed ends of a returning original Section remain distinct
        # incident owners; a source ID alone would incorrectly skip its peer.
        ports = tuple(_port_footprint(Key('ordinary_owned_prefix', 2*a.source.id+a.end), a, station,
            current.radial, context.network.settings.spacing, tolerance, partitions)
            for a, station in zip(current.approaches, current.outward_stations))
        missing = []
        for i, (a, station) in enumerate(zip(current.approaches, current.outward_stations)):
            profile = a.point(station, np.array((a.offset-.5*a.width, a.offset+.5*a.width)))
            if any(not segment_on_polygon_side(*profile, port.region.xy,
                    inside=False, tolerance=1e-7) for j, port in enumerate(ports) if i != j):
                missing.append(i)
        if not missing:
            from .ordinary_source_width_fit import fit
            fitted=fit(context,site,current,original_failure,_canonical_candidate)
            if fitted is not None:
                prepared,proof,record=fitted
                return _record(prepared,site,proof,records,tested,True,record)
            # A cleared entrance and a finite width transition can both be
            # necessary. Reverting the entrance would clip the side street.
            if records:
                from .local_width_continuation import fit_steps
                for index, corner in enumerate(current.corners):
                    if corner.continuation and corner.radius_mode == 'STYLE':
                        fitted = yield from fit_steps(context, site, current, original_failure.section_id,
                            _canonical_candidate, corner_index=index)
                        if fitted is not None:
                            prepared, proof, transition = fitted
                            return record_width_transition(prepared, site, proof, records, tested, transition)
            raise original_failure
        if tested >= maximum_candidates:
            raise ValueError('Ordinary mouth clearance exceeds its finite candidate budget')
        outward, step = advance(planner, current.approaches, current.outward_stations,
            missing, ports, ports, partitions, (site,), max_candidates=maximum_candidates-tested)
        tested += max(row['candidates_tested'] for row in step)
        records.extend(step)
        boundaries = tuple(BoundaryCurve(corner, (outward[i], outward[(i+1) % len(outward)]),
            context.network.settings.curve_tolerance) for i, corner in enumerate(current.corners))
        rows = current.rows
        for boundary in boundaries:
            chord = context.network.settings.spacing
            corner = boundary.corner
            if corner.center is not None:
                radius = corner.radius;error = context.network.settings.curve_tolerance
                chord = min(chord, 2.*math.sqrt(max(1e-12, 2.*radius*error-error**2)))
            rows = max(rows, int(math.ceil(boundary.length/chord)))
            if boundary.continuation is not None:
                rows = max(rows, boundary.continuation.intervals(context.network.settings.spacing,
                            context.network.settings.curve_tolerance))
        rows = 4*int(math.ceil(rows/4.))
        mouths = {token: Mouth(float(value), 2*current.radial)
                  for a, station in zip(current.approaches, outward)
                  for token, value in a.path.trims(station).items()}
        current = replace(current, boundaries=boundaries, outward_stations=tuple(outward),
                          mouths=mouths, rows=rows)
        group = _group(site, current)
        canonical=False
        try:
            proof = certify_layout(context, group, ())
        except SourceStripContainmentError as failure:
            try:
                current,proof=_canonical_candidate(context,site,current)
                canonical=True
            except SourceStripContainmentError as failure:
                original_failure=failure
                continue
        return _record(current,site,proof,records,tested,canonical)
    raise ValueError('Ordinary complete-mouth clearance exhausts its finite pass budget')


def prepare(context, junctions, layouts, *, targets=None, prepared_groups=(),source_singletons=None):
    return drain_steps(prepare_steps(context, junctions, layouts, targets=targets,
                                     prepared_groups=prepared_groups, source_singletons=source_singletons))


def prepare_steps(context, junctions, layouts, *, targets=None, prepared_groups=(),source_singletons=None):
    """Effective layout map after all raw neighbors exist, before ownership."""
    from .compound_preparation import CompoundExpansionRequired, DeferredCompoundFootprint
    effective = dict(layouts)
    search_layouts = dict(layouts)
    original_owners = {key: (key,) for key in layouts}
    assigned = set()
    for group in prepared_groups:
        inventory = _owner_mouths(context, group)
        for member in group.members:
            if member in assigned:
                raise ValueError('Ordinary mouth inventory repeats a prepared original owner')
            assigned.add(member)
            if member not in search_layouts:
                raise ValueError('Prepared ordinary-mouth neighbor lost its original physical site')
            search_layouts[member] = inventory
            original_owners[member] = tuple(group.members)
    planner = SimpleNamespace(context=context, layouts=search_layouts,
                              topology=RouteTopology(context.network))
    targets = set(layouts) if targets is None else set(targets)
    if not targets <= set(layouts):
        raise ValueError('Ordinary mouth targets need their original physical footprints')
    corrected = {}
    for site in sorted(junctions, key=lambda value: value.key):
        if site.key not in targets:
            continue
        yield
        original = effective[site.key]
        try:
            try:
                result = yield from correct_steps(planner, site.key, original)
                if result is not original:
                    from .mixed_bank_continuation import fit as fit_mixed_bank
                    result = fit_mixed_bank(context, site.key, result)
            except SourceStripContainmentError as failure:
                if source_singletons is None:raise
                from .ordinary_singleton_sources import build_steps
                yield
                group=yield from build_steps(context,site.key,original,search_layouts,junctions,failure)
                source_singletons[site.key]=group
                result=group.physical_layout
        except CompoundExpansionRequired as expansion:
            neighbors = tuple(sorted({member for neighbor in expansion.neighbors
                                      for member in original_owners[neighbor]}))
            effective[site.key] = DeferredCompoundFootprint(site.key, tuple(original.approaches),
                neighbors, 'Ordinary complete mouth requires original opposite-owner closure')
            continue
        if result is not original:
            effective[site.key] = result
            search_layouts[site.key] = result
            corrected[site.key] = result
    return effective, corrected


def singleton_groups(corrected, groups):
    covered = {site for group in groups for site in group.members}
    return tuple(_group(site, layout) for site, layout in sorted(corrected.items()) if site not in covered)


def finalize(context, junctions, layouts, groups):
    return drain_steps(finalize_steps(context, junctions, layouts, groups))


def finalize_steps(context, junctions, layouts, groups):
    """Audit eventual standalone sites after all established ownership succeeds."""
    from .compound_preparation import DeferredCompoundPlanner, DeferredCompoundFootprint
    groups = tuple(groups)
    occupied = {site for group in groups for site in group.members}
    targets = set(layouts)-occupied
    source_singletons={}
    effective, corrected = yield from prepare_steps(context, junctions, layouts, targets=targets,
                                                   prepared_groups=groups, source_singletons=source_singletons)
    # Preserve the complete original source operands of these accepted pools;
    # their physical layout is a port inventory, not a replacement area polygon.
    groups=(*groups,*source_singletons.values())
    pending = {site for site in targets if isinstance(effective[site], DeferredCompoundFootprint)}
    if pending:
        by_site = {site: group for group in groups for site in group.members}
        touched = {by_site[neighbor].members: by_site[neighbor] for site in pending
                   for neighbor in effective[site].neighbors if neighbor in by_site}
        for group in touched.values():
            if _family(group) != 'finite_scope_spur':
                raise PreparedOwnerClearanceRequired(sorted(pending), group, _family(group))
        if touched:
            from .finite_scope_compounds import compose
            groups = compose(context, junctions, effective, groups)
            occupied = {site for group in groups for site in group.members}
            pending -= occupied
        if pending:
            protected = {site: group for group in groups for site in group.members}
            inventory = dict(effective)
            for group in groups:
                mouths = _owner_mouths(context, group)
                for site in group.members:
                    inventory[site] = mouths
            class OriginalClosure(DeferredCompoundPlanner):
                def exterior_steps(self, members, **kwargs):
                    for member in members:
                        if member in protected:
                            group = protected[member]
                            raise PreparedOwnerClearanceRequired(members, group, _family(group))
                    return (yield from super().exterior_steps(members, **kwargs))
            groups = yield from OriginalClosure(context, junctions, inventory).build_steps(groups)
    groups = (*groups, *singleton_groups(corrected, groups))
    groups = finish_groups(groups, corrected)
    members = [site for group in groups for site in group.members]
    if len(members) != len(set(members)):
        raise ValueError('Ordinary mouth closure repeats an original physical owner')
    return groups


def finish_groups(groups, corrected):
    """Preserve provenance when a corrected site joins a larger source pool."""
    from copy import copy
    result = []
    for group in groups:
        records = tuple(corrected[key].ordinary_complete_mouth_clearance
                        for key in group.members if key in corrected)
        if records and group.physical_layout is not None:
            layout = copy(group.physical_layout)
            layout.ordinary_complete_mouth_clearance = dict(policy=ORDINARY_MOUTH_POLICY,
                original_members=[key.text() for key in group.members if key in corrected],
                corrections=records)
            group = replace(group, physical_layout=layout)
        result.append(group)
    return tuple(result)


def signature(group):
    """Exact effective-layout inputs for projections omitting the PLAN digest."""
    records = []
    for source in group.footprints:
        record = getattr(source.layout, 'ordinary_complete_mouth_clearance', None)
        if record:
            records.append((source.region.owner.text(), record))
    physical = getattr(group.physical_layout, 'ordinary_complete_mouth_clearance', None)
    if physical:
        records.append(('physical ordinary mouths', physical))
    return (ORDINARY_MOUTH_POLICY, tuple(records)) if records else None
