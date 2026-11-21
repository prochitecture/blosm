"""Failed-only normal charts for one admitted original junction source union.

The physical source union, complete profiles and ownership remain inputs. Only
the finite computational normals change; this does not infer a graph cycle or
consume a Section. Existing direct-quad and energy correction code is reused.
"""
from dataclasses import dataclass
import hashlib
import math
import numpy as np

from .model import End, Intersection, Key
from .routes import Route
from .physical_hole_cuts import PhysicalHoleCuts
from .source_identity import matches_current_section
from .source_height_chart_admission import admission_steps
from .updates import drain_steps


ORDINARY_PREFIX_CHART_POLICY = '3-single-junction-prefix-admission-after-ownership'


class OrdinarySourcePrefixUnavailable(ValueError):
    """This family is not a complete ordinary original-root source union."""


def _singleton_record(group):
    """Read the same site's correction before and after ownership finalization."""
    record = getattr(group.physical_layout, 'ordinary_complete_mouth_clearance', {})
    if not isinstance(record, dict) or len(group.members) != 1:
        return {}
    if 'corrections' not in record:
        return record
    members = [group.members[0].text()]
    rows = record['corrections']
    if not isinstance(rows, (tuple, list)) or len(rows) != 1 or record.get('original_members') != members:
        return {}
    row = rows[0]
    if not isinstance(row, dict) or row.get('original_members') != members or row.get('original_site') != members[0]:
        return {}
    return row


def eligible(group):
    """A narrow cheap gate; actual graph/profile/range admission is separate."""
    layout = group.physical_layout
    record = _singleton_record(group)
    return (len(group.members) == 1 and group.members[0].kind == 'intersection'
            and not group.claims and layout is not None and isinstance(record, dict)
            and record.get('correction') == 'FINITE_ORIGINAL_SINGLE_JUNCTION_SOURCE_UNION'
            and record.get('source_semantics') == 'one_original_junction_and_its_owned_prefixes'
            and bool(getattr(layout, 'source_cycle_full_width_proof', None)))


def _ranges(context, proof):
    if (not isinstance(proof, dict) or not proof.get('sections')
            or not math.isfinite(proof.get('maximum_error', float('inf')))
            or not 0. < proof['maximum_error'] <= .002):
        raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals need a complete finite source-width proof')
    rows = {}
    for count, row in enumerate(proof['sections'], 1):
        sid = row['section']; first, last = map(float, row['interval'])
        path = context[Key('path', sid)]
        if (count > 20000 or not matches_current_section(context,path.section)
                or not all(map(math.isfinite, (first, last)))
                or not 0. <= first < last <= path.length):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals changed a finite original source range')
        rows.setdefault(sid, []).append((first, last))
    result = []
    for sid, values in sorted(rows.items()):
        values = tuple(sorted(values))
        if any(a[1] > b[0] for a, b in zip(values, values[1:])):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals require disjoint ranges, never an unowned hull')
        result.append((sid, values))
    return tuple(result)


def _physical_token(group):
    layout = group.physical_layout
    digest = hashlib.sha256()
    for loop in group.loops:
        xy = np.asarray(loop.xy)
        digest.update(str((xy.dtype.str, xy.shape, loop.hole)).encode())
        digest.update(xy.tobytes())
    corners = tuple((c.first.source.id, c.first.end, c.second.source.id, c.second.end,
                     c.radius, c.requested_radius, c.radius_mode) for c in layout.corners)
    mouths = tuple(sorted((sid, end, mouth.station) for (sid, end), mouth in layout.mouths.items()))
    return digest.hexdigest(), corners, mouths, tuple(layout.outward_stations)


def _inventory(context, group):
    if not eligible(group):
        raise OrdinarySourcePrefixUnavailable('Ordinary prefix charts require an admitted original single-junction source union')
    layout = group.physical_layout; site = group.members[0]
    original = context.network.records().get(site)
    if type(original) is not Intersection:
        raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals require one actual original Junction')
    record = _singleton_record(group)
    if (record.get('original_site') != site.text()
            or record.get('original_members') != [site.text()]
            or record.get('consumed_sections') != []
            or record.get('original_paths_profiles_incidence_preserved') is not True
            or record.get('complete_source_grid') is not True
            or record.get('strict_authored_radii_preserved') is not True):
        raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals lost their original-root-union admission')
    approaches = tuple(layout.approaches)
    ends = tuple(End(a.source.id, a.end) for a in approaches)
    if (len(ends) != len(set(ends)) or set(ends) != set(original.ends)
            or len(approaches) != len(layout.outward_stations)
            or getattr(layout, 'partial_bundle_claims', ())
            or getattr(layout, 'finite_scope_ports', ())):
        raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals changed the original raw End inventory')
    sources = set()
    for approach in approaches:
        source = context.network.section_index[approach.source.id]
        if (approach.source is not source or not isinstance(approach.path, Route)
                or approach.profile.half or abs(approach.offset) > 1e-8
                or abs(approach.width - source.width) > 1e-8
                or approach.profile != context[Key('profile', source.id)].oriented(approach.end)
                or getattr(approach.path, 'absorbed_sites', ())
                or getattr(approach.path, 'mouth_piece', None) is not None
                or getattr(approach.path, 'deferred_profile_reserve', None) is not None):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals need centred complete original profiles')
        for piece in approach.path.pieces:
            sid = piece.section.id
            if (not matches_current_section(context,piece.section)
                    or piece.path is not context[Key('path', sid)]):
                raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals changed an original Path identity')
            sources.add(sid)
    if any(sources.intersection(bundle.sections) for bundle in context.network.bundles):
        raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals cannot infer lateral Bundle ownership')
    return site, tuple(sorted(ends)), sources


@dataclass(frozen=True)
class OrdinarySourcePrefixScope:
    network_uuid: str
    site: object
    ends: tuple
    ranges: tuple
    dependencies: tuple
    physical_token: tuple
    source_objects: tuple
    width_proof: object

    @classmethod
    def build(cls, context, group):
        from .cycle_dimensional_retry import certify_layout
        site, ends, sources = _inventory(context, group)
        proof = certify_layout(context, group, ())
        ranges = _ranges(context, proof)
        if ranges != _ranges(context, group.physical_layout.source_cycle_full_width_proof):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals rejected a stale owned-range proof')
        if not set(dict(ranges)) <= sources:
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix normals added an unowned original source')
        dependencies = tuple(sorted({site, *(Key(kind, sid) for sid in sources for kind in ('path', 'profile'))}))
        objects = tuple((key, context[key]) for key in dependencies if key.kind in ('path', 'profile'))
        return cls(context.network.uuid, site, ends, ranges, dependencies, _physical_token(group), objects, proof)

    def validate(self, context, group):
        site, ends, _ = _inventory(context, group)
        if (context.network.uuid != self.network_uuid or site != self.site or ends != self.ends
                or _physical_token(group) != self.physical_token
                or _ranges(context, group.physical_layout.source_cycle_full_width_proof) != self.ranges
                or any(context[key] is not value for key, value in self.source_objects)):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix scope is stale or belongs to another original graph')

    def candidates(self):
        ranked = sorted(self.ranges, key=lambda row: (-sum(last-first for first, last in row[1]), row[0]))[:3]
        axes = tuple((sid,) for sid, _ in ranked)
        return tuple(dict.fromkeys((*axes, tuple(sid for sid, _ in ranked))))


def _complete_profile_signature(context, group):
    """Cheap declared-input closure; construction independently reruns proof."""
    if not eligible(group):
        return None
    site, ends, sources = _inventory(context, group)
    ranges = _ranges(context, group.physical_layout.source_cycle_full_width_proof)
    dependencies = tuple(sorted({site, *(Key(kind, sid) for sid in sources for kind in ('path', 'profile'))}))
    return (ORDINARY_PREFIX_CHART_POLICY, context.network.uuid, site.text(), tuple((end.section_id, end.end) for end in ends),
            ranges, _physical_token(group), tuple((key.text(), context.digest(key)) for key in dependencies))


class OrdinarySourcePrefixCuts(PhysicalHoleCuts):
    def __init__(self, scope, axes, *, cadence_multiplier=1, **kwargs):
        if (type(scope) is not OrdinarySourcePrefixScope or not axes
                or tuple(axes) not in scope.candidates()):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix cuts need a bounded typed owned-axis candidate')
        super().__init__(source_ranges=dict(scope.ranges), cadence_multiplier=cadence_multiplier,
                         maximum_unpaired_span_widths=8., **kwargs)
        self.scope, self.axes = scope, tuple(axes)

    def _source_inventory(self, group, context):
        self.scope.validate(context, group)
        # All pavement intervals on a computational normal are included; this
        # flag grants no physical source ownership or longitudinal hull.
        return set(self.axes), set(self.axes)


def _complete_profile_steps(parent, group, context, previous):
    from .cycle_dimensional_retry import certify_layout
    from .partition import CompoundPartitioner
    from .source_cycle_atlas import SourceCycleEnergyChart
    from .four_side_charts import quality
    from .mesh import publication_coordinates
    scope = OrdinarySourcePrefixScope.build(context, group)

    class PrefixChart(SourceCycleEnergyChart):
        def _ordinary_source_prefix_steps(self, group, context, previous):
            raise OrdinarySourcePrefixUnavailable('Ordinary prefix charts cannot recursively retry themselves')

    failures = []
    for cadence in (1, 2):
        for axes in scope.candidates():
            policy = OrdinarySourcePrefixCuts(scope, axes, cadence_multiplier=cadence,
                spacing=parent.spacing, source_corridors=True, allow_filled=True,
                computational_extensions=True, existing_vertices=True)
            chart = PrefixChart(parent.patterns, spacing=parent.spacing, max_cells=parent.max_cells,
                quality_iterations=min(100, parent.quality_iterations), four_sided=True, recover_grass=True,
                partitioner=CompoundPartitioner(policy=policy, require_visibility=False))
            try:
                result = yield from chart.build_steps(group, context)
                if result.partition.group.members != group.members or result.partition.group.claims != group.claims:
                    raise ValueError('Ordinary prefix normals changed original physical ownership')
                if _physical_token(result.partition.group) != scope.physical_token:
                    raise ValueError('Ordinary prefix normals changed an exact physical boundary, mouth or radius')
                final = certify_layout(context, result.partition.group, ())
                if _ranges(context, final) != scope.ranges:
                    raise ValueError('Ordinary prefix normals changed original owned station intervals')
                meshes = (result.mesh, *(hole.mesh for hole in result.holes.values()))
                if sum(len(mesh.quads) for mesh in meshes) > parent.max_cells:
                    raise ValueError('Ordinary prefix atlas exceeds the unchanged compound cell budget')
                checks = []
                for mesh in meshes:
                    mesh.validate(); physical = quality(mesh)
                    published = mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))
                    published.validate(); stored = quality(published)
                    if (max(physical['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
                            or min(physical['minimum_corner_sine'], stored['minimum_corner_sine']) < .03):
                        raise ValueError('Ordinary prefix atlas fails unchanged physical or stored quad quality')
                    checks.append(dict(material=mesh.material, physical=physical, stored=stored))
                height_admission = None
                if parent.source_height_admission is not None:
                    height_admission = yield from admission_steps(parent.source_height_admission, result.partition.group, result)
            except ValueError as failure:
                failures.append(dict(axes=axes, cadence_multiplier=cadence, reason=str(failure)))
                continue
            result.metadata = dict(result.metadata, physical_hole_partition=True, source_corridors=True,
                source_width_visibility=final, ordinary_source_prefix_atlas=dict(
                    policy=ORDINARY_PREFIX_CHART_POLICY, original_site=scope.site.text(),
                    original_ends=[(end.section_id, end.end) for end in scope.ends],
                    source_ranges=dict(scope.ranges), selected_axes=axes, cadence_multiplier=cadence,
                    dependencies=[key.text() for key in scope.dependencies],
                    physical_boundary_changed=False, original_paths_profiles_preserved=True,
                    original_width_proof=scope.width_proof, final_width_proof=final,
                    maximum_axis_sets=4, maximum_attempts=8, actual_attempts=len(failures)+1,
                    maximum_bands_per_original_section=64, maximum_computational_cuts=policy.maximum_cuts,
                    unpaired_span_reports=policy.source_cut_reports, computational_extensions=policy.extension_reports,
                    stored_checks=checks, previous_chart_failure=previous, previous_prefix_failures=failures,
                    energy_recoveries=chart.energy_recoveries, cell_budget=parent.max_cells,
                    **(dict(source_height_admission=height_admission) if parent.source_height_admission is not None else {})))
            return result
    raise ValueError('No qualified finite ordinary-source-prefix atlas: '+str(failures))

from .model import fingerprint as _native_fingerprint
from .ordinary_mouth_clearance import _authenticate_native_half_profiles
from .profiles import memberships as _native_memberships
from .routes import route_key as _native_route_key
from .source_bank_reconciliation import owned_ranges as _native_owned_ranges

NATIVE_ORDINARY_PREFIX_CHART_POLICY = '2-clipped-original-circular-prefixes'


class NativePrefixAuthorityError(RuntimeError):
    """Actual source identity, ownership or width cannot be retried as a chart."""


def native_eligible(group):
    """Recognize claimed mixed unions before validating their authority.

    A missing width proof must fail admission, rather than making this claimed
    family disappear into the optional not-applicable path.
    """
    layout = group.physical_layout
    record = _singleton_record(group)
    return (len(group.members) == 1 and group.members[0].kind == 'intersection'
            and layout is not None and isinstance(record, dict)
            and record.get('correction') == 'FINITE_ORIGINAL_SINGLE_JUNCTION_SOURCE_UNION'
            and any(a.profile.half for a in layout.approaches))


def _native_path_token(path):
    # Canonical Path identity remains required. Capture its evaluated geometry
    # too, so mutation of an existing Path's arrays cannot reuse a scoped chart.
    arrays = []
    for name, value in sorted(vars(path).items()):
        if isinstance(value, np.ndarray):
            arrays.append((name, value.dtype.str, value.shape,
                           hashlib.sha256(value.tobytes()).hexdigest()))
    return (float(path.length), bool(path.exact), tuple(arrays))


def _native_prefix_inventory(context, group):
    if not native_eligible(group):
        raise NativePrefixAuthorityError('Native prefix charts require an admitted mixed native-profile root union')
    layout = group.physical_layout
    site = group.members[0]
    original = context.network.records().get(site)
    record = _singleton_record(group)
    ends = tuple(sorted(End(a.source.id, a.end) for a in layout.approaches))
    if (group.claims or type(original) is not Intersection or set(ends) != set(original.ends)
            or len(ends) != len(set(ends)) or len(ends) != len(layout.outward_stations)
            or record.get('original_site') != site.text()
            or record.get('original_members') != [site.text()]
            or record.get('source_semantics') != 'one_original_junction_and_its_owned_prefixes'
            or record.get('consumed_sections') != []
            or record.get('original_paths_profiles_incidence_preserved') is not True
            or record.get('complete_source_grid') is not True
            or record.get('strict_authored_radii_preserved') is not True
            or getattr(layout, 'partial_bundle_claims', ())
            or getattr(layout, 'finite_scope_ports', ())):
        raise NativePrefixAuthorityError('Native prefix charts changed actual original Ends or source-union admission')
    try:
        if not _authenticate_native_half_profiles(context, layout):
            raise ValueError('Prepared profiles are not canonical native inputs')
        ranges = _native_owned_ranges(context, group)
        if ranges != dict(_ranges(context, getattr(layout, 'source_cycle_full_width_proof', None))):
            raise ValueError('Actual finite source ranges differ from the admitted width proof')
        ownership = _native_memberships(context.network)
    except ValueError as failure:
        raise NativePrefixAuthorityError('Native prefix charts rejected changed native profile/range authority: '+str(failure)) from failure
    sources = {piece.section.id for arm in layout.approaches for piece in arm.path.pieces}
    if not set(ranges) <= sources:
        raise NativePrefixAuthorityError('Native prefix charts acquired an unowned original source')
    dependencies = {site, Key('settings', 0)}
    dependencies.update(_native_route_key(end) for end in ends)
    profiles = []
    objects = []
    tokens = []
    for sid in sorted(sources):
        section = context.network.section_index[sid]
        dependencies.add(section.key)
        for kind in ('path', 'profile'):
            key = Key(kind, sid)
            dependencies.add(key)
            objects.append((key, context[key]))
        member = ownership.get(sid)
        role = None
        if member is not None:
            bundle, lateral_role, forward = member
            dependencies.add(bundle.key)
            dependencies.update(Key('section', visit.section_id)
                                for path in bundle.members for visit in path.path)
            role = (bundle.key.text(), lateral_role, forward)
        profile = context[Key('profile', sid)]
        profiles.append((sid, role, profile.rise, profile.datum, profile.half))
        tokens.append((sid, _native_fingerprint(section), _native_path_token(context[Key('path', sid)]),
                       _native_fingerprint(profile)))
    objects.extend((_native_route_key(end), context[_native_route_key(end)]) for end in ends)
    # Only genuine unbundled FULL sources may select computational normals.
    # Native half roles remain height inputs and never become an axis mask.
    axes = tuple((sid,) for sid in sorted((sid for sid in ranges
        if sid not in ownership and context[Key('profile', sid)].half == 0),
        key=lambda sid: (-sum(b-a for a,b in ranges[sid]), sid))[:3])
    records = context.network.records()
    source_tokens = tuple((key.text(), _native_fingerprint(records[key]))
                          for key in sorted(dependencies) if key in records)
    return (site, ends, tuple(sorted(ranges.items())), tuple(sorted(dependencies)),
            tuple(objects), tuple(tokens), tuple(profiles), axes, source_tokens)


@dataclass(frozen=True)
class NativeOrdinarySourcePrefixScope:
    network_uuid: str
    site: object
    ends: tuple
    ranges: tuple
    dependencies: tuple
    source_objects: tuple
    source_tokens: tuple
    profile_roles: tuple
    axes: tuple
    input_records: tuple
    physical_token: tuple
    width_proof: object

    @classmethod
    def build(cls, context, group):
        from .cycle_dimensional_retry import certify_layout
        values = _native_prefix_inventory(context, group)
        try:
            proof = certify_layout(context, group, ())
            if _ranges(context, proof) != values[2]:
                raise ValueError('Fresh proof changed actual finite owned source ranges')
        except ValueError as failure:
            raise NativePrefixAuthorityError('Native prefix charts require fresh complete source width: '+str(failure)) from failure
        return cls(context.network.uuid, *values, _physical_token(group), proof)

    def validate(self, context, group):
        values = _native_prefix_inventory(context, group)
        previous = (self.site, self.ends, self.ranges, self.dependencies, self.source_objects,
                    self.source_tokens, self.profile_roles, self.axes, self.input_records)
        # Source objects must be canonical; numerical equality cannot alias a
        # different evaluated Path or Route into an already admitted scope.
        if (context.network.uuid != self.network_uuid
                or any(context[key] is not value for key, value in self.source_objects)
                or values[:4] != previous[:4] or values[5:] != previous[5:]
                or _physical_token(group) != self.physical_token):
            raise NativePrefixAuthorityError('Native prefix scope is stale or its actual source geometry changed')

    def candidates(self):
        return self.axes


class NativeOrdinarySourcePrefixCuts(PhysicalHoleCuts):
    def __init__(self, scope, axes, *, cadence_multiplier=1, **kwargs):
        if type(scope) is not NativeOrdinarySourcePrefixScope or tuple(axes) not in scope.candidates():
            raise NativePrefixAuthorityError('Native prefix cuts require an authenticated actual unbundled axis')
        self.scope, self.axes = scope, tuple(axes)
        super().__init__(source_ranges=dict(scope.ranges), cadence_multiplier=cadence_multiplier,
                         maximum_unpaired_span_widths=8., **kwargs)

    def _source_inventory(self, group, context):
        self.scope.validate(context, group)
        return set(self.axes), set()


def signature(context, group):
    if not native_eligible(group):
        return _complete_profile_signature(context, group)
    values = _native_prefix_inventory(context, group)
    site, ends, ranges, dependencies, _, tokens, profiles, axes, records = values
    return (NATIVE_ORDINARY_PREFIX_CHART_POLICY, context.network.uuid, site.text(),
            tuple((end.section_id, end.end) for end in ends), ranges, axes,
            _physical_token(group), tokens, profiles, records,
            tuple((key.text(), context.digest(key)) for key in dependencies))


def _native_prefix_steps(parent, group, context, previous):
    from .cycle_dimensional_retry import certify_layout
    from .partition import CompoundPartitioner, QuadChartGradingError
    from .source_cycle_atlas import SourceCycleEnergyChart
    from .quad_charts import (QuadChartCellBudgetError, QuadChartLayerBudgetError,
                              QuadChartCoreError, QuadChartQualityError, VisibilityKernelError)
    from .decomposition import VisibilityCutFamilyError
    from .physical_hole_cuts import PhysicalHoleCutError
    from .paired_interval_atlas import _checks
    scope = NativeOrdinarySourcePrefixScope.build(context, group)
    if not scope.candidates():
        return None
    failures = []
    computational = (QuadChartCellBudgetError, QuadChartLayerBudgetError, QuadChartCoreError,
                     QuadChartQualityError, QuadChartGradingError, VisibilityKernelError,
                     VisibilityCutFamilyError, PhysicalHoleCutError)

    class NativePrefixChart(SourceCycleEnergyChart):
        def _ordinary_source_prefix_steps(self, group, context, previous):
            raise NativePrefixAuthorityError('Native prefix charts cannot recursively retry themselves')

    for cadence in (1, 2):
        for axes in scope.candidates():
            policy = NativeOrdinarySourcePrefixCuts(scope, axes, cadence_multiplier=cadence,
                spacing=parent.spacing, source_corridors=True, hole_bridges=True, allow_filled=True,
                computational_extensions=True, existing_vertices=True)
            chart = NativePrefixChart(parent.patterns, spacing=parent.spacing, max_cells=parent.max_cells,
                quality_iterations=min(100, parent.quality_iterations), four_sided=True, recover_grass=True,
                partitioner=CompoundPartitioner(policy=policy, require_visibility=False))
            try:
                result = yield from chart.build_steps(group, context)
            except computational as failure:
                failures.append(dict(axes=axes, cadence_multiplier=cadence,
                                     failure_type=type(failure).__name__, reason=str(failure)))
                continue
            scope.validate(context, result.partition.group)
            if result.partition.group.members != group.members or result.partition.group.claims != group.claims:
                raise NativePrefixAuthorityError('Native prefix chart changed actual physical ownership')
            try:
                final = certify_layout(context, result.partition.group, ())
                if _ranges(context, final) != scope.ranges:
                    raise ValueError('Chart changed actual finite source ranges')
            except ValueError as failure:
                raise NativePrefixAuthorityError('Native prefix chart lost actual complete source width: '+str(failure)) from failure
            # These are the existing mandatory physical/stored mesh checks.
            # An untyped error is fatal, rather than a silently retried source
            # authority or admission failure. Only typed chart refusals retry.
            checks = _checks(parent, result)
            height = None
            if parent.source_height_admission is not None:
                try:
                    height = yield from admission_steps(parent.source_height_admission, result.partition.group, result)
                except ValueError as failure:
                    raise NativePrefixAuthorityError('Native prefix source-height admission refused: '+str(failure)) from failure
            result.metadata = dict(result.metadata, ordinary_source_prefix_atlas=dict(
                policy=NATIVE_ORDINARY_PREFIX_CHART_POLICY, original_site=scope.site.text(),
                original_ends=[(end.section_id, end.end) for end in scope.ends],
                source_ranges=dict(scope.ranges), selected_axes=axes, cadence_multiplier=cadence,
                dependencies=[key.text() for key in scope.dependencies], native_profile_roles=scope.profile_roles,
                actual_unbundled_axes_only=True, inferred_paired_mask=False, physical_claims_invented=False,
                physical_boundary_changed=False, original_paths_profiles_preserved=True,
                original_width_proof=scope.width_proof, final_width_proof=final,
                maximum_axis_sets=3, maximum_attempts=6, actual_attempts=len(failures)+1,
                maximum_bands_per_original_section=64, maximum_computational_cuts=policy.maximum_cuts,
                unpaired_span_reports=policy.source_cut_reports, computational_extensions=policy.extension_reports,
                stored_checks=checks, previous_chart_failure=previous, previous_prefix_failures=failures,
                energy_recoveries=chart.energy_recoveries, cell_budget=parent.max_cells,
                **(dict(source_height_admission=height) if height is not None else {})))
            return result
    raise QuadChartGradingError('No qualified authenticated native-prefix chart within existing bounds',
                               details=dict(candidates=failures, maximum_attempts=6))


def build(parent, group, context, previous):
    return drain_steps(build_steps(parent, group, context, previous))


def build_steps(parent, group, context, previous):
    if not native_eligible(group):
        return (yield from _complete_profile_steps(parent, group, context, previous))
    from .partition import QuadChartGradingError
    try:
        return (yield from _native_prefix_steps(parent, group, context, previous))
    except QuadChartGradingError as failure:
        attempts = (failure.details or {}).get('candidates', ())
        if not any(row['failure_type'] == 'QuadChartLayerBudgetError' for row in attempts):
            raise
        from .clipped_corner_sampling import resample
        candidate = resample(group, context.network.settings)
        if candidate is None:
            raise
        result = yield from _native_prefix_steps(parent, candidate, context, str(failure))
        if result is not None:
            result.metadata = dict(result.metadata, clipped_circular_fillets_resampled=True)
        return result
