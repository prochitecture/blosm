"""Retry a rejected ordinary junction chart along native divided-road approaches."""
from dataclasses import dataclass
import hashlib
from types import SimpleNamespace

import numpy as np

from .compounds import rim_label
from .cycle_dimensional_retry import certify_layout
from .junction_bank_transitions import original_corner
from .model import End, Intersection, Key, fingerprint
from .ordinary_mouth_clearance import _authenticate_native_half_profiles
from .ordinary_source_prefix_atlas import _native_path_token, _physical_token
from .physical_hole_cuts import PhysicalHoleCuts
from .source_bank_reconciliation import owned_ranges
from .source_cycle_atlas import SourceCycleEnergyChart


POLICY = '1-authenticated-ordinary-all-half-native-normal-geometry'


class HalfChartSourceError(RuntimeError):
    """The chart inputs no longer correspond to the original physical junction."""


def applicable(group):
    layout = group.physical_layout
    return (layout is not None and len(group.members) == 1 and not group.claims
            and len(group.footprints) == len(group.loops) == 1 and not group.loops[0].hole
            and all(arm.profile.half for arm in layout.approaches)
            and bool(getattr(layout, 'ordinary_complete_mouth_clearance', {}).get('corrections')))


@dataclass
class NativeHalfJunction:
    site: Key
    ends: tuple
    ranges: dict
    paths: tuple
    axes: tuple
    width_proof: dict

    @classmethod
    def read(cls, context, group):
        from .ownership import PLAN
        if not applicable(group):
            raise HalfChartSourceError('Expected one ordinary junction with original half crowns')
        site, = group.members
        source, = group.footprints
        layout = group.physical_layout
        plan = context.products.get(PLAN)
        owner = plan.by_site.get(site) if plan is not None else None
        if (owner is None or owner.members != group.members or owner.claims != group.claims
                or owner.physical_layout is not layout or len(owner.footprints) != 1
                or owner.footprints[0] is not source):
            raise HalfChartSourceError('Junction geometry no longer belongs to its current ownership plan')
        node = context.network.records().get(site)
        ends = tuple(End(a.source.id, a.end) for a in layout.approaches)
        if (type(node) is not Intersection or not 3 <= len(ends) <= 8 or len(set(ends)) != len(ends)
                or set(ends) != set(node.ends) or source.region.owner != site or source.region.boundary_only
                or any(getattr(source.layout, attr) is not getattr(layout, attr) for attr in
                       ('approaches', 'boundaries', 'corners', 'outward_stations', 'mouths', 'guide_weights'))
                or any(getattr(layout, attr, ()) for attr in ('finite_scope_ports', 'partial_bundle_claims',
                                                              'source_cycle_proofs'))):
            raise HalfChartSourceError('Expected the original ordinary approaches and corner geometry')
        try:
            native = _authenticate_native_half_profiles(context, layout)
        except ValueError as error:
            raise HalfChartSourceError('Native Bundle profile or orientation changed') from error
        record = layout.ordinary_complete_mouth_clearance
        if (not native or record.get('original_members') != [site.text()]
                or any(row.get('original_site') != site.text() or not row.get('complete_source_width')
                       or row.get('original_paths_profiles_incidence_preserved') is not True
                       or row.get('complete_source_grid') is not True for row in record['corrections'])):
            raise HalfChartSourceError('Expected a width-certified correction of the original mouths')
        scope = SimpleNamespace(context=context, group=group)
        if not all(original_corner(scope, node, corner) for corner in layout.corners):
            raise HalfChartSourceError('Corner radius or original approach changed')
        paths = tuple(sorted({piece.section.id for arm in layout.approaches for piece in arm.path.pieces}))
        ranges = owned_ranges(context, group)
        proof = certify_layout(context, group, ())
        if (not ranges or not ranges.keys() <= set(paths)
                or ranges != {row['section']: (tuple(row['interval']),) for row in proof['sections']}):
            raise HalfChartSourceError('Owned station ranges disagree with the complete-width proof')
        axes = sorted(ranges, key=lambda sid: (-sum(b-a for a, b in ranges[sid]), sid))[:3]
        return cls(site, ends, ranges, paths, tuple(axes), proof)

    def signature(self, context):
        from .profiles import memberships
        roles = memberships(context.network)
        keys = {self.site, Key('settings', 0)}
        for sid in self.paths:
            keys.update(Key(kind, sid) for kind in ('path', 'profile', 'section'))
            if sid in roles:
                bundle = roles[sid][0]
                keys.add(bundle.key)
                keys.update(Key('section', t.section_id) for m in bundle.members for t in m.path)
        return (POLICY, tuple((key.text(), context.digest(key)) for key in sorted(keys)),
                tuple(sorted(self.ranges.items())), self.axes)


def _source_geometry(group):
    areas = tuple((s.region.owner, s.region.xy.dtype.str, s.region.xy.shape,
                   hashlib.sha256(s.region.xy.tobytes()).digest(), s.segments) for s in group.footprints)
    corners = tuple((c.first.source.id, c.first.end, c.second.source.id, c.second.end, tuple(c.stations),
                     None if c.center is None else tuple(c.center), c.radius, c.requested_radius,
                     c.radius_mode, c.start_angle, c.sweep) for c in group.physical_layout.corners)
    return areas, corners


def _input_value(key, value):
    if key.kind == 'path':
        return fingerprint(value.section), _native_path_token(value)
    return fingerprint(value)


class NativeHalfChartScope:
    def __init__(self, context, group, source, axis):
        from .sampling import AnalyticalBoundarySampler
        if axis not in source.axes:
            raise HalfChartSourceError('The requested cut axis must be an owned approach')
        self.source, self.axis = source, axis
        self.uuid = context.network.uuid
        self.settings = fingerprint(context.network.settings)
        sampler = AnalyticalBoundarySampler(tolerance=context.network.settings.curve_tolerance,
                                            spacing=context.network.settings.spacing)
        self.rim = _physical_token(sampler.build(group))
        self.geometry = _source_geometry(group)
        keys = (Key(kind, sid) for sid in source.paths for kind in ('path', 'profile'))
        self.inputs = tuple((key, context[key], _input_value(key, context[key])) for key in keys)

    def validate(self, context, group):
        current = NativeHalfJunction.read(context, group)
        if (context.network.uuid != self.uuid or fingerprint(context.network.settings) != self.settings
                or current.site != self.source.site or current.ends != self.source.ends
                or current.ranges != self.source.ranges or self.axis not in current.axes
                or self.rim != _physical_token(group) or self.geometry != _source_geometry(group)
                or any(context[k] is not obj or token != _input_value(k, obj) for k, obj, token in self.inputs)):
            raise HalfChartSourceError('Original carriageway geometry changed during chart construction')
        return current.width_proof


class OrdinaryHalfNormalCuts(PhysicalHoleCuts):
    def __init__(self, scope, **kwargs):
        self.scope = scope
        super().__init__(source_ranges=dict(scope.source.ranges), source_corridors=True, hole_bridges=True,
                         computational_extensions=True, existing_vertices=True, maximum_unpaired_span_widths=8.,
                         **kwargs)

    def _source_inventory(self, group, context):
        self.scope.validate(context, group)
        return {self.scope.axis}, set()

    def cuts(self, group, context):
        self.scope.validate(context, group)
        sources = {s.region.owner: s for s in group.footprints}
        labels = [[rim_label(span, sources) for span in loop.spans] for loop in group.loops]
        points = group.loops[0].xy
        center = points.mean(axis=0)
        _, axes = np.linalg.eigh((points-center).T @ (points-center))
        self.extension_reports, self.endpoint_reuse_reports = [], []
        return self._cuts(group, context, (), center, axes[:, -1], labels)


class OrdinaryHalfNormalChart(SourceCycleEnergyChart):
    def __init__(self, *args, scope, context, **kwargs):
        self.scope, self.context = scope, context
        super().__init__(*args, **kwargs)

    def _repairable_leaf(self, group):
        self.scope.validate(self.context, group)
        return True


def recover(parent, group, context, previous):
    from .partition import CompoundPartitioner, QuadChartGradingError
    from .quad_charts import QuadChartCellBudgetError, QuadChartLayerBudgetError, QuadChartCoreError, VisibilityKernelError
    from .decomposition import VisibilityCutFamilyError
    from .physical_hole_cuts import PhysicalHoleCutError
    from .four_side_charts import quality
    from .mesh import publication_coordinates
    from .source_rim_normalization import _profiles

    source = NativeHalfJunction.read(context, group)
    profiles = _profiles(group)
    failures = []
    for cadence in (1, 2):
        for axis in source.axes:
            scope = NativeHalfChartScope(context, group, source, axis)
            cuts = OrdinaryHalfNormalCuts(scope, spacing=parent.spacing, cadence_multiplier=cadence)
            chart = OrdinaryHalfNormalChart(parent.patterns, scope=scope, context=context, spacing=parent.spacing,
                max_cells=parent.max_cells, quality_iterations=min(100, parent.quality_iterations), four_sided=True,
                source_height_admission=parent.source_height_admission,
                partitioner=CompoundPartitioner(cuts, require_visibility=False))
            try:
                result = chart.build(group, context)
            except (QuadChartGradingError, QuadChartCellBudgetError, QuadChartLayerBudgetError, QuadChartCoreError,
                    VisibilityKernelError, VisibilityCutFamilyError, PhysicalHoleCutError) as error:
                failures.append(dict(source=axis, cadence=cadence, error=str(error), error_type=type(error).__name__))
                continue
            proof = scope.validate(context, result.partition.group)
            after = _profiles(result.partition.group)
            if profiles.keys() != after.keys() or any(not np.array_equal(a, b) for k in profiles
                                                      for a, b in zip(profiles[k], after[k], strict=True)):
                raise HalfChartSourceError('Normal cuts changed an original profile')
            physical = quality(result.mesh)
            stored = quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz, 1., physical=True)))
            if (len(result.mesh.quads) > parent.max_cells
                    or max(physical['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
                    or min(physical['minimum_corner_sine'], stored['minimum_corner_sine']) < .03):
                raise HalfChartSourceError('Normal cuts exceed the physical or stored quad quality limits')
            result.metadata['ordinary_native_half_normal_geometry'] = dict(
                policy=POLICY, original_site=source.site.text(), actual_axis=axis, cadence=cadence,
                owned_native_ranges=source.ranges, paired_selector=[], actual_profiles_preserved=True,
                physical_curves_radii_mouths_claims_unchanged=True, before_source_width=source.width_proof,
                final_source_width=proof, source_cut_reports=cuts.source_cut_reports, previous_chart_failure=str(previous),
                previous_candidate_failures=failures, maximum_axes=3, maximum_cadences=2, maximum_bands=64,
                maximum_cuts=128, maximum_unpaired_span_widths=8., cell_budget=parent.max_cells)
            return result
    raise previous
