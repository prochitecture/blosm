"""Original curb sides, rather than hole topology, determine boundary Z."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.source_curb_conditions import SourceCurbConditions, SourceCurbConditionUnavailable
from tested_python_streets.compounds import CompoundGroup, RimSpan, RimSegment, SampledFootprint
from tested_python_streets.corner_patch import _PreparedCurve
from tested_python_streets.profiles import Crown
from tested_python_streets.regions import Region
from tested_python_streets.source_rim_normalization import normalize_exterior, ResampledSourceCurb
from tested_python_streets.surface_fields import BoundaryJetField
from test_python_street_source_cycles import triangle, prepare
from test_python_street_conditions import inputs
from _python_street_fixtures import compound_cross_network
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.updates import BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator
from tested_python_streets.compound_dimensional_corners import prepare as prepare_dimension
from tested_python_streets.corner_radii import CornerRadiusPolicy
from test_python_street_dimensional_corners import fixture as dimensional_fixture


def paired_triangle():
    network = triangle(2.)
    bundle = model.Bundle(1, (
        model.BundleMember((model.Traversal(3, False),)),
        model.BundleMember((model.Traversal(1), model.Traversal(2)))))
    return replace(network, bundles=(bundle,))


def data(context, group, *, explicit=True):
    layouts = {source.region.owner: source.layout for source in group.footprints}
    kwargs = {'context': context} if explicit else {}
    return [BoundaryConditions().build(group, loop, layouts, **kwargs) for loop in group.loops]


class SourceCurbConditionTests(unittest.TestCase):
    def assert_continuous(self, boundaries):
        for boundary in boundaries:
            for i, rule in enumerate(boundary.rules):
                height, gradient = rule.sample(np.linspace(0., 1., 11))
                self.assertTrue(np.isfinite(gradient).all())
                self.assertAlmostEqual(float(np.asarray(height).flat[-1]),
                    float(boundary.rules[(i+1) % len(boundary.rules)].sample(0.)[0]), places=7)
        self.assertLessEqual(BoundaryJetField(boundaries).metadata['maximum_sampled_fit_error'], 1e-8)

    def test_explicit_full_source_curbs_keep_an_ordinary_closed_block_at_zero(self):
        for mode in (0, 3):
            context, group = prepare(triangle(2., mode))
            before = context.network.to_dict()
            xy = [loop.xy.copy() for loop in group.loops]
            legacy = data(context, group, explicit=False)
            actual = data(context, group)
            self.assertEqual(context.network.to_dict(), before)
            for loop, points in zip(group.loops, xy):
                np.testing.assert_array_equal(loop.xy, points)
            self.assertEqual(len(actual), 2)
            np.testing.assert_array_equal(legacy[1].heights, .08)
            np.testing.assert_array_equal(actual[1].heights, 0.)
            self.assertGreater(float(np.linalg.norm(actual[1].gradients, axis=1).max()), 0.)
            self.assert_continuous(actual)
            for old, new in zip(legacy, actual):
                self.assertEqual(old.ports, new.ports)
                for port in old.ports:
                    np.testing.assert_array_equal(new.heights[list(port.vertices)], old.heights[list(port.vertices)])
                    np.testing.assert_array_equal(new.gradients[list(port.vertices)], old.gradients[list(port.vertices)])

    def test_real_inner_half_profile_plateau_is_not_lowered_with_the_ordinary_block(self):
        context, group = prepare(paired_triangle())
        boundaries = data(context, group)
        hole = next(boundary for boundary in boundaries if boundary.hole)
        np.testing.assert_array_equal(hole.heights, .08)
        np.testing.assert_array_equal(hole.gradients, 0.)
        self.assert_continuous(boundaries)

    def test_equal_original_corner_values_and_jets_are_anchors_after_contact_clipping(self):
        context, group = prepare(triangle(2.))
        source = next(source for source in group.footprints if source.region.owner.kind == 'cycle_corner')
        index = next(i for i, curve in enumerate(source.layout.boundaries)
                     if any(label.kind == 'curb' and label.index == i for label in source.segments))
        corner = source.layout.corners[0]
        self.assertEqual(corner.first.profile.value(-1.), corner.second.profile.value(1.))
        resolver = SourceCurbConditions(context, group)
        rule = resolver.source_rule(source.region.owner, index, .25, .75)
        self.assertTrue(rule.source_anchor)
        height, gradient = rule.sample(np.linspace(0., 1., 13))
        np.testing.assert_array_equal(height, 0.)
        expected = 2.*abs(corner.first.profile.derivative(-1.))/corner.first.width
        np.testing.assert_allclose(np.linalg.norm(gradient, axis=1), expected, rtol=0., atol=1e-12)

    def test_inconsistent_adjacent_hard_sources_are_rejected_without_height_priority(self):
        for paired in (False, True):
            context, group = prepare(paired_triangle() if paired else triangle(2.))
            changed = dict(context.products)
            profile = changed[model.Key('profile', 1)]
            changed[model.Key('profile', 1)] = Crown(.12, profile.datum, profile.half)
            stale = BuildContext(context.network, changed, context.patterns)
            with self.assertRaisesRegex(ValueError, 'Adjacent original source anchors disagree|true native profile/frame'):
                data(stale, group)

    def test_ordinary_original_boundary_curves_use_their_explicit_source_lead_domains(self):
        network = compound_cross_network(model)
        group, layouts = inputs(network)
        products = UpdateCoordinator((PathProvider(), ProfileProvider()), PatternLibrary()).prepare(network).products
        context = BuildContext(network, products, PatternLibrary())
        old = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
        actual = [BoundaryConditions().build(group, loop, layouts, context=context) for loop in group.loops]
        for before, after in zip(old, actual):
            self.assertEqual(before.ports, after.ports)
            for port in before.ports:
                np.testing.assert_array_equal(after.heights[list(port.vertices)], before.heights[list(port.vertices)])
                np.testing.assert_array_equal(after.gradients[list(port.vertices)], before.gradients[list(port.vertices)])
        self.assert_continuous(actual)

    def test_normalized_exterior_resolves_only_its_explicit_original_span_ledger(self):
        context, group = prepare(triangle(2.))
        normalized = normalize_exterior(group)
        original_hole = group.loops[1]
        self.assertIs(normalized.loops[1], original_hole)
        for a, b in zip(normalized.area_footprints, group.area_footprints):
            self.assertIs(a, b)
        conditions = data(context, normalized)
        np.testing.assert_array_equal(conditions[1].heights, 0.)
        self.assert_continuous(conditions)
        carrier = normalized.footprints[0]
        curve = carrier.layout.boundaries[0]
        saved = curve.source_spans
        curve.source_spans = saved[:-1]
        try:
            with self.assertRaisesRegex(ValueError, 'point/span order'):
                data(context, normalized)
        finally:
            curve.source_spans = saved

    def test_hole_only_carrier_uses_constant_original_sources_without_a_fictitious_mouth(self):
        context, group = prepare(triangle(2.))
        source = next(source for source in group.footprints if source.region.owner.kind == 'cycle_corner')
        index = next(label.index for label in source.segments if label.kind == 'curb')
        curve = source.layout.boundaries[index]
        parameters = curve.stations/curve.length
        spans = tuple(RimSpan(source.region.owner, RimSegment('curb', index, float(a), float(b)))
                      for a, b in zip(parameters, parameters[1:]))
        retained = ResampledSourceCurb(curve.points, spans, ((source.region.owner, index, curve),))
        owner = model.Key('hole_condition_carrier', 1)
        carrier = SampledFootprint(replace(source.region, owner=owner, boundary_only=True), (),
            SimpleNamespace(approaches=(), corners=source.layout.corners, boundaries=(retained,)))
        prepared = replace(group, footprints=(carrier, *group.footprints))
        resolver = SourceCurbConditions(context, prepared)
        for first, last in zip(retained.initial_parameters, retained.initial_parameters[1:]):
            rule = resolver.source_rule(owner, 0, float(first), float(last))
            height, gradient = rule.sample(np.linspace(0., 1., 9))
            np.testing.assert_array_equal(height, 0.)
            self.assertTrue(np.isfinite(gradient).all())

    def test_default_call_stays_byte_exact_and_explicit_policy_rejects_missing_source(self):
        context, group = prepare(triangle(2.))
        first, second = data(context, group, explicit=False), data(context, group, explicit=False)
        for a, b in zip(first, second):
            np.testing.assert_array_equal(a.heights, b.heights)
            np.testing.assert_array_equal(a.gradients, b.gradients)
            self.assertEqual(a.rules, b.rules)
        resolver = SourceCurbConditions(context, group)
        source = next(source for source in group.footprints if source.region.owner.kind == 'cycle_strip')
        changed = replace(source, region=replace(source.region, owner=model.Key('unknown_source', 1)))
        unsupported = replace(group, footprints=(changed,))
        with self.assertRaises(SourceCurbConditionUnavailable):
            SourceCurbConditions(context, unsupported).source_rule(changed.region.owner, 0, 0., .1)

    def test_dimensional_core_and_original_leads_keep_exact_declared_piece_maps(self):
        for inner in (False, True):
            a, b = dimensional_fixture(reverse=True, routed=True, angle=.713, translation=(1200., -1800.))
            if inner:
                a.profile, b.profile = Crown(.08, 0., -1), Crown(.08, 0., 1)
            corner = CornerRadiusPolicy.solve(a, b, 2., 'STYLE')
            source = prepare_dimension(model.Key('dimensional_curb', 1), corner,
                model.Settings(curve_tolerance=.00025), original_inner_half_profiles=inner).source
            points = source.region.xy.copy()
            group = CompoundGroup((), (), (source,), ())
            resolver = SourceCurbConditions({}, group)
            anchored, core = 0, 0
            for label in source.segments:
                if label.kind != 'curb':
                    continue
                rule = resolver.source_rule(source.region.owner, label.index, label.first, label.last)
                height, gradient = rule.sample(np.linspace(0., 1., 9))
                np.testing.assert_array_equal(height, .08 if inner else 0.)
                self.assertTrue(np.isfinite(gradient).all())
                if getattr(rule, 'corner_fragment', False):
                    core += 1
                else:
                    anchored += 1
                    station = np.interp(np.linspace(rule.first, rule.last, 9)*rule.curve.length,
                        rule.curve.cumulative, rule.curve.parameters)
                    expected = np.asarray([rule.curve.approach.physical_gradient(float(s), rule.q)
                                           for s in station])
                    np.testing.assert_array_equal(gradient, expected)
            self.assertGreater(core, 0)
            self.assertGreater(anchored, 0)
            np.testing.assert_array_equal(source.region.xy, points)

    def test_reversed_source_storage_and_world_transform_preserve_real_curb_sides(self):
        angle = .713
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        for paired in (False, True):
            network = paired_triangle() if paired else triangle(2.)
            sections = []
            for section in network.sections:
                points = np.asarray(section.points).copy()
                points[:, :2] = points[:, :2]@rotation.T+(1200., -1800.)
                sections.append(replace(section, points=tuple(map(tuple, points[::-1]))))
            sites = tuple(replace(site, ends=tuple(model.End(e.section_id, 1-e.end) for e in site.ends))
                          for site in network.intersections)
            bundles = tuple(replace(bundle, members=tuple(replace(member,
                path=tuple(model.Traversal(v.section_id, not v.forward) for v in member.path))
                for member in bundle.members)) for bundle in network.bundles)
            context, group = prepare(replace(network, sections=tuple(sections), intersections=sites, bundles=bundles))
            boundaries = data(context, group)
            hole = next(boundary for boundary in boundaries if boundary.hole)
            np.testing.assert_array_equal(hole.heights, .08 if paired else 0.)
            self.assert_continuous(boundaries)

    def test_clipped_corner_intersection_uses_surviving_anchors_not_hidden_zero_contacts(self):
        # Two genuine unequal circular arcs intersect before their zero-curb
        # contacts. Their independent original height polynomials disagree at
        # that cut. The retained chain joins the two exposed plateau anchors.
        center0, center1 = np.array((-1., 0.)), np.array((1., 0.))
        r0, r1 = 2., 2.5
        x = -1.+(r0*r0-r1*r1+4.)/4.
        point = np.array((x, np.sqrt(r0*r0-(x+1.)**2)))
        sources, retained, spans = [], [], []
        for i, (center, radius, heights) in enumerate(((center0, r0, (.08, 0.)), (center1, r1, (0., .08)))):
            angle = np.arctan2(*(point-center)[::-1])
            angles = np.sort(np.r_[np.linspace(0., np.pi, 65), angle])[::-1]
            points = center+radius*np.c_[np.cos(angles), np.sin(angles)]
            at = int(np.flatnonzero(angles == angle)[0])
            points[at] = point
            curve = _PreparedCurve(points)
            first = SimpleNamespace(profile=Crown(.08, 0., -1 if heights[0] else 0), width=4.)
            second = SimpleNamespace(profile=Crown(.08, 0., 1 if heights[1] else 0), width=4.)
            for approach in (first, second):
                approach.physical_gradient = lambda station, q, profile=approach.profile: np.array((0., profile.derivative(q)/2.))
            corner = SimpleNamespace(first=first, second=second, center=center, radius=radius,
                                     start_angle=np.pi, sweep=-np.pi, stations=(0., 0.))
            owner = model.Key('cycle_corner', i)
            source = SampledFootprint(Region(owner, points[::-1]), (), SimpleNamespace(corners=(corner,), boundaries=(curve,)))
            sources.append(source)
            parameters = curve.stations/curve.length
            begin, end = (0, at) if i == 0 else (at, len(points)-1)
            retained.extend(points[begin:end])
            spans.extend(RimSpan(owner, RimSegment('curb', 0, float(a), float(b)))
                         for a, b in zip(parameters[begin:end], parameters[begin+1:end+1]))
        retained.append(sources[1].layout.boundaries[0].points[-1])
        group = CompoundGroup((), (), tuple(sources), ())
        resolver = SourceCurbConditions({}, group)
        last0 = next(i for i, span in enumerate(spans) if span.owner == sources[1].region.owner)-1
        a, b = spans[last0].label, spans[last0+1].label
        old_a = resolver.source_rule(sources[0].region.owner, 0, a.first, a.last).sample(1.)[0]
        old_b = resolver.source_rule(sources[1].region.owner, 0, b.first, b.last).sample(0.)[0]
        self.assertGreater(abs(float(old_a-old_b)), .005)
        anchor = SimpleNamespace(sample=lambda t: (.08, np.zeros(2)))
        before = [source.layout.boundaries[0].points.copy() for source in sources]
        with self.assertRaises(SourceCurbConditionUnavailable):
            resolver.chain_rules(spans, np.asarray(retained), False)
        rules = resolver.chain_rules(spans, np.asarray(retained), False, before=anchor, after=anchor)
        for rule in rules:
            height, gradient = rule.sample(np.linspace(0., 1., 7))
            np.testing.assert_array_equal(height, .08)
            np.testing.assert_array_equal(gradient, 0.)
        for source, points in zip(sources, before):
            np.testing.assert_array_equal(source.layout.boundaries[0].points, points)


if __name__ == '__main__':
    unittest.main()
