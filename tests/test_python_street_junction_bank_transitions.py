"""Finite original-bank height influence, with exact real profile planes."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.updates import UpdateCoordinator, BuildContext
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.source_support import ClippedApproachSupport
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets.compound_cycles import source_cycle_group
from tested_python_streets.conditions import BoundaryConditions, ProfileEdge
from tested_python_streets.source_curb_conditions import apply, SourceCurbConditions, SourceAnchorConflict
from tested_python_streets.junction_bank_transitions import (
    JunctionBankTransitionUnavailable, certify_body_domains, gradient_maximum)
from tested_python_streets.surface_fields import BoundaryJetField
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder


def fixture(*, reverse=False, angle=0., translation=(0., 0.), datum=0.):
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    controls = (((0, 0, 0), (-30, -4, 0)), ((0, 0, 0), (30, 8, 0)),
                ((0, 0, 0), (-10, 30, 0)), ((0, -8, 0), (30, 0, 0)))
    sections = []
    for sid, points in enumerate(controls, 1):
        xyz = np.array(points, float)
        xyz[:, :2] = xyz[:, :2]@rotation.T+translation
        if reverse:
            xyz = xyz[::-1]
        sections.append(model.Section(sid, sid, tuple(map(tuple, xyz)), width=5., interpolation=0, datum=datum))
    end = int(reverse)
    node = model.Intersection(1, tuple(model.End(sid, end) for sid in (1, 2, 3)),
                              radius=3., radius_modes=('STYLE',)*3)
    bundle = model.Bundle(1, (model.BundleMember((model.Traversal(2, not reverse),)),
                             model.BundleMember((model.Traversal(4, not reverse),))))
    network = model.Network('finite-inner-bank', tuple(sections), (node,), bundles=(bundle,))
    patterns = PatternLibrary()
    products = dict(UpdateCoordinator((PathProvider(), ProfileProvider()), patterns).prepare(network).products)
    context = BuildContext(network, products, patterns)
    ClippedApproachSupport._approaches(context, node)
    layouts = {node.key: IntersectionProvider.footprint(context, node)}
    planner = DeferredCompoundPlanner(context, network.intersections, layouts)
    group = source_cycle_group(planner, (node.key,), set(), set(), {}, ())
    baseline = BoundaryConditions().build(group, group.loops[0],
        {fp.region.owner: fp.layout for fp in group.footprints})
    return context, group, baseline


class JunctionBankTransitionTests(unittest.TestCase):
    def test_consumed_bank_uses_both_original_junctions_in_either_storage_direction(self):
        from test_python_street_compounds import footprint_context
        from tested_python_streets.compounds import RimSegment, RimSpan, SharedClaim
        from tested_python_streets.finite_free_caps import PreparedCapCurb
        from tested_python_streets.source_curb_coordinates import CurbCoordinates
        from tested_python_streets.source_curb_conditions import _NativeCoordinateRule
        from tested_python_streets.junction_bank_transitions import _bank_span
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                end = int(reverse)
                roads = (((0,0),(50,0),7.), ((0,0),(0,40),14.), ((0,0),(0,-40),3.),
                         ((50,0),(50,40),14.), ((50,0),(50,-40),3.), ((0,-10),(50,-10),7.))
                sections = tuple(model.Section(i,i,tuple((x,y,0.) for x,y in
                    ((b,a) if i==0 and reverse else (a,b))),width=width, interpolation=0)
                    for i,(a,b,width) in enumerate(roads))
                nodes = tuple(model.Intersection(i, (model.End(0,tip),model.End(n,0),model.End(n+1,0)),
                    radius=2., radius_modes=('STYLE',)*3) for i,tip,n in ((1,end,1),(2,1-end,3)))
                bundle = model.Bundle(1, (model.BundleMember((model.Traversal(0,not reverse),)),
                                          model.BundleMember((model.Traversal(5,True),))))
                ctx = footprint_context(model.Network('claimed-bank',sections,nodes,(bundle,)))
                sources = {}
                for node in nodes:
                    for corner in IntersectionProvider.footprint(ctx,node).corners:
                        key = model.Key('cycle_corner',len(sources))
                        sources[key] = SimpleNamespace(region=SimpleNamespace(owner=key),
                            layout=SimpleNamespace(original_junction=node.key,corners=(corner,)))
                path = ctx[model.Key('path',0)]
                side = -1 if reverse else 1
                sd = path.length-np.array((4.,5.,6.)) if reverse else np.array((4.,5.,6.))
                offsets = np.full(3, side*3.5)
                curve = PreparedCapCurb(path.point(sd,offsets))
                curve.source_coordinates = CurbCoordinates(path,sd,offsets,np.full(3,side))
                key = model.Key('cycle_strip',0)
                sources[key] = SimpleNamespace(region=SimpleNamespace(owner=key),
                                               layout=SimpleNamespace(boundaries=(curve,)))
                endpoints = tuple(n.key for n in (nodes[::-1] if reverse else nodes))
                claim = SharedClaim(0,*endpoints,path.length,0.)
                group = SimpleNamespace(claims=(claim,), members=endpoints)
                resolver = SimpleNamespace(group=group, context=ctx, sources=sources)
                label = RimSegment('curb',0,0.,1.)
                span = RimSpan(key,label)
                rule = _NativeCoordinateRule(ctx,curve,0.,1.)
                domain, at, knots = _bank_span(resolver,span,label,rule,nodes[0].key)
                np.testing.assert_array_equal(at,(4.,6.))
                self.assertEqual(knots,(True,True))
                self.assertEqual(domain.far_junction,nodes[1].key)
                self.assertIsNone(domain.mouth_station)
                group.claims = (replace(claim,first_trim=path.length-1.),)
                with self.assertRaises(JunctionBankTransitionUnavailable):
                    _bank_span(resolver,span,label,rule,nodes[0].key)

    def test_real_native_prefix_rotated_and_reversed_keeps_profiles_zero_edges_and_geometry(self):
        for reverse, angle in ((False, 0.), (True, .713)):
            context, group, baseline = fixture(reverse=reverse, angle=angle, translation=(1200., -1800.))
            original = context.network.to_dict()
            xy = group.loops[0].xy.copy()
            radii = [(c.radius, c.radius_mode, c.stations) for fp in group.footprints
                     for c in getattr(fp.layout, 'corners', ())]
            actual = apply(context, group, group.loops[0], baseline, maximum_slope=.1)
            domain, = actual.source_transitions
            self.assertLessEqual(domain['continuous_gradient_bound'], .1)
            self.assertLess(domain['outward_interval'][1], domain['opposite_contact'])
            self.assertLess(domain['opposite_contact'], domain['unchanged_mouth'])
            self.assertEqual(domain['source_end'], int(reverse))
            self.assertEqual(context.network.to_dict(), original)
            np.testing.assert_array_equal(group.loops[0].xy, xy)
            self.assertEqual(radii, [(c.radius, c.radius_mode, c.stations) for fp in group.footprints
                                     for c in getattr(fp.layout, 'corners', ())])
            self.assertLessEqual(BoundaryJetField((actual,)).metadata['maximum_sampled_fit_error'], 1e-8)
            for i, rule in enumerate(actual.rules):
                self.assertAlmostEqual(float(rule.sample(1.)[0]),
                                       float(actual.rules[(i+1) % len(actual.rules)].sample(0.)[0]), places=8)
                self.assertLessEqual(gradient_maximum(rule), .1+1e-12)
                if isinstance(rule, ProfileEdge):
                    np.testing.assert_array_equal(rule.sample(np.linspace(0., 1., 9))[0],
                                                  baseline.rules[i].sample(np.linspace(0., 1., 9))[0])
            for port in baseline.ports:
                np.testing.assert_array_equal(actual.heights[list(port.vertices)], baseline.heights[list(port.vertices)])
                np.testing.assert_array_equal(actual.gradients[list(port.vertices)], baseline.gradients[list(port.vertices)])
            resolver = SourceCurbConditions(context, group)
            from tested_python_streets.compounds import rim_label
            for i, span in enumerate(group.loops[0].spans):
                label = rim_label(span, resolver.sources)
                if label.kind != 'curb':
                    continue
                source = resolver.source_rule(span.owner, label.index, label.first, label.last)
                if not getattr(source, 'corner_fragment', False) and getattr(source, 'z0', None) == 0.:
                    np.testing.assert_array_equal(actual.rules[i].sample(np.linspace(0., 1., 9))[0], 0.)

    def test_old_success_is_byte_exact_and_higher_real_limit_chooses_shorter_native_domain(self):
        context, group, baseline = fixture()
        old = apply(context, group, group.loops[0], baseline)
        unchanged = apply(context, group, group.loops[0], baseline, maximum_slope=1.)
        np.testing.assert_array_equal(old.heights, unchanged.heights)
        np.testing.assert_array_equal(old.gradients, unchanged.gradients)
        for first, second in zip(old.rules, unchanged.rules):
            self.assertIs(type(first), type(second))
            self.assertEqual(vars(first).keys(), vars(second).keys())
            for key, value in vars(first).items():
                if isinstance(value, np.ndarray):
                    np.testing.assert_array_equal(value, vars(second)[key])
                else:
                    self.assertEqual(value, vars(second)[key])
        self.assertFalse(getattr(unchanged, 'source_transitions', ()))
        low = apply(context, group, group.loops[0], baseline, maximum_slope=.1).source_transitions[0]
        high = apply(context, group, group.loops[0], baseline, maximum_slope=.2).source_transitions[0]
        self.assertLess(high['outward_interval'][1], low['outward_interval'][1])

    def test_missing_coordinates_and_generated_opposite_corner_cannot_supply_authority(self):
        context, group, baseline = fixture()
        for remove_coordinates in (True, False):
            footprints = []
            for fp in group.footprints:
                layout = copy(fp.layout)
                if remove_coordinates and fp.region.owner == model.Key('cycle_port', 4):
                    curves = list(layout.boundaries)
                    curves[0] = copy(curves[0]); curves[0].source_coordinates = None
                    layout.boundaries = tuple(curves)
                elif not remove_coordinates and fp.region.owner == model.Key('cycle_corner', 0):
                    layout.source_bank_fillet = SimpleNamespace(kind='generated')
                footprints.append(replace(fp, layout=layout))
            changed = replace(group, footprints=tuple(footprints))
            with self.assertRaises(JunctionBankTransitionUnavailable):
                apply(context, changed, changed.loops[0], baseline, maximum_slope=.1)

    def test_forged_source_coordinates_stay_a_fatal_source_error(self):
        context, group, baseline = fixture()
        footprints = []
        from tested_python_streets.source_curb_coordinates import CurbCoordinates
        for fp in group.footprints:
            if fp.region.owner == model.Key('cycle_port', 4):
                layout = copy(fp.layout); curves = list(layout.boundaries); curves[0] = copy(curves[0])
                old = curves[0].source_coordinates
                stations = old.stations.copy(); stations[1] += .001
                curves[0].source_coordinates = CurbCoordinates(old.path, stations, old.offsets, old.sides,
                                                              approach=old.approach)
                layout.boundaries = tuple(curves); fp = replace(fp, layout=layout)
            footprints.append(fp)
        with self.assertRaisesRegex(ValueError, 'physical vertices'):
            apply(context, replace(group, footprints=tuple(footprints)), group.loops[0], baseline, maximum_slope=.1)

    def test_fixed_crossfall_invalid_limits_and_real_body_overlap_reject(self):
        context, group, baseline = fixture()
        for value in (0., -1., np.nan, np.inf, True):
            with self.assertRaises(ValueError):
                apply(context, group, group.loops[0], baseline, maximum_slope=value)
        with self.assertRaisesRegex(JunctionBankTransitionUnavailable, 'crossfall'):
            apply(context, group, group.loops[0], baseline, maximum_slope=.01)
        actual = apply(context, group, group.loops[0], baseline, maximum_slope=.1)
        row = actual.source_transitions[0]
        guides = SimpleNamespace(ledger={'complete_supported_scope': True,
            'source_influence': [dict(section=row['section'], owned=row['native_interval'], core_ranges=[])]})
        self.assertEqual(certify_body_domains(guides, (actual,)), actual.source_transitions)
        guides.ledger['source_influence'][0]['core_ranges'] = [row['native_interval']]
        with self.assertRaisesRegex(JunctionBankTransitionUnavailable, 'hard original body'):
            certify_body_domains(guides, (actual,))
        guides.ledger['source_influence'] = []
        with self.assertRaisesRegex(JunctionBankTransitionUnavailable, 'complete original source-body'):
            certify_body_domains(guides, (actual,))

    def test_builder_passes_configured_limit_from_immutable_real_profile_planes(self):
        context, group, baseline = fixture()
        incoming = max(float(np.linalg.norm(baseline.gradients[list(port.vertices)], axis=1).max())
                       for port in baseline.ports)
        class Captured(Exception): pass
        seen = []
        def capture(c, g, loop, data, *, maximum_slope=None):
            self.assertIs(c, context); self.assertIs(g, group)
            seen.append(maximum_slope)
            raise Captured()
        chart = SimpleNamespace(partition=SimpleNamespace(group=group))
        def chart_steps(*args, **kwargs):
            yield
            return chart
        with patch('tested_python_streets.compound_surfaces.CompoundQuadChart.build_steps', side_effect=chart_steps), \
             patch('tested_python_streets.source_curb_conditions.apply', side_effect=capture):
            with self.assertRaises(Captured):
                CompoundSurfaceBuilder(context.patterns, source_heights=True,
                    maximum_slope=.37, slope_margin=7.).build(group, context)
        self.assertEqual(seen, [max(.37, 7.*incoming)])

    def test_unrelated_source_failure_is_never_treated_as_a_transition_candidate(self):
        context, group, baseline = fixture()
        with patch.object(SourceCurbConditions, '_chain_rules', side_effect=ValueError('unrelated source failure')):
            with self.assertRaisesRegex(ValueError, '^unrelated source failure$'):
                apply(context, group, group.loops[0], baseline, maximum_slope=.1)

    def test_scope_and_full_profile_bank_cannot_be_relabelled_as_inner_transition(self):
        context, group, baseline = fixture()
        endpoint = model.Key('profile', 2)
        old = context.products[endpoint]
        context.products[endpoint] = replace(old, half=0)
        try:
            with self.assertRaises((JunctionBankTransitionUnavailable, SourceAnchorConflict)):
                apply(context, group, group.loops[0], baseline, maximum_slope=.1)
        finally:
            context.products[endpoint] = old

        section = context.network.section_index[2]
        object.__setattr__(section, 'scope_ends', (0,))
        try:
            with self.assertRaises(JunctionBankTransitionUnavailable):
                apply(context, group, group.loops[0], baseline, maximum_slope=.1)
        finally:
            object.__setattr__(section, 'scope_ends', ())

    def test_available_original_bank_is_finite_and_common_nonzero_datum_is_preserved(self):
        context, group, baseline = fixture(datum=.7)
        with self.assertRaisesRegex(JunctionBankTransitionUnavailable, 'finite inner-bank knot'):
            apply(context, group, group.loops[0], baseline, maximum_slope=.064001)
        actual = apply(context, group, group.loops[0], baseline, maximum_slope=.1)
        self.assertEqual(actual.datum, .7)
        self.assertEqual(actual.ceiling, .7+.08)
        self.assertGreaterEqual(float(actual.heights.min()), .7)
        self.assertLessEqual(float(actual.heights.max()), .7+.08)

    def test_original_strict_override_cannot_reuse_stale_style_geometry(self):
        context, group, baseline = fixture()
        site = context.network.intersections[0]
        context.network = replace(context.network, intersections=(replace(site, radius_modes=('STRICT',)*3),))
        with self.assertRaisesRegex(JunctionBankTransitionUnavailable, 'actual original corner'):
            apply(context, group, group.loops[0], baseline, maximum_slope=.1)


if __name__ == '__main__': unittest.main()
