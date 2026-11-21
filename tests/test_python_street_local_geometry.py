"""Local transition and clipped-circle regressions in physical metres."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from test_python_street_continuations import through_junction
from tested_python_streets.clipped_corner_sampling import ClippedCircularCurb
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.local_width_continuation import LocalTransition, TransitionHeight, fit, fit_with_wider_approach_lead
from tested_python_streets.model import Connection, End, Key
from tested_python_streets.profiles import Crown
from tested_python_streets.sampling import CertifiedCurbSampler


class LocalGeometryTests(unittest.TestCase):
    def test_width_transition_uses_the_entrance_already_cleared_for_a_side_street(self):
        from types import SimpleNamespace
        from tested_python_streets import model
        from tested_python_streets.ordinary_mouth_clearance import correct
        from tested_python_streets.routes import RouteTopology

        # Original arm directions and widths of the second Bergamo Bundle fork.
        root = (-850.3736572265625, -824.6301879882812)
        controls = (((-859.437744140625,-828.6809692382812),root),
                    (root,(-840.5625,-819.7556762695312)),
                    ((-846.4586181640625,-869.1251220703125),(-848.5791625976562,-844.4564208984375),
                     (-850.1173706054688,-828.0365600585938),root),
                    ((-865.,-818.),(-845.,-808.)))
        sections = tuple(model.Section(i, street, tuple((x,y,0.) for x,y in points),
            width=width, interpolation=3) for i,(points,street,width) in enumerate(
                zip(controls,(10,10,11,12),(3.5,7.,2.,3.5))))
        site = model.Intersection(4,(End(0,1),End(1,0),End(2,1)),2.,(2.,2.,1.),('STYLE',)*3)
        bundle = model.Bundle(8,(model.BundleMember((model.Traversal(3,True),)),
            model.BundleMember((model.Traversal(0,True),model.Traversal(1,True)))))
        network = model.Network('cleared-width-transition',sections,(site,),(bundle,),
                                settings=model.Settings(cell_size=.7,exaggeration=1.))
        ctx = footprint_context(network)
        layout = IntersectionProvider.footprint(ctx,site)
        planner = SimpleNamespace(context=ctx,layouts={site.key:layout},topology=RouteTopology(network))
        result = correct(planner,site.key,layout)
        record = result.ordinary_complete_mouth_clearance
        advance, = record['records']
        self.assertEqual((advance['section'],advance['end']),(1,0))
        self.assertGreater(advance['station'], advance['previous_station'])
        self.assertFalse(record['original_mouth_planes_preserved'])
        self.assertEqual(result.mouths[1,0].station, advance['station'])
        transition = record['local_width_transition']
        self.assertLessEqual(transition['curvature_bound'], 1./transition['radius'])
        self.assertEqual({row['section'] for row in record['complete_source_width']['sections']}, {0,1,2})
        for i, boundary in enumerate(result.boundaries):
            if i != transition['corner']:
                self.assertIs(boundary.corner, layout.corners[i])

    def test_fork_keeps_the_wide_bank_when_a_circle_touches_behind_the_narrow_arm(self):
        from tested_python_streets.local_width_continuation import fit_steps
        from tested_python_streets.model import Intersection, Network, Section
        from tested_python_streets.ordinary_mouth_clearance import _canonical_candidate
        from tested_python_streets.source_strip_containment import SourceStripContainmentError
        from tested_python_streets.updates import drain_steps

        arms = (
            ((0.,0.), (3.38,.78), (7.32,1.07), (10.59,.67), (13.16,.17), (15.50,-1.32),
             (17.79,-3.82), (19.33,-7.22), (22.08,-18.92)),
            ((0.,0.), (2.41,2.92), (6.12,5.42), (10.29,7.18), (19.06,9.75), (34.03,12.92), (47.05,13.05)),
            ((0.,0.), (-2.13,-3.33), (-3.03,-6.59), (-3.11,-9.86), (-.65,-40.62)),
        )
        sections = tuple(Section(i, i, tuple((x,y,0.) for x,y in arm), width=width, interpolation=3)
                         for i,(arm,width) in enumerate(zip(arms,(3.5,3.5,7.))))
        site = Intersection(7, tuple(End(s.id,0) for s in sections), 2., (2.,)*3, ('STYLE',)*3)
        ctx = footprint_context(Network('unequal-fork', sections, (site,)))
        layout = IntersectionProvider.footprint(ctx,site)
        with self.assertRaises(SourceStripContainmentError) as failure:
            _canonical_candidate(ctx,site.key,layout)
        wide = sections[-1].id
        self.assertEqual(failure.exception.section_id, wide)
        candidate, _, evidence = drain_steps(fit_steps(ctx, site.key, layout, wide, _canonical_candidate,
                                                       extrapolated=True))
        i = evidence['corner']
        self.assertLess(min(layout.corners[i].stations), 0.)
        self.assertGreaterEqual(min(candidate.corners[i].stations), 0.)
        self.assertEqual(layout.mouths, candidate.mouths)
        self.assertLessEqual(evidence['curvature_bound'], 1./layout.corners[i].radius)
        for j in range(len(layout.corners)):
            if i != j:
                self.assertIs(layout.corners[j], candidate.corners[j])
        strict = tuple(replace(c, radius_mode='STRICT') for c in layout.corners)
        self.assertIsNone(drain_steps(fit_steps(ctx, site.key, replace(layout, corners=strict),
                                                wide, _canonical_candidate, extrapolated=True)))
        # A Bundle crest reaching an outer bank needs both original endpoint
        # heights and gradients. Opposing crests alone establish no passage.
        original = layout.corners[i]
        for arm in (original.first, original.second):
            arm.profile = Crown(.08, 0., -1)
        mixed, _, _ = drain_steps(fit_steps(ctx, site.key, layout, wide, _canonical_candidate,
                                             corner_index=i, extrapolated=True))
        core, = [piece.curve for piece in mixed.boundaries[i].pieces if isinstance(piece, LocalTransition)]
        height, gradient = TransitionHeight(core,0.,1.).sample(np.array((0.,1.)))
        np.testing.assert_allclose(height,(.08,0.),atol=1e-12)
        for at,(arm,side) in enumerate(((original.first,-1),(original.second,1))):
            np.testing.assert_allclose(gradient[at],arm.physical_gradient(core.mouths[at],side),atol=1e-7)
        original.second.profile = Crown(.08, 0., 1)
        self.assertIsNone(drain_steps(fit_steps(ctx,site.key,layout,wide,_canonical_candidate,
                                                corner_index=i,extrapolated=True)))

        # Both banks can also be the same crest along one declared member.
        from tested_python_streets.model import Bundle, BundleMember, Traversal
        first, second = original.first.source.id, original.second.source.id
        members = (BundleMember((Traversal(first, False), Traversal(second, True))),
                   BundleMember((Traversal(99, True),)))
        roads = tuple(replace(s, street_id=37) if s.id in (first, second) else s for s in sections)
        peer = Section(99, 99, ((50., 50., 0.), (80., 50., 0.)), width=4.)
        ctx = footprint_context(Network('same-crest', roads + (peer,), (site,), (Bundle(8, members),)))
        layout = IntersectionProvider.footprint(ctx, site)
        candidate, _, evidence = drain_steps(fit_steps(ctx,site.key,layout,wide,_canonical_candidate,
                                                       corner_index=i,extrapolated=True))
        core, = [p.curve for p in candidate.boundaries[i].pieces if isinstance(p,LocalTransition)]
        height, gradient = TransitionHeight(core,0.,1.).sample(np.linspace(0.,1.,101))
        np.testing.assert_array_equal(height, np.full(101, core.corner.first.profile.rise))
        np.testing.assert_array_equal(gradient, np.zeros((101,2)))
        self.assertEqual(candidate.mouths,layout.mouths)
        self.assertLessEqual(evidence['curvature_bound'],1./core.corner.radius)

    def test_width_fit_moves_only_the_wider_unopposed_mouth(self):
        network = through_junction((10., 4., 6.))
        site, = network.intersections
        site = replace(site, radii=(site.radius,)*3, radius_modes=('STYLE',)*3)
        ctx = footprint_context(replace(network, intersections=(site,)))
        layout = IntersectionProvider.footprint(ctx, site)
        i, = [i for i,c in enumerate(layout.corners) if c.continuation]
        c = layout.corners[i]
        arms = c.first, c.second
        wide = int(arms[1].width > arms[0].width)
        slot = (i+wide) % 3
        stations = list(layout.outward_stations)
        stations[slot], stations[(i+1-wide)%3] = .5, 28.
        layout = replace(layout, outward_stations=tuple(stations))
        sid = arms[1-wide].source.id
        result = fit_with_wider_approach_lead(ctx, site.key, layout, sid, lambda _, __, candidate: (candidate, {}))
        candidate, _, record = result
        self.assertGreater(candidate.outward_stations[slot], stations[slot])
        self.assertEqual([a for j,a in enumerate(candidate.outward_stations) if j != slot],
                         [a for j,a in enumerate(stations) if j != slot])
        self.assertFalse(record['original_mouth_planes_preserved'])
        self.assertEqual(record['advanced_section'], arms[wide].source.id)
        strict = list(layout.corners)
        strict[i] = replace(c, radius_mode='STRICT')
        self.assertIsNone(fit_with_wider_approach_lead(ctx, site.key, replace(layout, corners=tuple(strict)),
                                                       sid, None))
        other_end = End(arms[wide].source.id, 1-arms[wide].end)
        neighbor = Connection(2, (other_end, End(sid, 1-arms[1-wide].end)))
        connected = footprint_context(replace(ctx.network, connections=(neighbor,)))
        self.assertIsNone(fit_with_wider_approach_lead(connected, site.key, layout, sid, None))

    def test_local_transition_accepts_wide_bank_on_either_side(self):
        for widths in ((10., 4., 6.), (6., 4., 10.)):
            with self.subTest(widths=widths):
                network = through_junction(widths)
                ctx = footprint_context(network)
                site, = network.intersections
                layout = IntersectionProvider.footprint(ctx, site)
                i, = [i for i,c in enumerate(layout.corners) if c.continuation]
                corner = layout.corners[i]
                arms = corner.first, corner.second
                wide = int(arms[1].width > arms[0].width)
                arms[wide].profile = Crown(.08, 0., (-1, 1)[wide])
                stations = list(layout.outward_stations)
                stations[i] = .5 if wide == 0 else 28.
                stations[(i+1) % len(stations)] = .5 if wide == 1 else 28.
                layout = replace(layout, outward_stations=tuple(stations))
                candidate, _, report = fit(ctx, site.key, layout, arms[wide].source.id,
                                           lambda ctx, site, layout: (layout, None))
                self.assertEqual(candidate.outward_stations, tuple(stations))
                self.assertEqual(report['contacts'][wide], 0.)
                local, = [p for p in candidate.boundaries[i].pieces if isinstance(p, LocalTransition)]
                rule = TransitionHeight(local.curve, 0., 1.)
                heights, gradients = rule.sample(np.array([0., 1.]))
                self.assertAlmostEqual(heights[wide], .08)
                self.assertAlmostEqual(heights[1-wide], 0.)
                for k, arm in enumerate(arms):
                    expected = arm.physical_gradient(report['contacts'][k], (-1.,1.)[k])
                    np.testing.assert_allclose(gradients[k], expected, rtol=0., atol=1e-7)
                u = np.linspace(.02, .98, 51)
                numerical = (rule.sample(u+1e-5)[0] - rule.sample(u-1e-5)[0])/2e-5
                analytical = np.sum(rule.sample(u)[1]*local.curve.continuation.derivative(u), axis=1)
                np.testing.assert_allclose(numerical, analytical, rtol=0., atol=1e-8)

    def test_clipped_circle_keeps_original_angular_accuracy_and_avoids_last_sliver(self):
        network = through_junction()
        site, = network.intersections
        layout = IntersectionProvider.footprint(footprint_context(network), site)
        corner = next(c for c in layout.corners if c.center is not None)
        owner = Key('cycle_corner', 1)
        patch = corner_patch(owner, corner, .5, .001)
        prepared, = [c for c in patch.layout.boundaries if hasattr(c, 'original_circle_interval')]
        circle = ClippedCircularCurb(owner, 0, prepared, corner)
        phase = prepared.original_circle_interval.circle_parameters
        stop = circle.source_parameter(.99*phase[-2] + .01*phase[-1])
        parameters = circle.certified_parameters(CertifiedCurbSampler(.01, 1., 1000), circle.circle_bounds[0], stop)
        self.assertEqual(parameters[-1], stop)
        gaps = np.diff(circle.circle_parameter(parameters))
        self.assertLessEqual(max(gaps), max(np.diff(phase))+1e-12)
        self.assertLess(max(gaps)/min(gaps), 1.01)


if __name__ == '__main__':
    unittest.main()
