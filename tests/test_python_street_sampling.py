from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, short_junction_pair
from test_python_street_partition import inputs
from tested_python_streets import model
from tested_python_streets.compounds import RimSpan, rim_label
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.sampling import AnalyticalBoundarySampler


class PythonCompoundSamplingTests(unittest.TestCase):
    def test_every_retained_curb_is_on_its_original_analytical_curve(self):
        group, context = inputs(compound_cross_network(model))
        sampled = AnalyticalBoundarySampler().build(group)
        layouts = {p.region.owner: p.layout for p in group.footprints}
        sources = {p.region.owner: p for p in group.footprints}
        self.assertEqual(sampled.members, group.members)
        self.assertEqual(sampled.claims, group.claims)
        self.assertEqual(len(sampled.loops), 2)
        for loop in sampled.loops:
            for i, span in enumerate(loop.spans):
                label = rim_label(span, sources)
                if label.kind != 'curb':
                    continue
                curve = layouts[span.owner].boundaries[label.index]
                np.testing.assert_allclose(loop.xy[i], curve(label.first), atol=1e-8, rtol=0.)
                np.testing.assert_allclose(loop.xy[(i + 1) % len(loop.xy)], curve(label.last), atol=1e-8, rtol=0.)
                t = np.linspace(0., 1., 101)
                actual = curve(label.first + t * (label.last - label.first))
                chord = loop.xy[i] + t[:, None] * (loop.xy[(i + 1) % len(loop.xy)] - loop.xy[i])
                self.assertLessEqual(float(np.linalg.norm(actual - chord, axis=1).max()), .005 + 1e-8)

    def test_all_external_profile_samples_and_conditions_remain_authoritative(self):
        group, context = inputs(compound_cross_network(model))
        sampled = AnalyticalBoundarySampler().build(group)
        layouts = {p.region.owner: p.layout for p in group.footprints}
        first = BoundaryConditions().build(group, group.loops[0], layouts)
        second = BoundaryConditions().build(sampled, sampled.loops[0], layouts)
        self.assertEqual(len(first.ports), len(second.ports))
        for a in first.ports:
            b = next(p for p in second.ports if (p.section_id, p.end) == (a.section_id, a.end))
            self.assertEqual(a.source_u, b.source_u)
            self.assertEqual(a.station, b.station)
            np.testing.assert_array_equal(first.xy[list(a.vertices)], second.xy[list(b.vertices)])
            np.testing.assert_array_equal(first.heights[list(a.vertices)], second.heights[list(b.vertices)])
            np.testing.assert_array_equal(first.gradients[list(a.vertices)], second.gradients[list(b.vertices)])

    def test_incidental_collinear_samples_do_not_survive_as_tiny_edges(self):
        group, context = inputs(short_junction_pair(model, length=6.127))
        first = AnalyticalBoundarySampler().build(group)
        # Inject an incidental clipping vertex into a curb chord. Its position
        # does not introduce a new analytical curve or semantic ownership event.
        loop = first.loops[0]
        at = next(i for i, span in enumerate(loop.spans) if span.label.kind == 'curb')
        span = loop.spans[at]
        t = span.label.first + .001 * (span.label.last - span.label.first)
        curve = next(p.layout for p in group.footprints if p.region.owner == span.owner).boundaries[span.label.index]
        xy = np.insert(loop.xy, at + 1, curve(t), axis=0)
        spans = (*loop.spans[:at], RimSpan(span.owner, replace(span.label, last=t)),
                 RimSpan(span.owner, replace(span.label, first=t)), *loop.spans[at + 1:])
        split = replace(first, loops=(replace(loop, xy=xy, spans=spans),))
        restored = AnalyticalBoundarySampler().build(split)
        self.assertEqual(len(restored.loops[0].xy), len(first.loops[0].xy))
        np.testing.assert_allclose(restored.loops[0].xy, first.loops[0].xy, atol=1e-10, rtol=0.)

    def test_sampling_is_deterministic_and_does_not_mutate_old_footprints(self):
        group, context = inputs(compound_cross_network(model))
        original = [loop.xy.copy() for loop in group.loops]
        sampler = AnalyticalBoundarySampler()
        first, second = sampler.build(group), sampler.build(group)
        for a, b, old, expected in zip(first.loops, second.loops, group.loops, original):
            np.testing.assert_array_equal(a.xy, b.xy)
            self.assertEqual(a.spans, b.spans)
            np.testing.assert_array_equal(old.xy, expected)
            self.assertFalse(a.xy.flags.writeable)

    def test_missing_source_geometry_and_unbounded_requests_are_diagnosed(self):
        group, context = inputs(short_junction_pair(model))
        missing = replace(group, footprints=tuple(replace(p, layout=None) for p in group.footprints))
        with self.assertRaisesRegex(ValueError, 'source footprint layouts'):
            AnalyticalBoundarySampler().build(missing)
        with self.assertRaisesRegex(ValueError, 'sample budget'):
            AnalyticalBoundarySampler(max_samples=4).build(group)
        for settings in ({'tolerance': 0.}, {'tolerance': float('nan')}, {'max_samples': 3}):
            with self.assertRaises(ValueError):
                AnalyticalBoundarySampler(**settings)


if __name__ == '__main__':
    unittest.main()
