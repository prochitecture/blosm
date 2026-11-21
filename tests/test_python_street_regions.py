import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Key
from tested_python_streets.regions import BoundaryIncidenceError, Region, RegionUnion


def box(identity, x0, y0, x1, y1):
    return Region(Key('intersection', identity), np.array(((x0, y0), (x1, y0), (x1, y1), (x0, y1))))


class PythonRegionUnionTests(unittest.TestCase):
    def test_indexed_union_preserves_owner_ties_and_clipped_edge_parameters(self):
        from tested_python_streets.prepared_source_classifier import SourceProofClassifier
        from tested_python_streets.regions import _classify
        square = box(1, 0., 0., 3., 2.).xy
        fractions = np.arange(40)/40.
        rim = (square[:, None] + fractions[None, :, None]
               * (np.roll(square, -1, axis=0)-square)[:, None]).reshape(-1, 2)
        regions = (Region(Key('intersection', 1), rim), Region(Key('intersection', 2), rim),
                   box(3, 1.01, .5, 4., 1.52))
        with patch.object(SourceProofClassifier, '__call__',
                          lambda self, points, xy, epsilon: _classify(points, xy, epsilon)):
            reference, = RegionUnion().build(regions)
        actual, = RegionUnion().build(reversed(regions))
        np.testing.assert_array_equal(actual.xy, reference.xy)
        self.assertEqual(actual.spans, reference.spans)
        self.assertEqual({s.owner.id for s in actual.spans}, {1, 3})
        self.assertTrue(any(0. < s.first < 1. or 0. < s.last < 1. for s in actual.spans))

    def assert_provenance(self, regions, loops):
        source = {r.owner: r for r in regions}
        for loop in loops:
            self.assertEqual(len(loop.xy), len(loop.spans))
            self.assertFalse(loop.xy.flags.writeable)
            for i, span in enumerate(loop.spans):
                region = source[span.owner]
                a, b = region.xy[span.segment], region.xy[(span.segment + 1) % len(region.xy)]
                np.testing.assert_allclose(loop.xy[i], a + span.first * (b - a), atol=1.1e-7, rtol=0.)
                np.testing.assert_allclose(loop.xy[(i + 1) % len(loop.xy)], a + span.last * (b - a), atol=1.1e-7, rtol=0.)

    def test_crossing_edges_split_and_keep_only_union_boundary(self):
        regions = (box(1, 0, 0, 3, 2), box(2, 2, 1, 4, 3))
        loops = RegionUnion().build(regions)
        self.assertEqual(len(loops), 1)
        self.assertAlmostEqual(loops[0].area, 9.)
        self.assertEqual({span.owner for span in loops[0].spans}, {r.owner for r in regions})
        self.assert_provenance(regions, loops)

    def test_partial_collinear_overlap_and_shared_edge_contact_merge(self):
        for regions, area in (((box(1, 0, 0, 2, 2), box(2, 1, 0, 3, 1)), 5.),
                              ((box(1, 0, 0, 2, 2), box(2, 2, 0, 4, 2)), 8.)):
            with self.subTest(area=area):
                loops = RegionUnion().build(regions)
                self.assertEqual(len(loops), 1)
                self.assertAlmostEqual(loops[0].area, area)
                self.assert_provenance(regions, loops)

    def test_containment_duplicate_regions_and_disconnected_components(self):
        for regions, areas in (((box(1, 0, 0, 4, 4), box(2, 1, 1, 2, 2)), [16.]),
                               ((box(1, 0, 0, 2, 2), box(2, 0, 0, 2, 2)), [4.]),
                               ((box(1, 0, 0, 1, 1), box(2, 2, 0, 3, 1)), [1., 1.])):
            with self.subTest(areas=areas):
                loops = RegionUnion().build(regions)
                self.assertEqual(sorted(loop.area for loop in loops), areas)
                self.assert_provenance(regions, loops)

    def test_ring_union_preserves_clockwise_hole_for_separator_ownership(self):
        regions = (box(1, 0, 0, 4, 1), box(2, 3, 1, 4, 3), box(3, 0, 3, 4, 4), box(4, 0, 1, 1, 3))
        loops = RegionUnion().build(regions)
        self.assertEqual([loop.hole for loop in loops], [False, True])
        self.assertEqual([loop.area for loop in loops], [16., -4.])
        self.assert_provenance(regions, loops)

    def test_near_coincident_differently_sampled_boundary_uses_first_owner(self):
        # A straight curb and two chords describe the same boundary within the
        # numerical welding distance. Retaining both introduces a graph branch
        # at the extra sample even though there is no physical point contact.
        straight = np.array(((0., 0.), (2., 0.), (2., 1.), (0., 1.)))
        sampled = np.array(((1., 0.), (1.5, 5e-8), (3., 0.), (3., 1.), (1., 1.)))
        for angle in (0., .713):
            c, s = np.cos(angle), np.sin(angle)
            rotation = np.array(((c, -s), (s, c)))
            for offset in (np.zeros(2), np.array((1e6, -2e6))):
                for identities in ((261, 1211), (1211, 261)):
                    with self.subTest(angle=angle, offset=tuple(offset), owners=identities):
                        regions = tuple(Region(Key('intersection', owner), xy @ rotation.T + offset)
                                        for owner, xy in zip(identities, (straight, sampled)))
                        loop, = RegionUnion().build(regions)
                        self.assertAlmostEqual(loop.area, 3., places=6)
                        self.assertEqual({span.owner for span in loop.spans}, {r.owner for r in regions})
                        self.assert_provenance(regions, (loop,))
                        # All of the overlapping bottom interval belongs to
                        # the lower owner, independent of input enumeration.
                        local = (loop.xy - offset) @ rotation
                        following = np.roll(local, -1, axis=0)
                        middle = (local + following) * .5
                        overlap = ((middle[:, 0] > 1.) & (middle[:, 0] < 2.)
                                   & (abs(middle[:, 1]) < 1e-7))
                        self.assertTrue(np.any(overlap))
                        self.assertEqual({loop.spans[i].owner.id for i in np.flatnonzero(overlap)}, {261})
                        repeated, = RegionUnion().build(regions[::-1])
                        np.testing.assert_array_equal(loop.xy, repeated.xy)
                        self.assertEqual(loop.spans, repeated.spans)

    def test_enumeration_and_large_coordinate_translation_preserve_geometry(self):
        regions = (box(1, 0, 0, 3, 2), box(2, 2, 1, 4, 3))
        first = RegionUnion().build(regions)
        reversed_order = RegionUnion().build(regions[::-1])
        np.testing.assert_array_equal(first[0].xy, reversed_order[0].xy)
        offset = np.array((1_000_000., -2_000_000.))
        shifted = tuple(Region(r.owner, r.xy + offset) for r in regions)
        actual = RegionUnion().build(shifted)
        np.testing.assert_allclose(first[0].xy, actual[0].xy - offset, rtol=0., atol=1e-8)
        self.assert_provenance(shifted, actual)

    def test_shallow_crossing_is_not_a_coincident_owner_tie_at_its_midpoint(self):
        # The post-crossing midpoint is 86 nm from the other edge, while
        # its endpoint is 175 nm away. A midpoint-only owner tie retains
        # both outgoing boundaries and produces a spurious graph branch.
        first = np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))
        second = np.array(((.025, 5.63e-7), (1.036, -2.02e-7),
                           (1.036, 1.), (.025, 1.)))
        for identities in ((1141, 1142), (1142, 1141)):
            original = tuple(Region(Key('intersection', owner), xy)
                             for owner, xy in zip(identities, (first, second)))
            for angle, offset in ((0., np.zeros(2)), (.713, np.array((1000., -1100.)))):
                c, s = np.cos(angle), np.sin(angle)
                rotation = np.array(((c, -s), (s, c)))
                regions = tuple(Region(r.owner, r.xy @ rotation.T + offset) for r in original)
                expected, = RegionUnion(epsilon=1e-9).build(regions)
                loop, = RegionUnion().build(regions)
                self.assertAlmostEqual(loop.area, expected.area, places=10)
                self.assert_provenance(regions, (loop,))
                distances = np.linalg.norm(loop.xy[:, None]-expected.xy, axis=-1)
                self.assertLess(float(distances.min(axis=0).max()), 1e-8)
                self.assertLess(float(distances.min(axis=1).max()), 1e-8)
                repeated, = RegionUnion().build(regions[::-1])
                np.testing.assert_array_equal(loop.xy, repeated.xy)
                self.assertEqual(loop.spans, repeated.spans)

    def test_point_only_contact_and_nonsimple_region_are_diagnosed(self):
        with self.assertRaisesRegex(BoundaryIncidenceError, 'point-only'):
            RegionUnion().build((box(1, 0, 0, 1, 1), box(2, 1, 1, 2, 2)))
        crossing = Region(Key('intersection', 1), np.array(((0., 0.), (3., 0.), (1., 3.), (1., -1.), (3., 3.), (0., 3.))))
        with self.assertRaisesRegex(ValueError, 'itself'):
            RegionUnion().build((crossing,))

    def test_segment_and_intersection_budgets_are_checked_before_dense_work(self):
        regions = (box(1, 0, 0, 3, 2), box(2, 2, 1, 4, 3))
        with self.assertRaisesRegex(ValueError, 'segment budget'):
            RegionUnion(max_segments=7).build(regions)
        with self.assertRaisesRegex(ValueError, 'intersection-test budget'):
            RegionUnion(max_pair_tests=2).build(regions)
        self.assertEqual(RegionUnion().build(()), ())


if __name__ == '__main__':
    unittest.main()
