import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Key
from tested_python_streets.regions import Region, RegionUnion


def box(identity, x0, y0, x1, y1):
    return Region(Key('intersection', identity), np.array(((x0, y0), (x1, y0), (x1, y1), (x0, y1))))


class PythonRegionUnionTests(unittest.TestCase):
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

    def test_point_only_contact_and_nonsimple_region_are_diagnosed(self):
        with self.assertRaisesRegex(ValueError, 'point-only'):
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
