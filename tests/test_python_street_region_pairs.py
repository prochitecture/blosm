import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Key
from tested_python_streets.region_pairs import PairBudget, SegmentPairIndex
from tested_python_streets.regions import Region, RegionUnion
from test_python_street_regions import box


class ForcedIndexedUnion(RegionUnion):
    """Compare the new narrow path with the retained exhaustive reference."""
    def _split_pair(self, first, second, offsets, cuts, *, same=False):
        self._split_pair_indexed(first, second, offsets, cuts,
                                 PairBudget(self.max_pair_tests), same=same)


class PythonRegionPairTests(unittest.TestCase):
    def assert_same_union(self, regions):
        expected, actual = RegionUnion().build(regions), ForcedIndexedUnion().build(regions)
        self.assertEqual(len(expected), len(actual))
        for first, second in zip(expected, actual):
            self.assertEqual(first.xy.tobytes(), second.xy.tobytes())
            self.assertEqual(first.spans, second.spans)

    def test_indexed_splitting_retains_exact_union_buffers_and_provenance(self):
        cases = ((box(1, 0, 0, 3, 2), box(2, 2, 1, 4, 3)),
                 (box(1, 0, 0, 2, 2), box(2, 1, 0, 3, 1)),
                 (box(1, 0, 0, 2, 2), box(2, 2, 0, 4, 2)),
                 (box(1, 0, 0, 4, 4), box(2, 1, 1, 2, 2)),
                 (box(1, 0, 0, 2, 2), box(2, 0, 0, 2, 2)),
                 (box(1, 0, 0, 4, 1), box(2, 3, 1, 4, 3),
                  box(3, 0, 3, 4, 4), box(4, 0, 1, 1, 3)))
        for regions in cases:
            for offset in (np.zeros(2), np.array((1e6, -2e6))):
                with self.subTest(regions=len(regions), offset=tuple(offset)):
                    self.assert_same_union(tuple(Region(r.owner, r.xy + offset) for r in regions))

    def test_rotated_unequal_sampled_loops_match_exhaustive_arithmetic(self):
        for seed in range(8):
            random = np.random.default_rng(seed)
            regions = []
            for index in range(3):
                angle = np.sort(random.uniform(0., 2. * np.pi, 19 + index * 7))
                xy = np.column_stack((3. * np.cos(angle), 2. * np.sin(angle)))
                xy += (index * 1.7, index * .4)
                regions.append(Region(Key('intersection', index), xy))
            self.assert_same_union(regions)

    def test_candidate_boxes_are_complete_including_near_epsilon_and_translation(self):
        random = np.random.default_rng(412)
        for offset in (np.zeros(2), np.array((1e9, -1e9))):
            low = random.uniform(-20., 20., (220, 2)) + offset
            high = low + random.uniform(0., 6., (220, 2))
            # Explicit contacts and near contacts on both sides of epsilon.
            low[1], high[1] = high[0] + .9e-7, high[0] + 2.
            index = SegmentPairIndex(low, high, 1e-7)
            for same in (False, True):
                batches = tuple(index.batches(low, high, PairBudget(1_000_000), same=same))
                actual = [(int(i), int(j)) for rows, columns in batches for i, j in zip(rows, columns)]
                nearby = np.all(high[:, None] + 1e-7 >= low, axis=-1)
                nearby &= np.all(high + 1e-7 >= low[:, None], axis=-1)
                if same:
                    rows, columns = np.indices(nearby.shape)
                    nearby &= (columns > rows + 1) & ~((rows == 0) & (columns == len(low) - 1))
                self.assertEqual(actual, list(zip(*np.nonzero(nearby))))
                self.assertTrue(all(len(rows) <= SegmentPairIndex.BATCH_SIZE for rows, _ in batches))

    def test_long_edges_and_large_queries_do_not_expand_unbounded_grids(self):
        low = np.zeros((512, 2))
        high = np.full_like(low, 100.)
        index = SegmentPairIndex(low, high, 1e-7)
        self.assertEqual(len(index.long_edges), len(low))
        self.assertFalse(index.buckets)
        batches = tuple(index.batches(low[:1], high[:1], PairBudget(513)))
        self.assertEqual(sum(len(rows) for rows, _ in batches), 512)
        with self.assertRaisesRegex(ValueError, 'intersection-test budget'):
            tuple(index.batches(low[:1], high[:1], PairBudget(512)))

    def test_dense_sampled_curve_uses_spatial_fallback_under_unchanged_limit(self):
        angle = np.linspace(0., 2. * np.pi, 4100, endpoint=False)
        region = Region(Key('intersection', 1), np.column_stack((10. * np.cos(angle), 10. * np.sin(angle))))
        union = RegionUnion()
        self.assertEqual(union.max_pair_tests, 8_000_000)
        self.assertGreater(len(angle) * (len(angle) - 1) // 2, union.max_pair_tests)
        loops = union.build((region,))
        self.assertEqual(len(loops), 1)
        self.assertEqual(len(loops[0].xy), len(angle))
        self.assertAlmostEqual(loops[0].area, region.xy.shape[0] * 50. * np.sin(2. * np.pi / len(angle)))
        self.assertTrue(all(span.first == 0. and span.last == 1. for span in loops[0].spans))

    def test_repeated_bucket_entries_are_charged_before_candidate_set_work(self):
        low, high = np.zeros((20, 2)), np.ones((20, 2))
        index = SegmentPairIndex(low, high, 1e-7)
        self.assertFalse(index.long_edges)
        with self.assertRaisesRegex(ValueError, 'intersection-test budget'):
            tuple(index.batches(low[:1], high[:1], PairBudget(21)))
        actual = tuple(index.batches(low[:1], high[:1], PairBudget(10_000)))
        self.assertEqual(sum(len(rows) for rows, _ in actual), 20)

    def test_contacts_and_self_intersections_keep_original_diagnostics(self):
        for union in (RegionUnion, ForcedIndexedUnion):
            with self.assertRaisesRegex(ValueError, 'point-only'):
                union().build((box(1, 0, 0, 1, 1), box(2, 1, 1, 2, 2)))
            crossing = Region(Key('intersection', 1), np.array(((0., 0.), (3., 0.), (1., 3.),
                                                               (1., -1.), (3., 3.), (0., 3.))))
            with self.assertRaisesRegex(ValueError, 'itself'):
                union().build((crossing,))


if __name__ == '__main__':
    unittest.main()
