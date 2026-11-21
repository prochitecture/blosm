"""Shared-coordinate branches must not choose a Bundle through road by order."""
import importlib.util
from itertools import permutations
from pathlib import Path
from types import SimpleNamespace
import unittest

spec = importlib.util.spec_from_file_location('bundle_continuations',
    Path(__file__).resolve().parents[1] / 'way' / 'bundle_continuations.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class Street:
    def __init__(self, source_id, a, b):
        self.src, self.dst = a, b
        self.head = self.tail = SimpleNamespace(sourceWayIds=frozenset((source_id,)))


def record(street, point):
    return dict(street=street, end=point)


class BundleMergePairingTests(unittest.TestCase):
    def setUp(self):
        self.a, self.b = (0., 0.), (0., 10.)
        self.before = [Street('main', (-50., 0.), self.a),
                       Street('opposite', self.b, (-50., 10.))]
        self.after = [Street('main', self.a, (50., 0.)),
                      Street('opposite', (50., 10.), self.b)]
        self.branch = Street('bus_lane', self.a, (50., -8.))

    def test_original_way_continuation_is_independent_of_all_container_orders(self):
        left = [record(s, p) for s, p in zip(self.before, (self.a, self.b))]
        right = [record(s, p) for s, p in zip(self.after, (self.a, self.b))]
        right.append(record(self.branch, self.a))
        for lhs in permutations(left):
            for rhs in permutations(right):
                for groups in ((lhs, rhs), (rhs, lhs)):
                    with self.subTest(groups=groups):
                        result = module.bundle_merge_pairs(dict(enumerate(groups)))
                        self.assertEqual(result, ((self.a, self.before[0], self.after[0]),
                                                  (self.b, self.after[1], self.before[1])))
                        # Selection does not remove or mutate either original branch.
                        self.assertIn(record(self.branch, self.a), right)
                        self.assertEqual(self.branch.src, self.a)

    def test_uncertain_original_continuation_does_not_pick_an_arbitrary_branch(self):
        left = [record(s, p) for s, p in zip(self.before, (self.a, self.b))]
        right = [record(s, p) for s, p in zip(self.after, (self.a, self.b))]
        right.append(record(self.branch, self.a))
        for source in ('main', 'unknown'):
            self.branch.head.sourceWayIds = frozenset((source,))
            self.after[0].head.sourceWayIds = frozenset((source,))
            self.assertIsNone(module.bundle_merge_pairs({0: left, 1: right}))

    def test_two_simple_members_do_not_require_matching_osm_way_ids(self):
        self.after[0].head.sourceWayIds = frozenset(('different_source',))
        groups = {i: [record(s, p) for s, p in zip(streets, (self.a, self.b))]
                  for i, streets in enumerate((self.before, self.after))}
        self.assertEqual(len(module.bundle_merge_pairs(groups)), 2)

    def test_reversed_or_mismatched_incidence_does_not_create_a_continuation(self):
        self.after[0].src, self.after[0].dst = self.after[0].dst, self.after[0].src
        groups = {i: [record(s, p) for s, p in zip(streets, (self.a, self.b))]
                  for i, streets in enumerate((self.before, self.after))}
        self.assertIsNone(module.bundle_merge_pairs(groups))
        groups[1][0]['end'] = (0., .001)
        self.assertIsNone(module.bundle_merge_pairs(groups))


if __name__ == '__main__':
    unittest.main()
