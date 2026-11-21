"""A median-side support can join only consecutive native crest banks."""
import unittest
from types import SimpleNamespace

import _python_street_package
from tested_python_streets.bundle_crest_roots import crest_bundle
from tested_python_streets.model import Bundle, BundleMember, Traversal
from tested_python_streets.profiles import Crown


def arm(section, end, half, rise=.08):
    return SimpleNamespace(source=SimpleNamespace(id=section), end=end, profile=Crown(rise, 0., half))


class BundleCrestTests(unittest.TestCase):
    def setUp(self):
        self.bundle = Bundle(7, (BundleMember((Traversal(1), Traversal(2, False), Traversal(3))),
                                 BundleMember((Traversal(4),))))
        self.network = SimpleNamespace(bundles=(self.bundle,))

    def check_pair(self, first, second):
        return crest_bundle(self.network, SimpleNamespace(first=first, second=second))

    def test_consecutive_reversed_traversal(self):
        self.assertEqual(self.check_pair(arm(1, 1, -1), arm(2, 1, 1)), (self.bundle, 0, .08))
        self.assertEqual(self.check_pair(arm(3, 0, -1), arm(2, 0, 1)), (self.bundle, 0, .08))

    def test_unrelated_profiles_or_incidence(self):
        first = arm(1, 1, -1)
        for second in (arm(2, 0, 1), arm(3, 0, 1), arm(4, 0, 1), arm(2, 1, 0),
                       arm(2, 1, -1), arm(2, 1, 1, .07)):
            with self.subTest(section=second.source.id, end=second.end, profile=second.profile):
                self.assertIsNone(self.check_pair(first, second))


if __name__ == '__main__':
    unittest.main()
