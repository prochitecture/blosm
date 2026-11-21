"""Source-backed Bundle eligibility and continuation, in isolated Blender."""
from itertools import permutations
from types import SimpleNamespace
import unittest

from test_street_preparation import source, graph
from blosm.way.item import Intersection, Section
from blosm.way.item.bundle import (bundleCenterline, Bundle, canBeMerged,
                                  mergeBundles, StreetGroup, removeSplittingStreets,
                                  _longitudinalGroupLink)


class BundleSourcePreparationTests(unittest.TestCase):
    def test_group_overlap_keeps_longitudinal_original_connector(self):
        # Both ends intersect a second heuristic group, but this is one
        # carriageway's longitudinal source segment, not a median crossing.
        for reverse in (False, True):
            for order in permutations(range(4)):
                with self.subTest(reverse=reverse, order=order):
                    points = ((0., 10.), (40., 10.), (48., 10.), (100., 11.))
                    if reverse:
                        points = points[::-1]
                    upper = [source(a, b) for a, b in zip(points, points[1:])]
                    lower = source((100., 0.), (0., 0.))
                    streets = upper+[lower]
                    group = StreetGroup([streets[i] for i in order])
                    waymap=graph(streets)
                    for _, node in waymap.iterNodes(Intersection):
                        node.processIntersection()
                    links = [(s.pred, s.succ) for s in streets]
                    controls = [tuple(s.head.centerline) for s in streets]
                    marker = group.id+1000
                    intersections = {p: [group.id, marker] for p in points[1:3]}
                    changed, _, removed = removeSplittingStreets(None, group.id, group, intersections)
                    self.assertFalse(changed)
                    self.assertFalse(removed)
                    self.assertEqual(links, [(s.pred, s.succ) for s in streets])
                    self.assertEqual(controls, [tuple(s.head.centerline) for s in streets])
                    self.assertEqual(set(group), set(streets))

    def test_longitudinal_link_requires_unique_directed_source_incidence(self):
        link=source((0., 0.), (8., 0.))
        before=source((-40., 0.), (0., 0.))
        after=source((8., 0.), (40., 0.))
        self.assertTrue(_longitudinalGroupLink(link, [before, link, after]))
        # A one-ended candidate can continue a source toward the import frame.
        self.assertTrue(_longitudinalGroupLink(link, [before, link]))
        curved=source((8., 0.), (40., 10.))
        self.assertTrue(_longitudinalGroupLink(link, [before, link, curved]))
        turn=source((8., 0.), (40., 32.))
        self.assertFalse(_longitudinalGroupLink(link, [before, link, turn]))
        transverse=source((0., -10.), (0., 0.))
        self.assertFalse(_longitudinalGroupLink(link, [link, transverse]))
        near=source((-40., 0.), (-.001, 0.))
        self.assertFalse(_longitudinalGroupLink(link, [near, link]))
        wrong_direction=source((40., 0.), (8., 0.))
        self.assertFalse(_longitudinalGroupLink(link, [before, link, wrong_direction]))
        ambiguous=source((8., 0.), (40., 1.))
        self.assertFalse(_longitudinalGroupLink(link, [before, link, after, ambiguous]))

    def test_open_roundabout_arc_stays_an_ordinary_source_street(self):
        for kind in ('roundabout', 'circular'):
            street = source((0., 0.), (10., 5.))
            street.head.tags['junction'] = kind
            controls = tuple(street.head.centerline)
            self.assertIsNone(bundleCenterline(street))
            self.assertIsInstance(street.head, Section)
            self.assertEqual(tuple(street.head.centerline), controls)
            self.assertIsNone(street.bundle)
        ordinary = source((0., 0.), (10., 5.))
        self.assertIsNotNone(bundleCenterline(ordinary))

    def test_merged_bundle_keeps_main_source_continuation_and_all_branch_geometry(self):
        for order in permutations(range(3)):
            with self.subTest(order=order):
                a, b = (0., 0.), (0., 10.)
                first = [source((-50., 0.), a), source(b, (-50., 10.))]
                second = [source(a, (50., 0.)), source((50., 10.), b)]
                bus = source(a, (50., -8.))
                upper_branch = source(b, (0., 35.))
                all_streets = first + second + [bus, upper_branch]
                for sid, street in zip(('main', 'opposite', 'main', 'opposite', 'bus'),
                                       first + second + [bus]):
                    street.head.sourceWayIds = frozenset((sid,))
                sections = [s.head for s in all_streets]
                controls = [tuple(s.centerline) for s in sections]
                waymap = graph(all_streets)
                nodes = [node for _, node in waymap.iterNodes(Intersection)]
                for node in nodes:
                    node.processIntersection()
                bundles = (Bundle(), Bundle())
                members = (first, [second[0], second[1], bus])
                for bundle, streets in zip(bundles, members):
                    for street in streets:
                        street.bundle = bundle
                rows = [[dict(street=s, end=s.src if s in (first[1], second[0], bus) else s.dst,
                              bundle=bundle, type='tail' if i == 0 else 'head')
                         for s in streets] for i, (bundle, streets) in enumerate(zip(bundles, members))]
                involved = {bundles[0]: rows[0], bundles[1]: [rows[1][i] for i in order]}
                manager = SimpleNamespace(streets={s.id: s for s in all_streets},
                    bundles={s.id: s for s in bundles}, minorIntersections={},
                    majorIntersections={n.location: n for n in nodes if n.startConnector})
                self.assertTrue(canBeMerged(manager, involved))
                mergeBundles(manager, involved)
                result, = manager.bundles.values()
                through = next(s for s in result.streetsHead if sections[0] in list(s.iterItems()))
                self.assertIn(sections[2], list(through.iterItems()))
                self.assertNotIn(bus.head, list(through.iterItems()))
                self.assertIs(manager.streets[bus.id], bus)
                self.assertIs(bus.pred.intersection, manager.minorIntersections[bus.src])
                self.assertEqual([tuple(s.centerline) for s in sections], controls)
                self.assertEqual({s.id for street in manager.streets.values() for s in street.iterItems()
                                  if isinstance(s, Section)}, {s.id for s in sections})


if __name__ == '__main__':
    unittest.main()
