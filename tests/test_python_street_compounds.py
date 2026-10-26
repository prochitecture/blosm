import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, short_junction_pair
from check_python_street_footprints import ApproachInputs
from tested_python_streets import model
from tested_python_streets.compounds import CompoundPlanner
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import BuildContext, UpdateCoordinator


def footprint_context(network):
    patterns = PatternLibrary()
    prepared = UpdateCoordinator((PathProvider(), ProfileProvider(), ApproachInputs()), patterns).prepare(network)
    context = BuildContext(network, prepared.products, patterns)
    return context


def footprints(network):
    context = footprint_context(network)
    return {site.key: IntersectionProvider.footprint(context, site) for site in network.intersections}


class PythonCompoundPlanTests(unittest.TestCase):
    def test_short_link_has_one_combined_outline_without_erasing_source_sites(self):
        network = short_junction_pair(model)
        layouts = footprints(network)
        groups = CompoundPlanner().build(layouts)
        self.assertEqual(len(groups), 1)
        group = groups[0]
        self.assertEqual(group.members, (model.Key('intersection', 1), model.Key('intersection', 2)))
        self.assertEqual([c.section_id for c in group.claims], [1])
        self.assertEqual(len(group.loops), 1)
        self.assertFalse(group.loops[0].hole)
        self.assertEqual(len(network.intersections), 2)
        self.assertEqual(len(network.sections), 5)
        for span in group.loops[0].spans:
            sampled = next(p for p in group.footprints if p.region.owner == span.owner)
            self.assertIn(sampled.segments[span.segment].kind, ('mouth', 'curb'))

    def test_longer_link_keeps_independent_junctions(self):
        self.assertEqual(CompoundPlanner().build(footprints(short_junction_pair(model, length=20.))), ())

    def test_divided_cross_retains_the_separator_hole_and_physical_corner_arcs(self):
        network = compound_cross_network(model)
        groups = CompoundPlanner().build(footprints(network))
        self.assertEqual(len(groups), 1)
        group = groups[0]
        self.assertEqual(len(group.members), 4)
        self.assertEqual(len(group.claims), 4)
        self.assertEqual([loop.hole for loop in group.loops], [False, True])
        hole = group.loops[1]
        # New arrangement vertices lie on the sampled chords. Preserve the
        # intended analytical radius to the stated physical approximation
        # tolerance, rather than demanding an exact circle from a polygon.
        np.testing.assert_allclose(np.linalg.norm(hole.xy, axis=1), 2., rtol=0., atol=network.settings.curve_tolerance)
        self.assertLess(abs(hole.area + np.pi * 4.), .12)
        for span in hole.spans:
            sampled = next(p for p in group.footprints if p.region.owner == span.owner)
            self.assertEqual(sampled.segments[span.segment].kind, 'curb')


if __name__ == '__main__':
    unittest.main()
