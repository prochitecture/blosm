"""Actual source regions decide whether a closed incidence cycle has a hole."""
import unittest
import numpy as np

import _python_street_package
from tested_python_streets import model
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets.compounds import rim_label
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.source_support import ClippedApproachSupport
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder


def triangle(width,mode=0):
    a,b,c=(0.,0.,0.),(12.,0.,0.),(6.,10.392304845,0.)
    paths=((a,b),(b,c),(c,a),(a,(-15.,-10.,0.)),
           (b,(27.,-10.,0.)),(c,(6.,30.,0.)))
    sections=tuple(model.Section(i,i,points,width=width,interpolation=mode)
                   for i,points in enumerate(paths,1))
    sites=(model.Intersection(1,(model.End(1,0),model.End(3,1),model.End(4,0))),
           model.Intersection(2,(model.End(1,1),model.End(2,0),model.End(5,0))),
           model.Intersection(3,(model.End(2,1),model.End(3,0),model.End(6,0))))
    return model.Network('source-cycle',sections,sites)


def make_planner(network):
    patterns=PatternLibrary()
    products=dict(UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(network).products)
    context=BuildContext(network,products,patterns)
    layouts={}
    for site in network.intersections:
        ClippedApproachSupport._approaches(context,site)
        layouts[site.key]=IntersectionProvider.footprint(context,site)
    return DeferredCompoundPlanner(context,network.intersections,layouts)


def prepare(network):
    planner=make_planner(network)
    group=planner.exterior(tuple(site.key for site in network.intersections))
    return planner.context,group


class SourceCycleTests(unittest.TestCase):
    def test_interrupted_exterior_restarts_with_exact_source_spans_and_island(self):
        network = triangle(2.)
        authored = network.to_dict()
        planner = make_planner(network)
        members = tuple(site.key for site in network.intersections)
        interrupted = planner.exterior_steps(members)
        for _ in range(12):
            self.assertIsNone(next(interrupted))
        interrupted.close()
        restarted = planner.exterior(members)
        _, fresh = prepare(network)
        self.assertEqual(network.to_dict(), authored)
        self.assertEqual(restarted.claims, fresh.claims)
        self.assertEqual(len(restarted.loops), 2)
        self.assertEqual(len(restarted.loops), len(fresh.loops))
        for actual, expected in zip(restarted.loops, fresh.loops):
            self.assertEqual(actual.spans, expected.spans)
            np.testing.assert_array_equal(actual.xy, expected.xy)

    def test_physical_width_decides_hole_without_filling_the_graph_cycle(self):
        for width,holes in ((2.,1),(8.,0)):
            with self.subTest(width=width):
                network=triangle(width)
                before=network.to_dict()
                context,group=prepare(network)
                self.assertEqual(network.to_dict(),before)
                self.assertEqual(sum(loop.hole for loop in group.loops),holes)
                self.assertEqual({claim.section_id for claim in group.claims},{1,2,3})
                sources={source.region.owner:source for source in group.footprints}
                self.assertFalse(any(rim_label(span,sources).kind=='internal'
                    for loop in group.loops for span in loop.spans))
                sampled=AnalyticalBoundarySampler(spacing=network.settings.spacing,
                    tolerance=network.settings.curve_tolerance).build(group)
                data=[BoundaryConditions().build(sampled,loop,
                    {p.region.owner:p.layout for p in sampled.footprints}) for loop in sampled.loops]
                self.assertEqual(sum(len(boundary.ports) for boundary in data),3)
                for boundary in data:
                    if boundary.hole:
                        np.testing.assert_array_equal(boundary.heights,.08)

    def test_filled_source_pool_preserves_rim_and_checks_drains_and_joins(self):
        context,group=prepare(triangle(8.))
        surface=CompoundSurfaceBuilder(context.patterns).build(group,context)
        self.assertFalse(surface.holes)
        self.assertLessEqual(surface.mesh.validate()['maximum_edge_ratio'],8.)
        reconciled=surface.metadata['source_cycle_reconciliation']
        self.assertTrue(reconciled['exact_boundary'])
        self.assertEqual(reconciled['joins']['isolated_interior_minima'],0)
        self.assertLessEqual(reconciled['joins']['maximum_entrance_longitudinal_slope'],.03)


if __name__=='__main__':unittest.main()
