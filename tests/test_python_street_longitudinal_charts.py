"""Direct quad charts for long curved fork regions; physical rims stay fixed."""
from types import SimpleNamespace
from unittest.mock import patch
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.compounds import CompoundGroup, RimSegment, RimSpan
from tested_python_streets.regions import BoundaryLoop
from tested_python_streets.model import Key
from tested_python_streets.longitudinal_cuts import LongitudinalCuts, even_sided_partition
from tested_python_streets.reflex_cuts import ReflexCurbCuts
from tested_python_streets.partition import (CompoundPartitioner, CompoundQuadChart,
    BoundaryRef, CutPolicy, JunctionCut)
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.four_side_charts import quality
from _python_street_compound_source_fixture import compound_source_fixture


def compound_context(network):
    """Evaluate physical preparation without unrelated surface providers."""
    from tested_python_streets.backend import PythonMeshBackend
    from tested_python_streets.ownership import PLAN
    from tested_python_streets.updates import BuildContext, DependencyIndex
    backend = PythonMeshBackend(network.settings)
    graph = DependencyIndex(tuple(task for provider in backend.providers
                                  for task in provider.tasks(network)))
    wanted, pending = set(), [PLAN]
    while pending:
        key = pending.pop()
        if key in wanted or key not in graph.tasks:
            continue
        wanted.add(key)
        pending.extend(graph.tasks[key].dependencies)
    context = BuildContext(network, {}, backend.patterns)
    for key in graph.order:
        if key in wanted:
            context.products[key] = graph.tasks[key].build(context)
    group = context[PLAN].by_site[Key('intersection', 794)]
    for junction in network.intersections:
        if junction.key in group.members:
            context.products[junction.layout_key] = group.physical_layout
    return context, group


def fork_region(angle=0., *, curved=False):
    points, labels = [], []
    owner = Key('junction', 10)
    def side(values, kind, index):
        points.extend(values[:-1])
        n = len(values)-1
        labels.extend(RimSpan(owner,RimSegment(kind,index,i/n,(i+1)/n)) for i in range(n))
    length = 100.
    side(np.column_stack((np.linspace(0.,length,251),np.zeros(251))), 'curb', 0)
    side(np.column_stack((np.full(17,length),np.linspace(0.,4.,17))), 'mouth', 0)
    t = np.linspace(-.5*math.pi,-1.5*math.pi,23)
    side(np.column_stack((length+np.cos(t),5.+np.sin(t))), 'curb', 1)
    side(np.column_stack((np.full(17,length),np.linspace(6.,10.,17))), 'mouth', 1)
    side(np.column_stack((np.linspace(length,0.,251),np.full(251,10.))), 'curb', 2)
    side(np.column_stack((np.zeros(17),np.linspace(10.,0.,17))), 'mouth', 2)
    points = np.asarray(points)
    if curved:
        points[:,1] += .0007*points[:,0]**2
    rotation = np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    points = points @ rotation.T + (1200.,-1800.)
    return CompoundGroup((owner,),(),(),(BoundaryLoop(points,tuple(labels)),))


def compact_return_region():
    points, labels = [], []
    owner = Key('junction', 11)
    def side(values, kind, index):
        points.extend(values[:-1])
        n = len(values)-1
        labels.extend(RimSpan(owner,RimSegment(kind,index,i/n,(i+1)/n)) for i in range(n))
    side(np.column_stack((np.linspace(0.,8.,17),np.zeros(17))), 'mouth',0)
    t=np.linspace(math.pi,0.,25)
    side(np.column_stack((10.+2.*np.cos(t),2.*np.sin(t))), 'curb',0)
    side(np.column_stack((np.linspace(12.,20.,17),np.zeros(17))), 'mouth',1)
    side(np.column_stack((np.full(41,20.),np.linspace(0.,20.,41))), 'curb',1)
    side(np.column_stack((np.linspace(20.,0.,17),np.full(17,20.))), 'mouth',2)
    side(np.column_stack((np.zeros(17),np.linspace(20.,12.,17))), 'curb',2)
    t=np.linspace(.5*math.pi,-.5*math.pi,25)
    side(np.column_stack((2.*np.cos(t),10.+2.*np.sin(t))), 'curb',3)
    side(np.column_stack((np.zeros(17),np.linspace(8.,0.,17))), 'curb',4)
    return CompoundGroup((owner,),(),(),(BoundaryLoop(np.asarray(points),tuple(labels)),))


class LongitudinalQuadTests(unittest.TestCase):
    def test_direct_charts_preserve_curved_fork_profiles_and_physical_outline(self):
        for angle in (0., .7, 2.9):
            with self.subTest(angle=angle):
                group = fork_region(angle,curved=True)
                context = SimpleNamespace(network=SimpleNamespace(settings=SimpleNamespace(curve_tolerance=.005)))
                chart = CompoundQuadChart(PatternLibrary(),spacing=.7,
                    partitioner=CompoundPartitioner(LongitudinalCuts(.7),require_visibility=False),
                    quality_iterations=20,four_sided=True)
                # This fixture already supplies the sampled analytical rim.
                with patch.object(AnalyticalBoundarySampler,'build',lambda self,value:value):
                    result = chart._build(group,context,require_quality=True)
                self.assertLess(result.partition.area_error,1e-7)
                self.assertLess(len(result.mesh.quads),10000)
                published=result.mesh.with_positions(publication_coordinates(result.mesh.xyz,1.,physical=True))
                report=quality(published)
                self.assertLessEqual(report['maximum_edge_ratio'],8.)
                self.assertGreaterEqual(report['minimum_corner_sine'],.03)
                original=group.loops[0]
                boundary=result.mesh.xyz[result.mesh.seams[0],:2]
                for point in original.xy:
                    self.assertTrue(np.any(np.all(boundary==point,axis=1)))
                for ref in result.boundary_refs.values():
                    if 1e-10<ref.fraction<1.-1e-10:
                        self.assertEqual(original.spans[ref.edge].label.kind,'curb')

    def test_fractional_cut_endpoint_resolves_external_loop_parity_on_curb(self):
        owner = Key('junction', 10)
        points = np.array(((0.,0.),(4.,0.),(4.,4.),(0.,4.)))
        labels = tuple(RimSpan(owner,RimSegment('curb',i,0.,1.)) for i in range(4))
        group = CompoundGroup((owner,),(),(),(BoundaryLoop(points,labels),))
        class FractionalCut(CutPolicy):
            def cuts(self, group, context):
                return (JunctionCut(None,owner,owner,BoundaryRef(0,0,.25),BoundaryRef(0,2,0.)),)
        partition = CompoundPartitioner(FractionalCut()).build(group,None)
        self.assertEqual(len(partition.boundary_loops[0]) % 2,0)
        self.assertEqual(len(partition.blocks),2)
        self.assertLess(partition.area_error,1e-10)
        self.assertTrue(all(len(rim) % 2 == 0 for rim in partition.sample(.7).rims))

    def test_compact_return_caps_use_conforming_reflex_charts(self):
        group=compact_return_region()
        self.assertTrue(ReflexCurbCuts.eligible(group))
        context=SimpleNamespace(network=SimpleNamespace(settings=SimpleNamespace(curve_tolerance=.005)))
        chart=CompoundQuadChart(PatternLibrary(),spacing=.7,
            partitioner=CompoundPartitioner(ReflexCurbCuts(),require_visibility=False),
            quality_iterations=100,four_sided=True)
        with patch.object(AnalyticalBoundarySampler,'build',lambda self,value:value):
            result=chart._build(group,context,require_quality=True)
        self.assertGreater(len(result.partition.blocks),1)
        self.assertLess(result.partition.area_error,1e-8)
        self.assertLessEqual(quality(result.mesh)['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(quality(result.mesh)['minimum_corner_sine'],.03)
        for ref in result.boundary_refs.values():
            if 1e-10<ref.fraction<1.-1e-10:
                self.assertEqual(group.loops[0].spans[ref.edge].label.kind,'curb')
        with self.assertRaisesRegex(ValueError,'graph-node budget'):
            even_sided_partition(result.partition,max_nodes=len(result.partition.xy)-1)

    def test_partition_eligibility_and_budget_are_explicit(self):
        group=fork_region()
        for values in ({'minimum_aspect':1.9},{'minimum_aspect':float('nan')},{'all_noses':1}):
            with self.assertRaises(ValueError):LongitudinalCuts(.7,**values)
        with self.assertRaisesRegex(ValueError,'cut budget'):
            CompoundPartitioner(LongitudinalCuts(.7,maximum_cuts=2)).build(group,None)
        with self.assertRaisesRegex(ValueError,'elongated corridor'):
            CompoundPartitioner(LongitudinalCuts(20.)).build(group,None)

    def test_original_compound_clears_stored_shape_height_and_profile_guards(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.conditions import BoundaryConditions
        for angle, shift in ((0., (0., 0.)), (.47, (-1000., 2000.))):
            with self.subTest(angle=angle):
                network = compound_source_fixture(angle, shift)
                original = network.to_dict()
                context, group = compound_context(network)
                surface = CompoundSurfaceBuilder(context.patterns).build(group, context)
                self.assertEqual(network.to_dict(), original)
                self.assertTrue(surface.chart.metadata['reflex_curb_partition'])
                self.assertLess(len(surface.mesh.quads), 5000)
                self.assertLess(surface.chart.partition.area_error, 1e-7)
                stored = surface.mesh.with_positions(publication_coordinates(
                    surface.mesh.xyz, network.settings.exaggeration, physical=True))
                self.assertLessEqual(stored.validate()['maximum_edge_ratio'], 8.)
                self.assertGreaterEqual(float(stored.xyz[:, 2].min()), 0.)
                self.assertLessEqual(float(stored.xyz[:, 2].max()), .080001)
                layouts = {p.region.owner: p.layout for p in surface.chart.partition.group.footprints}
                conditions = [BoundaryConditions().build(surface.chart.partition.group, loop, layouts)
                              for loop in surface.chart.partition.group.loops]
                for vertex, ref in surface.boundary_refs.items():
                    expected = conditions[ref.loop].edge_samples(ref.edge, ref.fraction)
                    np.testing.assert_allclose(surface.mesh.xyz[vertex],
                        (*expected[0], expected[1]), atol=1e-8, rtol=0.)

    def test_reflex_cut_limits_and_profile_sample_parity_are_explicit(self):
        group = compact_return_region()
        with self.assertRaisesRegex(ValueError, 'cut budget'):
            CompoundPartitioner(ReflexCurbCuts(maximum_cuts=1),
                require_visibility=False).build(group, None)
        for value in (0., math.pi, float('nan')):
            with self.assertRaises(ValueError):
                ReflexCurbCuts(minimum_turn=value)
        context = SimpleNamespace(network=SimpleNamespace(settings=SimpleNamespace(curve_tolerance=.005)))
        policy = ReflexCurbCuts()
        cuts = policy.cuts(group, context)
        self.assertTrue(cuts)
        for cut in cuts:
            self.assertEqual(cut.first.fraction, 0.)
            self.assertEqual(cut.last.fraction, 0.)
            self.assertTrue(0 <= cut.last.edge < len(group.loops[0].xy))
        # A target on an authored 16-interval profile uses an existing even
        # sample. It cannot add a new subdivision to an independent road seam.
        points = group.loops[0].xy
        labels = [span.label for span in group.loops[0].spans]
        self.assertEqual(policy.profile_endpoint(points, labels, 3, .2) % 2, 0)


if __name__=='__main__':unittest.main()
