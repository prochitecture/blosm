"""Generated cut contacts may reuse a nearby vertex without changing a rim."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_physical_holes import holed_region
from test_python_street_partition import inputs
from _python_street_fixtures import compound_cross_network
from tested_python_streets import model
from tested_python_streets.compounds import rim_label
from tested_python_streets.partition import BoundaryRef, CompoundPartitioner, CompoundQuadChart
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.physical_hole_cuts import PhysicalHoleCuts
from tested_python_streets.sampling import AnalyticalBoundarySampler


class CutEndpointReuseTests(unittest.TestCase):
    def test_only_strict_near_end_generated_contacts_reuse_the_chord_vertex(self):
        group=holed_region(((-2.,-1.,2.,1.),))
        before=[loop.xy.copy() for loop in group.loops]
        policy=PhysicalHoleCuts(near_endpoint_fraction=.125)
        for fraction,edge in ((.1,0),(.9,1)):
            self.assertEqual(policy.endpoint(group,(0.,0,0,fraction)),BoundaryRef(0,edge,0.))
        for fraction in (0.,.125,.5,.875,1.):
            self.assertEqual(policy.endpoint(group,(0.,0,0,fraction)),BoundaryRef(0,0,fraction))
        self.assertEqual(len(policy.endpoint_reuse_reports),2)
        for saved,loop in zip(before,group.loops):np.testing.assert_array_equal(saved,loop.xy)
        for row in policy.endpoint_reuse_reports:
            self.assertLess(min(row['fraction'],1.-row['fraction']),.125)
            self.assertGreater(row['displacement'],0.)

    def test_default_zero_preserves_old_cut_references_and_positive_mesh_byte_exact(self):
        group=holed_region(((-2.,-1.,2.,1.),))
        context=SimpleNamespace(network=SimpleNamespace(settings=SimpleNamespace(curve_tolerance=.005)))
        outputs=[]
        for policy in (PhysicalHoleCuts(),PhysicalHoleCuts(near_endpoint_fraction=0.)):
            chart=CompoundQuadChart(PatternLibrary(),four_sided=True,
                partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
            # This fixture already supplies its complete literal polygon rim.
            with patch.object(AnalyticalBoundarySampler,'build',lambda self,value:value):
                outputs.append(chart.build(group,context))
            self.assertFalse(policy.endpoint_reuse_reports)
        first,last=outputs
        self.assertEqual(first.partition.cuts,last.partition.cuts)
        np.testing.assert_array_equal(first.mesh.xyz,last.mesh.xyz)
        np.testing.assert_array_equal(first.mesh.quads,last.mesh.quads)
        self.assertEqual(first.metadata,last.metadata)

    def test_original_profiles_and_all_rim_vertices_survive_transformed_partitions(self):
        group,context=inputs(compound_cross_network(model))
        group=replace(group,physical_layout=object())
        sources={footprint.region.owner:footprint for footprint in group.footprints}
        for angle,offset in ((0.,(0.,0.)),(.713,(1200.,-1800.))):
            c,s=np.cos(angle),np.sin(angle);rotation=np.array(((c,-s),(s,c)))
            from tested_python_streets.regions import BoundaryLoop
            transformed=replace(group,loops=tuple(BoundaryLoop(loop.xy@rotation.T+offset,loop.spans)
                for loop in group.loops))
            before=[loop.xy.copy() for loop in transformed.loops]
            policy=PhysicalHoleCuts(near_endpoint_fraction=.125)
            result=CompoundPartitioner(policy=policy,require_visibility=False).build(transformed,context)
            self.assertLess(result.area_error,1e-7)
            for li,loop in enumerate(transformed.loops):
                boundary=result.xy[list(result.boundary_loops[li])]
                for i,point in enumerate(loop.xy):
                    self.assertTrue(np.any(np.all(boundary==point,axis=1)))
                    if rim_label(loop.spans[i],sources).kind=='mouth':
                        np.testing.assert_array_equal(point,before[li][i])
                np.testing.assert_array_equal(loop.xy,before[li])

    def test_option_validation_and_original_work_limits_remain_mandatory(self):
        for value in (True,False,-.001,.1250001,float('inf'),float('nan'),'0.1',None):
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'near-end'):
                PhysicalHoleCuts(near_endpoint_fraction=value)
        with self.assertRaisesRegex(ValueError,'near-end'):
            PhysicalHoleCuts(existing_vertices=True,near_endpoint_fraction=.125)
        group=holed_region(((-4.,-1.,-2.,1.),(2.,-1.,4.,1.)))
        with self.assertRaisesRegex(ValueError,'hole budget'):
            PhysicalHoleCuts(near_endpoint_fraction=.125,maximum_holes=1).cuts(group,None)
        with self.assertRaisesRegex(ValueError,'cut budget'):
            PhysicalHoleCuts(near_endpoint_fraction=.125,maximum_cuts=2).cuts(group,None)
        with self.assertRaisesRegex(ValueError,'graph-node budget'):
            CompoundPartitioner(policy=PhysicalHoleCuts(near_endpoint_fraction=.125),max_nodes=4).build(group,None)


if __name__=='__main__':unittest.main()
