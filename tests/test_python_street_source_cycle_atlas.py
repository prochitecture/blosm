"""Late source charts preserve successful seeds and exact physical guides."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.source_cycle_atlas import SourceCycleEnergyChart
from tested_python_streets.source_cycle_atlas import build
from tested_python_streets.four_side_charts import quality


def grid(displaced=False):
    xyz=np.array([(float(x),float(y),0.) for y in range(3) for x in range(3)])
    if displaced:xyz[4,:2]=(.01,1.)
    quads=np.array(((0,1,4,3),(1,2,5,4),(3,4,7,6),(4,5,8,7)),dtype=np.int32)
    return Mesh(xyz,quads,seams={0:np.array((0,1,2,5,8,7,6,3),dtype=np.int32)})


class SourceCycleAtlasTests(unittest.TestCase):
    def atlas_dispatch(self, accept):
        group=SimpleNamespace(claims=[SimpleNamespace(section_id=sid) for sid in (1,2,3,4)])
        context={}
        from tested_python_streets.model import Key
        for sid in (1,2,3,4):context[Key('path',sid)]=SimpleNamespace(length=float(sid))
        parent=SimpleNamespace(patterns=PatternLibrary(),spacing=.7,max_cells=50000,quality_iterations=100)
        calls=[]
        class Chart:
            def __init__(self,*args,**kwargs):
                self.policy=kwargs['partitioner'].policy;self.energy_recoveries=[]
                self.max_cells=kwargs['max_cells']
            def build(self,value,ctx):
                token=(self.policy.section_id,self.policy.near_endpoint_fraction)
                calls.append(token)
                self_outer.assertEqual(self.max_cells,50000)
                self_outer.assertEqual(self.policy.maximum_cuts,128)
                if not accept(token):raise ValueError('unqualified bounded chart')
                return SimpleNamespace(mesh=grid(),holes={},partition=SimpleNamespace(group=value),metadata={})
        self_outer=self
        with patch('tested_python_streets.source_cycle_atlas.eligible',return_value=True),\
             patch('tested_python_streets.source_cycle_atlas.certify_layout',return_value={'complete':True}),\
             patch('tested_python_streets.source_cycle_atlas.SourceCycleEnergyChart',Chart):
            try:result=build(parent,group,context,[])
            except ValueError:return calls,None
        return calls,result

    def test_all_original_corridors_precede_the_bounded_endpoint_retry(self):
        calls,result=self.atlas_dispatch(lambda token:token==(4,.125))
        self.assertEqual(calls,[(4,0.),(3,0.),(2,0.),(4,.125)])
        details=result.metadata['original_source_cycle_atlas']
        self.assertEqual(details['maximum_corridor_attempts'],6)
        self.assertEqual(details['maximum_original_corridors'],3)
        self.assertEqual(details['actual_atlas_attempts'],4)
        self.assertEqual(details['near_endpoint_fraction'],.125)

    def test_default_success_keeps_original_result_metadata_and_skips_retry(self):
        calls,result=self.atlas_dispatch(lambda token:token==(2,0.))
        self.assertEqual(calls,[(4,0.),(3,0.),(2,0.)])
        details=result.metadata['original_source_cycle_atlas']
        self.assertEqual(details['maximum_corridor_attempts'],3)
        for key in ('maximum_endpoint_policies','actual_atlas_attempts','near_endpoint_fraction','endpoint_reuse_reports'):
            self.assertNotIn(key,details)
        np.testing.assert_array_equal(result.mesh.xyz,grid().xyz)
        np.testing.assert_array_equal(result.mesh.quads,grid().quads)

    def test_total_failed_family_is_exactly_six_candidates_over_three_sources(self):
        calls,result=self.atlas_dispatch(lambda token:False)
        self.assertIsNone(result)
        self.assertEqual(calls,[(4,0.),(3,0.),(2,0.),(4,.125),(3,.125),(2,.125)])

    def test_already_qualified_chart_keeps_original_seed_and_skips_late_energy(self):
        mesh=grid();before=mesh.xyz.copy();chart=SourceCycleEnergyChart(PatternLibrary())
        result,details=chart._correct(mesh,retry_rotations=True)
        self.assertIs(result,mesh)
        np.testing.assert_array_equal(result.xyz,before)
        self.assertFalse(chart.energy_recoveries)
        self.assertEqual(details['stop'],'already_qualified')

    def test_attachment_releases_interior_guides_but_keeps_every_physical_point_and_id(self):
        mesh=grid(True);mesh.validate();chart=SourceCycleEnergyChart(PatternLibrary())
        ids=mesh.seams[0].copy();boundary=mesh.xyz[ids].copy();refs={int(i):('original',int(i)) for i in ids}
        result,newrefs,guides,details=chart.recover_separator_attachment(mesh,refs,np.arange(9),
            dict(qualified=False),[])
        self.assertIs(newrefs,refs)
        np.testing.assert_array_equal(result.seams[0],ids)
        np.testing.assert_array_equal(result.xyz[ids],boundary)
        self.assertGreater(np.linalg.norm(result.xyz[4]-mesh.xyz[4]),.1)
        self.assertLessEqual(quality(result)['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(quality(result)['minimum_corner_sine'],.03)
        self.assertTrue(details['source_attachment_global_correction']['interior_guides_released'])


if __name__=='__main__':unittest.main()
