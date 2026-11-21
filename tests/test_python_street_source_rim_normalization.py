"""Derived curb projection preserves profiles, lineage and layout metadata."""
from types import SimpleNamespace
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Key,Section
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach,JunctionLayout
from tested_python_streets.finite_free_caps import PreparedCapCurb
from tested_python_streets.compounds import CompoundGroup,sample_footprint,rim_label
from tested_python_streets.regions import RegionUnion
from tested_python_streets.source_rim_normalization import normalize,normalize_exterior,ResampledSourceCurb,_profiles


def fixture():
    a=Approach(Path(Section(1,1,((10.,0.,0.),(30.,0.,0.)),width=4.)),0,Crown(.08,0.))
    b=Approach(Path(Section(2,2,((0.,0.,0.),(-20.,0.,0.)),width=4.)),0,Crown(.08,0.))
    curves=(PreparedCapCurb(np.array(((10.,2.),(5.,2.0000005),(0.,2.)))),
            PreparedCapCurb(np.array(((0.,-2.),(5.,-2.0000005),(10.,-2.)))))
    layout=JunctionLayout((a,b),(),curves,(0.,0.),{},10,2,np.zeros(2),[0.,0.],mesh_strategy='boundary_chart')
    layout.dynamic_proof=SimpleNamespace(tokens=('original',1));layout.source_cycle_terminal_ports=[dict(bundle=7)]
    source=sample_footprint(Key('intersection',10),layout)
    return CompoundGroup((Key('intersection',10),),(),(source,),RegionUnion().build((source.region,)),layout)


class SourceRimNormalizationTests(unittest.TestCase):
    def test_profiles_are_byte_exact_and_all_dynamic_proofs_and_original_operands_survive(self):
        before=fixture();after=normalize(before)
        self.assertIs(after.physical_layout.dynamic_proof,before.physical_layout.dynamic_proof)
        self.assertIs(after.physical_layout.source_cycle_terminal_ports,before.physical_layout.source_cycle_terminal_ports)
        self.assertIs(after.footprints[1],before.footprints[0])
        self.assertIs(normalize(after),after)
        self.assertLess(len(after.loops[0].xy),len(before.loops[0].xy))
        def mouth_points(g):
            inventory={s.region.owner:s for s in g.footprints};loop=g.loops[0]
            return np.asarray([loop.xy[i] for i,s in enumerate(loop.spans) if rim_label(s,inventory).kind=='mouth'])
        np.testing.assert_array_equal(mouth_points(before),mouth_points(after))
        for curve,report in zip(after.physical_layout.boundaries,after.physical_layout.source_rim_normalization):
            self.assertLessEqual(report['maximum_projection_error'],1e-6)
            self.assertGreater(len(curve.original_polyline),len(curve.points))
            self.assertTrue(curve.source_curves)
        with self.assertRaises(ValueError):normalize(after,maximum_projection=1e-7)

    def test_projection_budget_and_sharp_original_corners_are_preserved(self):
        points=np.array(((0.,0.),(1.,0.),(1.,1.),(2.,1.)))
        curve=ResampledSourceCurb(points,(),(),1e-6)
        np.testing.assert_array_equal(curve.points,points)
        for invalid in (0.,float('nan'),.001):
            with self.assertRaises(ValueError):ResampledSourceCurb(points,(),(),invalid)
            with self.assertRaises(ValueError):normalize(fixture(),maximum_projection=invalid)

    def test_exterior_carrier_never_becomes_hole_filling_area(self):
        from test_python_street_source_cycles import triangle,prepare
        from tested_python_streets.regions import _classify
        _,before=prepare(triangle(2.))
        original_union=RegionUnion().build(source.region for source in before.area_footprints)
        profiles=_profiles(before)
        after=normalize_exterior(before)
        self.assertIs(after.loops[1],before.loops[1])
        np.testing.assert_array_equal(after.loops[1].xy,before.loops[1].xy)
        self.assertTrue(after.footprints[0].region.boundary_only)
        self.assertIs(normalize_exterior(after),after)
        legacy=normalize(replace(before,loops=(before.loops[0],)))
        with self.assertRaisesRegex(ValueError,'boundary-only area provenance'):
            normalize_exterior(replace(legacy,loops=(*legacy.loops,*before.loops[1:])))
        with self.assertRaisesRegex(ValueError,'boundary-only'):
            RegionUnion().build(source.region for source in after.footprints)
        current_union=RegionUnion().build(source.region for source in after.area_footprints)
        self.assertEqual(len(original_union),len(current_union))
        for a,b in zip(original_union,current_union):
            np.testing.assert_array_equal(a.xy,b.xy)
            self.assertEqual(a.spans,b.spans)
        probe=before.loops[1].xy.mean(axis=0)[None,:]
        # A source-area query is the original operands, not the carrier disk.
        self.assertFalse(any(_classify(probe,source.region.xy,1e-8)[0][0]
                             for source in after.area_footprints))
        self.assertTrue(_classify(probe,after.footprints[0].region.xy,1e-8)[0][0])
        current=_profiles(after)
        self.assertEqual(set(profiles),set(current))
        for key in profiles:
            for a,b in zip(profiles[key],current[key]):np.testing.assert_array_equal(a,b)


if __name__=='__main__':unittest.main()
