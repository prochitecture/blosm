"""Exact region classification/provenance parity after an AABB rejection."""
import unittest
import numpy as np

import _python_street_package
from tested_python_streets.model import Key
from tested_python_streets.regions import Region,RegionUnion,_classify


class FullClassificationUnion(RegionUnion):
    def _classification_candidates(self,points,bounds):
        return np.ones(len(points),bool)


class RegionClassificationBoundsTests(unittest.TestCase):
    def compare(self,regions):
        optimized=RegionUnion().build(regions)
        reference=FullClassificationUnion().build(regions)
        self.assertEqual(len(optimized),len(reference))
        for a,b in zip(optimized,reference):
            np.testing.assert_array_equal(a.xy,b.xy)
            self.assertEqual(a.spans,b.spans)

    def test_partial_coincidence_holes_disjoint_regions_and_owner_ties_match_exactly(self):
        coordinates=[((0.,0.),(8.,0.),(8.,2.),(0.,2.)),
            ((0.,6.),(8.,6.),(8.,8.),(0.,8.)),
            ((0.,0.),(2.,0.),(2.,8.),(0.,8.)),
            ((6.,0.),(8.,0.),(8.,8.),(6.,8.)),
            ((0.,0.),(4.,0.),(8.,0.),(8.,2.),(0.,2.)),
            ((30.,0.),(34.,0.),(34.,4.),(30.,4.))]
        for angle in (0.,.37):
            rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
            regions=tuple(Region(Key('test',i),np.asarray(p)@rotation.T+(10.,-30.)) for i,p in enumerate(coordinates))
            self.compare(regions)
            self.compare(regions[::-1])

    def test_epsilon_coincident_boundaries_keep_the_original_tie(self):
        for delta in (0.,.25e-7,-.25e-7):
            first=Region(Key('test',0),np.array(((0.,0.),(8.,0.),(8.,4.),(0.,4.))))
            second=Region(Key('test',1),np.array(((4.,delta),(12.,delta),(12.,4.),(4.,4.))))
            self.compare((second,first))

    def test_broad_phase_overselects_narrow_arithmetic_at_large_coordinate_range(self):
        polygon=np.array(((-1e12,-1.),(1000.,0.),(-1e12,1.)))
        points=np.array(((1000.+2e-7,0.),))
        inside,boundary=_classify(points,polygon,1e-7)
        self.assertTrue(bool(inside[0] or boundary[0]>=0))
        self.assertTrue(RegionUnion()._classification_candidates(points,
            (polygon.min(axis=0),polygon.max(axis=0)))[0])


if __name__=='__main__':unittest.main()
