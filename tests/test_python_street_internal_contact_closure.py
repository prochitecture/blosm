"""A redundant internal contact lobe cannot remove or invent physical curbs."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section
from tested_python_streets.paths import Path
from tested_python_streets.intersections import Approach
from tested_python_streets.profiles import Crown
from tested_python_streets.compound_displaced_contacts import SourceContactWalk
from tested_python_streets.internal_contact_closure import normalize


class InternalContactClosureTests(unittest.TestCase):
    @staticmethod
    def fixture(width=10.):
        path=Path(Section(1,1,((-10.,0.,0.),(10.,0.,0.)),width=width,interpolation=0))
        approach=Approach(path,0,Crown(.08,0.))
        corner=SimpleNamespace(source_contact_walks=(SourceContactWalk(approach,approach),),stations=(20.,))
        return np.array(((0.,0.),(2.,2.),(0.,2.),(2.,0.))),['internal']*3+['curb'],corner

    def test_only_the_wholly_source_contained_internal_lobe_is_removed(self):
        xy,kinds,corner=self.fixture()
        before=xy.copy()
        result,labels,proof=normalize(xy,kinds,corner,.3,.00025)
        self.assertTrue(proof['continuous_internal_containment'])
        self.assertEqual(proof['section'],1)
        self.assertAlmostEqual(proof['discarded_area'],1.)
        actual=[(result[i],result[(i+1)%len(result)]) for i,kind in enumerate(labels) if kind=='curb']
        self.assertEqual(len(actual),1)
        self.assertEqual({tuple(point) for point in actual[0]},
                         {tuple(point) for point in xy[[3,0]]})
        np.testing.assert_array_equal(xy,before)

    def test_a_curb_lobe_or_uncovered_internal_lobe_still_rejects(self):
        xy,kinds,corner=self.fixture()
        with self.assertRaisesRegex(ValueError,'physical source curb'):
            normalize(xy,['internal','curb','internal','curb'],corner,.3,.00025)
        xy,kinds,corner=self.fixture(width=.2)
        with self.assertRaisesRegex(ValueError,'not covered by an actual finite source strip'):
            normalize(xy,kinds,corner,.3,.00025)

    def test_pair_allocation_and_ambiguous_non_crossing_closures_remain_failures(self):
        xy,kinds,corner=self.fixture()
        with self.assertRaisesRegex(ValueError,'finite pair budget'):
            normalize(xy,kinds,corner,.3,.00025,maximum_pairs=0)
        with self.assertRaisesRegex(ValueError,'single proper'):
            normalize(xy[[0,3,1,2]],kinds,corner,.3,.00025)


if __name__=='__main__':unittest.main()
