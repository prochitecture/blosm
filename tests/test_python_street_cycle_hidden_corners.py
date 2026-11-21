"""An unrelated median hole cannot hide a failed original internal fillet."""
from types import SimpleNamespace
import unittest
import numpy as np

import _python_street_package
from tested_python_streets.compound_cycle_hidden import (certify_hidden_corners,
    certify_partially_hidden_corners, partial_corner_candidate)
from tested_python_streets.intersections import Mouth
from tested_python_streets.model import Section,End
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.regions import BoundaryLoop


def fixture(reverse=False):
    paths=[Path(Section(1,1,((0.,0.,0.),(10.,0.,0.)),width=2.,interpolation=0)),
           Path(Section(2,2,((0.,0.,0.),(0.,10.,0.)),width=2.,interpolation=0))]
    if reverse:
        paths=[Path(Section(p.section.id,p.section.street_id,p.section.points[::-1],
                            width=2.,interpolation=0)) for p in paths]
    ends=(End(1,int(reverse)),End(2,int(reverse)))
    approaches={e:SimpleNamespace(path=Route((Piece(p,not reverse),))) for e,p in zip(ends,paths)}
    partitions={p.section.id:np.array((0.,5.,10.)) for p in paths}
    outer=BoundaryLoop(np.array(((-5.,-5.),(20.,-5.),(20.,20.),(-5.,20.))),())
    return ends,approaches,partitions,outer


class HiddenCycleCornerTests(unittest.TestCase):
    def test_partial_hidden_corner_uses_only_its_actual_owned_external_prefix(self):
        for reverse in (False,True):
            ends,approaches,partitions,outer=fixture(reverse)
            partitions[2]=np.array((0.,3.,5.,7.,10.))
            mouth = 7. if reverse else 3.
            # The complete external road goes outside this owner. Only its
            # real mouth-bounded prefix belongs to the compound.
            outer=BoundaryLoop(np.array(((-5.,-5.),(20.,-5.),(20.,4.),(-5.,4.))),())
            records=certify_partially_hidden_corners(((*ends,'failed STYLE circle'),),approaches,
                partitions,(outer,),{1},{(2,int(reverse)):Mouth(mouth,8)})
            self.assertEqual(records[0]['exposed_physical_boundary_intervals'],0)
            self.assertTrue(records[0]['external_remainder_unchanged'])
            row=next(row for row in records[0]['source_curbs'] if row['section']==2)
            self.assertEqual(row['original_interval'],[7.,10.] if reverse else [0.,3.])
            self.assertEqual(row['scope'],'actual_owned_source_prefix')
            self.assertFalse(partial_corner_candidate(*ends,tuple(approaches[e] for e in ends),set(),{1}))
            self.assertFalse(partial_corner_candidate(*ends,tuple(approaches[e] for e in ends),set(ends),{1,2}))

    def test_partial_hidden_curb_must_not_coincide_with_outer_boundary_or_a_hole(self):
        ends,approaches,partitions,outer=fixture()
        partitions[2]=np.array((0.,3.,5.,10.));mouths={(2,0):Mouth(3.,8)}
        # All tested curbs are contained, but a positive part of Source1's
        # facing curb remains exposed on this actual exterior.
        exterior=BoundaryLoop(np.array(((-5.,-5.),(20.,-5.),(20.,1.),(2.,1.),(2.,4.),(-5.,4.))),())
        with self.assertRaisesRegex(ValueError,'physical boundary retains'):
            certify_partially_hidden_corners(((*ends,'local failure'),),approaches,
                partitions,(exterior,),{1},mouths)
        hole=BoundaryLoop(np.array(((4.,2.),(6.,2.),(6.,1.),(4.,1.))),())
        with self.assertRaisesRegex(ValueError,'physical boundary retains'):
            certify_partially_hidden_corners(((*ends,'local failure'),),approaches,
                partitions,(outer,hole),{1},mouths)
        tiny=BoundaryLoop(np.array(((3.1,1.2),(3.2,1.2),(3.2,.8),(3.1,.8))),())
        with self.assertRaisesRegex(ValueError,'continuous owned facing-curb'):
            certify_partially_hidden_corners(((*ends,'local failure'),),approaches,
                partitions,(outer,tiny),{1},mouths)

    def test_partial_hidden_corner_cannot_skip_complete_consumed_curbs_or_invent_a_mouth(self):
        ends,approaches,partitions,outer=fixture()
        with self.assertRaisesRegex(ValueError,'actual positive owned source prefix'):
            certify_partially_hidden_corners(((*ends,'local failure'),),approaches,
                partitions,(outer,),{1},{})
        clipped=BoundaryLoop(np.array(((-5.,-5.),(9.,-5.),(9.,4.),(-5.,4.))),())
        with self.assertRaisesRegex(ValueError,'continuous owned facing-curb'):
            certify_partially_hidden_corners(((*ends,'local failure'),),approaches,
                partitions,(clipped,),{1},{(2,0):Mouth(3.,8)})

    def test_remote_hole_preserves_complete_original_facing_curbs_in_both_directions(self):
        for reverse in (False,True):
            ends,approaches,partitions,outer=fixture(reverse)
            hole=BoundaryLoop(np.array(((15.,15.),(15.,18.),(18.,18.),(18.,15.))),())
            record,=certify_hidden_corners(((*ends,'actual local solve failure'),),approaches,
                                           partitions,(outer,hole),{1,2})
            self.assertTrue(record['continuous_pavement_closure'])
            self.assertEqual(record['exposed_hole_intervals'],0)
            self.assertEqual([c['original_side'] for c in record['source_curbs']],
                             [1,-1] if reverse else [-1,1])

    def test_hole_facing_curb_exposure_and_interior_crossing_remain_errors(self):
        for top,bottom in ((2.,1.),(1.5,.5)):
            ends,approaches,partitions,outer=fixture()
            hole=BoundaryLoop(np.array(((4.,top),(6.,top),(6.,bottom),(4.,bottom))),())
            with self.assertRaisesRegex(ValueError,'source curb|facing-curb coverage'):
                certify_hidden_corners(((*ends,'local failure'),),approaches,partitions,(outer,hole),{1,2})

    def test_missing_consumed_source_and_outside_strip_cannot_be_certified(self):
        ends,approaches,partitions,outer=fixture()
        with self.assertRaisesRegex(ValueError,'complete original source'):
            certify_hidden_corners(((*ends,'local failure'),),approaches,partitions,(outer,),{1})
        clipped=BoundaryLoop(np.array(((-5.,-5.),(4.,-5.),(4.,20.),(-5.,20.))),())
        with self.assertRaisesRegex(ValueError,'facing-curb coverage'):
            certify_hidden_corners(((*ends,'local failure'),),approaches,partitions,(clipped,),{1,2})


if __name__=='__main__':unittest.main()
