"""An already-permitted style fit must preserve original finite roots."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section
from tested_python_streets.paths import Path
from tested_python_streets.intersections import Approach
from tested_python_streets.profiles import Crown
from tested_python_streets.corner_radii import CornerRadiusPolicy,CornerContactError
from tested_python_streets.corner_source_support import complete_source_roots,raise_source_root_refusal


class CornerSourceSupportTests(unittest.TestCase):
    def test_optional_refusal_reports_the_actual_finite_facing_root_instead_of_a_reach_failure(self):
        import math
        controls=(((0.,0.,0.),(-60.,0.,0.)),
            ((0.,0.,0.),(-60.*math.cos(.4),-60.*math.sin(.4),0.)),
            ((0.,0.,0.),(60.,0.,0.)))
        sources=tuple(Section(i,i,points,width=width,interpolation=0)
            for i,(points,width) in enumerate(zip(controls,(10.5,7.,7.)),1))
        arms=tuple(Approach(Path(source),0,Crown(.08,0.)) for source in sources)
        corner=CornerRadiusPolicy.solve(*arms[:2],2.,'STRICT')
        self.assertTrue(all(0.<station<arm.length for arm,station in zip(arms,corner.stations)))
        records=[];before=tuple(source.points for source in sources)
        result=complete_source_roots(corner,.7,.00025,support_approaches=arms,
            audit_all_contacts=True,audit_facing_halves=True,on_refusal=records.append)
        self.assertIs(result,False)
        record,=records
        self.assertEqual((record['section'],record['end'],record['side']),(1,0,-1))
        self.assertEqual(record['reason'],'EXPOSED_ORIGINAL_FACING_ROOT_HALF')
        self.assertGreater(record['exposed_length'],1.)
        self.assertEqual(record['interval'][0],0.);self.assertLess(record['interval'][1],.5)
        with self.assertRaisesRegex(CornerContactError,'Original Source 1/end0.*LEFT root interval'):
            complete_source_roots(corner,.7,.00025,support_approaches=arms,
                audit_all_contacts=True,audit_facing_halves=True,on_refusal=raise_source_root_refusal)
        self.assertEqual(tuple(source.points for source in sources),before)

    def test_optional_refusal_diagnostics_leave_a_successful_candidate_unchanged(self):
        from test_python_street_corner_patch import approaches
        first,last=approaches('positive')
        corner=CornerRadiusPolicy.solve(first,last,4.,'STRICT')
        points=corner.center.copy();stations=corner.stations;records=[]
        plain=complete_source_roots(corner,.5,.00025,audit_all_contacts=True,audit_facing_halves=True)
        observed=complete_source_roots(corner,.5,.00025,audit_all_contacts=True,audit_facing_halves=True,
            on_refusal=records.append)
        self.assertIs(plain,True);self.assertIs(observed,True);self.assertEqual(records,[])
        np.testing.assert_array_equal(corner.center,points);self.assertEqual(corner.stations,stations)

    def test_forced_audit_distinguishes_facing_sector_from_another_mixed_corner_sector(self):
        from test_python_street_corner_patch import approaches as actual_approaches
        for case in ('negative','mixed'):
            for routed in (False,True):
                first,last=actual_approaches(case,routed=routed)
                corner=CornerRadiusPolicy.solve(first,last,4.,'STRICT')
                self.assertLess(min(corner.stations),0.)
                before=(first.source.points,last.source.points)
                self.assertEqual(complete_source_roots(corner,.5,.00025,
                    support_approaches=(first,last),audit_all_contacts=True),case=='negative')
                self.assertTrue(complete_source_roots(corner,.5,.00025,
                    support_approaches=(first,last),audit_all_contacts=True,audit_facing_halves=True))
                self.assertEqual((first.source.points,last.source.points),before)

    def test_actual1198_keeps_original_radius_when_its_third_same_node_road_covers_the_roots(self):
        controls=(((1642.1929931640625,986.462646484375),(1658.0205078125,956.9894409179688),
                   (1664.4317626953125,947.7738647460938)),
                  ((1664.4317626953125,947.7738647460938),(1671.4512939453125,932.781005859375)),
                  ((1664.4317626953125,947.7738647460938),(1663.8153076171875,955.5772094726562)))
        sources=tuple(Section(i,i,tuple((x,y,0.) for x,y in points),width=5.) for i,points in enumerate(controls))
        approaches=tuple(Approach(Path(source),end,Crown(.08,0.)) for source,end in zip(sources,(1,0,0)))
        first,last,support=approaches
        preferred=CornerRadiusPolicy.solve(first,last,2.,'STRICT')
        self.assertFalse(complete_source_roots(preferred,.7,.00025))
        self.assertTrue(complete_source_roots(preferred,.7,.00025,support_approaches=approaches))
        accepted=CornerRadiusPolicy.solve(first,last,2.,'STRICT',
            accept=lambda corner:complete_source_roots(corner,.7,.00025,support_approaches=approaches))
        self.assertEqual(accepted.radius,2.)
        elsewhere=Approach(Path(Section(10,10,((0.,0.,0.),(0.,20.,0.)),width=5.)),0,Crown(.08,0.))
        with self.assertRaisesRegex(ValueError,'same original physical junction'):
            complete_source_roots(preferred,.7,.00025,support_approaches=(elsewhere,))

    def test_actual1211_fits_within_style_range_and_keeps_explicit_radius_strict(self):
        controls=(((956.9209594726562,-477.24456787109375),(903.8540649414062,-485.6000061035156),
                   (882.2139892578125,-489.27606201171875),(876.3306274414062,-496.05609130859375)),
                  ((1106.570556640625,-375.6564025878906),(1069.6148681640625,-412.8089904785156),
                   (1044.6640625,-434.241455078125),(1023.8369750976562,-449.0052490234375),
                   (1001.3356323242188,-461.5205383300781),(981.447021484375,-470.0279541015625),
                   (956.9209594726562,-477.24456787109375)))
        for angle,offset in ((0.,(0.,0.)),(.713,(1200.,-1800.))):
            rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
            sources=tuple(Section(i,i,tuple((x,y,0.) for x,y in np.array(points)@rotation.T+offset),width=5.)
                          for i,points in enumerate(controls))
            approaches=tuple(Approach(Path(source),end,Crown(.08,0.)) for source,end in zip(sources,(0,1)))
            before=tuple(source.points for source in sources)
            preferred=CornerRadiusPolicy.solve(*approaches,3.,'STRICT')
            self.assertGreater(min(preferred.stations),0.)
            self.assertFalse(complete_source_roots(preferred,.7,.00025))
            self.assertFalse(complete_source_roots(preferred,.7,.00025,
                audit_all_contacts=True,audit_facing_halves=True))
            from types import SimpleNamespace
            from tested_python_streets.model import Network
            from tested_python_streets.compound_corner_contacts import retry
            context=SimpleNamespace(network=Network('facing-exposure-regression',sources,()))
            retried=retry(context,*approaches,3.,'STYLE',set(),set(),approaches)
            self.assertIsNotNone(retried)
            self.assertNotEqual(retried.radius,3.)
            accept=lambda corner:complete_source_roots(corner,.7,.00025)
            fitted=CornerRadiusPolicy.solve(*approaches,3.,'STYLE',accept=accept)
            self.assertAlmostEqual(fitted.radius,2.4)
            self.assertLess(max(fitted.stations),0.)
            self.assertEqual(fitted.radius_attempts,(3.,2.4000000000000004))
            self.assertEqual(tuple(source.points for source in sources),before)
            with self.assertRaises(CornerContactError):
                CornerRadiusPolicy.solve(*approaches,3.,'STRICT',accept=accept)
            with self.assertRaisesRegex(ValueError,'requires forced'):
                complete_source_roots(preferred,.7,.00025,audit_facing_halves=True)


if __name__=='__main__':unittest.main()
