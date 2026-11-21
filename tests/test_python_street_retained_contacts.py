"""Only a failed declared original circular contact enters source-union recovery."""
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Settings,Key
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerContactError
from tested_python_streets.compound_curbs import (_local_corners,RetainedCornerContactError,
    _composite,RetainedSourceRouteError)


class RetainedContactTests(unittest.TestCase):
    def test_actual_internal_1649_failure_keeps_original_incidence_radius_and_mode(self):
        controls={821:(((-573.083740234375,-752.0703735351562),(-572.3983764648438,-744.0220336914062),(-573.2526245117188,-727.6023559570312)),0),
                  825:(((-578.0994873046875,-728.1251220703125),(-575.5103759765625,-744.856689453125),(-573.083740234375,-752.0703735351562)),1)}
        for reverse,angle in ((False,0.),(True,.713)):
            rot=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
            approaches=[];ends=[]
            for sid,(xy,end) in controls.items():
                points=np.asarray(xy)@rot.T+(1200.,-1800.)
                if reverse:points=points[::-1]
                actual_end=1-end if reverse else end
                source=Section(sid,sid,tuple((x,y,0.) for x,y in points),width=7.)
                approaches.append(Approach(Path(source),actual_end,Crown(.08,0.)))
                ends.append(End(sid,actual_end))
            for mode in ('STYLE','STRICT'):
                site=Intersection(1649,tuple(ends),2.,radius_modes=(mode,mode))
                planner=SimpleNamespace(sites={site.key:site},layouts={site.key:SimpleNamespace(approaches=approaches)},
                    context=SimpleNamespace(network=SimpleNamespace(settings=Settings())))
                before=[a.source.points for a in approaches]
                with self.subTest(reverse=reverse,mode=mode),self.assertRaises(RetainedCornerContactError) as caught:
                    _local_corners(planner,(tuple(ends),))
                error=caught.exception
                self.assertIsInstance(error.__cause__,CornerContactError)
                self.assertEqual((error.junction,error.first,error.last),(site.key,*ends))
                self.assertEqual((error.radius,error.mode),(2.,mode))
                self.assertEqual(before,[a.source.points for a in approaches])

    def test_undeclared_cross_site_pair_is_not_a_recoverable_local_contact(self):
        a=Approach(Path(Section(1,1,((0.,0.,0.),(20.,0.,0.)),width=7.)),0,Crown(.08,0.))
        b=Approach(Path(Section(2,2,((0.,0.,0.),(0.,20.,0.)),width=7.)),0,Crown(.08,0.))
        planner=SimpleNamespace(layouts={Key('intersection',1):SimpleNamespace(approaches=(a,)),
                                         Key('intersection',2):SimpleNamespace(approaches=(b,))})
        with self.assertRaisesRegex(ValueError,'one original junction') as caught:
            _local_corners(planner,((End(1,0),End(2,0)),))
        self.assertNotIsInstance(caught.exception,RetainedCornerContactError)

    def test_equal_length_does_not_make_distinct_original_curbs_one_reversed_route(self):
        a=Approach(Path(Section(1,1,((0.,0.,0.),(30.,0.,0.)),width=4.)),0,Crown(.08,0.))
        b=Approach(Path(Section(2,2,((30.,.01,0.),(0.,.01,0.)),width=4.)),0,Crown(.08,0.))
        first=SimpleNamespace(first=a,second=a,stations=(0.,0.),center=None)
        last=SimpleNamespace(first=b,second=b,stations=(0.,0.),center=None)
        original=(a.source.points,b.source.points)
        with self.assertRaisesRegex(RetainedSourceRouteError,'opposite original incidence'):
            _composite(None,(first,last),(0.,0.))
        self.assertEqual((a.source.points,b.source.points),original)


if __name__=='__main__':unittest.main()
