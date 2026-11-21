"""An actual four-arm outer root keeps its original circles and full sources."""
from dataclasses import replace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Settings,Key
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.cyclic_root_support import prepare


def fixture(reverse=False,angle=0.,translation=(0.,0.)):
    controls=((838,0,7.,1,((-1066.207763671875,-759.0186767578125),(-1051.8868408203125,-768.6168212890625))),
      (839,0,2.,0,((-1066.207763671875,-759.0186767578125),(-1059.2529296875,-760.2109985351562),(-1041.7257080078125,-764.57763671875),(-1040.71435546875,-764.8338623046875),(-1008.3982543945312,-772.8764038085938))),
      (298,1,3.5,0,((-1074.3209228515625,-751.8260498046875),(-1069.427978515625,-755.2666625976562),(-1066.207763671875,-759.0186767578125))),
      (842,1,3.5,0,((-1085.3521728515625,-749.2081298828125),(-1066.207763671875,-759.0186767578125))))
    rot=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    arms={}
    for sid,end,width,half,points in controls:
        xy=np.array(points)@rot.T+translation
        if reverse:xy=xy[::-1];end=1-end
        source=Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=width,rise=.07999999821186066)
        path=Route((Piece(Path(source),end==0),))
        arms[sid]=Approach(path,end,Crown(source.rise,0.,half if end==0 else -half))
    ends=tuple(End(a.source.id,a.end) for a in arms.values())
    site=Intersection(1683,ends,radii=tuple(1. if e.section_id==839 else 2. for e in ends),radius_modes=('STYLE',)*4)
    donor=CornerRadiusPolicy.solve(arms[839],arms[298],1.,'STYLE')
    support=(CornerRadiusPolicy.solve(arms[842],arms[838],2.,'STYLE'),
             CornerRadiusPolicy.solve(arms[838],arms[839],2.,'STYLE'))
    return site,arms,donor,support


class CyclicRootSupportTests(unittest.TestCase):
    def test_actual_four_arm_curved_sources_reverse_rotate_and_keep_complete_root_curve(self):
        for reverse,angle,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            site,arms,donor,support=fixture(reverse,angle,translation)
            before={sid:(a.source,a.profile) for sid,a in arms.items()}
            circles=[(c.center.copy(),c.stations,c.radius) for c in support]
            record=prepare(Key('patch',1),donor,Settings(curve_tolerance=.00025),site=site,
                original_approaches=tuple(arms.values()),support_corners=support,required_roots=((arms[838],-1),))
            self.assertEqual(record.report['original_degree'],4)
            self.assertEqual(record.report['preferred_minimum_radius'],1.)
            self.assertLessEqual(record.report['curvature_bound'],1.+1e-8)
            self.assertEqual(record.report['maximum_outward_envelope_m'],.5)
            self.assertLessEqual(record.report['outward_envelope_m'],.5)
            self.assertLessEqual(record.report['reach_steps'],16)
            self.assertTrue(record.report['complete_new_curb_exposed'])
            self.assertTrue(record.report['all_original_root_caps_geometrically_internal'])
            for c,(center,stations,radius) in zip(support,circles):
                np.testing.assert_array_equal(c.center,center)
                self.assertEqual((c.stations,c.radius),(stations,radius))
            for index,a in enumerate((donor.first,donor.second)):
                side=-1. if index==0 else 1.
                np.testing.assert_allclose(record.curve(float(index)),a.point(record.stations[index],side*.5*a.width),atol=1e-10)
                self.assertLess(record.stations[index],a.length)
            self.assertEqual(before,{sid:(a.source,a.profile) for sid,a in arms.items()})
            self.assertEqual(float(arms[838].profile.value(1.)),.07999999821186066)

    def test_authored_strict_plateau_and_unrelated_incidence_are_not_root_support(self):
        site,arms,donor,support=fixture();settings=Settings(curve_tolerance=.00025)
        arguments=dict(site=site,original_approaches=tuple(arms.values()),support_corners=support,
                       required_roots=((arms[838],-1),))
        with self.assertRaisesRegex(ValueError,'STYLE'):
            prepare(Key('patch',1),replace(donor,radius_mode='STRICT'),settings,**arguments)
        with self.assertRaisesRegex(ValueError,'zero-height outer'):
            prepare(Key('patch',1),donor,settings,**dict(arguments,required_roots=((arms[838],1),)))
        changed=replace(site,ends=(*site.ends[:-1],End(999,0)))
        with self.assertRaisesRegex(ValueError,'declared original incidence'):
            prepare(Key('patch',1),donor,settings,**dict(arguments,site=changed))
        with self.assertRaisesRegex(ValueError,'distinct actual local'):
            prepare(Key('patch',1),donor,settings,**dict(arguments,support_corners=(*support,support[0])))


if __name__=='__main__':unittest.main()
