"""An actual wider third root uses finite, globally exposed source curbs."""
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
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.regions import Region
from tested_python_streets.continuations import ContinuationReachError
from tested_python_streets.cyclic_root_support import prepare


def fixture(reverse=False,angle=0.,translation=(0.,0.)):
    controls={
      635:(0,2.,((1468.292724609375,-89.6155776977539),(1483.294677734375,-78.32487487792969))),
      636:(0,2.,((1468.292724609375,-89.6155776977539),(1545.934814453125,-89.03250885009766),(1567.3800048828125,-99.75933074951172))),
      1025:(1,2.,((1448.73974609375,-104.31349182128906),(1468.292724609375,-89.6155776977539))),
      1038:(1,7.,((1469.2406005859375,-133.72018432617188),(1468.611572265625,-104.68817901611328),(1468.292724609375,-89.6155776977539))),
      595:(1,5.,((1180.124755859375,-319.2843322753906),(1303.6412353515625,-194.608642578125),
          (1346.2811279296875,-151.57626342773438),(1375.0208740234375,-124.932373046875),
          (1405.2393798828125,-106.49235534667969),(1429.59228515625,-95.08871459960938),
          (1457.1048583984375,-84.5304183959961),(1483.294677734375,-78.32487487792969)))}
    matrix=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    arms={};paths={}
    for sid,(end,width,points) in controls.items():
        xy=(np.asarray(points)-np.array((1468.,-89.)))@matrix.T+translation
        if reverse:xy=xy[::-1];end=1-end
        source=Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=width,interpolation=3)
        path=Path(source);paths[sid]=path
        profile=Crown(.08,0.,(-1 if not reverse else 1) if sid==595 else 0)
        arms[sid]=Approach(Route((Piece(path,end==0),)),end,profile)
    originals=tuple(arms[sid] for sid in (636,635,1025,1038))
    ends=tuple(End(a.source.id,a.end) for a in originals)
    site=Intersection(2075,ends,radii=(1.,1.,1.,2.),radius_modes=('STYLE',)*4)
    donor=CornerRadiusPolicy.solve(arms[635],arms[1025],1.,'STYLE')
    support=(CornerRadiusPolicy.solve(arms[636],arms[635],1.,'STYLE'),
             CornerRadiusPolicy.solve(arms[1025],arms[1038],1.,'STYLE'),
             CornerRadiusPolicy.solve(arms[1038],arms[636],2.,'STYLE'))
    settings=Settings(curve_tolerance=.00025)
    far635=Approach(Route((Piece(paths[635],reverse),)),int(not reverse),Crown(.08,0.))
    remote=CornerRadiusPolicy.solve(arms[595],far635,1.92,'STYLE')
    remote_source=corner_patch(Key('other_original_corner',1),remote,settings.spacing,settings.curve_tolerance)
    at=paths[635].samples(0.,paths[635].length,.5,tolerance=.00025)
    strip=np.vstack((paths[635].point(at,1.),paths[635].point(at,-1.)[::-1]))
    from tested_python_streets.regions import signed_area
    if signed_area(strip)<0.:strip=strip[::-1]
    other=(remote_source.region,Region(Key('complete_original_strip',635),strip))
    return site,arms,donor,support,originals,other,settings


class CyclicProjectionRecoveryTests(unittest.TestCase):
    def test_original_wide_root_and_remote_fillet_require_one_finite_asymmetric_curb(self):
        for reverse,angle,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            with self.subTest(reverse=reverse):
                site,arms,donor,support,originals,other,settings=fixture(reverse,angle,translation)
                args=dict(site=site,original_approaches=originals,support_corners=support,
                    required_roots=((arms[1038],1),))
                with self.assertRaises(ContinuationReachError):prepare(Key('curb',1),donor,settings,**args)
                unchanged=[(c.radius,c.stations,c.center.copy()) for c in support]
                record=prepare(Key('curb',1),donor,settings,**args,source_regions=other,consumed={635})
                proof=record.report['cyclic_projection_recovery']
                self.assertEqual(record.report['preferred_minimum_radius'],1.)
                self.assertLessEqual(record.report['curvature_bound'],1.+1e-8)
                self.assertEqual(proof['maximum_outward_bound'],2.5)
                self.assertLessEqual(record.report['outward_envelope_m'],2.5)
                self.assertGreater(record.report['outward_envelope_m'],.5)
                self.assertTrue(proof['global_whole_curve_exposure'])
                self.assertTrue(record.report['all_original_root_caps_geometrically_internal'])
                self.assertTrue(any(t.get('qualified') and t['held_end']==0 for t in proof['trials']))
                self.assertLessEqual(proof['bounded_families_per_radius'],6)
                for c,(radius,stations,center) in zip(support,unchanged):
                    self.assertEqual((c.radius,c.stations),(radius,stations));np.testing.assert_array_equal(c.center,center)

    def test_no_unproven_held_source_or_strict_donor_can_enter_the_new_family(self):
        site,arms,donor,support,originals,other,settings=fixture()
        args=dict(site=site,original_approaches=originals,support_corners=support,
            required_roots=((arms[1038],1),),source_regions=other)
        with self.assertRaises(ContinuationReachError):prepare(Key('curb',1),donor,settings,**args)
        with self.assertRaisesRegex(ValueError,'STYLE'):
            prepare(Key('curb',1),replace(donor,radius_mode='STRICT'),settings,**args,consumed={635})


if __name__=='__main__':unittest.main()
