"""A finite STYLE curb retains both original roots without shrinking radius."""
from dataclasses import replace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compound_root_envelopes import prepare, prepare_adjacent, _complete_roots
from tested_python_streets.compound_dimensional_corners import prepare as prepare_dimensional
from tested_python_streets.compound_dimensional_corners import DimensionalLeadExhausted
from tested_python_streets.compound_root_envelopes import prepare_dimensional as prepare_finite_width
from tested_python_streets.compound_cycles import _footprint
from tested_python_streets.compounds import rim_label
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.intersections import Approach, ApproachExtentError
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.regions import RegionUnion
from tested_python_streets.routes import Piece, Route
from tested_python_streets.sampling import CertifiedCurbSampler


def actual274(angle=0.,translation=(0.,0.),reverse=False,routed=False):
    controls=((138,0,((-1555.914306640625,57.77984619140625),(-1563.561279296875,87.13671112060547),
                      (-1575.1964111328125,125.09959411621094))),
              (81,1,((-1323.4530029296875,-236.1239471435547),(-1336.9072265625,-218.6884307861328),
                     (-1352.7115478515625,-205.6829071044922),(-1396.557861328125,-185.19044494628906),
                     (-1409.2049560546875,-177.6623992919922),(-1422.23291015625,-168.33087158203125),
                     (-1439.2734375,-153.57711791992188),(-1448.4193115234375,-144.02378845214844),
                     (-1456.063720703125,-133.3909912109375),(-1542.7215576171875,25.88353157043457),
                     (-1555.914306640625,57.77984619140625))))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    result=[]
    for sid,end,points in controls:
        xy=np.asarray(points)@rotation.T+translation
        if reverse:xy=xy[::-1];end=1-end
        source=model.Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=3.5)
        path=Path(source)
        if routed:path=Route((Piece(path,end==0),))
        result.append(Approach(path,end,Crown(.08,0.)))
    return tuple(result)


def exposed_construction(corner,source,settings):
    sources=[source]
    for i,a in enumerate((corner.first,corner.second)):
        at=np.unique(np.r_[np.linspace(0.,8.,161),max(0.,corner.stations[i])])
        points=np.vstack((a.point(at,.5*a.width),a.point(at,-.5*a.width)[::-1]))
        sources.append(_footprint(model.Key('root_strip_proof',i),points,internal_edges=(2*len(at)-1,)))
    loops=RegionUnion().build(s.region for s in sources)
    inventory={s.region.owner:s for s in sources}
    return sum(float(np.linalg.norm(loop.xy[(i+1)%len(loop.xy)]-loop.xy[i]))
               for loop in loops for i,span in enumerate(loop.spans) if rim_label(span,inventory).kind=='internal')


def actual2799(angle=0., translation=(0., 0.), reverse=False):
    controls=((1339,0,3.5,((-695.6260375976562,-1054.3682861328125),
        (-696.6661376953125,-1032.29345703125),(-696.066650390625,-1028.60888671875),
        (-694.6581420898438,-1025.92626953125),(-691.9501953125,-1022.0637817382812),
        (-689.6469116210938,-1018.7689819335938))),
        (1340,0,7.,((-695.6260375976562,-1054.3682861328125),(-704.042724609375,-1041.365234375),
        (-711.5264282226562,-1034.1285400390625),(-718.4659423828125,-1030.665771484375),
        (-728.9066162109375,-1028.148681640625))),
        (1400,1,7.,((-691.9347534179688,-1093.2303466796875),(-695.6260375976562,-1054.3682861328125))))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    approaches={}
    for sid,end,width,points in controls:
        xy=np.asarray(points)@rotation.T+translation
        if reverse:xy=xy[::-1];end=1-end
        source=model.Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=width,interpolation=3)
        half=(-1 if not reverse else 1) if sid==1340 else 0
        approaches[sid]=Approach(Path(source),end,Crown(.08,0.,half))
    a,b,required=(approaches[sid] for sid in (1400,1339,1340))
    original=CornerRadiusPolicy.solve(a,b,2.,'STYLE',math.radians(15.))
    proxy=replace(original,center=None,continuation=True,stations=(0.,0.))
    settings=model.Settings(curve_tolerance=.00025)
    dimensional=prepare_dimensional(model.Key('original_dimensional',1),proxy,settings)
    corner=replace(proxy,stations=dimensional.stations);corner.dimensional_corner=dimensional
    support=(CornerRadiusPolicy.solve(b,required,2.,'STYLE'),CornerRadiusPolicy.solve(required,a,2.,'STYLE'))
    return corner,required,support,settings


class RootEnvelopeTests(unittest.TestCase):
    def test_typed_finite_width_root_spends_only_real_leads_without_reducing_radius(self):
        settings=model.Settings(curve_tolerance=.00025)
        for mode,reverse,angle,translation in (('STYLE',False,0.,(0.,0.)),('STRICT',True,.713,(1200.,-1800.))):
            matrix=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
            approaches=[]
            for sid,width,finish in ((1,10.8,(4.796,0.)),(2,14.,(-40.,0.))):
                points=np.array(((0.,0.),finish))@matrix.T+translation
                if reverse:points=points[::-1]
                source=model.Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width)
                path=Route((Piece(Path(source),not reverse),))
                approaches.append(Approach(path,int(reverse),Crown(.08,0.)))
            corner=CornerRadiusPolicy.solve(*approaches,4.,mode)
            original=tuple((a.source.points,a.width,a.profile) for a in approaches)
            with self.assertRaises(DimensionalLeadExhausted) as caught:
                prepare_dimensional(model.Key('finite_width',1),corner,settings)
            record=prepare_finite_width(model.Key('finite_width',1),corner,settings,caught.exception)
            self.assertEqual(record.corner.radius,4.)
            self.assertEqual(record.corner.radius_mode,mode)
            self.assertLessEqual(record.report['curvature_bound'],.25*(1.+1e-8))
            self.assertEqual(record.stations[0],caught.exception.last_finite_stations[0])
            self.assertTrue(_complete_roots(record.corner,record.source,settings,facing_halves=True))
            self.assertTrue(record.report['complete_original_facing_root_halves'])
            self.assertFalse(record.report['complete_original_root_caps'])
            opposite=CornerRadiusPolicy.solve(*approaches[::-1],4.,mode)
            with self.assertRaises(DimensionalLeadExhausted) as other:
                prepare_dimensional(model.Key('finite_width',2),opposite,settings)
            support=prepare_finite_width(model.Key('finite_width',2),opposite,settings,other.exception)
            self.assertTrue(_complete_roots(record.corner,record.source,settings,support_regions=(support.source.region,)))
            self.assertEqual(original,tuple((a.source.points,a.width,a.profile) for a in approaches))
            for i,a in enumerate(approaches):self.assertLess(record.stations[i],a.length)
            with self.assertRaisesRegex(ValueError,'typed narrow-lead'):
                prepare_finite_width(model.Key('finite_width',1),corner,settings,ValueError('unrelated'))

    def test_actual2799_supports_adjacent_half_profile_root_without_changing_its_acute_circles(self):
        for angle,translation,reverse in ((0.,(0.,0.),False),(.713,(1200.,-1800.),True)):
            with self.subTest(reverse=reverse):
                corner,required,support,settings=actual2799(angle,translation,reverse)
                original=corner.dimensional_corner
                sources=tuple(corner_patch(model.Key('fixed_circle',i),c,settings.spacing,settings.curve_tolerance).region
                              for i,c in enumerate(support))
                self.assertFalse(_complete_roots(corner,original.source,settings,
                                                 extra_roots=(required,),support_regions=sources))
                controls=tuple(a.source.points for a in (corner.first,corner.second,required))
                circles=tuple((c.radius,c.stations,c.center.copy()) for c in support)
                record=prepare_adjacent(model.Key('root_support',1),corner,settings,
                                        required_roots=(required,),support_corners=support)
                self.assertTrue(_complete_roots(record.corner,record.source,settings,
                                                extra_roots=(required,),support_regions=sources))
                self.assertEqual(record.stations[0],original.stations[0])
                self.assertGreater(record.stations[1],original.stations[1])
                self.assertEqual(record.report['preferred_minimum_radius'],2.)
                self.assertLessEqual(record.report['curvature_bound'],.5*(1.+1e-8))
                self.assertGreater(record.report['outward_envelope_m'],0.)
                self.assertLessEqual(record.report['outward_envelope_m'],record.report['maximum_outward_envelope_m'])
                self.assertEqual(record.report['required_original_roots'][0]['profile_half'],-1)
                self.assertFalse(hasattr(record.corner,'dimensional_corner'))
                self.assertEqual(controls,tuple(a.source.points for a in (corner.first,corner.second,required)))
                for c,(radius,stations,center) in zip(support,circles):
                    self.assertEqual(c.radius,radius);self.assertEqual(c.stations,stations)
                    np.testing.assert_array_equal(c.center,center)
                for t,a,station,distance,sign in ((0.,corner.first,record.stations[0],-.5*corner.first.width,-1.),
                                                (1.,corner.second,record.stations[1],.5*corner.second.width,1.)):
                    np.testing.assert_allclose(record.curve(t),a.point(station,distance),atol=1e-10)
                    delta=record.curve(min(1.,t+1e-6))-record.curve(max(0.,t-1e-6))
                    expected=sign*a.jet(station,distance)[1]
                    np.testing.assert_allclose(delta/np.linalg.norm(delta),expected/np.linalg.norm(expected),atol=1e-6)
                fine=CertifiedCurbSampler(settings.curve_tolerance,settings.spacing,8192)
                coarse=CertifiedCurbSampler(.02,settings.spacing,8192)
                np.testing.assert_array_equal(fine.sample(record.curve,0.,1.),coarse.sample(record.curve,0.,1.))

    def test_adjacent_root_support_requires_original_local_sources_and_certified_style_transition(self):
        corner,required,support,settings=actual2799()
        with self.assertRaisesRegex(ValueError,'STYLE'):
            prepare_adjacent(model.Key('root_support',1),replace(corner,radius_mode='STRICT'),settings,
                             required_roots=(required,),support_corners=support)
        with self.assertRaisesRegex(ValueError,'already certified'):
            prepare_adjacent(model.Key('root_support',1),replace(corner),settings,
                             required_roots=(required,),support_corners=support)
        for roots,circles in (((),support),((corner.first,),support),((required,),())):
            with self.assertRaises(ValueError):
                prepare_adjacent(model.Key('root_support',1),corner,settings,required_roots=roots,support_corners=circles)
        moved=replace(required.source,points=tuple((x+10.,y,z) for x,y,z in required.source.points))
        remote=Approach(Path(moved),required.end,Crown(.08,0.,-1))
        with self.assertRaisesRegex(ValueError,'same-node original'):
            prepare_adjacent(model.Key('root_support',1),corner,settings,
                             required_roots=(remote,),support_corners=support)

    def test_actual274_reversed_routed_and_transformed_curbs_keep_both_roots_and_tangent_jets(self):
        settings=model.Settings(curve_tolerance=.00025)
        for angle,translation,reverse,routed in ((0.,(0.,0.),False,False),(.713,(1200.,-1800.),True,True)):
            with self.subTest(reverse=reverse,routed=routed):
                a,b=actual274(angle,translation,reverse,routed);before=(a.source.points,b.source.points)
                original=CornerRadiusPolicy.solve(a,b,10.,'STYLE')
                self.assertGreater(exposed_construction(original,
                    corner_patch(model.Key('circle',1),original,settings.spacing,settings.curve_tolerance),settings),.03)
                record=prepare(model.Key('envelope',1),original,settings)
                self.assertEqual(record.report['preferred_minimum_radius'],10.)
                self.assertLessEqual(record.report['curvature_bound'],(1.+1e-8)/10.)
                self.assertLessEqual(record.report['reach_steps'],16)
                self.assertGreater(record.report['outward_envelope_m'],0.)
                self.assertLessEqual(record.report['outward_envelope_m'],record.report['maximum_outward_envelope_m'])
                self.assertLess(exposed_construction(record.corner,record.source,settings),1e-7)
                np.testing.assert_allclose(record.curve(0.),a.point(record.stations[0],-.5*a.width),atol=1e-10)
                np.testing.assert_allclose(record.curve(1.),b.point(record.stations[1],.5*b.width),atol=1e-10)
                for t,approach,station,distance,sign in ((0.,a,record.stations[0],-.5*a.width,-1.),
                                                       (1.,b,record.stations[1],.5*b.width,1.)):
                    h=1e-6;lo=max(0.,t-h);hi=min(1.,t+h)
                    velocity=record.curve(hi)-record.curve(lo)
                    expected=sign*approach.jet(station,distance)[1]
                    np.testing.assert_allclose(velocity/np.linalg.norm(velocity),expected/np.linalg.norm(expected),atol=1e-6)
                fine=CertifiedCurbSampler(settings.curve_tolerance,settings.spacing,8192)
                coarse=CertifiedCurbSampler(.02,settings.spacing,8192)
                np.testing.assert_array_equal(fine.sample(record.curve,0.,1.),coarse.sample(record.curve,0.,1.))
                self.assertEqual(before,(a.source.points,b.source.points))

    def test_strict_derived_and_exhausted_finite_profiles_remain_rejected(self):
        a,b=actual274();corner=CornerRadiusPolicy.solve(a,b,10.,'STYLE');settings=model.Settings()
        with self.assertRaisesRegex(ValueError,'STYLE'):
            prepare(model.Key('envelope',1),replace(corner,radius_mode='STRICT'),settings)
        with self.assertRaises(ApproachExtentError):
            prepare(model.Key('envelope',1),corner,settings,(a.length,b.length))
        for attribute in ('half','width','offset','absorbed_sites','mouth_piece','deferred_profile_reserve'):
            a,b=actual274();corner=CornerRadiusPolicy.solve(a,b,10.,'STYLE')
            if attribute=='half':a.profile=replace(a.profile,half=1)
            elif attribute=='width':a.width+=.1
            elif attribute=='offset':a.offset=.1
            else:setattr(a.path,attribute,(model.Key('intersection',1),))
            with self.subTest(attribute=attribute),self.assertRaisesRegex(ValueError,'unmodified original'):
                prepare(model.Key('envelope',1),corner,settings)


if __name__=='__main__':unittest.main()
