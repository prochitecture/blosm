"""Finite original bank fillets preserve sides, orientation and authored radii."""
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Network,Key
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.source_curb_rounding import FiniteSourceBank,fit_banks,patch
from tested_python_streets.corner_radii import CornerContactError


def fixture(reverse=False,angle=0.,shift=(0.,0.)):
    # Two positive banks cross with a reflex turn in the physical union rim.
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    controls=(np.array(((-20.,0.),(20.,0.))),np.array(((-20.,20.),(20.,-20.))))
    sections=[];products={}
    for sid,(points,width) in enumerate(zip(controls,(5.,2.)),1):
        if reverse:points=points[::-1]
        points=points@rotation.T+shift
        section=Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width,interpolation=0)
        sections.append(section);products[Key('path',sid)]=Path(section)
        products[Key('profile',sid)]=Crown(.08,0.)
    context=BuildContext(Network('native-bank-rounding',tuple(sections)),products,PatternLibrary())
    side=-1 if reverse else 1
    banks=tuple(FiniteSourceBank(products[Key('path',sid)],side,((0.,products[Key('path',sid)].length),)) for sid in (1,2))
    # Native offset banks intersect analytically, independent of the solver.
    guess=np.array((2.5-math.sqrt(2.),-2.5))@rotation.T+shift
    return context,banks,guess


class SourceBankRoundingTests(unittest.TestCase):
    def test_finite_bank_orientation_keeps_the_exposed_circle_under_reversal_and_rigid_motion(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1800.,-2500.))):
            context,banks,guess=fixture(reverse,angle,shift)
            before=context.network.to_dict()
            result=fit_banks(context,*banks,1.,'STRICT',guess)
            self.assertEqual(result.corner.radius,1.)
            self.assertEqual(result.corner.radius_attempts,(1.,))
            self.assertEqual(context.network.to_dict(),before)
            self.assertIs(result.banks[0].path,banks[0].path)
            for arm,station,bank,cut,side in zip((result.corner.first,result.corner.second),
                    result.corner.stations,banks,result.intersection,(-1.,1.)):
                native=lambda s:s if arm.end==0 else arm.length-s
                self.assertTrue(bank.contains(native(station),native(cut)))
                point=arm.point(station,arm.offset+.5*side*arm.width)
                self.assertAlmostEqual(float(np.linalg.norm(point-result.corner.center)),1.,places=8)
                tangent=arm.jet(station,arm.offset+.5*side*arm.width)[1]
                self.assertAlmostEqual(float((point-result.corner.center)@tangent),0.,places=8)

    def test_disjoint_support_does_not_become_a_hull_and_style_cannot_leave_the_ranges(self):
        context,banks,guess=fixture()
        bank=FiniteSourceBank(banks[0].path,1,((0.,19.),(21.,40.)))
        self.assertFalse(bank.contains(18.,22.))
        tiny=FiniteSourceBank(banks[0].path,1,((0.,.01),))
        for mode in ('STRICT','STYLE'):
            with self.assertRaises(CornerContactError):fit_banks(context,tiny,banks[1],1.,mode,guess)

    def test_false_identity_invalid_ranges_and_untyped_or_nonfinite_inputs_are_rejected(self):
        context,banks,guess=fixture()
        other,_,_=fixture()
        with self.assertRaisesRegex(ValueError,'Path identity'):fit_banks(other,*banks,1.,'STRICT',guess)
        for side,ranges in ((True,((0.,1.),)),(0,((0.,1.),)),(1,((2.,1.),)),
                            (1,((0.,2.),(1.,3.))),(1,((0.,math.inf),))):
            with self.assertRaises(ValueError):FiniteSourceBank(banks[0].path,side,ranges)
        for radius,mode,point in ((0.,'STYLE',guess),(math.nan,'STYLE',guess),(1.,'UNKNOWN',guess),
                                  (1.,'STYLE',(math.nan,0.))):
            with self.assertRaises(ValueError):fit_banks(context,*banks,radius,mode,point)

    def test_preferred_radius_fits_the_real_support_but_strict_radius_does_not_shrink(self):
        context,banks,guess=fixture()
        original=fit_banks(context,*banks,1.,'STRICT',guess)
        narrowed=[]
        for bank,arm,cut in zip(banks,(original.corner.first,original.corner.second),original.intersection):
            station=cut if arm.end==0 else arm.length-cut
            narrowed.append(FiniteSourceBank(bank.path,bank.side,((station-.36,station+.36),)))
        with self.assertRaises(CornerContactError):fit_banks(context,*narrowed,1.,'STRICT',guess)
        fitted=fit_banks(context,*narrowed,1.,'STYLE',guess)
        self.assertEqual(fitted.corner.radius,.8)
        self.assertEqual(fitted.corner.requested_radius,1.)
        self.assertEqual(fitted.corner.radius_attempts,(1.,.8))

    def test_nearly_tangent_circle_uses_the_exact_native_contact_endpoints(self):
        from dataclasses import replace
        sections=(Section(1,1,((-1022.,-600.,0.),(-1122.,-600.,0.)),width=10.8,interpolation=3),
            Section(2,1,((-1122.,-600.,0.),(-1222.,-599.4,0.)),width=10.8,interpolation=3))
        products={Key('path',s.id):Path(s,.005) for s in sections}
        products.update({Key('profile',s.id):Crown(.08,0.) for s in sections})
        context=BuildContext(Network('near-tangent-native-contact',sections),products,PatternLibrary())
        banks=tuple(FiniteSourceBank(products[Key('path',s.id)],1,((0.,products[Key('path',s.id)].length),)) for s in sections)
        fitted=fit_banks(context,*banks,4.,'STRICT',np.array((-1122.,-594.6)))
        partitions=tuple(np.unique((0.,contact,crossing,arm.length)) for arm,contact,crossing in zip(
            (fitted.corner.first,fitted.corner.second),fitted.corner.stations,fitted.intersection))
        source=patch(Key('native_rounding',0),fitted.corner,fitted.intersection,.7,.00025,partitions=partitions)
        index=next(segment.index for segment in source.segments if segment.kind=='curb')
        curve=source.layout.boundaries[index]
        actual=[arm.point(at,arm.offset+side*.5*arm.width) for arm,at,side in zip(
            (fitted.corner.first,fitted.corner.second),fitted.corner.stations,(-1,1))]
        self.assertEqual({tuple(point) for point in curve.points[[0,-1]]},{tuple(point) for point in actual})
        with self.assertRaisesRegex(ValueError,'representable native contact agreement'):
            patch(Key('native_rounding',1),replace(fitted.corner,center=fitted.corner.center+(.000001,0.)),
                  fitted.intersection,.7,.00025,partitions=partitions)


if __name__=='__main__':unittest.main()
