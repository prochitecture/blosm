"""Native zero seams and original negative auxiliary leads retain exact physical jets."""
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Network,Key
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.source_offset_rules import SourceOffsetRule
from tested_python_streets.surface_fields import _compile_rule,_RulePiece


def fixture(*,reverse=False,angle=0.,shift=(0.,0.),interpolation=3,q=-1.):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    xy=np.array(((0.,0.),(10.,1.),(20.,0.)))@rotation.T+shift
    if reverse:xy=xy[::-1]
    section=Section(1,1,tuple((x,y,0.) for x,y in xy),width=5.4,interpolation=interpolation)
    path=Path(section,.005);profile=Crown(.08,0.)
    context=BuildContext(Network('native-auxiliary-source',(section,)),
        {Key('path',1):path,Key('profile',1):profile},PatternLibrary())
    approach=Approach(Route((Piece(path,not reverse),)),int(reverse),profile)
    parameters=np.array((.0901108495,0.,-.0293911389))
    offset=approach.offset+.5*approach.width*q
    points=approach.point(parameters,offset)
    cumulative=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
    lead=SimpleNamespace(approach=approach,parameters=parameters,cumulative=cumulative,
        offset=offset,length=float(cumulative[-1]))
    return SourceOffsetRule(context,lead,0.,1.)


class SourceOffsetAuxiliaryTests(unittest.TestCase):
    def test_auxiliary_covector_is_exact_on_reversed_and_rigid_source_frames(self):
        for interpolation in (0,3):
            for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1347.,-1800.))):
                for q in (-1.,1.):
                    with self.subTest(interpolation=interpolation,reverse=reverse,q=q):
                        rule=fixture(reverse=reverse,angle=angle,shift=shift,interpolation=interpolation,q=q)
                        seam,=rule.polynomial_breaks()
                        self.assertEqual(rule._station_anchors[seam],0.)
                        height,gradient=rule.sample_piece(np.array((0.,.2,.5,.8,1.)),seam,1.)
                        _,tangent,normal=rule.curve.approach.frame(0.)
                        basis=np.stack((tangent,.5*rule.curve.approach.width*normal))
                        expected=np.linalg.solve(basis,(0.,rule.curve.approach.profile.derivative(q)))
                        np.testing.assert_array_equal(gradient,np.tile(expected,(5,1)))
                        np.testing.assert_allclose(gradient@tangent,0.,atol=1e-15,rtol=0.)
                        np.testing.assert_allclose(gradient@(.5*rule.curve.approach.width*normal),
                            rule.curve.approach.profile.derivative(q),atol=1e-15,rtol=0.)
                        np.testing.assert_array_equal(height,np.zeros(5))
                        self.assertLessEqual(_compile_rule(_RulePiece(rule,seam,1.),1e-8)[1],1e-8)

    def test_positive_native_branch_and_true_zero_endpoint_compile_without_curve_changes(self):
        for reverse in (False,True):
            rule=fixture(reverse=reverse,angle=.713,shift=(1347.,-1800.))
            values=np.linspace(0.,1.,101)
            original=np.interp(values*rule.curve.length,rule.curve.cumulative,rule.curve.parameters)
            seam,=rule.polynomial_breaks()
            self.assertEqual(rule._station_anchors[seam],0.)
            _,gradient=rule.sample_piece(np.array((0.,.5,1.)),0.,seam)
            station=rule._stations(np.array((0.,.5*seam,seam)));station[-1]=0.
            _,tangent,normal=rule.curve.approach.frame(station)
            np.testing.assert_allclose(np.sum(gradient*tangent,axis=1),0.,atol=1e-14,rtol=0.)
            np.testing.assert_allclose(np.sum(gradient*(.5*rule.curve.approach.width*normal),axis=1),
                rule.curve.approach.profile.derivative(rule.q),atol=1e-14,rtol=0.)
            self.assertLessEqual(_compile_rule(_RulePiece(rule,0.,seam),1e-8)[1],1e-8)
            np.testing.assert_array_equal(rule._stations(values),original)

    def test_an_undeclared_interval_crossing_zero_is_not_clamped(self):
        rule=fixture()
        with self.assertRaisesRegex(ValueError,'undeclared original zero'):
            rule._gradient(np.array((-.01,.02)),(-.01,.02))


if __name__=='__main__':unittest.main()
