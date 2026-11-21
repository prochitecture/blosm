"""Declared carrier/source maps preserve one-sided original jets and physical frames."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Key
from tested_python_streets.paths import right
from tested_python_streets.source_rim_normalization import ResampledSourceCurb
from tested_python_streets.source_curb_conditions import _MappedRule
from tested_python_streets.surface_fields import _compile_rule,_RulePiece


class ConstantJet:
    def __init__(self,gradient):self.gradient=np.asarray(gradient,dtype=float)
    def sample(self,t):
        shape=np.asarray(t).shape
        return np.zeros(shape),np.broadcast_to(self.gradient,shape+(2,)).copy()


class NativePieces:
    """A native one-sided API permits only its full declared branch bounds."""
    def polynomial_breaks(self):return (.4,)
    def sample_piece(self,t,first,last):
        if (first,last) not in ((0.,.4),(.4,1.)):raise ValueError('Unproven native branch')
        t=np.asarray(t)
        if np.any(t<0.) or np.any(t>1.):raise ValueError('Native samples leave their branch')
        return np.zeros(t.shape),np.broadcast_to((0.,2. if first==0. else 5.),t.shape+(2,)).copy()
    def sample(self,t):
        t=np.asarray(t)
        return np.zeros(t.shape),np.stack((np.zeros(t.shape),np.where(t<.4,2.,5.)),axis=-1)


def fixture(points,rules,first=0.,last=1.):
    owner=Key('normalized_source_rim',1)
    curve=ResampledSourceCurb(np.asarray(points),tuple(object() for _ in rules),())
    resolver=SimpleNamespace(mapped={},sources={owner:SimpleNamespace(layout=SimpleNamespace(approaches=()))},
        chain_rules=lambda *args,**kwargs:tuple(rules))
    return _MappedRule(resolver,curve,first,last,owner,0)


class MappedSourceRuleTests(unittest.TestCase):
    def test_one_ulp_before_a_source_knot_keeps_its_left_branch(self):
        start=np.nextafter(.75,0.)
        self.assertEqual(.5*(start+.75),.75)
        rule=fixture(((0.,0.),(1.5,0.),(2.,0.)),
                     (ConstantJet((0.,1.)),ConstantJet((0.,2.))),start,1.)
        first,=rule.polynomial_breaks()
        for lo,hi,expected in ((0.,first,1.),(first,1.,2.)):
            _,gradient=rule.sample_piece(np.asarray((0.,1.)),lo,hi)
            np.testing.assert_array_equal(gradient,np.array(((0.,expected),)*2))
        np.testing.assert_array_equal(rule.sample(0.)[1],(0.,1.))

    def test_real_span_map_knots_compile_without_changing_carrier_geometry(self):
        rule=fixture(((0.,0.),(1.,2e-7),(2.,0.),(3.,0.)),
                     (ConstantJet((0.,1.)),ConstantJet((0.,2.)),ConstantJet((0.,3.))))
        original=rule.curve.points.copy();breaks=rule.polynomial_breaks()
        self.assertEqual(len(breaks),2)
        for first,last in zip((0.,*breaks),(*breaks,1.)):
            _,error,_=_compile_rule(_RulePiece(rule,first,last),1e-8)
            self.assertLessEqual(error,1e-8)
        np.testing.assert_array_equal(rule.curve.points,original)
        self.assertEqual(rule.curve.maximum_projection_error,2e-7)

    def test_clipped_native_piece_keeps_the_declared_one_sided_parent_branch(self):
        rule=fixture(((0.,0.),(3.,0.)),(NativePieces(),),.1,.9)
        breakpoint,=rule.polynomial_breaks()
        self.assertAlmostEqual(breakpoint,.375)
        for first,last,expected in ((0.,breakpoint,2.),(breakpoint,1.,5.)):
            height,gradient=rule.sample_piece(np.array((0.,.5,1.)),first,last)
            np.testing.assert_array_equal(height,np.zeros(3))
            np.testing.assert_array_equal(gradient,np.array(((0.,expected),)*3))
            self.assertLessEqual(_compile_rule(_RulePiece(rule,first,last),1e-8)[1],1e-8)
        with self.assertRaisesRegex(ValueError,'consecutive'):rule.sample_piece(.5,0.,1.)

    def test_published_original_span_end_uses_the_inside_jet(self):
        whole=fixture(((0.,0.),(1.,0.),(2.,0.),(3.,0.)),
                      (ConstantJet((0.,1.)),ConstantJet((0.,2.)),ConstantJet((0.,3.))))
        rule=fixture(whole.curve.original_polyline,whole.rules,1./3.,2./3.)
        self.assertEqual(rule.polynomial_breaks(),())
        _,gradient=rule.sample(np.array((0.,.25,.5,.75,1.)))
        np.testing.assert_array_equal(gradient,np.array(((0.,2.),)*5))
        self.assertLessEqual(_compile_rule(rule,1e-8)[1],1e-8)

    def test_a_single_old_piece_retains_the_original_sample_arithmetic(self):
        rule=fixture(((1200.,-1800.),(1203.,-1796.)),(ConstantJet((.2,.3)),),.13,.91)
        values=np.array((0.,.1,.3,.7,1.))
        height,gradient=rule.sample(values)
        tangent=rule.curve.original_polyline[1]-rule.curve.original_polyline[0]
        tangent/=np.linalg.norm(tangent)
        old=np.tile((.2,.3),(len(values),1))
        along,inward=old@tangent,old@-right(tangent)
        expected=along[:,None]*rule.arc_scale*rule.tangent-inward[:,None]*right(rule.tangent)
        self.assertEqual(rule.polynomial_breaks(),())
        np.testing.assert_array_equal(height,np.zeros(len(values)))
        np.testing.assert_array_equal(gradient,expected)

    def test_a_projected_carrier_piece_uses_its_actual_map_jacobian(self):
        rule=fixture(((0.,0.),(1.,2e-7),(2.,0.)),
                     (ConstantJet((.7,.2)),ConstantJet((.4,.9))))
        first,last=0.,rule.polynomial_breaks()[0]
        _,gradient=rule.sample_piece(np.array((0.,.5,1.)),first,last)
        old_tangent=(rule.curve.original_polyline[1]-rule.curve.original_polyline[0])/rule.old_stations[1]
        physical_tangent=rule.curve.points[1]-rule.curve.points[0]
        physical_tangent/=np.linalg.norm(physical_tangent)
        scale=(rule.old_at[1]-rule.old_at[0])/rule.curve.length
        native=np.array((.7,.2))
        expected=(native@old_tangent)*scale*physical_tangent-(native@-right(old_tangent))*right(physical_tangent)
        np.testing.assert_allclose(gradient,np.tile(expected,(3,1)),atol=1e-15,rtol=0.)


if __name__=='__main__':unittest.main()
