"""Declared source pieces change jet evaluation, never physical edge weights."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_surfaces import boundary_data
from tested_python_streets.conditions import PlateauEdge
from tested_python_streets.surface_fields import BoundaryJetField


class Kink:
    def sample(self,t):
        at=np.asarray(t)
        return np.full(at.shape,.04),np.stack((np.zeros(at.shape),.02+.04*abs(at-.371)),axis=-1)


class DeclaredKink(Kink):
    def polynomial_breaks(self):return (.371,)


class Linear:
    def sample(self,t):
        at=np.asarray(t)
        return np.full(at.shape,.04),np.stack((np.zeros(at.shape),.02+.04*at),axis=-1)


class SplitLinear(Linear):
    def polynomial_breaks(self):return (.123,.4567)


class OneSided:
    def polynomial_breaks(self):return (.5,)
    def sample(self,t):
        at=np.asarray(t)
        return np.full(at.shape,.04),np.stack((np.zeros(at.shape),np.where(at<.5,.02,.07)),axis=-1)
    def sample_piece(self,t,first,last):
        at=np.asarray(t);value=.02 if .5*(first+last)<.5 else .07
        return np.full(at.shape,.04),np.stack((np.zeros(at.shape),np.full(at.shape,value)),axis=-1)


def boundary(rule):
    return replace(boundary_data(),heights=np.full(4,.04),gradients=np.zeros((4,2)),
                   rules=(rule,PlateauEdge(.04),PlateauEdge(.04),PlateauEdge(.04)))


class BoundaryJetPieceTests(unittest.TestCase):
    def test_declared_knot_meets_original_fit_limit_and_keeps_original_edge_weights(self):
        with self.assertRaisesRegex(ValueError,'different approximation branch'):
            BoundaryJetField((boundary(Kink()),))
        data=boundary(DeclaredKink());field=BoundaryJetField((data,))
        self.assertLessEqual(field.metadata['maximum_sampled_fit_error'],1e-8)
        self.assertEqual(field.metadata['declared_polynomial_pieces']['prediction_branches'],5)
        np.testing.assert_array_equal(field.origins,data.xy)
        np.testing.assert_array_equal(field.edges,np.roll(data.xy,-1,axis=0)-data.xy)
        points=np.array(((.2,.1),(.6,.3),(.37,.02),(.8,.4)))
        x,y=points.T;gradient=.02+.04*abs(x-.371)
        # Four original unit edges retain their inverse-fourth-distance weights.
        edge_weight_fraction=1./(1.+(y/(1.-x))**4+(y/(1.-y))**4+(y/x)**4)
        expected=.04+.06*np.tanh(y*gradient/.06)*edge_weight_fraction
        np.testing.assert_allclose(field(points),expected,atol=2e-13,rtol=0.)

    def test_redundant_declared_breaks_preserve_the_original_field_and_batching(self):
        points=np.random.default_rng(1729).uniform(.01,.99,(1001,2))
        original=BoundaryJetField((boundary(Linear()),))
        split=BoundaryJetField((boundary(SplitLinear()),))
        small=BoundaryJetField((boundary(SplitLinear()),),scratch_bytes=4096)
        np.testing.assert_allclose(split(points),original(points),atol=5e-16,rtol=0.)
        np.testing.assert_allclose(small(points),split(points),atol=5e-16,rtol=0.)
        self.assertEqual(split.batch_rows,original.batch_rows)
        for _,bounds,coefficients in split._pieces:
            self.assertFalse(bounds.flags.writeable);self.assertFalse(coefficients.flags.writeable)

    def test_true_one_sided_source_jets_are_not_averaged_at_a_declared_break(self):
        field=BoundaryJetField((boundary(OneSided()),))
        points=np.array(((.499,.1),(.5,.1),(.501,.1)))
        x,y=points.T;gradient=np.where(x<.5,.02,.07)
        weight=1./(1.+(y/(1.-x))**4+(y/(1.-y))**4+(y/x)**4)
        expected=.04+.06*np.tanh(y*gradient/.06)*weight
        np.testing.assert_allclose(field(points),expected,atol=2e-13,rtol=0.)

    def test_invalid_or_excessive_declarations_refuse_before_any_fit(self):
        for values in ([.5],(0.,), (1.,), (float('nan'),), (.7,.3), (.5,.5),
                       tuple(np.linspace(.0001,.9999,19999))):
            rule=DeclaredKink();rule.polynomial_breaks=lambda:values
            with patch.object(rule,'sample',side_effect=AssertionError('Budget must be checked before fitting')):
                with self.assertRaises(ValueError):BoundaryJetField((boundary(rule),))
        rule=DeclaredKink();rule.polynomial_breaks=False
        with self.assertRaises(ValueError):BoundaryJetField((boundary(rule),))

    def test_explicit_pieces_cannot_bypass_residual_degree_or_branch_guards(self):
        class Unresolved(DeclaredKink):
            def sample(self,t):
                at=np.asarray(t)
                return np.full(at.shape,.04),np.stack((np.zeros(at.shape),np.sin(1200.*at)),axis=-1)
        with self.assertRaisesRegex(ValueError,'different approximation branch'):
            BoundaryJetField((boundary(Unresolved()),))


if __name__=='__main__':unittest.main()
