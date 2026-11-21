"""Actual append ancestry separates original arcs from auxiliary curb leads."""
from copy import copy
from dataclasses import replace,FrozenInstanceError
import unittest
import numpy as np
import _python_street_package
from test_python_street_corner_patch import approaches
from tested_python_streets.model import Key
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.corner_patch import corner_patch,_PreparedCurve
from tested_python_streets.source_circular_intervals import declare


class SourceCircularIntervalTests(unittest.TestCase):
    def test_real_positive_negative_and_mixed_producer_phases_retain_original_leads(self):
        for case in ('positive','negative','mixed'):
            for routed in (False,True):
                for angle,shift in ((0.,(0.,0.)),(.713,(1200.,-1800.))):
                    a,b=approaches(case,angle,shift,routed)
                    corner=CornerRadiusPolicy.solve(a,b,4.,'STRICT')
                    source=corner_patch(Key('source_circle',1),corner,.5,.00025)
                    records=[(curve,curve.original_circle_interval) for curve in source.layout.boundaries
                        if hasattr(curve,'original_circle_interval')]
                    self.assertEqual(len(records),1)
                    curve,record=records[0];self.assertTrue(record.validate(curve,corner))
                    self.assertEqual(len(record.auxiliary_leads),sum(value<0. for value in corner.stations))
                    self.assertEqual(record.maximum_chord_error,.00025)
                    phase=np.linspace(0.,1.,31)
                    np.testing.assert_allclose(record.circle_parameter(record.source_parameter(phase)),phase,atol=2e-15,rtol=0.)
                    for lead in record.auxiliary_leads:
                        self.assertTrue(lead.validate(curve,corner))
                        self.assertTrue(all(at<=0. for at in lead.outward_stations))
                    if case=='positive':self.assertEqual(record.prepared_bounds,(0.,1.))
                    else:self.assertNotEqual(record.prepared_bounds,(0.,1.))

    def test_changed_vertices_parameters_corner_or_lead_identity_are_rejected(self):
        a,b=approaches('mixed',routed=True);corner=CornerRadiusPolicy.solve(a,b,4.,'STRICT')
        source=corner_patch(Key('source_circle',1),corner,.5,.00025)
        curve=next(curve for curve in source.layout.boundaries if hasattr(curve,'original_circle_interval'))
        record=curve.original_circle_interval
        with self.assertRaises(FrozenInstanceError):record.maximum_chord_error=.1
        moved=copy(curve);moved.points=curve.points.copy()
        moved.points[record.vertex_indices[1],0]=np.nextafter(moved.points[record.vertex_indices[1],0],np.inf)
        with self.assertRaisesRegex(ValueError,'source knots or vertices'):record.validate(moved,corner)
        with self.assertRaisesRegex(ValueError,'Corner identity'):record.validate(curve,replace(corner))
        fake=replace(record,circle_parameters=tuple(np.asarray(record.circle_parameters)**2))
        with self.assertRaisesRegex(ValueError,'original analytical arc'):fake.validate(curve,corner)
        lead,=record.auxiliary_leads
        fake=replace(lead,original_approach=a if lead.original_approach is b else b)
        with self.assertRaisesRegex(ValueError,'source/vertex ancestry'):fake.validate(curve,corner)

    def test_finite_piece_bounds_reverse_maps_and_legacy_no_inference(self):
        a,b=approaches('mixed');corner=CornerRadiusPolicy.solve(a,b,4.,'STRICT')
        source=corner_patch(Key('source_circle',1),corner,.5,.00025)
        curve=next(curve for curve in source.layout.boundaries if hasattr(curve,'original_circle_interval'))
        record=curve.original_circle_interval
        for at in (np.nextafter(record.prepared_bounds[0],-np.inf),np.nextafter(record.prepared_bounds[1],np.inf)):
            with self.assertRaisesRegex(ValueError,'finite interval'):record.circle_parameter(at)
        legacy=_PreparedCurve(curve.points.copy());self.assertFalse(hasattr(legacy,'original_circle_interval'))
        # Explicit producer reversal preserves the same real arc and controls;
        # the metadata maps circle phase backwards through prepared distance.
        reversed_curve=_PreparedCurve(curve.points[::-1].copy())
        phases=np.full(len(curve.points),np.nan);phases[np.asarray(record.vertex_indices)]=record.circle_parameters
        tags=[None]*len(curve.points)
        reversed_record=declare(corner,reversed_curve,phases[::-1],tags,.00025)
        self.assertTrue(reversed_record.validate(reversed_curve,corner))
        self.assertGreater(reversed_record.circle_parameters[0],reversed_record.circle_parameters[-1])
        phase=np.linspace(0.,1.,31)
        np.testing.assert_allclose(reversed_record.circle_parameter(reversed_record.source_parameter(phase)),phase,atol=2e-15,rtol=0.)


if __name__=='__main__':unittest.main()
