"""Independent geometric and allocation contracts for regular corner-preserving rows."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.polyline_sampling import ProtectedPolylineSampler
from tested_python_streets.bundle_angular_curves import polyline_distance


class PythonProtectedSamplingTests(unittest.TestCase):
    def test_exact_corners_and_semantic_partitions_remain_mandatory(self):
        points = np.array(((0.,0.),(1.,0.),(3.,0.),(3.,2.),(3.,3.),(5.,3.)))
        sampled,owners,report = ProtectedPolylineSampler(spacing=.7,tolerance=.002).sample(points,(2,4))
        for point in points[[0,2,4,5]]:
            self.assertTrue(np.any(np.all(sampled==point,axis=1)))
        self.assertEqual(len(owners),len(sampled)-1)
        self.assertEqual(set(owners),{0,2,4})
        self.assertLessEqual(float(polyline_distance(sampled,points).max()),1e-12)
        self.assertLessEqual(float(polyline_distance(points,sampled).max()),1e-12)
        self.assertLessEqual(report['maximum_station_error'],1e-12)

    def test_dense_curved_reference_and_rigid_transform_stay_within_budget(self):
        t = np.linspace(0.,1.,1501)
        points = np.c_[8.*t,.2*np.sin(4.*np.pi*t)]
        angle = .713
        rotation = np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        for original in (points,points@rotation.T+(1234.,-5678.)):
            sampled,_,report = ProtectedPolylineSampler(spacing=.7,tolerance=.002).sample(original,(500,1000))
            self.assertLessEqual(float(polyline_distance(original,sampled).max()),.002+1e-10)
            self.assertLessEqual(float(polyline_distance(sampled,original).max()),.002+1e-10)
            self.assertLessEqual(report['maximum_station_error'],.002)
            self.assertLessEqual(len(sampled),4097)

    def test_bad_inputs_and_excessive_sampling_fail_before_large_allocation(self):
        for spacing,tolerance in ((0.,.001),(.7,0.),(float('inf'),.001)):
            with self.assertRaises(ValueError):
                ProtectedPolylineSampler(spacing=spacing,tolerance=tolerance)
        sampler = ProtectedPolylineSampler(spacing=.1,tolerance=.001,max_intervals=4)
        with self.assertRaisesRegex(ValueError,'budget'):
            sampler.sample(np.array(((0.,0.),(100.,0.))))
        with self.assertRaises(ValueError):
            sampler.sample(np.array(((0.,0.),(0.,0.))))
