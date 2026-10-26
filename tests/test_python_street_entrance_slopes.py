"""Longitudinal entrance limits are distinct from crown and total-slope bounds."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.entrance_slopes import EntranceSlopeGuard


class PythonEntranceSlopeTests(unittest.TestCase):
    def test_large_transverse_grade_and_speed_variation_are_not_longitudinal_grade(self):
        xy = np.array(((0.,0.),(2.,0.),(2.,1.),(0.,1.)))
        rotation = np.array(((np.cos(.63),-np.sin(.63)),(np.sin(.63),np.cos(.63))))
        world = xy @ rotation.T + np.array((1e6,-2e6))
        physical = .2*xy[:,1]
        mesh = Mesh(np.column_stack((world, np.full(4,37.))), np.array(((0,1,2,3),)))
        velocity = np.array((.5,1.,3.,4.))[:,None]*rotation[:,0]
        value = EntranceSlopeGuard().validate(mesh,(0,3),velocity,physical)
        self.assertLess(value,1e-8)

    def test_held_profiles_and_height_ceiling_do_not_hide_a_first_row_ramp(self):
        pattern = PatternLibrary().grid(2,8)
        u,v = pattern.uv.T
        x,y = .5*u,4.*v
        profile = .08*(2*v-v*v)
        heights = profile+.6*x*(1.-2*x)*v*(1.-v)**2
        mesh = Mesh(np.column_stack((x,y,heights)),pattern.quads)
        np.testing.assert_array_equal(heights[(u==0.)|(u==1.)],profile[(u==0.)|(u==1.)])
        self.assertGreaterEqual(heights.min(),0.)
        self.assertLessEqual(heights.max(),.08)
        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'],8.)
        self.assertLess(mesh.slope_statistics()['maximum'],.1)
        velocity = np.tile((1.,0.),(len(x),1))
        with self.assertRaisesRegex(ValueError,'Entrance longitudinal.*exceeds 3%'):
            EntranceSlopeGuard().validate(mesh,np.flatnonzero(u==0.),velocity,heights)

    def test_an_interior_direction_peak_is_not_missed_by_corner_projections(self):
        xy = np.array(((0.,0.),(2.,0.),(2.,1.),(0.,1.)))
        heights = .04*xy[:,0]
        mesh = Mesh(np.column_stack((xy,heights)),np.array(((0,1,2,3),)))
        velocity = np.array(((.5,np.sqrt(.75)),(.5,np.sqrt(.75)),
                             (.5,-np.sqrt(.75)),(.5,-np.sqrt(.75))))
        # Each endpoint projects to only 2%; the midpoint's direction gives 4%.
        self.assertAlmostEqual(float(np.array((.04,0.))@velocity[0]),.02)
        self.assertAlmostEqual(float(np.array((.04,0.))@velocity[3]),.02)
        guard = EntranceSlopeGuard()
        self.assertAlmostEqual(guard.measure(mesh,(0,3),velocity,heights),.04)
        with self.assertRaisesRegex(ValueError,'exceeds 3%'):
            guard.validate(mesh,(0,3),velocity,heights)

    def test_zero_inspection_height_cannot_bypass_the_physical_limit(self):
        xy = np.array(((0.,0.),(2.,0.),(2.,1.),(0.,1.)))
        mesh = Mesh(np.column_stack((xy,np.full(4,37.))),np.array(((0,1,2,3),)))
        velocity = np.tile((1.,0.),(4,1))
        guard = EntranceSlopeGuard()
        self.assertAlmostEqual(guard.validate(mesh,(0,3),velocity,.025*xy[:,0]),.025)
        with self.assertRaisesRegex(ValueError,'exceeds 3%'):
            guard.validate(mesh,(0,3),velocity,.04*xy[:,0])

    def test_conservative_bound_covers_dense_bilinear_edge_samples(self):
        random = np.random.default_rng(941)
        t = np.linspace(0.,1.,1001)
        for _ in range(32):
            xy = np.array(((0.,0.),(2.2,.1),(2.,1.1),(-.1,1.)))+random.uniform(-.04,.04,(4,2))
            heights = random.uniform(0.,.12,4)
            mesh = Mesh(np.column_stack((xy,heights)),np.array(((0,1,2,3),)))
            mesh.validate()
            angles = random.uniform(-1.2,1.2,2)
            speeds = np.exp(random.uniform(-2.,2.,2))
            first,last = speeds[:,None]*np.column_stack((np.cos(angles),np.sin(angles)))
            velocity = np.array((first,first,last,last))
            p = np.column_stack((xy,heights))
            du = (1.-t[:,None])*(p[1]-p[0])+t[:,None]*(p[2]-p[3])
            dv = p[3]-p[0]
            determinant = du[:,0]*dv[1]-du[:,1]*dv[0]
            gradient = np.column_stack(((du[:,2]*dv[1]-dv[2]*du[:,1])/determinant,
                                        (du[:,0]*dv[2]-dv[0]*du[:,2])/determinant))
            direction = (1.-t[:,None])*first+t[:,None]*last
            direction /= np.linalg.norm(direction,axis=1)[:,None]
            sampled = float(abs(np.sum(gradient*direction,axis=1)).max())
            bound = EntranceSlopeGuard().measure(mesh,(0,3),velocity,heights)
            self.assertGreaterEqual(bound+1e-12,sampled)


if __name__=='__main__':unittest.main(verbosity=2)
