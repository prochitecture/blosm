"""Numerical navigation invariants, independent of bpy or a display server."""
import importlib.util
import math
from pathlib import Path
import random
import unittest

spec=importlib.util.spec_from_file_location('panorama_navigation',Path(__file__).resolve().parents[1]/'viewport/navigation.py')
nav=importlib.util.module_from_spec(spec)
spec.loader.exec_module(nav)


class NavigationTests(unittest.TestCase):
    def upright(self,control):
        rows=control.matrix()
        right,up,back=map(tuple,zip(*rows))
        self.assertGreater(up[2],0)
        self.assertEqual(right[2],0)
        for axis in (right,up,back):self.assertAlmostEqual(sum(x*x for x in axis),1,places=12)
        for a,b in ((right,up),(up,back),(back,right)):
            self.assertAlmostEqual(sum(x*y for x,y in zip(a,b)),0,places=12)
        cross=(right[1]*up[2]-right[2]*up[1],right[2]*up[0]-right[0]*up[2],right[0]*up[1]-right[1]*up[0])
        self.assertAlmostEqual(sum(x*y for x,y in zip(cross,back)),1,places=12)

    def test_opening_direction_and_pole_fallback(self):
        control=nav.PanoramaControls((0,math.sqrt(.75),.5))
        self.assertAlmostEqual(control.yaw,0)
        self.assertAlmostEqual(control.pitch,math.pi/6)
        self.upright(control)
        for sign in (-1,1):
            control=nav.PanoramaControls((0,0,sign),(0,-1,0))
            self.assertAlmostEqual(control.yaw,math.pi/2)
            self.assertAlmostEqual(control.pitch,sign*nav.PITCH_LIMIT)
            self.upright(control)
        with self.assertRaises(ValueError):nav.PanoramaControls((0,0,0))

    def test_grabbing_direction_and_immediate_reversal_at_limit(self):
        control=nav.PanoramaControls();control.resize(1000,1000)
        control.drag(10,10)
        self.assertLess(control.yaw,0)
        self.assertLess(control.pitch,0)
        control.drag(0,-1e8)
        self.assertEqual(control.pitch,nav.PITCH_LIMIT)
        control.drag(0,1)
        self.assertLess(control.pitch,nav.PITCH_LIMIT)
        control.drag(0,1e8)
        self.assertEqual(control.pitch,-nav.PITCH_LIMIT)
        control.drag(0,-1)
        self.assertGreater(control.pitch,-nav.PITCH_LIMIT)

    def test_ten_thousand_mixed_gestures_never_flip_or_roll(self):
        rng=random.Random(1441)
        control=nav.PanoramaControls((.4,.3,.5));control.resize(1500,800)
        for _ in range(10000):
            control.drag(rng.uniform(-3000,3000),rng.uniform(-1800,1800))
            self.assertLessEqual(abs(control.pitch),nav.PITCH_LIMIT)
            self.assertGreaterEqual(control.yaw,-math.pi)
            self.assertLessEqual(control.yaw,math.pi)
            self.upright(control)

    def test_full_turns_have_no_drift(self):
        control=nav.PanoramaControls();control.resize(1000,1000)
        for _ in range(100):control.drag(1000*2*math.pi/control.fov,0)
        self.assertAlmostEqual(control.yaw,0,places=10)
        self.assertEqual(control.pitch,0)

    def test_drag_tracks_resolution_and_field_of_view(self):
        first,second=nav.PanoramaControls(),nav.PanoramaControls()
        first.resize(1200,800);second.resize(2400,1600)
        first.drag(100,50);second.drag(200,100)
        self.assertAlmostEqual(first.yaw,second.yaw)
        self.assertAlmostEqual(first.pitch,second.pitch)
        normal,zoomed=nav.PanoramaControls(),nav.PanoramaControls()
        normal.resize(800,800);zoomed.resize(800,800)
        zoomed.wheel(1000);zoomed.advance(10)
        normal.drag(10,10);zoomed.drag(10,10)
        self.assertAlmostEqual(zoomed.yaw,normal.yaw/2)
        self.assertAlmostEqual(zoomed.pitch,normal.pitch/2)

    def test_wheel_is_smooth_monotone_and_rate_independent(self):
        first,second=nav.PanoramaControls(),nav.PanoramaControls()
        opening=first.fov
        for control in (first,second):control.wheel(4)
        self.assertEqual(first.fov,opening)
        previous=first.fov
        for _ in range(6):
            first.advance(1/30)
            self.assertLess(first.fov,previous)
            self.assertGreaterEqual(first.fov,first.target_fov)
            previous=first.fov
        for _ in range(24):second.advance(1/120)
        self.assertAlmostEqual(first.fov,second.fov,places=12)
        first.advance(10)
        self.assertFalse(first.animating)
        self.assertEqual(first.fov,first.target_fov)
        first.wheel(-1)
        self.assertGreater(first.target_fov,first.fov)

    def test_projection_bounds_across_viewport_shapes(self):
        for width,height in ((1000,1000),(1600,900),(400,1200),(5000,50),(0,0)):
            control=nav.PanoramaControls();control.resize(width,height)
            for steps in (-10000,10000):
                control.wheel(steps);control.advance(10)
                self.assertGreaterEqual(control.fov,control.minimum_fov)
                self.assertLessEqual(control.fov,control.maximum_fov)
                self.assertLessEqual(control.fov,nav.MAX_FOV)
                hfov=2*math.atan(math.tan(control.fov/2)*control.width/control.height)
                self.assertLessEqual(hfov,nav.MAX_HORIZONTAL_FOV+1e-12)
            control.resize(20,2000)
            self.assertLessEqual(control.fov,control.maximum_fov)

    def test_stopping_zoom_preserves_displayed_value(self):
        control=nav.PanoramaControls();control.wheel(5);control.advance(.02)
        current=control.fov
        control.stop_zoom()
        self.assertFalse(control.animating)
        self.assertFalse(control.advance(1))
        self.assertEqual(control.fov,current)


if __name__=='__main__':unittest.main()
