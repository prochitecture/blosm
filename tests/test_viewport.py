"""Pure layout/input and reference spherical-projection checks."""
import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    result=importlib.util.module_from_spec(spec)
    sys.modules[name]=result
    spec.loader.exec_module(result)
    return result
ui=module('viewport_ui_test','viewport/ui.py')
pano=module('panorama_math_test','mapillary/panorama_math.py')


class ViewportTests(unittest.TestCase):
    def test_direction_landmarks(self):
        for direction,expected in [((0,0,-1),(.5,.5)),((1,0,0),(.75,.5)),((-1,0,0),(.25,.5)),
                                   ((0,1,0),(.5,1)),((0,-1,0),(.5,0))]:
            for actual,value in zip(pano.image_uv(direction),expected):self.assertAlmostEqual(actual,value)

    def test_exif_uses_project_then_pixel_rotation(self):
        # Independently construct basic camera direction from a raw SfM bearing.
        raw=(.3,.4,.5)
        raw_uv=(.5+math.atan2(raw[0],raw[2])/(2*math.pi),.5+math.atan2(raw[1],math.hypot(raw[0],raw[2]))/math.pi)
        for orientation,angle in [(1,0),(3,math.pi),(6,math.pi/2),(8,3*math.pi/2)]:
            c,s=math.cos(angle),math.sin(angle)
            x,y,z=raw
            local=(c*x-s*y,-(s*x+c*y),-z)
            u,v=raw_uv
            expected={1:(u,1-v),3:(1-u,v),6:(1-v,1-u),8:(v,u)}[orientation]
            for actual,value in zip(pano.image_uv(local,orientation),expected):self.assertAlmostEqual(actual,value)

    def test_seam_and_layout_validation(self):
        self.assertLess(pano.image_uv((1e-8,0,1))[0],1)
        self.assertGreaterEqual(pano.image_uv((-1e-8,0,1))[0],0)
        self.assertTrue(pano.valid_image_size(1024,2048,6))
        self.assertFalse(pano.valid_image_size(2048,1024,6))
        with self.assertRaises(ValueError):pano.image_uv((0,0,-1),2)

    def test_responsive_layout(self):
        for width,height in [(1200,800),(500,300),(140,100),(20,20)]:
            toolbar,card=ui.shell_layout(ui.Rect(0,0,width,height))
            self.assertLessEqual(toolbar.x+toolbar.width,width)
            self.assertLessEqual(toolbar.y+toolbar.height,height)
            if width<850:self.assertIsNone(card)

    def test_capture_focus_and_composition(self):
        values=[]
        button=ui.Button('button','Run',lambda:values.append('click'))
        button.rect=ui.Rect(0,0,50,30)
        root=ui.Container('root',[button]);root.rect=ui.Rect(0,0,100,100)
        router=ui.Router(root)
        router.handle(ui.InputEvent('LEFTMOUSE','PRESS',10,10))
        router.handle(ui.InputEvent('LEFTMOUSE','RELEASE',90,90))
        self.assertEqual(values,[])
        router.handle(ui.InputEvent('TAB','PRESS'))
        router.handle(ui.InputEvent('RET','PRESS'))
        self.assertEqual(values,['click'])
        button.visible=False;router.reset_capture()
        self.assertIsNone(router.focus)
        self.assertIsNone(root.hit(10,10))

    def test_slider_outside_drag_and_resize(self):
        values=[.5]
        slider=ui.Slider('opacity','Opacity',lambda:values[-1],values.append)
        slider.rect=ui.Rect(10,10,100,30)
        root=ui.Container('root',[slider]);root.rect=ui.Rect(0,0,200,100)
        router=ui.Router(root)
        router.handle(ui.InputEvent('LEFTMOUSE','PRESS',50,20))
        router.handle(ui.InputEvent('MOUSEMOVE','NOTHING',300,20))
        self.assertEqual(values[-1],1.)
        router.reset_capture()
        router.handle(ui.InputEvent('MOUSEMOVE','NOTHING',-50,20))
        self.assertEqual(values[-1],1.)


if __name__=='__main__': unittest.main()
