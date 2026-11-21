"""Reverse Route round trips retain exact native source mouth geometry."""
from copy import copy
import unittest
from unittest import mock
import numpy as np
import _python_street_package
from tested_python_streets.model import Section
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.ordinary_singleton_sources import _same_mouth_representation


def fixture():
    first=np.array((-94.82969665527344,-95.42214965820312))
    last=np.array((-89.46419525146484,-69.71858215332031))
    direction=(last-first)/np.linalg.norm(last-first)
    controls=[(first,last)];last=first
    for length in (20.,10.,30.):
        first=last-length*direction;controls.append((first,last));last=first
    pieces=[]
    for sid,points in enumerate(controls,1):
        section=Section(sid,sid,tuple((x,y,0.) for x,y in points),width=5.4,interpolation=3)
        pieces.append(Piece(Path(section),False))
    route=Route(pieces);arm=Approach(route,1,Crown(.08,0.))
    s=4.480662966425596
    t=pieces[0].path.length-(pieces[0].path.length-s)
    return arm,s,t


class MouthRepresentationTests(unittest.TestCase):
    def test_reverse_multi_piece_roundtrip_keeps_every_frame_and_jacobian(self):
        a,s,t=fixture();b=copy(a)
        self.assertEqual(len(a.path.pieces),4);self.assertNotEqual(s,t)
        self.assertEqual(a.path.trims(s),b.path.trims(t))
        self.assertTrue(_same_mouth_representation(a,s,b,t))
        for x,y in zip(a.frame(s),b.frame(t)):np.testing.assert_array_equal(x,y)
        q=np.linspace(-1.,1.,65);distance=a.offset+.5*a.width*q
        for x,y in zip(a.jet(s,distance),b.jet(t,distance)):np.testing.assert_array_equal(x,y)
        np.testing.assert_array_equal(a.physical_gradient(s,q),b.physical_gradient(t,q))
        # The original exact-parameter branch still does no additional sampling.
        with mock.patch.object(a,'frame',side_effect=AssertionError('old exact path')):
            self.assertTrue(_same_mouth_representation(a,s,a,s))

    def test_real_station_identity_orientation_and_profile_changes_are_rejected(self):
        a,s,t=fixture();b=copy(a)
        self.assertFalse(_same_mouth_representation(a,s,b,t+.001))
        self.assertFalse(_same_mouth_representation(a,s,b,a.length+1e-12))
        b.path=Route(a.path.pieces)
        self.assertFalse(_same_mouth_representation(a,s,b,t))
        b=copy(a);b.profile=Crown(a.profile.rise,a.profile.datum,a.profile.half)
        self.assertEqual(b.profile,a.profile);self.assertIsNot(b.profile,a.profile)
        self.assertFalse(_same_mouth_representation(a,s,b,t))
        b=copy(a);b.end=0
        self.assertFalse(_same_mouth_representation(a,s,b,t))
        b=copy(a);b.width+=.01
        self.assertFalse(_same_mouth_representation(a,s,b,t))

    def test_non_equal_center_or_derivative_frame_and_custom_gradient_are_rejected(self):
        for at_derivative in (False,True):
            a,s,t=fixture();b=copy(a);original=a.path.frame
            witness=min(t+1e-4,a.length) if at_derivative else t
            def changed(station):
                p,u,n=original(station)
                if float(station)==witness:p=p+np.array((1e-6,0.))
                return p,u,n
            with mock.patch.object(a.path,'frame',side_effect=changed):
                self.assertFalse(_same_mouth_representation(a,s,b,t))
        a,s,t=fixture();b=copy(a)
        with mock.patch.object(b,'physical_gradient',return_value=np.zeros((5,2))):
            self.assertFalse(_same_mouth_representation(a,s,b,t))


if __name__=='__main__':unittest.main()
