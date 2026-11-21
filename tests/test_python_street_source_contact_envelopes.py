"""Original curvature, contact aliases and finite source-phase evaluation."""
from dataclasses import replace
import math
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Key
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.continuations import ContinuationReachError
from tested_python_streets.source_contact_envelopes import (SourceJetContinuation,
    source_curvature_vector,contact_patch,phase_parameters)
from tested_python_streets.compound_displaced_contacts import retry,patch
from test_python_street_displaced_source_contacts import fixture,context


class SourceContactEnvelopeTests(unittest.TestCase):
    def test_original_offset_curvature_vector_survives_reverse_rotation_and_translation(self):
        points=np.array(((0.,0.,0.),(10.,1.,0.),(20.,8.,0.),(35.,8.,0.)))
        path=Path(Section(1,1,tuple(map(tuple,points)),width=4.))
        forward=Approach(Route((Piece(path,True),)),0,Crown(.08,0.))
        backward=Approach(Route((Piece(path,False),)),1,Crown(.08,0.))
        station=.4*path.length;distance=-2.
        expected=source_curvature_vector(forward,station,distance)
        np.testing.assert_allclose(source_curvature_vector(backward,path.length-station,-distance),expected,atol=1e-12)
        # Compare to the geometric normal component of the actual offset jet.
        step=1e-2;p=forward.point(station,distance)
        velocity=(forward.point(station+step,distance)-forward.point(station-step,distance))/(2.*step)
        acceleration=(forward.point(station+step,distance)-2.*p+forward.point(station-step,distance))/step**2
        tangent=velocity/np.linalg.norm(velocity)
        numerical=(acceleration-tangent*(acceleration@tangent))/(velocity@velocity)
        np.testing.assert_allclose(expected,numerical,rtol=1e-5,atol=1e-7)
        angle=.71;rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
        transformed=points.copy();transformed[:,:2]=transformed[:,:2]@rotation.T+(1200.,-800.)
        other=Path(replace(path.section,points=tuple(map(tuple,transformed[::-1]))))
        reversed=Approach(Route((Piece(other,False),)),1,Crown(.08,0.))
        np.testing.assert_allclose(source_curvature_vector(reversed,station,distance),expected@rotation.T,rtol=1e-9,atol=1e-10)
        with self.assertRaises(ValueError):source_curvature_vector(forward,path.length+.1,distance)

    def test_hermite_retains_endpoint_geometry_and_certifies_the_whole_original_radius(self):
        ka,kb=np.array((0.,.005)),np.array((0.,-.005))
        curve=SourceJetContinuation((0.,0.),(30.,0.),(1.,0.),(1.,0.),10.,(ka,kb))
        np.testing.assert_array_equal(curve(np.array((0.,1.))),((0.,0.),(30.,0.)))
        for index,expected in ((0,ka),(-1,kb)):
            velocity=curve.first_derivative[index];acceleration=curve.second_derivative[index]
            tangent=velocity/np.linalg.norm(velocity)
            vector=(acceleration-tangent*(acceleration@tangent))/(velocity@velocity)
            np.testing.assert_allclose(tangent,(1.,0.),atol=1e-12)
            np.testing.assert_allclose(vector,expected,atol=1e-12)
        self.assertLessEqual(curve.curvature_bound,.1*(1.+1e-8))
        self.assertFalse(curve.controls.flags.writeable)
        with self.assertRaises(ContinuationReachError):
            SourceJetContinuation((0.,0.),(30.,0.),(1.,0.),(1.,0.),1000.,(ka,kb))

    def test_explicit_source_patch_and_phase_cuts_retain_original_aliases_and_geometry(self):
        n=fixture();before=n.to_dict();c,a,b=context(n)
        corner=retry(c,a,b,10.,'STYLE',{Key('intersection',10),Key('intersection',11)},{1})
        original=patch(Key('proof',1),corner,.7,.00025)
        index=next(label.index for label in original.segments if label.kind=='curb')
        actual_arc=original.layout.boundaries[index]
        settings=SimpleNamespace(spacing=.7,curve_tolerance=.00025)
        prepared=contact_patch(Key('proof',2),corner,actual_arc,corner.stations,settings,corner.source_contact_walks)
        values=phase_parameters(prepared,corner,corner.stations)
        self.assertEqual(values,phase_parameters(prepared,corner,corner.stations))
        self.assertEqual(set(values),{2,4})
        for approach,station in zip((corner.first,corner.second),corner.stations):
            for (sid,end),contact in approach.path.trims(station).items():
                self.assertIn(contact,values[sid])
                self.assertLessEqual(len(values[sid]),18)
                self.assertTrue(all(0.<=s<=c[Key('path',sid)].length for s in values[sid]))
        self.assertEqual(n.to_dict(),before)
        with self.assertRaises(ValueError):phase_parameters(prepared,corner,corner.stations,maximum_edges=17)


if __name__=='__main__':unittest.main()
