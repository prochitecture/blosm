"""Coincident lead contacts use source identity and finite station domains."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Network,Key
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.routes import Piece,Route
from tested_python_streets.intersections import Approach
from tested_python_streets.coincident_source_curbs import contact


class LeadCurve:
    def __init__(self,approach,offset):
        length=approach.length
        self.lead=(approach,np.array((0.,length)),np.array((0.,length)),offset)
        self.leads=(self.lead,self.lead)
        self.length=length+1.
        self.lengths=np.array((0.,length,length+1.,length+1.))
    def __call__(self,t):
        a,parameters,cumulative,offset=self.lead
        return a.point(np.interp(np.asarray(t)*self.length,cumulative,parameters),offset)


def fixture(second_source=1,second_offset=-.5):
    n=Network('coincident_curbs',(Section(1,1,((0.,0.,0.),(20.,0.,0.))),
        Section(2,2,((0.,0.,0.),(20.,0.,0.)))))
    pat=PatternLibrary();p=UpdateCoordinator((PathProvider(),ProfileProvider()),pat).prepare(n).products
    c=BuildContext(n,p,pat)
    a=Approach(Route((Piece(c[Key('path',1)],True),)),0,c[Key('profile',1)])
    b=Approach(Route((Piece(c[Key('path',second_source)],False),)),1,c[Key('profile',second_source)])
    return LeadCurve(a,.5),LeadCurve(b,second_offset)


class CoincidentCurbTests(unittest.TestCase):
    def test_reverse_incidence_finds_a_common_original_source_station(self):
        curves=fixture();initial=np.array((.4,.6));lower=np.array((.2,.2));upper=np.array((.8,.8))
        result=contact(curves,initial,lower,upper)
        self.assertIsNotNone(result)
        self.assertTrue(np.all(result>=lower));self.assertTrue(np.all(result<=upper))
        np.testing.assert_allclose(curves[0](result[0]),curves[1](result[1]),atol=1e-12,rtol=0.)

    def test_equal_coordinates_other_offsets_and_disjoint_domains_do_not_connect(self):
        initial=np.array((.4,.6));lower=np.array((.2,.2));upper=np.array((.8,.8))
        self.assertIsNone(contact(fixture(second_source=2),initial,lower,upper))
        self.assertIsNone(contact(fixture(second_offset=-1.),initial,lower,upper))
        self.assertIsNone(contact(fixture(),initial,np.array((.05,.05)),np.array((.1,.1))))


if __name__=='__main__':unittest.main()
