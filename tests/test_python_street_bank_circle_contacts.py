"""Analytical tangency on both orientations of a Section, with known contacts."""
import math
from types import SimpleNamespace
import unittest
import numpy as np

import _python_street_package
from tested_python_streets.bank_circle_contacts import OriginalCornerArc
from tested_python_streets.bank_circle_preparation import VisibleCorner, fit
from tested_python_streets.finite_free_caps import PreparedCapCurb
from tested_python_streets.model import Key, Section
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.source_bank_reconciliation import SourceBankReconciliationError
from tested_python_streets.corner_radii import CornerContactError
from tested_python_streets.curb_contact import CircularCurbContact
from tested_python_streets.intersections import Approach
from tested_python_streets.local_bank_contacts import _fit_on_intervals
from tested_python_streets.source_bank_reconciliation import NativeBankWitness, fit_switch


class BankCircleContactTests(unittest.TestCase):
    def test_original_circle_and_finite_auxiliary_lead(self):
        from tested_python_streets.arc_lead_joins import fit as fit_lead
        origin=np.array((1520.,-172.))
        circle=SimpleNamespace(center=origin,radius=8.,requested_radius=8.,radius_mode='STYLE',
                               start_angle=1.,sweep=-2.)
        arc=SimpleNamespace(declaration=SimpleNamespace(original_corner=circle),exposed=(.2,.8))
        for circle_first in (True,False):
            sign=1 if circle_first else -1
            tangent=np.array((-sign*math.sqrt(3)/2,-.5))
            contact=origin+np.array((7.,-sign*math.sqrt(3)))
            lead=SimpleNamespace(reference_points=(contact-2*tangent,contact+2*tangent),
                                 prepared_parameters=(.3,1.))
            bank=SimpleNamespace(lead=lead,exposed=(.3,1.),declaration=SimpleNamespace(
                original_corner=SimpleNamespace(radius=2.,requested_radius=2.,radius_mode='STYLE')))
            join,phase,parameter=fit_lead(arc,bank,circle_first)
            np.testing.assert_allclose(join.center,origin+[6.,0.],rtol=0,atol=1e-10)
            self.assertAlmostEqual(join.radius,2.)
            self.assertAlmostEqual(join.sweep,-math.pi/3)
            self.assertAlmostEqual(phase,.5)
            self.assertAlmostEqual(parameter,.65)
            bank.exposed=(.3001,.3002)
            self.assertIsNone(fit_lead(arc,bank,circle_first))

    def visible(self, reverse, angle, origin, phases=(.1,.9)):
        rotation = np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
        def xy(points):
            return np.asarray(points) @ rotation.T + origin
        controls = xy([(-2.,.5),(2.,.5)])
        if reverse:
            controls = controls[::-1]
        section = Section(6,6,tuple((x,y,0.) for x,y in controls),width=1.,interpolation=0)
        path = Path(section)
        corner = SimpleNamespace(center=xy((-2.,-.5)), radius=2.5,
                                 start_angle=angle+math.pi/4, sweep=-math.pi/2)
        theta = corner.start_angle+corner.sweep*np.linspace(0.,1.,17)
        points = corner.center+corner.radius*np.column_stack((np.cos(theta),np.sin(theta)))
        arc = OriginalCornerArc(Key('cycle_corner',8),0,PreparedCapCurb(points),corner)
        visible = VisibleCorner(None,0,arc,phases,SimpleNamespace(path=path,side=-1 if reverse else 1),
                                (1.,3.),(1.,3.))
        return {Key('profile',6):Crown(.08,0.)}, visible, xy

    def test_contacts_on_original_bank_and_circle_after_reversal_and_rotation(self):
        proposal = SimpleNamespace(radius=.5,mode='STYLE')
        for reverse,angle,origin in ((False,0.,(0.,0.)),(True,.65,(900.,-580.))):
            with self.subTest(reverse=reverse):
                context,visible,xy = self.visible(reverse,angle,origin)
                radius,solved = fit(context,visible,proposal)
                self.assertEqual(radius,.5)
                self.assertAlmostEqual(solved['poly_contact'],.5)
                self.assertAlmostEqual(solved['bank_outward_contact'],2.)
                self.assertAlmostEqual(solved['sweep'],-math.pi/2)
                for key,point in [('first',(0.,0.)),('last',(.5,-.5)),('center',(0.,-.5))]:
                    np.testing.assert_allclose(solved[key],xy(point),atol=1e-9,rtol=0.)

    def test_known_tangency_outside_exposed_arc_is_rejected(self):
        context,visible,_ = self.visible(False,0.,(0.,0.),phases=(.55,.95))
        with self.assertRaises(SourceBankReconciliationError):
            fit(context,visible,SimpleNamespace(radius=.5,mode='STYLE'))

    def test_remote_endpoint_direction_does_not_select_the_local_turn(self):
        sections = (Section(1,1,((0.,0.,0.),(30.,-23.,0.),(34.,-26.,0.)),width=2.,interpolation=3),
                    Section(2,2,((0.,0.,0.),(18.,-25.,0.),(32.,-35.,0.)),width=5.,interpolation=3))
        context = {Key('path',s.id):Path(s) for s in sections}
        context.update({Key('profile',s.id):Crown(.08,0.) for s in sections})
        # Use the actual analytical crossing of the two offset banks.
        arms = [Approach(context[Key('path',s.id)],end,context[Key('profile',s.id)])
                for s,end in zip(sections,(0,1))]
        crossing,point,_ = CircularCurbContact(*arms,0.,1.).solve(
            (5.,arms[1].length-5.),np.zeros(2),np.array([arm.length for arm in arms]))
        stations = crossing[0], arms[1].length-crossing[1]
        banks = [NativeBankWitness(arm.path,-1,station,((0.,29.),),(),Key('cycle_strip',s.id),0,'')
                 for arm,station,s in zip(arms,stations,sections)]
        directions = [bank.side*bank.path.frame(bank.station)[1] for bank in banks]
        a,b = directions
        turn = math.atan2(a[0]*b[1]-a[1]*b[0],float(a@b))
        proposal = SimpleNamespace(first=banks[0],last=banks[1],point=point,radius=3.,mode='STRICT')
        class Context(dict):
            network = SimpleNamespace(settings=SimpleNamespace(curve_tolerance=.01))
        context = Context(context)
        with self.assertRaises(CornerContactError):
            fit_switch(context,proposal)
        domains = ((29.,stations[0]),(stations[1],.3))
        result = _fit_on_intervals(context,proposal,domains,turn)
        self.assertEqual(result.corner.radius,3.)
        self.assertEqual(result.corner.radius_attempts,(3.,))
        self.assertLess(result.corner.sweep,0.)
        for arm,sd,cut,side in zip(arms,result.corner.stations,result.intersection,(-1,1)):
            p = arm.point(sd,arm.offset+.5*side*arm.width)
            self.assertAlmostEqual(np.linalg.norm(p-result.corner.center),3.,places=8)
            self.assertLessEqual(np.linalg.norm(arm.point(cut,arm.offset+.5*side*arm.width)-point),1e-8)


if __name__ == '__main__':
    unittest.main()
