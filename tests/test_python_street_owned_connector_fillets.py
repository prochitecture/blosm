"""Logical source sectors and finite actual contact profiles stay distinct."""
from dataclasses import replace
import unittest
import _python_street_package
from _python_street_owned_connector_fixture import fixture
from tested_python_streets.source_owned_connector_fillets import (
    discover,required,insert,bind_contacts,signature)
from tested_python_streets.generated_junction_curbs import is_generated_junction_area


class OwnedConnectorTests(unittest.TestCase):
    def test_acute_turn_can_contact_the_next_section_beyond_a_consumed_section(self):
        from tested_python_streets.model import End, Intersection, Network, Section, Settings
        from tested_python_streets.backend import PythonMeshBackend
        from tested_python_streets.updates import ProductReady
        from tested_python_streets.successor_bank_contacts import validate_exposure

        origin, fork, far = (0.,0.), (7.5,-5.01), (36.65,-19.8)
        controls = (
            (3.5, (origin,(3.37,-5.66),(7.06,-9.71),(11.75,-12.99),(18.94,-15.97),(28.41,-18.28),far)),
            (2., (origin,(4.71,-2.83),fork)),
            (2., (fork,(11.78,-7.2),(15.83,-8.75),(20.82,-9.78),(27.19,-9.69))),
            (2., ((9.37,2.38),(8.43,-1.73),fork)),
            (7., ((-14.08,20.86),origin)),
            (7., (far,(116.68,-5.49))),
            (3.5, ((6.73,-28.49),(12.81,-24.85),(22.42,-22.06),far)))
        sections = tuple(Section(i,i,tuple((x,y,0.) for x,y in points),width=width)
                         for i,(width,points) in enumerate(controls,1))
        sites = (
            Intersection(1,(End(5,1),End(1,0),End(2,0)),2.,(2.,2.,1.),('STYLE',)*3),
            Intersection(2,(End(4,1),End(2,1),End(3,0)),1.,(1.,)*3,('STYLE',)*3),
            Intersection(3,(End(6,0),End(1,1),End(7,1)),2.,(2.,)*3,('STYLE',)*3))
        network = Network('short-fork',sections,sites,settings=Settings(exaggeration=1.))
        steps = PythonMeshBackend(network.settings).prepare_steps(network,emit_products=True)
        try:
            plan = next(event.product for event in steps
                        if isinstance(event,ProductReady) and event.key.kind == 'compound_plan')
        finally:
            steps.close()
        group = plan.by_site[sites[0].key]
        record, = group.physical_layout.source_owned_connector_fillets
        self.assertEqual(tuple(c.path.section.id for c in record.contacts),(1,3))
        self.assertAlmostEqual(record.corner.radius,.64)
        validate_exposure(group)

    def test_real_connector_reserves_actual_wide_native_contact_without_changing_logical_corner(self):
        for reverse,angle in ((False,0.),(True,.7)):
            planner,group,grids=fixture(reverse=reverse,angle=angle)
            found=discover(planner,group);self.assertEqual(len(found),1)
            proposal,logical,walk,*_=found[0]
            self.assertEqual(proposal.first.path.section.id,11)
            self.assertEqual(proposal.last.path.section.id,14)
            self.assertNotEqual(walk.original_origin.site,walk.contact_origin.site)
            self.assertEqual(walk.original_origin.endpoint.section_id,12)
            self.assertEqual(walk.contact_origin.endpoint.section_id,14)
            record=required(planner,group,());self.assertEqual(record.corner.radius,2.)
            self.assertEqual(tuple(v.path.section.width for v in record.contacts),(5.4,10.8))
            self.assertFalse(any(v.owned_intervals for v in record.contacts))
            self.assertEqual(tuple(p.section.id for p in record.connector_walk.pieces),(12,))
            sources=[];stations={};old=group.physical_layout.corners
            insert(planner,group.members,{12,13},sources,stations,(record,),.7,.00025)
            self.assertTrue(is_generated_junction_area(sources[0]))
            self.assertIs(group.physical_layout.corners,old)
            self.assertNotIn(record.corner,old)
            self.assertEqual(len(stations),2)

    def test_unconsumed_or_scoped_connector_cannot_authorize_a_successor(self):
        planner,group,_=fixture(consume=False)
        self.assertFalse(discover(planner,group));self.assertIsNone(required(planner,group,()))
        # A scope handoff cannot simultaneously be a real junction incidence.
        with self.assertRaisesRegex(ValueError,'scope boundary'):
            fixture(scope=True)

    def test_through_site_strict_radius_participates_in_the_original_policy(self):
        planner,group,_=fixture(through_radius=4.,through_mode='STRICT')
        self.assertFalse(discover(planner,group))
        planner,group,_=fixture(through_radius=2.,through_mode='STRICT')
        record=required(planner,group,());self.assertIsNotNone(record)
        self.assertEqual(record.proposal.mode,'STRICT');self.assertEqual(record.corner.radius,2.)

    def test_reservation_cannot_publish_against_its_old_short_mouths(self):
        planner,group,_=fixture();record=required(planner,group,())
        sources=[];insert(planner,group.members,{12,13},sources,{},(record,),.7,.00025)
        with self.assertRaisesRegex(ValueError,'final finite native ownership'):
            bind_contacts(planner.context,replace(group,footprints=(*group.footprints,*sources)))
        with self.assertRaisesRegex(ValueError,'failed to remove'):
            required(planner,group,(record,))

    def test_unchanged_groups_have_no_new_dependency_signature(self):
        planner,group,_=fixture()
        self.assertIsNone(signature(group))
        self.assertIs(bind_contacts(planner.context,group),group)


if __name__=='__main__':unittest.main()
