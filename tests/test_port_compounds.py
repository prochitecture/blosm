"""Compound topology must address spans rather than whole-Street IDs."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
import unittest

from test_physical_ports import fixture, section, m as ports

spec = importlib.util.spec_from_file_location('tested_port_compounds',
    Path(__file__).resolve().parents[1]/'way/street_network/port_compounds.py')
m = importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)


class PortCompoundTests(unittest.TestCase):
    def setUp(self):
        self.snapshot=fixture()
        self.plan=ports.PreparedPhysicalPorts(('residential',)).prepare(
            self.snapshot, physical_ids=(1001,1002), interior_ids=(21,22))
        self.reader=m.PreparedPortCompound()

    def test_physical_to_interior_does_not_consume_remaining_same_street(self):
        before=deepcopy(self.plan)
        compound=self.reader.prepare(self.plan,(1,10,10))
        self.assertEqual(compound.site_ids,(1001,21))
        self.assertEqual(len(compound.external_ports),4)
        remaining=[p for p in compound.external_ports if p.endpoint.street_id==1]
        self.assertEqual([p.endpoint for p in remaining],[ports.SectionEnd(1,11,0)])
        self.assertEqual({p.endpoint for p in compound.internal_ports},
                         {ports.SectionEnd(1,10,0),ports.SectionEnd(1,10,1)})
        self.assertEqual([b[0] for b in compound.bridges],[1,3])
        self.assertEqual(self.plan,before)

    def test_two_internal_sites_retain_two_external_spans_of_same_street(self):
        compound=self.reader.prepare(self.plan,(1,11,12))
        same=[p.endpoint for p in compound.external_ports if p.endpoint.street_id==1]
        self.assertEqual(set(same),{ports.SectionEnd(1,10,1),ports.SectionEnd(1,13,0)})
        self.assertEqual(len(compound.external_ports),4)

    def test_unselected_owner_and_missing_span_reject(self):
        with self.assertRaises(ValueError):self.reader.prepare(self.plan,(1,999,999))
        removed=replace(self.plan,sites=tuple(s for s in self.plan.sites if s.junction_id!=21))
        with self.assertRaisesRegex(ValueError,'Both original'):self.reader.prepare(removed,(1,10,10))

    def test_changed_positions_and_grouping_preserve_plan(self):
        expected=self.reader.prepare(self.plan,(1,11,12))
        self.snapshot.streets[1]['items'][0]['points']=((-40.,0.),(-20.,1.3),(0.,0.))
        self.snapshot.streets[1]['object']='One complete Street'
        plan=ports.PreparedPhysicalPorts(('residential',)).prepare(
            self.snapshot,physical_ids=(1001,1002),interior_ids=(21,22))
        self.assertEqual(expected,self.reader.prepare(plan,(1,11,12)))

    def test_second_shared_span_is_not_silently_absorbed(self):
        source=self.snapshot.streets[1]
        for item in source['items']:
            if item.get('id') in (21,22) and item['kind']=='Intersection':
                item['attachments'].append(dict(street_id=8,section_id=80,
                    end=0 if item['id']==21 else 1,connector_id=800+item['id'],category='residential'))
        self.snapshot.streets[8]=dict(id=8,status='SUPPORTED',bundle_id=None,endpoints=(21,22),
            items=[section(80,((0.,0.),(20.,-15.),(40.,0.)))])
        plan=ports.PreparedPhysicalPorts(('residential',)).prepare(self.snapshot,interior_ids=(21,22))
        with self.assertRaisesRegex(ValueError,'second shared span'):
            self.reader.prepare(plan,(1,11,12))

    def test_site_storage_order_does_not_reorder_boundary(self):
        expected=self.reader.prepare(self.plan,(1,11,12))
        reordered=replace(self.plan,sites=tuple(reversed(self.plan.sites)))
        self.assertEqual(self.reader.prepare(reordered,(1,11,12)),expected)


if __name__=='__main__':unittest.main()
