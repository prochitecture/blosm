"""Canonical source endpoints and solver views retain their distinct roles."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
import unittest

from test_physical_ports import fixture, m as ports
from test_physical_ports import scoped_fixture

spec=importlib.util.spec_from_file_location('tested_port_source_layout',
    Path(__file__).resolve().parents[1]/'way/street_network/port_source_layout.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)


class PortSourceLayoutTests(unittest.TestCase):
    def setUp(self):
        self.source=fixture()
        self.plan=ports.PreparedPhysicalPorts(('residential',)).prepare(self.source,interior_ids=(21,22))

    def test_excluded_section_controls_remain_in_the_same_editable_source(self):
        source=scoped_fixture()
        plan=ports.PreparedPhysicalPorts(('residential',),section_scopes=True).prepare(source,physical_ids=(1001,))
        layout=m.PortSourceLayout.prepare(source,plan)
        self.assertIn(40,layout.section_ids);self.assertIn(41,layout.section_ids)
        used=[r.vertex for r in layout.controls]+[r.alias for r in layout.controls if r.alias>=0]
        self.assertEqual(sorted(used),list(range(layout.source_vertex_count)))
        self.assertEqual(sum(e.section_id==41 for e in layout.edges),1)
        source.streets[4]['items'][-1]['category']='residential'
        with self.assertRaises(ValueError):m.PortSourceLayout.prepare(source,plan)

    def test_complete_sections_keep_both_editable_endpoint_roles(self):
        before=deepcopy(self.source);layout=m.PortSourceLayout.prepare(self.source,self.plan)
        self.assertEqual(self.source,before)
        self.assertEqual(layout.section_ids,(10,11,12,13,20,30))
        self.assertEqual((layout.source_vertex_count,layout.source_edge_count),(14,8))
        self.assertEqual(len(layout.controls),11)
        self.assertEqual([(r.vertex,r.alias) for r in layout.controls if r.alias>=0],[(2,3),(4,5),(6,7)])
        # Each authored point is used exactly once, including both sides of
        # each Section seam. Different Streets are never welded by location.
        used=[r.vertex for r in layout.controls]+[r.alias for r in layout.controls if r.alias>=0]
        self.assertEqual(sorted(used),list(range(14)))
        self.assertEqual([(e.first,e.last) for e in layout.edges[:6]],[(0,1),(1,2),(3,4),(5,6),(7,8),(8,9)])

    def test_object_grouping_and_geometry_do_not_change_private_identity(self):
        original=m.PortSourceLayout.prepare(self.source,self.plan)
        for street in self.source.streets.values():street['object_group']='one combined object'
        self.source.streets[1]['items'][0]['points']=((-40.,0.),(-20.,2.),(0.,0.))
        self.assertEqual(original,m.PortSourceLayout.prepare(self.source,self.plan))

    def test_stale_section_order_and_counts_reject(self):
        self.source.streets[1]['items'][0]['points']+=((1.,1.),)
        with self.assertRaisesRegex(ValueError,'topology changed'):
            m.PortSourceLayout.prepare(self.source,self.plan)
        self.source=fixture()
        self.source.streets[1]['items'][0]['id']=999
        with self.assertRaisesRegex(ValueError,'Section order changed'):
            m.PortSourceLayout.prepare(self.source,self.plan)

    def test_same_total_count_cannot_hide_changed_section_membership(self):
        spans=list(self.plan.spans)
        at=next(i for i,s in enumerate(spans) if s.street_id==1 and len(s.section_ids)==2)
        spans[at]=replace(spans[at],section_ids=tuple(reversed(spans[at].section_ids)))
        with self.assertRaisesRegex(ValueError,'Section order changed'):
            m.PortSourceLayout.prepare(self.source,replace(self.plan,spans=tuple(spans)))

    def test_same_total_count_cannot_hide_changed_span_seam(self):
        items=self.source.streets[1]['items']
        items[0]['points']=items[0]['points'][:2]
        items[2]['points']=((0.,0.),(7.,.5),(15.,1.))
        with self.assertRaisesRegex(ValueError,'Section seam changed'):
            m.PortSourceLayout.prepare(self.source,self.plan)


if __name__=='__main__':unittest.main()
