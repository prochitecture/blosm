"""Original spans must survive multi-junction ownership and contraction."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
import unittest
from test_physical_ports import fixture, m as ports
from test_physical_ports import scoped_fixture
from test_port_compounds import m as compounds

spec=importlib.util.spec_from_file_location('tested_port_network',Path(__file__).resolve().parents[1]/'way/street_network/port_network.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)


class PortNetworkTests(unittest.TestCase):
    def setUp(self):
        self.reader=m.PreparedPortNetwork(compounds.PreparedPortCompound())
        self.plan=ports.PreparedPhysicalPorts(('residential',)).prepare(fixture(),physical_ids=(1001,1002),interior_ids=(21,22))

    def test_category_handoff_is_not_a_free_end_or_unresolved_road(self):
        plan=ports.PreparedPhysicalPorts(('residential',),section_scopes=True).prepare(
            scoped_fixture(),physical_ids=(1001,1002),interior_ids=(21,22))
        result=self.reader.prepare(plan)
        road,path=[r for r in result.routes if r.span.street_id==4]
        self.assertEqual((road.kind,path.kind),('CATEGORY_APPROACH','DEFERRED_CATEGORY'))
        self.assertFalse(result.unresolved)
        self.assertEqual(road.span.end.owner_kind,'CATEGORY_BOUNDARY')

    def test_complete_source_street_has_three_distinct_corridors(self):
        result=self.reader.prepare(self.plan)
        self.assertEqual([r.span.identity for r in result.routes if r.kind=='CORRIDOR'],[(1,10,10),(1,11,12),(1,13,13)])
        self.assertEqual(sum(r.kind=='APPROACH' for r in result.routes),6)
        self.assertFalse(result.unresolved)
        self.assertEqual(tuple(r.span for r in result.routes),self.plan.spans)

    def test_compound_keeps_both_other_spans_on_same_street(self):
        result=self.reader.prepare(self.plan,compound_spans=((1,11,12),))
        self.assertEqual(result.ordinary_site_ids,(1001,1002))
        self.assertEqual([r.span.identity for r in result.routes if r.kind=='COMPOUND'],[(1,11,12)])
        self.assertEqual([r.span.identity for r in result.routes if r.kind=='CORRIDOR'],[(1,10,10),(1,13,13)])

    def test_two_disjoint_compounds_leave_the_middle_corridor(self):
        result=self.reader.prepare(self.plan,compound_spans=((1,10,10),(1,13,13)))
        self.assertFalse(result.ordinary_site_ids)
        self.assertEqual([r.span.identity for r in result.routes if r.kind=='CORRIDOR'],[(1,11,12)])

    def test_unselected_end_is_never_made_free(self):
        plan=ports.PreparedPhysicalPorts(('residential',)).prepare(fixture(),interior_ids=(21,))
        result=self.reader.prepare(plan)
        kinds={r.span.identity:r.kind for r in result.routes}
        self.assertEqual(kinds[(1,10,10)],'REMOTE_CONNECTION')
        self.assertEqual(kinds[(1,11,12)],'REMOTE_CONNECTION')
        self.assertEqual(kinds[(1,13,13)],'UNSELECTED')
        self.assertEqual(kinds[(2,20,20)],'APPROACH')

    def test_overlapping_compounds_reject(self):
        for ids in (((1,10,10),(1,11,12)),((1,11,12),(1,11,12))):
            with self.assertRaisesRegex(ValueError,'overlap'):self.reader.prepare(self.plan,compound_spans=ids)

    def test_missing_selected_incidence_rejects(self):
        site=self.plan.sites[0]
        broken=replace(self.plan,sites=(replace(site,ports=site.ports[1:]),*self.plan.sites[1:]))
        with self.assertRaisesRegex(ValueError,'missing'):self.reader.prepare(broken)

    def test_conflicting_endpoint_owner_rejects(self):
        span=self.plan.spans[0]
        broken=replace(self.plan,spans=(replace(span,start=replace(span.start,junction_id=999)),*self.plan.spans[1:]))
        with self.assertRaisesRegex(ValueError,'different junction'):self.reader.prepare(broken)

    def logical_fixture(self):
        from test_port_clusters import fixture as logical_fixture, m as logical
        data=logical_fixture()
        plan=ports.PreparedPhysicalPorts(('residential',)).prepare(data,physical_ids=range(10,14))
        return plan,logical.PreparedLogicalPortClusters().prepare(data,plan,5000)

    def test_logical_surface_owns_internal_cycle_once(self):
        plan,cluster=self.logical_fixture()
        network=self.reader.prepare(plan,logical_clusters=(cluster,))
        self.assertFalse(network.ordinary_site_ids)
        self.assertFalse(network.compounds)
        self.assertEqual(network.logical_clusters,(cluster,))
        self.assertEqual(tuple(r.span for r in network.routes),plan.spans)
        self.assertEqual(sum(r.kind=='LOGICAL_CLUSTER' for r in network.routes),4)
        self.assertEqual(sum(r.kind=='APPROACH' for r in network.routes),8)
        self.assertFalse(network.unresolved)

    def test_logical_duplicate_owner_rejects(self):
        plan,cluster=self.logical_fixture()
        with self.assertRaisesRegex(ValueError,'overlap'):
            self.reader.prepare(plan,logical_clusters=(cluster,cluster))

    def test_logical_missing_or_changed_internal_span_rejects(self):
        plan,cluster=self.logical_fixture()
        for spans in (cluster.internal_spans[:-1],(replace(cluster.internal_spans[0],
            first_control=100),*cluster.internal_spans[1:])):
            with self.assertRaisesRegex(ValueError,'exact original internal span'):
                self.reader.prepare(plan,logical_clusters=(replace(cluster,internal_spans=spans),))

    def test_logical_missing_or_internal_external_port_rejects(self):
        plan,cluster=self.logical_fixture()
        for ring in (cluster.external_ring[:-1],(plan.sites[0].ports[1],*cluster.external_ring[1:])):
            with self.assertRaisesRegex(ValueError,'exact original external port'):
                self.reader.prepare(plan,logical_clusters=(replace(cluster,external_ring=ring),))


if __name__=='__main__':unittest.main()
