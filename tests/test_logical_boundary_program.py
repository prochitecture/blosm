"""A native chain program must preserve the complete original boundary walk."""
from dataclasses import replace
import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest

from test_shared_port_clusters import fixture, m, ports

spec=importlib.util.spec_from_file_location('tested_logical_boundary_program',
    Path(__file__).resolve().parents[1]/'way/gn/editable_input/curved_logical_chain_bindings.py')
program=importlib.util.module_from_spec(spec);spec.loader.exec_module(program)


class LogicalBoundaryProgramTests(unittest.TestCase):
    def setUp(self):
        data=fixture()
        self.plan=ports.PreparedPhysicalPorts(('residential',)).prepare(data,physical_ids=(10,11,12))
        self.cluster=m.PreparedLogicalPortClusters().prepare(data,self.plan,5000,
            allow_internal_connectors=True,allow_shared_owners=True)
        start=0;ranges=[]
        for site in self.plan.sites:
            ranges.append((site.junction_id,start,len(site.ports)));start+=len(site.ports)
        self.compiled=NS(junction_ranges=ranges)

    def test_three_corners_and_two_spans_remain_separate_steps(self):
        sites,links,steps=program.records(self.plan,self.compiled,self.cluster)
        self.assertEqual(len(sites),3);self.assertEqual(len(links),8)
        last=links[-1];self.assertEqual(last['count'],3)
        chain=steps[last['first']:last['first']+last['count']]
        self.assertEqual([self.plan.spans[r['span']].street_id for r in chain[:-1]],[86,70])
        self.assertEqual(chain[-1]['span'],-1)
        self.assertEqual([r['owner'] for r in chain],[7,7,7])
        self.assertEqual([r['footprint'] for r in chain],[1,1,0])
        self.assertEqual(sum(link['count'] for link in links),len(steps))

    def test_a_jump_to_a_different_original_corner_rejects(self):
        boundary=self.cluster.boundaries[-1]
        corners=(boundary.original_corners[0],boundary.original_corners[2],boundary.original_corners[1])
        cluster=replace(self.cluster,boundaries=(*self.cluster.boundaries[:-1],replace(boundary,original_corners=corners)))
        with self.assertRaisesRegex(ValueError,'original physical corners'):
            program.records(self.plan,self.compiled,cluster)

    def test_reversed_span_without_reversing_its_incidence_rejects(self):
        boundary=self.cluster.boundaries[-1]
        swapped=(boundary.traversed_spans[1],boundary.traversed_spans[0])
        cluster=replace(self.cluster,boundaries=(*self.cluster.boundaries[:-1],replace(boundary,traversed_spans=swapped)))
        with self.assertRaises(ValueError):program.records(self.plan,self.compiled,cluster)

    def test_missing_corner_or_broken_external_cycle_rejects(self):
        boundary=self.cluster.boundaries[-1]
        for changes in (dict(original_corners=boundary.original_corners[:-1]),
                        dict(end_port=boundary.start_port)):
            cluster=replace(self.cluster,boundaries=(*self.cluster.boundaries[:-1],replace(boundary,**changes)))
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                program.records(self.plan,self.compiled,cluster)

    def test_unrepresented_footprint_and_stale_site_binding_reject(self):
        cluster=replace(self.cluster,footprint_spans=(*self.cluster.footprint_spans,(999,9990,9990)))
        with self.assertRaisesRegex(ValueError,'footprint'):program.records(self.plan,self.compiled,cluster)
        ranges=self.compiled.junction_ranges.copy();jid,start,count=ranges[0];ranges[0]=(jid,start,count-1)
        with self.assertRaisesRegex(ValueError,'stale'):
            program.records(self.plan,NS(junction_ranges=ranges),self.cluster)


if __name__=='__main__':unittest.main()
