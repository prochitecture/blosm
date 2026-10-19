"""Boundary source metadata must preserve identity and leave input geometry intact."""
from copy import deepcopy
import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest

path=Path(__file__).resolve().parents[1]/'way/clipped_context.py'
spec=importlib.util.spec_from_file_location('clipped_context_tests',path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def segment(a,p,b,q):return NS(id1=a,v1=p,id2=b,v2=q)


class ClippedSourceContextTests(unittest.TestCase):
    def setUp(self):
        self.shared=(5.,-1.)
        self.road=NS(source_id='road',category='unclassified',element=NS(tags={'name':'Road'}),
            segments=[segment('a',(5.,5.),'shared',self.shared),segment('shared',self.shared,'b',(5.,-5.))])
        self.aisle=NS(source_id='aisle',category='service',element=NS(tags={'service':'parking_aisle'}),
            segments=[segment('c',(2.,5.),'shared',self.shared)])
        self.context=m.ClippedSourceContext((0.,10.,0.,10.))
        self.t=5./6.
        self.context.observe(self.road.segments[0],('road','a','shared',0.,self.t),((5.,5.),(5.,0.)))
        self.context.observe(self.aisle.segments[0],('aisle','c','shared',0.,self.t),((2.,5.),(4.5,0.)))

    def finalized(self):return self.context.finalize((self.road,self.aisle))

    def street(self,reverse=False):
        points=[(5.,5.,0.),(5.,0.,0.)];interval=(0.,self.t);endpoints=[7,None]
        if reverse:points.reverse();interval=interval[::-1];endpoints.reverse()
        return dict(id=100,endpoints=endpoints,items=[dict(kind='Section',id=10,points=points,
            source_segments=[('road','a','shared',*interval)])])

    def test_full_outside_leg_is_retained_without_changing_original_geometry(self):
        before=deepcopy([(s.id1,s.id2,s.v1,s.v2) for w in (self.road,self.aisle) for s in w.segments])
        report=self.finalized().report();node,=report['nodes']
        self.assertEqual(node['osm_node_id'],'shared')
        self.assertEqual(node['topology_kind'],'JUNCTION')
        self.assertEqual(len(node['incidences']),2)
        self.assertEqual({r['other_node_id'] for r in node['legs']},{'a','b','c'})
        self.assertEqual(before,[(s.id1,s.id2,s.v1,s.v2) for w in (self.road,self.aisle) for s in w.segments])

    def test_reversed_section_resolves_to_the_same_original_node(self):
        raw=self.finalized().report()
        for reverse in (False,True):
            report=m.resolve_open_ends(raw,[self.street(reverse)])
            end,=report['nodes'][0]['open_ends']
            self.assertEqual((end['street_id'],end['section_id'],end['end']),(100,10,0 if reverse else 1))
        self.assertNotIn('open_ends',raw['nodes'][0])

    def test_reversed_accepted_interval_keeps_correct_source_end(self):
        context=m.ClippedSourceContext((0.,10.,0.,10.))
        context.observe(self.road.segments[0],('road','a','shared',self.t,0.),((5.,0.),(5.,5.)))
        context.finalize((self.road,self.aisle))
        cut,=context.report()['nodes'][0]['incidences']
        self.assertEqual(cut['source_end'],1)

    def test_no_proximity_based_connection(self):
        other=deepcopy(self.aisle);other.segments[0].id2='different'
        context=m.ClippedSourceContext((0.,10.,0.,10.))
        context.observe(self.road.segments[0],('road','a','shared',0.,self.t),((5.,5.),(5.,0.)))
        context.observe(other.segments[0],('aisle','c','different',0.,self.t),((2.,5.),(4.5,0.)))
        context.finalize((self.road,other))
        self.assertEqual({n['osm_node_id'] for n in context.report()['nodes']},{'shared','different'})

    def test_moved_source_position_is_reported_not_rebound_silently(self):
        raw=self.finalized().report();street=self.street();street['items'][0]['points'][-1]=(6.,0.,0.)
        report=m.resolve_open_ends(raw,[street])
        self.assertEqual(report['nodes'][0]['open_ends'],[])
        self.assertEqual(report['issues'][0]['code'],'CLIPPED_SOURCE_POSITION_CHANGED')

    def test_reports_own_their_nested_values(self):
        context=self.finalized();report=context.report();report['nodes'][0]['legs'][0]['tags']['name']='changed'
        self.assertNotEqual(context.report(),report)

    def test_conflicting_osm_node_coordinates_reject(self):
        other=deepcopy(self.aisle);other.segments[0].v2=(5.,-2.)
        with self.assertRaises(ValueError):
            self.context.observe(other.segments[0],('aisle','c','shared',0.,self.t),((2.,5.),(4.5,0.)))

    def test_interior_split_does_not_create_clipped_context(self):
        context=m.ClippedSourceContext((0.,10.,0.,10.))
        context.observe(self.road.segments[0],('road','a','shared',.1,.5),((5.,4.4),(5.,2.)))
        context.finalize((self.road,));self.assertEqual(context.report()['nodes'],[])

    def test_source_record_must_match_original_segment(self):
        with self.assertRaises(ValueError):
            self.context.observe(self.road.segments[0],('road','wrong','shared',0.,self.t),((5.,5.),(5.,0.)))

    def test_capture_lifecycle_rejects_further_mutation(self):
        with self.assertRaises(ValueError):self.context.report()
        self.finalized()
        with self.assertRaises(ValueError):self.context.finalize(())
        with self.assertRaises(ValueError):self.context.observe(None,(),())

    def test_existing_junction_is_not_reassigned_to_outside_context(self):
        street=self.street();street['endpoints'][1]=99
        report=m.resolve_open_ends(self.finalized().report(),[street])
        self.assertEqual(report['nodes'][0]['open_ends'],[])


if __name__=='__main__':unittest.main()
