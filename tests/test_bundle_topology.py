"""Bundle addressing regressions; usable with ordinary Python, without Blender."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

# Load the Blender-independent module without importing the add-on's bpy entry.
name='tested_bundle_topology'
spec=importlib.util.spec_from_file_location(name,Path(__file__).resolve().parents[1]/'way/street_network/bundle_topology.py')
module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
PreparedBundleTopology=module.PreparedBundleTopology
BundleTopologyError=module.BundleTopologyError


def section(sid, points):
    return dict(kind='Section',id=sid,points=tuple(tuple(p)+(0.,) for p in points))


def sample(reverse=False, junction=True):
    first=section(10,((0.,0.),(10.,0.),(20.,0.)))
    second=section(11,((20.,0.),(30.,0.),(40.,0.)))
    opposite=((0.,10.),(20.,10.),(40.,10.))
    streets={1:dict(id=1,items=[first,second]),
             2:dict(id=2,items=[section(20,tuple(reversed(opposite)) if reverse else opposite)]),
             3:dict(id=3,items=[section(30,((20.,0.),(20.,-20.)))])}
    if junction:
        streets[1]['items'].insert(1,dict(kind='Intersection',id=8,location=(20.,0.),
            surface_treatment='NATIVE_JUNCTION',through=((10,1),(11,0)),attachments=[
                dict(street_id=3,section_id=30,end=0,connector_id=17,category='service',side='RIGHT')]))
    bundle=dict(id=7,head=(1,2),tail=(1,2),head_locations=((0.,0.),(0.,10.)),
                tail_locations=((40.,0.),(40.,10.)),owned_streets=(1,2))
    return bundle,streets


def chained(reverse=False):
    paths = {1: ((0., 10.), (20., 10.)), 3: ((20., 10.), (40., 10.)),
             2: ((0., 0.), (40., 0.))}
    streets = {sid: dict(id=sid, bundle_id=7,
        items=[section(sid*10, tuple(reversed(points)) if reverse else points)])
        for sid, points in paths.items()}
    bundle = dict(id=7, head=(1, 2), tail=(3, 2), owned_streets=(1, 2, 3),
        head_locations=((0., 10.), (0., 0.)), tail_locations=((40., 10.), (40., 0.)))
    ring = [dict(connector_id=sid, item_kind='Street', item_id=sid,
        endpoints=[dict(street_id=sid, terminal_id=sid*10, terminal_kind='Section',
                        end=(1-end if reverse else end))]) for sid, end in ((1, 1), (3, 0))]
    origins = dict(physical=[dict(junction_id=8, location=(20., 10.), ring=ring)])
    return bundle, streets, origins


class BundleTopologyTests(unittest.TestCase):
    def setUp(self):
        self.reader=PreparedBundleTopology(('residential','service'))

    def reject(self,bundle,streets,code):
        with self.assertRaises(BundleTopologyError) as cm:self.reader.read(bundle,streets)
        self.assertEqual(cm.exception.code,code)

    def test_internal_junction_partitions_only_its_own_member(self):
        bundle,streets=sample();original=deepcopy((bundle,streets))
        result=self.reader.read(bundle,streets)
        self.assertEqual((bundle,streets),original)
        self.assertEqual([m.street_id for m in result.members],[1,2])
        self.assertEqual([len(m.spans) for m in result.members],[2,1])
        self.assertEqual(result.members[0].spans[0].end.junction_id,8)
        self.assertEqual(result.members[0].spans[1].start.junction_id,8)
        self.assertEqual(result.junctions[0].attachments[0].connector_id,17)
        self.assertTrue(result.junctions[0].attachments[0].carriageway)
        self.assertEqual(result.junctions[0].attachments[0].side,'RIGHT')
        self.assertEqual(result.required_features,('INTERNAL_JUNCTIONS',))
        self.assertEqual({r.start.section_id for m in result.members for s in m.spans for r in s.sections},{10,11,20})
        json.dumps(result.report())

    def test_longitudinal_chain_retains_streets_sections_and_lateral_roles(self):
        for reverse in (False, True):
            bundle, streets, origins = chained(reverse)
            before = deepcopy((bundle, streets, origins))
            result = self.reader.read(bundle, streets, origins=origins)
            self.assertEqual(result.chains, ((1, 3), (2,)))
            self.assertEqual([m.street_id for m in result.members], [1, 3, 2])
            self.assertEqual([m.head_order for m in result.members], [0, 0, 1])
            self.assertTrue(all(m.forward != reverse and m.complete for m in result.members))
            a, b, c = result.members
            self.assertEqual(a.spans[0].start, module.SpanBoundary('HEAD'))
            self.assertEqual(a.spans[0].end, module.SpanBoundary(junction_id=8))
            self.assertEqual(b.spans[0].start, module.SpanBoundary(junction_id=8))
            self.assertEqual(b.spans[0].end, module.SpanBoundary('TAIL'))
            self.assertEqual(c.spans[0].start, module.SpanBoundary('HEAD'))
            self.assertEqual(result.chain_joins[0], module.ChainJoin(8,
                (1, 10, int(not reverse)), (3, 30, int(reverse))))
            self.assertEqual(result.required_features, ('LONGITUDINAL_CHAINS',))
            self.assertEqual((bundle, streets, origins), before)
            json.dumps(result.report())

    def test_chain_join_requires_original_physical_incidence_not_proximity(self):
        for change in ('missing', 'duplicate', 'terminal', 'location', 'branch', 'direction', 'owner', 'order', 'aggregate'):
            bundle, streets, origins = chained()
            site = origins['physical'][0]
            if change == 'missing': origins['physical'] = []
            elif change == 'duplicate': site['ring'].append(deepcopy(site['ring'][0]))
            elif change == 'terminal': site['ring'][0]['endpoints'][0]['terminal_id'] = 999
            elif change == 'location': site['location'] = (20., .01)
            elif change == 'branch':
                streets[4] = dict(id=4, bundle_id=7, items=[section(40, ((20., 10.), (40., 20.)))])
                bundle['owned_streets'] += (4,)
                site['ring'].append(dict(item_kind='Street', endpoints=[dict(street_id=4, terminal_id=40,
                    terminal_kind='Section', end=0)]))
            elif change == 'direction':
                streets[3]['items'][0]['points'] = tuple(reversed(streets[3]['items'][0]['points']))
                site['ring'][1]['endpoints'][0]['end'] = 1
            elif change == 'owner': streets[3]['bundle_id'] = 99
            elif change == 'aggregate': site['ring'][0]['item_kind'] = 'Bundle'
            else:
                bundle['tail'] = bundle['tail'][::-1]
                bundle['tail_locations'] = bundle['tail_locations'][::-1]
            with self.subTest(change=change), self.assertRaises(BundleTopologyError):
                self.reader.read(bundle, streets, origins=origins)

    def test_chain_mapping_apply_consumes_saved_original_sites(self):
        bundle, streets, origins = chained(); issues = []
        snapshot = SimpleNamespace(bundles={7: bundle}, streets=streets, junction_origins=origins,
                                   issue=lambda *args: issues.append(args))
        self.reader.apply(snapshot)
        self.assertEqual(bundle['topology']['chains'], ((1, 3), (2,)))
        self.assertFalse(issues)

    def test_reverse_authored_member_preserves_bundle_head_to_tail(self):
        bundle,streets=sample(reverse=True)
        result=self.reader.read(bundle,streets)
        a,b=result.members
        self.assertTrue(a.forward);self.assertFalse(b.forward)
        self.assertFalse(b.spans[0].sections[0].forward)
        self.assertEqual(b.spans[0].sections[0].start,module.SectionAnchor(20,1,1.))
        self.assertEqual(b.spans[0].sections[0].end,module.SectionAnchor(20,0,0.))

    def test_partial_extents_have_real_section_references(self):
        bundle,streets=sample(reverse=True)
        bundle['head_locations']=((5.,0.),(5.,10.))
        bundle['tail_locations']=((35.,0.),(35.,10.))
        result=self.reader.read(bundle,streets)
        self.assertEqual(result.required_features,('INTERNAL_JUNCTIONS','PARTIAL_EXTENTS'))
        self.assertFalse(any(m.complete for m in result.members))
        a=result.members[0]
        self.assertEqual(a.spans[0].sections[0].start,module.SectionAnchor(10,0,.5))
        self.assertEqual(a.spans[1].sections[0].end,module.SectionAnchor(11,1,.5))
        self.assertEqual(result.members[1].spans[0].sections[0].start,module.SectionAnchor(20,1,.75))

    def test_reversed_partial_extent_orders_spans_and_boundaries(self):
        bundle,streets=sample()
        bundle['head_locations']=((35.,0.),(35.,10.))
        bundle['tail_locations']=((5.,0.),(5.,10.))
        result=self.reader.read(bundle,streets)
        first=result.members[0]
        self.assertFalse(first.forward)
        self.assertEqual([s.sections[0].start.section_id for s in first.spans],[11,10])
        self.assertEqual(first.spans[0].start.bundle_end,'HEAD')
        self.assertEqual(first.spans[-1].end.bundle_end,'TAIL')
        self.assertFalse(any(r.forward for s in first.spans for r in s.sections))

    def test_junction_at_extent_end_does_not_make_an_empty_span(self):
        bundle,streets=sample();bundle['head_locations']=((20.,0.),(20.,10.))
        result=self.reader.read(bundle,streets)
        span,=result.members[0].spans
        self.assertEqual(span.start,module.SpanBoundary('HEAD',8))
        self.assertEqual(span.sections[0].start.section_id,11)

    def test_junction_outside_extent_does_not_cut_the_member(self):
        bundle,streets=sample();bundle['head_locations']=((25.,0.),(25.,10.))
        result=self.reader.read(bundle,streets)
        self.assertFalse(result.junctions)
        self.assertEqual(result.required_features,('PARTIAL_EXTENTS',))

    def test_shared_section_endpoint_is_one_incidence(self):
        bundle,streets=sample(junction=False)
        bundle['head_locations']=((20.,0.),(0.,10.))
        result=self.reader.read(bundle,streets)
        self.assertEqual(result.members[0].spans[0].sections[0].start,module.SectionAnchor(11,0,0.))

    def test_more_than_two_members_does_not_lose_the_middle_member(self):
        bundle,streets=sample(junction=False)
        streets[4]=dict(id=4,items=[section(40,((0.,20.),(40.,20.)))])
        bundle.update(head=(1,2,4),tail=(1,2,4),head_locations=((0.,0.),(0.,10.),(0.,20.)),
                      tail_locations=((40.,0.),(40.,10.),(40.,20.)),owned_streets=(1,2,4))
        result=self.reader.read(bundle,streets)
        self.assertEqual(len(result.members),3)
        self.assertEqual(result.required_features,('MULTIPLE_MEMBERS',))

    def test_owned_interior_street_is_not_an_invented_third_outer_member(self):
        bundle,streets=sample();bundle['owned_streets']=(1,2,3)
        result=self.reader.read(bundle,streets)
        self.assertEqual(len(result.members),2)
        self.assertEqual(result.interior_streets,(3,))
        self.assertIn('INTERIOR_STREETS',result.required_features)

    def test_changing_member_order_is_not_silently_sorted(self):
        bundle,streets=sample();bundle.update(tail=(2,1),tail_locations=((40.,10.),(40.,0.)))
        self.reject(bundle,streets,'MEMBER_ORDER_CHANGE')

    def test_duplicate_section_across_members_is_rejected(self):
        bundle,streets=sample();streets[2]['items'][0]['id']=10
        self.reject(bundle,streets,'DUPLICATE_SECTION')

    def test_one_member_is_not_a_bundle(self):
        bundle,streets=sample();bundle.update(head=(1,),tail=(1,),
            head_locations=((0.,0.),),tail_locations=((40.,0.),))
        self.reject(bundle,streets,'INSUFFICIENT_MEMBERS')

    def test_missing_membership_is_not_inferred_from_owned_streets(self):
        bundle,streets=sample();bundle['tail']=(1,3)
        self.reject(bundle,streets,'UNRESOLVED_MEMBERSHIP')

    def test_unlocated_and_ambiguous_endpoints_are_rejected(self):
        bundle,streets=sample();bundle['head_locations']=((1.,1.),(0.,10.))
        self.reject(bundle,streets,'UNLOCATED_ENDPOINT')
        bundle,streets=sample(junction=False)
        streets[1]['items']=[section(10,((0.,0.),(10.,0.),(10.,10.),(0.,0.),(40.,0.)))]
        self.reject(bundle,streets,'AMBIGUOUS_ENDPOINT')

    def test_duplicate_or_absent_section_reference_is_rejected(self):
        bundle,streets=sample();streets[1]['items'][1]['through']=((10,1),(99,0))
        self.reject(bundle,streets,'INVALID_THROUGH_REFERENCE')
        bundle,streets=sample();streets[1]['items'].append(deepcopy(streets[1]['items'][0]))
        self.reject(bundle,streets,'DUPLICATE_SECTION')

    def test_nearby_attachment_does_not_establish_connectivity(self):
        bundle,streets=sample();streets[3]['items'][0]=section(30,((20.,.01),(20.,-20.)))
        self.reject(bundle,streets,'ATTACHMENT_LOCATION_MISMATCH')

    def test_duplicate_attachment_and_junction_are_rejected(self):
        bundle,streets=sample();item=streets[1]['items'][1]
        item['attachments']*=2;self.reject(bundle,streets,'DUPLICATE_ATTACHMENT')
        bundle,streets=sample();streets[1]['items'].insert(2,deepcopy(streets[1]['items'][1]))
        self.reject(bundle,streets,'DUPLICATE_JUNCTION')

    def test_collapsed_disconnected_and_nonfinite_paths_are_rejected(self):
        for points,code in ((((20.,0.),(20.,0.)),'DEGENERATE_EDGE'),
                            (((21.,0.),(40.,0.)),'DISCONNECTED_CHAIN'),
                            (((20.,0.),(float('nan'),0.)),'INVALID_SECTION')):
            bundle,streets=sample();streets[1]['items'][-1]=section(11,points)
            self.reject(bundle,streets,code)

    def test_report_distinguishes_mapping_from_surface_support(self):
        bundle,streets=sample();bundle.update(status='UNSUPPORTED',reason='Surface patch not implemented')
        issues=[]
        snapshot=SimpleNamespace(bundles={7:bundle},streets=streets,issue=lambda *args:issues.append(args))
        self.reader.apply(snapshot)
        self.assertEqual(bundle['status'],'UNSUPPORTED')
        self.assertEqual(bundle['topology']['status'],'MAPPED')
        self.assertFalse(issues)
        bundle['head_locations']=((1.,1.),(0.,10.));self.reader.apply(snapshot)
        self.assertEqual(bundle['topology']['status'],'UNRESOLVED')
        self.assertEqual(issues[0][2],'BUNDLE_TOPOLOGY_UNLOCATED_ENDPOINT')

    def test_deferred_attachment_is_not_promoted_to_carriageway(self):
        bundle,streets=sample();streets[1]['items'][1]['attachments'][0]['category']='footway'
        result=self.reader.read(bundle,streets)
        self.assertFalse(result.junctions[0].attachments[0].carriageway)


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BundleTopologyTests))
    if not result.wasSuccessful():raise RuntimeError('Bundle topology regression failed')
