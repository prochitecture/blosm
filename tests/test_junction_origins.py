"""Physical-site identity survives logical Bundle aggregation without geometry edits."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

spec = importlib.util.spec_from_file_location('tested_junction_origins',
    Path(__file__).resolve().parents[1]/'way/junction_origins.py')
module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module
spec.loader.exec_module(module)


class Section(NS):
    pass


class Street(NS):
    pass


class Bundle(NS):
    pass


class Intersection(NS):
    def __iter__(self):
        return iter(self.connectors)


def street(sid, a, b):
    section = Section(id=10*sid, centerline=(a, b), width=6., offset=0.)
    return Street(id=sid, src=a, dst=b, head=section, tail=section, pred=None, succ=None)


def junction(jid, location, endpoints):
    return Intersection(id=jid, location=location,
        connectors=[NS(id=jid*10+i, item=s, leaving=leaving) for i,(s,leaving) in enumerate(endpoints)])


def sample(reverse=False):
    a,b,c=(-20.,0.),(0.,0.),(20.,0.)
    streets={1:street(1,a,(-20.,30.)), 2:street(2,b,(0.,30.)),
             3:street(3,b,(0.,-30.)), 4:street(4,c,(20.,-30.)),
             5:street(5,a,b), 6:street(6,b,c),
             7:street(7,a,(-20.,-30.)), 8:street(8,c,(20.,30.))}
    directions={sid:True for sid in streets}
    if reverse:
        for sid in (2,3,5):
            s=streets[sid]; s.src,s.dst=s.dst,s.src
            s.head.centerline=s.head.centerline[::-1]; directions[sid]=False
    def refs(items):
        return [(streets[sid], directions[sid] if outgoing else not directions[sid]) for sid,outgoing in items]
    physical=[junction(10,a,refs(((1,True),(5,True),(7,True)))),
              junction(20,b,refs(((2,True),(5,False),(6,True),(3,True)))),
              junction(30,c,refs(((4,True),(6,False),(8,True))))]
    left=Bundle(id=1,streetsHead=[streets[1],streets[2]],headLocs=[a,b],
                streetsTail=[streets[1],streets[2]],tailLocs=[(-20.,30.),(0.,30.)])
    right=Bundle(id=2,streetsHead=[streets[3],streets[4]],headLocs=[b,c],
                 streetsTail=[streets[3],streets[4]],tailLocs=[(0.,-30.),(20.,-30.)])
    logical=[junction(100,(-10.,0.),[(left,True)]+refs(((5,True),(5,False),(7,True),(6,True)))),
             junction(200,(10.,0.),[(right,True)]+refs(((6,False),(8,True))))]
    return physical, logical, streets


class JunctionOriginTests(unittest.TestCase):
    def setUp(self):
        self.recorder=module.JunctionOriginRecorder()

    def test_removed_splitter_keeps_exact_physical_source_without_changing_logical_inventory(self):
        nodes,_,sources=sample()
        physical=self.recorder.capture(nodes)
        splitter=sources[5]
        logical={sid:s for sid,s in sources.items() if sid!=5}
        before=deepcopy(sources)
        retained=self.recorder.retain_physical_sources(logical.values(),physical,(splitter,splitter))
        self.assertEqual(set(retained),set(sources))
        self.assertIs(retained[5],splitter)
        self.assertNotIn(5,logical)
        self.assertEqual(sources,before)
        ends=[(j.junction_id,e.end,e.terminal_id) for j in physical for c in j.ring
              for e in c.endpoints if e.street_id==5]
        self.assertEqual(ends,[(10,0,50),(20,1,50)])

    def test_unreferenced_removed_candidate_does_not_become_an_extra_physical_street(self):
        nodes,_,sources=sample()
        unrelated=street(90,sources[5].src,sources[5].dst)
        retained=self.recorder.retain_physical_sources(
            sources.values(),self.recorder.capture(nodes),(unrelated,))
        self.assertEqual(set(retained),set(sources))

    def test_physical_source_identity_collision_is_not_silently_overwritten(self):
        nodes,_,sources=sample()
        with self.assertRaisesRegex(ValueError,'different source objects'):
            self.recorder.retain_physical_sources(sources.values(),self.recorder.capture(nodes),
                                                 (deepcopy(sources[5]),))

    def test_missing_physical_source_is_not_inferred_from_a_coincident_candidate(self):
        nodes,_,sources=sample();physical=self.recorder.capture(nodes)
        wrong=street(90,sources[5].src,sources[5].dst)
        del sources[5]
        retained=self.recorder.retain_physical_sources(sources.values(),physical,(wrong,))
        self.assertNotIn(5,retained)
        self.assertNotIn(90,retained)

    def test_seams_preserve_each_authored_direction_and_current_section_owner(self):
        for ends in ((0,0),(1,1),(0,1),(1,0)):
            with self.subTest(ends=ends):
                sources = [street(1,(0.,0.),(-40.,0.)), street(2,(0.,0.),(40.,0.))]
                for s,end in zip(sources,ends):
                    if end:s.src,s.dst=s.dst,s.src;s.head.centerline=s.head.centerline[::-1]
                node=junction(20,(0.,0.),[])
                # Stale Street backrefs must not affect identity resolution.
                node.leaveWays=[NS(section=s.head,leaving=end==0,street=object()) for s,end in zip(sources,ends)]
                before=deepcopy(sources)
                seams=self.recorder.capture_seams((node,),sources)
                seam,=seams
                self.assertEqual([(e.street_id,e.end,e.terminal_id) for e in seam.endpoints],
                                 [(1,ends[0],10),(2,ends[1],20)])
                self.assertEqual(sources,before)
                free=self.recorder.capture_free_endpoints(sources,seams=seams)
                self.assertEqual({e.endpoint.key for e in free},{(1,1-ends[0]),(2,1-ends[1])})
                result=self.recorder.resolve(self.recorder.capture((node,)),(),free_endpoints=free,seams=seams)
                self.assertEqual(result.seams,seams)
                json.dumps(result.report())
                with self.assertRaises(FrozenInstanceError):seam.junction_id=999

    def test_seam_identity_does_not_use_equal_coordinates_or_section_ids(self):
        a,b=street(1,(0.,0.),(-40.,0.)),street(2,(0.,0.),(40.,0.))
        node=junction(20,(0.,0.),[])
        node.leaveWays=[NS(section=a.head,leaving=True),NS(section=deepcopy(b.head),leaving=True)]
        with self.assertRaisesRegex(ValueError,'no current terminal'):
            self.recorder.capture_seams((node,),(a,b))
        self.assertFalse(self.recorder.capture_seams((),(a,b)))
        node.leaveWays[1].section=b.head;b.src=(.001,0.)
        with self.assertRaisesRegex(ValueError,'node location'):
            self.recorder.capture_seams((node,),(a,b))

    def test_seam_duplicate_and_conflicting_ownership_rejects(self):
        a,b=street(1,(0.,0.),(-40.,0.)),street(2,(0.,0.),(40.,0.))
        node=junction(20,(0.,0.),[])
        node.leaveWays=[NS(section=s.head,leaving=True) for s in (a,b)]
        seams=self.recorder.capture_seams((node,),(a,b))
        with self.assertRaisesRegex(ValueError,'unique'):
            self.recorder.capture_seams((node,node),(a,b))
        with self.assertRaisesRegex(ValueError,'conflicting'):
            self.recorder.resolve((),(),seams=seams,free_endpoints=self.recorder.capture_free_endpoints((a,b)))
        physical=self.recorder.capture((junction(30,a.src,[(a,True)]),))
        with self.assertRaisesRegex(ValueError,'conflicting'):
            self.recorder.resolve(physical,(),seams=seams)
        node.leaveWays[1]=node.leaveWays[0]
        with self.assertRaisesRegex(ValueError,'duplicate ownership'):
            self.recorder.capture_seams((node,),(a,b))
        node.leaveWays[1]=NS(section=b.head,leaving=True);a.pred=object()
        with self.assertRaisesRegex(ValueError,'connector owner'):
            self.recorder.capture_seams((node,),(a,b))

    def test_seam_survives_later_logical_bundle_relinking(self):
        a,b=street(1,(0.,0.),(-40.,0.)),street(2,(0.,0.),(40.,0.))
        node=junction(20,(0.,0.),[])
        node.leaveWays=[NS(section=s.head,leaving=True) for s in (a,b)]
        seams=self.recorder.capture_seams((node,),(a,b))
        bundle=Bundle(id=7,streetsHead=[a],headLocs=[a.src])
        logical=junction(100,(5.,5.),[(bundle,True)])
        result=self.recorder.resolve(self.recorder.capture((node,)),(logical,),seams=seams)
        record,=result.logical[0].endpoints
        self.assertEqual((record.status,record.physical_junction_id,record.physical_connector_id),
                         ('ORIGINAL_SEAM',20,None))
        self.assertEqual(result.components()[0]['status'],'UNRESOLVED')

    def test_free_origin_survives_bundle_relinking_in_either_direction(self):
        for end in (0,1):
            with self.subTest(end=end):
                source=street(1,(0.,0.),(40.,0.))
                setattr(source,'succ' if end==0 else 'pred',object())
                free=self.recorder.capture_free_endpoints((source,))
                self.assertEqual(len(free),1)
                self.assertEqual(free[0].endpoint.end,end)
                original_location=source.src if end==0 else source.dst
                bundle=Bundle(id=7,streetsHead=[source],headLocs=[original_location],
                    streetsTail=[source],tailLocs=[original_location])
                logical=junction(10,(20.,10.),[(bundle,True)])
                setattr(source,'pred' if end==0 else 'succ',logical.connectors[0])
                result=self.recorder.resolve((),(logical,),free_endpoints=free)
                endpoint,=result.logical[0].endpoints
                self.assertEqual(endpoint.status,'ORIGINAL_FREE')
                self.assertIsNone(endpoint.physical_junction_id)
                self.assertEqual(result.components()[0]['status'],'UNRESOLVED')
                self.assertEqual(result.free_endpoints,free)
                source.src=(-100.,0.)
                self.assertEqual(free[0].location,original_location)
                with self.assertRaises(FrozenInstanceError):free[0].location=(1.,2.)
                json.dumps(result.report())

    def test_missing_connector_attribute_is_not_a_proven_free_end(self):
        source=street(1,(0.,0.),(40.,0.));del source.pred
        with self.assertRaises(AttributeError):self.recorder.capture_free_endpoints((source,))

    def test_free_origin_rejects_duplicates_and_reports_conflicts_or_changed_terminal(self):
        source=street(1,(0.,0.),(40.,0.));source.succ=object()
        free=self.recorder.capture_free_endpoints((source,))
        with self.assertRaises(ValueError):self.recorder.capture_free_endpoints((source,source))
        bundle=Bundle(id=7,streetsHead=[source],headLocs=[source.src])
        logical=junction(10,(20.,10.),[(bundle,True)])
        with self.assertRaises(ValueError):self.recorder.resolve((),(logical,),free_endpoints=free+free)
        physical=self.recorder.capture((junction(20,source.src,[(source,True)]),))
        conflict=self.recorder.resolve(physical,(logical,),free_endpoints=free)
        self.assertEqual(conflict.logical[0].endpoints[0].status,'CONFLICTING_FREE_ORIGIN')
        source.head=Section(id=999)
        changed=self.recorder.resolve((),(logical,),free_endpoints=free)
        self.assertEqual(changed.logical[0].endpoints[0].status,'CHANGED_TERMINAL')

    def test_shared_site_has_one_component_in_both_authored_directions(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                old,new,streets=sample(reverse); before=deepcopy(streets)
                result=self.recorder.resolve(self.recorder.capture(old),new)
                self.assertEqual(streets,before)
                component,=result.components()
                self.assertEqual(component['logical_junction_ids'],(100,200))
                self.assertEqual(component['physical_junction_ids'],(10,20,30))
                self.assertEqual(component['shared_physical_junction_ids'],(20,))
                self.assertEqual(component['status'],'MAPPED')
                self.assertEqual(len(result.logical[0].endpoints),6)
                self.assertEqual(len(result.logical[1].endpoints),4)
                json.dumps(result.report())

    def test_returning_street_retains_two_different_site_owners(self):
        old,new,_=sample()
        result=self.recorder.resolve(self.recorder.capture(old),new)
        refs=[e for e in result.logical[0].endpoints if e.endpoint.street_id==5]
        self.assertEqual([(e.endpoint.end,e.physical_junction_id) for e in refs],[(0,10),(1,20)])

    def test_snapshot_does_not_follow_mutated_source_objects(self):
        old,new,streets=sample(); physical=self.recorder.capture(old)
        report=repr(physical)
        streets[1].id=999; old[0].connectors.reverse(); old[0].location=(100.,100.)
        self.assertEqual(repr(physical),report)
        with self.assertRaises(FrozenInstanceError):physical[0].junction_id=8

    def test_original_ring_order_is_preserved(self):
        old,new,_=sample(); physical=self.recorder.capture(old)
        self.assertEqual([e.endpoints[0].street_id for e in physical[1].ring],[2,5,6,3])
        result=self.recorder.resolve(physical,new)
        self.assertEqual([c.connector_id for c in result.logical[0].junction.ring],
                         [c.id for c in new[0].connectors])

    def test_missing_or_duplicate_origin_is_explicit(self):
        for duplicate in (False,True):
            with self.subTest(duplicate=duplicate):
                old,new,_=sample()
                original=old[0].connectors[0]
                if duplicate:old[1].connectors.append(original)
                else:old[0].connectors.remove(original)
                result=self.recorder.resolve(self.recorder.capture(old),new)
                row=next(e for e in result.logical[0].endpoints if e.endpoint.street_id==1)
                self.assertEqual(row.status,'AMBIGUOUS_ORIGIN' if duplicate else 'MISSING_ORIGIN')
                self.assertEqual(result.components()[0]['status'],'UNRESOLVED')

    def test_same_coordinates_do_not_create_connectivity(self):
        old,new,_=sample(); physical=self.recorder.capture(old)
        first=module.JunctionOrigin(90,(0.,0.),())
        second=module.JunctionOrigin(91,(0.,0.),())
        ref=module.StreetEndpoint(1,0,10,'Section')
        logical=tuple(module.LogicalJunctionOrigins(j,(module.EndpointOrigin(1,ref,pid,1,'EXACT'),))
                      for j,pid in ((first,10),(second,30)))
        result=module.JunctionOrigins(physical,logical)
        self.assertEqual(len(result.components()),2)

    def test_changed_terminal_does_not_resolve_by_street_id_alone(self):
        old,new,streets=sample(); physical=self.recorder.capture(old)
        streets[1].head=Section(id=999,centerline=streets[1].head.centerline)
        result=self.recorder.resolve(physical,new)
        self.assertEqual(result.logical[0].endpoints[0].status,'CHANGED_TERMINAL')
        self.assertEqual(result.components()[0]['status'],'UNRESOLVED')

    def test_partial_bundle_endpoint_does_not_use_nearest_source_end(self):
        old,new,_=sample(); physical=self.recorder.capture(old)
        new[0].connectors[0].item.headLocs[0]=(-20.,1.)
        result=self.recorder.resolve(physical,new)
        row=result.logical[0].endpoints[0]
        self.assertEqual(row.endpoint.end,-1)
        self.assertEqual(row.status,'MISSING_ORIGIN')

    def test_invalid_member_count_and_duplicate_site_ids_reject(self):
        old,new,_=sample()
        with self.assertRaises(ValueError):self.recorder.capture(old+[old[0]])
        physical=self.recorder.capture(old)
        new[0].connectors[0].item.headLocs.pop()
        with self.assertRaises(ValueError):self.recorder.resolve(physical,new)

    def test_damaged_repeating_connector_iterator_is_bounded(self):
        old,_,_=sample();first=old[0].connectors[0]
        class BrokenIntersection(Intersection):
            def __iter__(self):
                while True:yield first
        bad=BrokenIntersection(id=999,location=(0.,0.))
        with self.assertRaisesRegex(ValueError,'repeats before closing'):
            self.recorder.capture((bad,))


if __name__=='__main__':
    unittest.main()
