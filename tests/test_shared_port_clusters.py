"""Shared logical rings retain identity and an arbitrarily long boundary walk."""
from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace as NS
import unittest

from test_port_clusters import m, ports, section


def fixture():
    locations={10:(0.,0.),11:(20.,0.),12:(40.,0.)}
    incidences={10:((100,1),(101,1),(102,0),(70,0)),
                11:((104,0),(86,0),(70,1),(103,1)),
                12:((106,0),(107,0),(86,1),(105,1))}
    endpoint_sites={(sid,end):pid for pid,ring in incidences.items() for sid,end in ring}
    bundles={300:(102,103),301:(105,104)}
    street_bundles={sid:bid for bid,members in bundles.items() for sid in members}
    streets={}
    for sid in sorted({sid for sid,_ in endpoint_sites}):
        owners=[endpoint_sites.get((sid,end)) for end in (0,1)]
        points=[locations[pid] if pid is not None else (float(sid),60.) for pid in owners]
        streets[sid]=dict(id=sid,status='SUPPORTED',endpoints=tuple(owners),
            bundle_id=street_bundles.get(sid),items=[section(sid*10,points)])
    def endpoint(sid,end):
        return dict(street_id=sid,terminal_id=sid*10,terminal_kind='Section',end=end)
    physical=[]; lookup={}
    for pid,ring in incidences.items():
        connectors=[]
        for at,(sid,end) in enumerate(ring):
            cid=pid*10+at
            connectors.append(dict(connector_id=cid,item_kind='Street',endpoints=[endpoint(sid,end)]))
            lookup[sid,end]=(pid,cid)
        physical.append(dict(junction_id=pid,location=locations[pid],ring=connectors))
    definitions={5000:((301,((105,1),(104,0))), (106,((106,0),)), (107,((107,0),)),
                       (86,((86,1),)), (86,((86,0),)), (70,((70,1),))),
                 5001:((300,((102,0),(103,1))), (100,((100,1),)), (101,((101,1),)), (70,((70,0),)))}
    logical=[]
    for jid,connectors in definitions.items():
        ring=[]; provenance=[]
        for at,(item,ends) in enumerate(connectors):
            # Deliberate connector-ID reuse across distinct logical rings.
            cid=1000+at
            ring.append(dict(connector_id=cid,item_id=item,item_kind='Bundle' if len(ends)==2 else 'Street',
                leaving=ends[0][1]==0,endpoints=[endpoint(*key) for key in ends]))
            for key in ends:
                pid,physical_cid=lookup[key]
                provenance.append(dict(logical_connector_id=cid,endpoint=endpoint(*key),status='EXACT',
                    physical_junction_id=pid,physical_connector_id=physical_cid))
        logical.append(dict(junction=dict(junction_id=jid,ring=ring),endpoints=provenance))
    bundle_rows={bid:dict(id=bid,status='SUPPORTED',topology=dict(status='MAPPED',members=[
        dict(street_id=sid,head_order=i,tail_order=i,complete=True,forward=True)
        for i,sid in enumerate(members)])) for bid,members in bundles.items()}
    return NS(streets=streets,bundles=bundle_rows,junction_origins=dict(physical=physical,logical=logical,
        components=[dict(logical_junction_ids=[5000,5001],physical_junction_ids=[10,11,12],
                         shared_physical_junction_ids=[11],status='MAPPED',unresolved=[])]))


class SharedLogicalClusterTests(unittest.TestCase):
    def setUp(self):
        self.snapshot=fixture()
        self.plan=ports.PreparedPhysicalPorts(('residential',)).prepare(self.snapshot,physical_ids=(10,11,12))
        self.reader=m.PreparedLogicalPortClusters()

    def result(self,jid=5000):
        return self.reader.prepare(self.snapshot,self.plan,jid,
            allow_shared_owners=True,allow_internal_connectors=True)

    def test_three_original_corners_and_two_internal_spans_are_retained(self):
        before=deepcopy(self.snapshot)
        result=self.result()
        self.assertEqual(result.owner_ids,(5000,5001))
        self.assertEqual(result.physical_junction_ids,(10,11,12))
        self.assertEqual([p.endpoint.street_id for p in result.external_ring],list(range(100,108)))
        chain=result.boundaries[-1]
        self.assertEqual([s[0][0] for s in chain.traversed_spans],[86,70])
        self.assertEqual([p.junction_id for p in chain.original_corners],[12,11,10])
        self.assertEqual(result.footprint_spans,((70,700,700),(86,860,860)))
        self.assertEqual(len(result.internal_groups),4)
        self.assertEqual(len({(g.logical_junction_id,g.connector_id) for g in result.groups}),6)
        self.assertEqual(self.snapshot,before)

    def test_either_original_owner_selects_the_same_complete_component(self):
        self.assertEqual(self.result(5000),self.result(5001))

    def test_no_implicit_admission_to_the_older_single_owner_branch(self):
        for shared,internal in ((False,True),(True,False),(False,False)):
            with self.subTest(shared=shared,internal=internal),self.assertRaises(ValueError):
                self.reader.prepare(self.snapshot,self.plan,5000,
                    allow_shared_owners=shared,allow_internal_connectors=internal)

    def test_input_order_and_cyclic_origins_do_not_change_ownership(self):
        expected=self.result()
        self.snapshot.junction_origins['logical'].reverse()
        self.snapshot.junction_origins['components'][0]['logical_junction_ids'].reverse()
        self.snapshot.junction_origins['components'][0]['physical_junction_ids'].reverse()
        for row in self.snapshot.junction_origins['logical']:
            ring=row['junction']['ring'];row['junction']['ring']=ring[2:]+ring[:2]
            row['endpoints'].reverse()
        # Internal memberships are canonicalized separately from circular rings.
        self.plan=replace(self.plan,sites=tuple(replace(s,ports=s.ports[1:]+s.ports[:1])
                           for s in reversed(self.plan.sites)),spans=tuple(reversed(self.plan.spans)))
        actual=self.result()
        self.assertEqual(actual.external_ring,expected.external_ring)
        self.assertEqual(actual.groups,expected.groups)
        self.assertEqual(actual.boundaries,expected.boundaries)
        self.assertEqual(set(actual.internal_groups),set(expected.internal_groups))

    def test_missing_shared_site_or_original_owner_rejects(self):
        for change in ('shared','logical','site','duplicate'):
            self.setUp()
            component=self.snapshot.junction_origins['components'][0]
            if change=='shared':component['shared_physical_junction_ids']=[]
            elif change=='logical':self.snapshot.junction_origins['logical'].pop()
            elif change=='site':component['physical_junction_ids'].remove(11)
            else:component['logical_junction_ids'].append(5000)
            with self.subTest(change=change),self.assertRaises(ValueError):self.result()

    def test_duplicate_or_wrong_original_incidence_rejects(self):
        for change in ('duplicate','connector','terminal','one_internal_end'):
            self.setUp()
            row=self.snapshot.junction_origins['logical'][0]
            if change=='duplicate':row['endpoints'].append(deepcopy(row['endpoints'][0]))
            elif change=='connector':row['endpoints'][0]['physical_connector_id']+=1
            elif change=='terminal':row['endpoints'][0]['endpoint']['terminal_id']+=1
            else:
                row['endpoints']=[r for r in row['endpoints'] if r['endpoint']['street_id']!=70]
                row['junction']['ring']=[r for r in row['junction']['ring'] if r['item_id']!=70]
            with self.subTest(change=change),self.assertRaises(ValueError):self.result()

    def test_incorrect_logical_cyclic_order_rejects(self):
        row=self.snapshot.junction_origins['logical'][0]
        row['junction']['ring'][1],row['junction']['ring'][2]=row['junction']['ring'][2],row['junction']['ring'][1]
        with self.assertRaises(ValueError):self.result()

    def test_joint_network_owns_each_original_span_and_site_once(self):
        from test_port_network import m as network, compounds
        cluster=self.result()
        result=network.PreparedPortNetwork(compounds.PreparedPortCompound()).prepare(self.plan,logical_clusters=(cluster,))
        self.assertFalse(result.ordinary_site_ids)
        self.assertFalse(result.unresolved)
        self.assertEqual(sum(r.kind=='LOGICAL_CLUSTER' for r in result.routes),2)
        self.assertEqual(sum(r.kind=='APPROACH' for r in result.routes),8)
        self.assertEqual(tuple(r.span for r in result.routes),self.plan.spans)
        for ids in ((),(5000,5000),(4999,5000),(-1,5000),(5000,),(5000,7000)):
            with self.subTest(ids=ids),self.assertRaises(ValueError):
                network.PreparedPortNetwork(compounds.PreparedPortCompound()).prepare(self.plan,
                    logical_clusters=(replace(cluster,logical_junction_ids=ids),))

    def test_no_original_logical_owner_can_be_hidden_by_component_declaration(self):
        # An extra logical ring touches the selected physical sites but is
        # deliberately omitted from the component. Its incidences are evidence
        # of a shared owner even if all necessary pavement ports were present.
        extra=deepcopy(self.snapshot.junction_origins['logical'][0])
        extra['junction']['junction_id']=7000
        self.snapshot.junction_origins['logical'].append(extra)
        with self.assertRaisesRegex(ValueError,'omits an original logical owner'):self.result()


if __name__=='__main__':unittest.main()
