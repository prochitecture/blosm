"""Transport complete source profiles through physical junction ownership."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

from test_junction_origins import sample, module as origins
from test_bounded_bundle_network import contract

spec=importlib.util.spec_from_file_location('tested_physical_junctions',
    Path(__file__).resolve().parents[1]/'way/street_network/physical_junctions.py')
module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)


def snapshot(reverse=False):
    before,after,streets=sample(reverse)
    recorder=origins.JunctionOriginRecorder()
    result=NS(junction_origins=recorder.resolve(recorder.capture(before),after).report(),streets={},bundles={})
    for sid,street in streets.items():
        result.streets[sid]=dict(id=sid,status='SUPPORTED',bundle_id=1 if sid<3 else 2 if sid<5 else None,
            items=[dict(id=street.head.id,points=street.head.centerline,width=6.,offset=0.,category='residential')])
    for bid,ids in ((1,(1,2)),(2,(3,4))):
        members=[dict(street_id=sid,complete=True,head_order=i,tail_order=i,
                      forward=not (reverse and sid in (2,3))) for i,sid in enumerate(ids)]
        result.bundles[bid]=dict(topology=dict(status='MAPPED',members=members))
    return result


def sections(s,sid):
    return tuple(contract.SectionInput(row['id'],row['points'],row['width'],row['offset'],row['category'])
                 for row in s.streets[sid]['items'])


class PhysicalJunctionTests(unittest.TestCase):
    def setUp(self):
        self.reader=module.PreparedPhysicalJunctions(contract,sections,{'residential':{'radius':2.}})

    def test_component_owns_each_site_once_without_source_changes(self):
        for reverse in (False,True):
            s=snapshot(reverse);old=deepcopy(s)
            result=self.reader.prepare(s,100)
            self.assertEqual(result,self.reader.prepare(s,200))
            self.assertEqual(s,old)
            self.assertEqual([j.source_id for j in result.junctions],[10,20,30])
            self.assertEqual(result.logical_junction_ids,(100,200))
            for j in result.junctions:
                for street in j.streets:
                    self.assertEqual(street.sections,sections(s,street.source_id))
                self.assertTrue(all(c.radius==2. for c in j.connectors))
                self.assertTrue(all(c.leaving==c.members[0].leaving for c in j.connectors))

    def test_half_profile_sign_follows_authored_end_and_lateral_role(self):
        normal=self.reader.prepare(snapshot(),100)
        reversed=self.reader.prepare(snapshot(True),100)
        a=[(j.source_id,[p.half_profile for p in j.profiles]) for j in normal.junctions]
        b=[(j.source_id,[p.half_profile for p in j.profiles]) for j in reversed.junctions]
        self.assertEqual(a,b)
        self.assertEqual(a,[(10,[1,0,0]),(20,[-1,0,0,1]),(30,[-1,0,0])])
        self.assertTrue(all(not p.separator for j in normal.junctions for p in j.profiles))

    def test_consecutive_streets_keep_the_same_bundle_half_profile(self):
        from test_bundle_topology import chained, PreparedBundleTopology
        for reverse in (False, True):
            bundle, streets, records = chained(reverse)
            bundle['topology'] = PreparedBundleTopology().read(bundle, streets, origins=records).report()
            s = NS(streets=streets, bundles={7: bundle})
            expected = -1 if reverse else 1
            for sid in (1, 3):
                self.assertEqual(self.reader._half(s, streets[sid], 0), (expected, 7))
                self.assertEqual(self.reader._half(s, streets[sid], 1), (-expected, 7))
            self.assertEqual(self.reader._half(s, streets[2], 0), (-expected, 7))
            bundle['topology']['chains'] = ((1, 3), (2, 3))
            with self.assertRaisesRegex(ValueError, 'distinct complete lateral roles'):
                self.reader._half(s, streets[1], 0)

    def test_unresolved_profiles_are_retained_only_by_explicit_policy(self):
        from test_physical_ports import fixture, m as ports
        s=fixture()
        s.bundles={7:dict(topology=dict(status='UNRESOLVED',reason='Three incident members'))}
        s.streets[2]['bundle_id']=7
        plan=ports.PreparedPhysicalPorts(('residential',)).prepare(s,interior_ids=(21,22))
        before=deepcopy(s)
        with self.assertRaisesRegex(ValueError,'Bundle 7.*Three incident members'):
            self.reader.prepare_port_designs(s,plan)
        result=self.reader.prepare_port_designs(s,plan,retain_unresolved=True)
        self.assertEqual(len(result),sum(len(site.ports) for site in plan.sites))
        self.assertEqual(s,before)
        rejected=[d for d in result if d.profile_error]
        self.assertTrue(rejected)
        self.assertTrue(all(d.street_id==2 and d.bundle_id==7 and d.half_profile is None for d in rejected))
        self.assertTrue(all(d.half_profile==0 for d in result if not d.profile_error))
        s.streets[2]['items'][0]['category']='undefined'
        with self.assertRaisesRegex(ValueError,'declared radius'):
            self.reader.prepare_port_designs(s,plan,retain_unresolved=True)

    def test_rejected_profile_policy_does_not_hide_missing_bundle(self):
        from test_physical_ports import fixture, m as ports
        s=fixture();s.bundles={};s.streets[2]['bundle_id']=7
        plan=ports.PreparedPhysicalPorts(('residential',)).prepare(s,interior_ids=(21,22))
        with self.assertRaisesRegex(ValueError,'absent Bundle'):
            self.reader.prepare_port_designs(s,plan,retain_unresolved=True)

    def test_original_ring_ids_and_complete_returning_source_survive(self):
        s=snapshot();result=self.reader.prepare(s,100)
        physical={r['junction_id']:r for r in s.junction_origins['physical']}
        ends=[]
        for j in result.junctions:
            self.assertEqual([c.source_id for c in j.connectors],[r['connector_id'] for r in physical[j.source_id]['ring']])
            for c in j.connectors:
                m=c.members[0];street=j.streets[m.street_index]
                if street.source_id==5:ends.append((j.source_id,m.leaving,street.sections[0].source_id))
        self.assertEqual(ends,[(10,True,50),(20,False,50)])

    def test_missing_unsupported_or_changed_sources_are_not_omitted(self):
        for change in ('missing','unsupported','terminal','location','partial'):
            s=snapshot()
            if change=='missing':del s.streets[5]
            elif change=='unsupported':s.streets[5]['status']='UNSUPPORTED'
            elif change=='terminal':s.streets[5]['items'][0]['id']=999
            elif change=='location':s.streets[5]['items'][0]['points']=((-21.,0.),(0.,0.))
            else:s.bundles[1]['topology']['members'][0]['complete']=False
            with self.subTest(change=change),self.assertRaises(ValueError):self.reader.prepare(s,100)

    def test_unresolved_origin_and_duplicate_surface_end_owner_reject(self):
        s=snapshot();s.junction_origins['components'][0]['status']='UNRESOLVED'
        with self.assertRaises(ValueError):self.reader.prepare(s,100)
        s=snapshot();site=s.junction_origins['physical'][0]
        site['ring']=site['ring']+(site['ring'][0],)
        with self.assertRaisesRegex(ValueError,'multiple physical'):self.reader.prepare(s,100)

    def test_no_radius_substitution(self):
        s=snapshot();s.streets[5]['items'][0]['category']='undefined'
        with self.assertRaisesRegex(ValueError,'declared radius'):self.reader.prepare(s,100)

    def test_declared_interior_street_keeps_full_profile_and_original_ownership(self):
        s=snapshot();expected=self.reader.prepare(s,100)
        s.streets[5]['bundle_id']=1
        s.bundles[1]['topology']['interior_streets']=[5]
        before=deepcopy(s)
        result=self.reader.prepare(s,100)
        self.assertEqual(result,expected)
        self.assertEqual(s,before)
        self.assertEqual(s.streets[5]['bundle_id'],1)
        # Both physical ends keep their existing ordinary connector role.
        self.assertEqual([p.half_profile for j in result.junctions
            for street,p in zip(j.streets,j.profiles) if street.source_id==5],[0,0])

    def test_undeclared_or_conflicting_interior_profile_roles_reject(self):
        for mutation in ('undeclared','duplicate','longitudinal','unmapped'):
            s=snapshot();s.streets[5]['bundle_id']=1
            topology=s.bundles[1]['topology']
            topology['interior_streets']=[5]
            if mutation=='undeclared':topology['interior_streets']=[]
            elif mutation=='duplicate':topology['interior_streets']=[5,5]
            elif mutation=='longitudinal':topology['members'][0]['street_id']=5
            else:topology['status']='UNRESOLVED'
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                self.reader.prepare(s,100)


if __name__=='__main__':unittest.main()
