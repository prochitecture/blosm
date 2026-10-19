"""Absorbed physical sources remain explicit without changing logical ownership."""
from copy import deepcopy
import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest

spec=importlib.util.spec_from_file_location('tested_physical_sources',
    Path(__file__).resolve().parents[1]/'way/street_network/physical_sources.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
coverage_spec=importlib.util.spec_from_file_location('tested_physical_coverage',
    Path(__file__).resolve().parents[1]/'way/street_network/coverage.py')
coverage=importlib.util.module_from_spec(coverage_spec);coverage_spec.loader.exec_module(coverage)


def fixture():
    row=dict(id=25,items=[dict(kind='Section',id=125,points=((0.,0.),(20.,0.)),
        width=5.,offset=0.,tags={})],status='SUPPORTED',reason=None,endpoints=(999,999),bundle_id=None)
    sites=[dict(junction_id=jid,ring=[dict(endpoints=[dict(street_id=25,end=end,
        terminal_kind='Section',terminal_id=125)])]) for end,jid in enumerate((10,20))]
    snapshot=NS(streets={},retained_physical_streets={25:row},bundles={},
                junction_origins=dict(physical=sites),issues=[],intersections={})
    snapshot.issue=lambda kind,sid,code,message:snapshot.issues.append((kind,sid,code,message))
    return snapshot


class PhysicalSourceTests(unittest.TestCase):
    def test_physical_ownership_replaces_obsolete_connectors_without_changing_source(self):
        s=fixture();items=deepcopy(s.retained_physical_streets[25]['items'])
        module.PhysicalStreetInventory().review(s)
        row=s.retained_physical_streets[25]
        self.assertEqual(row['items'],items)
        self.assertEqual(row['endpoints'],(10,20))
        self.assertEqual(row['endpoint_space'],'PHYSICAL_JUNCTION')
        self.assertEqual(row['coverage'],'PENDING_PHYSICAL_PROVIDER')
        self.assertEqual(row['status'],'SUPPORTED')
        self.assertFalse(s.streets)
        self.assertIs(module.PhysicalStreetInventory.view(s).streets[25],row)

    def test_missing_or_ambiguous_owner_does_not_become_a_free_end(self):
        for duplicate in (False,True):
            s=fixture();sites=s.junction_origins['physical']
            if duplicate:sites.append(deepcopy(sites[0]))
            else:sites.pop()
            module.PhysicalStreetInventory().review(s)
            self.assertEqual(s.retained_physical_streets[25]['status'],'UNSUPPORTED')
            self.assertIn('unique',s.retained_physical_streets[25]['reason'])

    def test_changed_terminal_or_reowned_section_rejects_without_changing_logical_record(self):
        for changed in (True,False):
            s=fixture()
            if changed:s.junction_origins['physical'][0]['ring'][0]['endpoints'][0]['terminal_id']=999
            else:
                s.streets[90]=deepcopy(s.retained_physical_streets[25]);s.streets[90]['id']=90
            logical=deepcopy(s.streets)
            module.PhysicalStreetInventory().review(s)
            self.assertEqual(s.retained_physical_streets[25]['status'],'UNSUPPORTED')
            self.assertEqual(s.streets,logical)

    def test_existing_source_rejection_is_preserved(self):
        s=fixture();s.retained_physical_streets[25].update(status='UNSUPPORTED',reason='Unsupported elevation')
        module.PhysicalStreetInventory().review(s)
        self.assertEqual(s.retained_physical_streets[25]['reason'],'Unsupported elevation')

    def test_duplicate_inventory_identity_rejects(self):
        s=fixture();s.streets[25]=s.retained_physical_streets[25]
        with self.assertRaisesRegex(ValueError,'active and retained'):
            module.PhysicalStreetInventory.view(s)

    def test_retained_length_is_pending_and_not_lost_from_coverage(self):
        s=fixture();module.PhysicalStreetInventory().review(s)
        report=coverage.CoverageReport().build(s)
        self.assertEqual(report['street_counts'],{'PENDING_PHYSICAL_PROVIDER':1})
        self.assertEqual(report['centerline_length_m'],{'PENDING_PHYSICAL_PROVIDER':20.})
        s.streets[25]=s.retained_physical_streets[25]
        with self.assertRaisesRegex(ValueError,'both active and retained'):
            coverage.CoverageReport().build(s)


class PhysicalBundleTests(unittest.TestCase):
    def prepare(self):
        from test_bounded_bundle_network import topology
        s=fixture();s.retained_physical_bundles={}
        s.retained_physical_streets[25]['bundle_id']=7
        other=deepcopy(s.retained_physical_streets[25]);other.update(id=26,bundle_id=7)
        other['items'][0].update(id=126,points=((20.,10.),(0.,10.)))
        s.retained_physical_streets[26]=other
        bundle=NS(id=7,streetsHead=[NS(id=25),NS(id=26)],streetsTail=[NS(id=25),NS(id=26)],
                  headLocs=[(0.,0.),(0.,10.)],tailLocs=[(20.,0.),(20.,10.)])
        return s,bundle,topology.PreparedBundleTopology(('residential',))

    def test_preserve_original_two_member_roles_and_pending_coverage(self):
        s,b,reader=self.prepare();before=deepcopy(s.retained_physical_streets)
        module.PhysicalBundleInventory().read(s,(b,),reader)
        self.assertEqual(before,s.retained_physical_streets);self.assertFalse(s.bundles)
        row=s.retained_physical_bundles[7]
        self.assertEqual(row['status'],'SUPPORTED')
        self.assertEqual([m['forward'] for m in row['topology']['members']],[True,False])
        self.assertEqual(row['coverage'],'PENDING_PHYSICAL_PROVIDER')
        self.assertIs(module.PhysicalStreetInventory.view(s).bundles[7],row)
        self.assertEqual(coverage.CoverageReport().build(s)['bundles'][0]['status'],'PENDING_PHYSICAL_PROVIDER')

    def test_missing_member_and_partial_extent_stay_explicit(self):
        for partial in (False,True):
            s,b,reader=self.prepare()
            if partial:b.headLocs[0]=(5.,0.)
            else:s.retained_physical_streets.pop(26)
            module.PhysicalBundleInventory().read(s,(b,),reader)
            self.assertEqual(s.retained_physical_bundles[7]['status'],'UNSUPPORTED')
            self.assertTrue(s.retained_physical_bundles[7]['reason'])

    def test_duplicate_bundle_ids_cannot_be_counted_or_merged(self):
        s,b,reader=self.prepare();module.PhysicalBundleInventory().read(s,(b,),reader)
        with self.assertRaisesRegex(ValueError,'duplicate'):
            module.PhysicalBundleInventory().read(s,(b,),reader)
        s.bundles[7]=s.retained_physical_bundles[7]
        with self.assertRaisesRegex(ValueError,'active and retained'):module.PhysicalStreetInventory.view(s)
        with self.assertRaisesRegex(ValueError,'both active and retained'):coverage.CoverageReport().build(s)

    def test_report_loader_keeps_retained_roles_and_rejects_duplicate_ids(self):
        s,b,reader=self.prepare();module.PhysicalBundleInventory().read(s,(b,),reader)
        report={key:list(getattr(s,key).values()) for key in
                ('streets','bundles','retained_physical_streets','retained_physical_bundles')}
        report['junction_origins']=s.junction_origins
        view=module.PhysicalStreetInventory.from_report(report)
        self.assertEqual(set(view.streets),{25,26});self.assertEqual(set(view.bundles),{7})
        report['retained_physical_bundles']*=2
        with self.assertRaisesRegex(ValueError,'Duplicate original IDs'):
            module.PhysicalStreetInventory.from_report(report)


if __name__=='__main__':unittest.main()
