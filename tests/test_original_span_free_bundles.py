"""Free Bundle integration keeps complete inputs and once-owned native output."""
from copy import deepcopy
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from test_bundle_topology import sample

package=ModuleType('_free_bundle_preparation_tests')
package.__path__=[str(Path(__file__).resolve().parents[1]/'way/street_network')]
sys.modules[package.__name__]=package
Planner=importlib.import_module(package.__name__+'.original_span_plan').OriginalSpanPlanner
Coverage=importlib.import_module(package.__name__+'.original_span_coverage').OriginalSurfaceCoverage
gn=ModuleType('_free_bundle_bindings_tests')
gn.__path__=[str(Path(__file__).resolve().parents[1]/'way/gn/editable_input')]
sys.modules[gn.__name__]=gn
rows=importlib.import_module(gn.__name__+'.free_bundle_bindings').path_rows


def fixture(reverse=False):
    bundle,streets=sample(reverse=reverse,junction=False);del streets[3]
    for street in streets.values():
        street.update(status='SUPPORTED',bundle_id=7,endpoints=(None,None))
        for section in street['items']:section.update(category='residential',width=5.,offset=0.)
    return SimpleNamespace(streets=streets,bundles={7:bundle},junction_origins={})


class FreeBundleTests(unittest.TestCase):
    def planner(self,enabled=True):
        return Planner(('residential',),free_bundles=enabled,section_scopes=True,boundary_runs=True)

    def test_complete_multisection_reverse_preserves_inputs_and_route_ownership(self):
        for reverse in (False,True):
            source=fixture(reverse);before=deepcopy(source)
            result=self.planner().prepare(source,standalone=True)
            self.assertFalse(result.issues);component,=result.components
            self.assertEqual(source,before)
            self.assertEqual(component.ports.source_streets,(1,2))
            self.assertEqual([r.kind for r in component.network.routes],['FREE_BUNDLE']*2)
            self.assertFalse(component.network.unresolved)
            bundle,=component.free_bundles
            self.assertEqual(bundle.members[1].forward,not reverse)
            table=rows(list(source.streets.values()),component.free_bundles,{10:(.12,0.),11:(.12,0.),20:(.12,0.)})
            self.assertEqual([(r['st_path_id'],r['st_path_start'],r['st_path_edge']) for r in table],
                             [(10,0,0),(11,3,2),(20,6,4)])
            self.assertEqual(table[0]['st_path_end_seam'],3)
            self.assertEqual(table[1]['st_path_start_seam'],2)
            self.assertEqual(table[0]['st_path_chain_reversed'],reverse)
            owners=Coverage.free_bundle_owners(tuple(tuple(r['st_path_'+key] for key in
                ('bundle_id','street_id','id','chain_index','chain_order','chain_owner')) for r in table))
            proof=Coverage(owners)
            missing=proof.evaluate([])
            self.assertEqual(missing['failed_streets'],[1,2]);self.assertEqual(missing['missing_bundles'],[7])
            self.assertTrue(proof.evaluate([(8,7)])['complete'])
            with self.assertRaises(ValueError):proof.evaluate([(5,0)])

    def test_old_asset_does_not_silently_emit_independent_carriageways(self):
        result=self.planner(False).prepare(fixture(),standalone=True)
        self.assertFalse(result.components);self.assertTrue(result.issues)

    def test_partial_attached_and_order_exchange_stay_explicit(self):
        for defect in ('partial','attached','order','internal_cut'):
            source=fixture()
            if defect=='partial':source.bundles[7]['head_locations']=((5.,0.),(5.,10.))
            elif defect=='attached':source.streets[1]['endpoints']=(9,None)
            elif defect=='order':
                source.bundles[7]['tail']=(2,1)
                source.bundles[7]['tail_locations']=((40.,10.),(40.,0.))
            else:
                source.junction_origins['selection_boundaries']=[dict(endpoint=dict(street_id=1,section_id=s,end=e),trim_m=1.)
                    for s,e in ((10,1),(11,0))]
            result=self.planner().prepare(source,standalone=True)
            self.assertFalse(any(c.free_bundles for c in result.components),defect)
            self.assertTrue(any(i['kind']=='BUNDLE' for i in result.issues),defect)

    def test_incompatible_profile_does_not_rewrite_rise(self):
        source=fixture();component,=self.planner().prepare(source,standalone=True).components
        with self.assertRaisesRegex(ValueError,'shared profile'):
            rows(list(source.streets.values()),component.free_bundles,{10:(.12,0.),11:(.12,0.),20:(.08,0.)})

    def test_whole_bundle_does_not_request_duplicate_separator_coverage(self):
        from_module=importlib.import_module(package.__name__+'.original_span_separators')
        source=fixture();component,=self.planner().prepare(source,standalone=True).components
        result=from_module.OriginalSpanSeparatorPlanner(('residential',)).prepare(source,component,free_ends=True)
        self.assertFalse(result.closed);self.assertFalse(result.free_ended);self.assertFalse(result.pending)

    def test_dimension_bindings_keep_the_original_seam_and_partner_chain(self):
        source=fixture();component,=self.planner().prepare(source,standalone=True).components
        transition=SimpleNamespace(before_section=10,after_section=11,source_id=77,
                                   minimum_half_length=1.5,maximum_edge_slope=.25)
        table=rows(list(source.streets.values()),component.free_bundles,
                   {10:(.12,0.),11:(.12,0.),20:(.12,0.)},transitions=(transition,))
        self.assertEqual(table[0]['st_path_end_transition'],77)
        self.assertEqual(table[1]['st_path_start_transition'],77)
        self.assertEqual(table[0]['st_path_end_transition_half'],1.5)
        self.assertEqual(table[0]['st_path_end_transition_partner'],table[1]['st_path_chain_index'])
        self.assertEqual(table[0]['st_path_end_seam'],table[1]['st_path_start'])

    def test_incomplete_duplicate_or_conflicting_coverage_rejects(self):
        valid=[(7,1,10,0,0,True),(7,1,11,0,1,True),(7,2,20,1,0,False)]
        for bad in (valid[:-1],valid+valid[:1],[(7,1,10,0,1,True)]+valid[1:],
                    valid[:2]+[(7,2,20,0,0,False)],valid[:2]+[(7,2,20,1,0,True)]):
            with self.subTest(rows=bad),self.assertRaises(ValueError):Coverage.free_bundle_owners(bad)


if __name__=='__main__':unittest.main()
