"""Import selection must retain original ownership and explicit frontiers."""
from copy import deepcopy
import importlib
from pathlib import Path
import sys
from types import ModuleType
import unittest
from test_physical_ports import fixture, mixed_fixture, seam_fixture

package=ModuleType('_original_span_import_tests')
package.__path__=[str(Path(__file__).resolve().parents[1]/'way/street_network')]
sys.modules[package.__name__]=package
m=importlib.import_module(package.__name__+'.original_span_plan')


class OriginalSpanImportTests(unittest.TestCase):
    def setUp(self):self.planner=m.OriginalSpanPlanner(('residential',))

    def test_category_seam_is_handoff_not_a_failed_patch_or_free_end(self):
        source=seam_fixture();source.streets[8]['items'][0]['category']='footway'
        source.streets[8].update(status='UNSUPPORTED',reason='Reserved component/category: footway')
        before=deepcopy(source)
        planner=m.OriginalSpanPlanner(('residential',),seam_connections=True,section_scopes=True,category_connections=True)
        current=planner.prepare(source)
        self.assertFalse(current.issues);component,=current.components
        self.assertFalse(component.frontiers);self.assertFalse(component.network.unresolved)
        self.assertNotIn(2000,{s.junction_id for s in component.ports.sites})
        self.assertEqual(next(r.kind for r in component.network.routes if r.span.street_id==4),'CATEGORY_APPROACH')
        self.assertEqual(source,before)
        legacy=m.OriginalSpanPlanner(('residential',),seam_connections=True,section_scopes=True).prepare(source)
        self.assertTrue(legacy.components[0].network.unresolved)

    def test_seam_provider_capability_keeps_complete_network_and_original_streets(self):
        source=seam_fixture();before=deepcopy(source)
        current=m.OriginalSpanPlanner(('residential',),seam_connections=True).prepare(source)
        self.assertFalse(current.issues);self.assertEqual(len(current.components),1)
        component,=current.components
        self.assertFalse(component.frontiers);self.assertFalse(component.network.unresolved)
        self.assertEqual({s.junction_id for s in component.ports.sites},{1001,1002,21,22,2000})
        self.assertEqual(source,before)
        self.assertEqual({r.span.street_id:r.kind for r in component.network.routes if r.span.street_id in (4,8)},
                         {4:'CORRIDOR',8:'APPROACH'})
        self.assertTrue(self.planner.prepare(source).components[0].network.unresolved)

    def test_shared_logical_owners_have_one_surface_and_keep_original_ids(self):
        from test_shared_port_clusters import fixture as shared_fixture
        source=shared_fixture();before=deepcopy(source)
        legacy=m.OriginalSpanPlanner(('residential',),internal_logical_ends=True).prepare(source)
        self.assertEqual({r['id'] for r in legacy.issues if r['kind']=='LOGICAL'},{5000,5001})
        planner=m.OriginalSpanPlanner(('residential',),internal_logical_ends=True,shared_logical_owners=True)
        current=planner.prepare(source)
        self.assertFalse(current.issues)
        self.assertEqual(len(current.components),1)
        network=current.components[0].network
        cluster,=network.logical_clusters
        self.assertEqual(cluster.owner_ids,(5000,5001))
        self.assertEqual(cluster.physical_junction_ids,(10,11,12))
        self.assertFalse(network.ordinary_site_ids)
        self.assertEqual(sum(r.kind=='LOGICAL_CLUSTER' for r in network.routes),2)
        self.assertEqual(source,before)

    def test_shared_owner_capability_requires_internal_footprint_support(self):
        with self.assertRaisesRegex(ValueError,'retained internal ends'):
            m.OriginalSpanPlanner(('residential',),shared_logical_owners=True)

    def test_internal_logical_ends_require_the_matching_asset_capability(self):
        from test_port_clusters import fixture as logical_fixture, retain_internal_connectors
        source=logical_fixture();retain_internal_connectors(source)
        before=deepcopy(source)
        legacy=self.planner.prepare(source)
        self.assertTrue(any(row.get('id')==5000 and 'external port' in row['reason'] for row in legacy.issues))
        current=m.OriginalSpanPlanner(('residential',),internal_logical_ends=True).prepare(source)
        self.assertFalse(current.issues)
        self.assertEqual(source,before)
        self.assertEqual(len(current.components),1)
        cluster,=current.components[0].network.logical_clusters
        self.assertEqual(cluster.footprint_spans,((1,10,10),))
        self.assertEqual(len(cluster.internal_groups),2)

    def test_preserved_free_origin_restores_approach_without_inventing_a_junction(self):
        source=fixture();source.streets[4]['endpoints']=(1001,999)
        before=self.planner.prepare(source)
        route=next(r for c in before.components for r in c.network.routes if r.span.street_id==4)
        self.assertEqual(route.kind,'REMOTE_CONNECTION')
        source.junction_origins['free_endpoints']=[dict(
            endpoint=dict(street_id=4,terminal_id=40,terminal_kind='Section',end=1),
            location=(-60.,-20.))]
        result=self.planner.prepare(source)
        route=next(r for c in result.components for r in c.network.routes if r.span.street_id==4)
        self.assertEqual(route.kind,'APPROACH')
        self.assertEqual(route.span.identity,(4,40,40))
        self.assertFalse(result.components[0].network.unresolved)
        self.assertEqual({s.junction_id for c in result.components for s in c.ports.sites},
                         {s.junction_id for c in before.components for s in c.ports.sites})

    def test_automatic_import_covers_all_original_sites_without_mutating_sources(self):
        source=fixture();before=deepcopy(source)
        result=self.planner.prepare(source)
        self.assertEqual(source,before)
        self.assertFalse(result.issues)
        self.assertEqual(len(result.components),1)
        component=result.components[0]
        self.assertEqual({s.junction_id for s in component.ports.sites},{21,22,1001,1002})
        self.assertFalse(component.network.unresolved)
        self.assertEqual(component.ports.source_streets,tuple(range(1,8)))

    def test_sections_of_one_street_never_get_duplicate_editable_owners(self):
        result=self.planner.prepare(fixture(),expand=False)
        self.assertEqual(len(result.components),1)
        self.assertFalse(result.components[0].network.unresolved)

    def test_automatic_mixed_ring_retains_deferred_identity_without_a_false_frontier(self):
        source=mixed_fixture();before=deepcopy(source)
        result=self.planner.prepare(source)
        self.assertEqual(source,before)
        self.assertFalse(result.issues)
        self.assertEqual(len(result.components),1)
        component=result.components[0]
        self.assertFalse(component.frontiers)
        self.assertFalse(component.network.unresolved)
        self.assertEqual(component.ports.source_streets,tuple(range(1,8)))
        deferred,=component.ports.deferred
        self.assertEqual((deferred.junction_id,deferred.connector_id,deferred.category),(1001,999,'footway'))
        self.assertEqual(deferred.endpoint.street_id,99)

    def test_legacy_asset_keeps_explicit_mixed_ring_frontier(self):
        planner=m.OriginalSpanPlanner(('residential',),defer_physical_categories=False)
        result=planner.prepare(mixed_fixture())
        self.assertNotIn(1001,{s.junction_id for c in result.components for s in c.ports.sites})
        self.assertTrue(any(row.get('id')==1001 and 'Street 99 is unsupported' in row['reason']
                            for row in result.issues))

    def test_explicit_selection_keeps_remote_connections_pending(self):
        result=self.planner.prepare(fixture(),physical_ids=(),interior_ids=(21,),expand=False)
        self.assertEqual(len(result.components),1)
        self.assertTrue(result.components[0].network.unresolved)
        self.assertEqual({r.span.street_id for r in result.components[0].network.routes},{1,2})

    def test_bad_explicit_seed_rejects_instead_of_changing_selection(self):
        with self.assertRaises(ValueError):self.planner.prepare(fixture(),physical_ids=(9999,))

    def test_automatic_unsupported_ring_preserves_reason_and_other_sites(self):
        source=fixture();source.streets[7].update(status='UNSUPPORTED',reason='Unsupported grade test')
        result=self.planner.prepare(source)
        self.assertTrue(any('Unsupported grade test' in row['reason'] for row in result.issues))
        self.assertEqual({s.junction_id for c in result.components for s in c.ports.sites},{21,22,1001})
        self.assertTrue(result.components[0].network.unresolved)

    def test_unknown_compound_does_not_silently_disappear(self):
        with self.assertRaisesRegex(ValueError,'absent'):
            self.planner.prepare(fixture(),compound_spans=((99,990,990),))

    def short_fixture(self, gap=8., bend=0.):
        source=fixture();items=source.streets[1]['items']
        for item in items:
            if item['id']==11 and item['kind']=='Section':item['points']=((0.,0.),(gap*.5,bend))
            if item['id']==12 and item['kind']=='Section':item['points']=((gap*.5,bend),(gap,0.))
            if item['id']==22 and item['kind']=='Intersection':item['location']=(gap,0.)
            if item['id']==13 and item['kind']=='Section':item['points']=((gap,0.),(60.,-.5),(80.,0.))
        source.streets[3]['items'][0]['points']=((gap,0.),(gap,-30.))
        return source

    def test_short_span_policy_selects_existing_compound_without_source_changes(self):
        source=self.short_fixture();before=deepcopy(source)
        planner=m.OriginalSpanPlanner(('residential',),category_properties={'residential':{'radius':2.}})
        result=planner.prepare(source)
        self.assertEqual(source,before)
        self.assertEqual([c.span_identity for c in result.components[0].network.compounds],[(1,11,12)])
        self.assertEqual(result.components[0].compound_decisions[0]['candidate_limit_m'],10.)

    def test_short_span_policy_uses_full_control_path_length_and_strict_limit(self):
        planner=m.OriginalSpanPlanner(('residential',),category_properties={'residential':{'radius':2.}})
        for source in (self.short_fixture(10.),self.short_fixture(8.,10.)):
            self.assertFalse(planner.prepare(source).components[0].network.compounds)

    def test_authored_terminal_radii_drive_the_same_short_span_policy(self):
        source=self.short_fixture();before=deepcopy(source);seen=[]
        def radius(inventory,end):
            self.assertIs(inventory,source)
            seen.append(end.endpoint)
            return 1.5 if end.endpoint.end==0 else 3.5
        planner=m.OriginalSpanPlanner(('residential',),radius_for_end=radius)
        component,=planner.prepare(source).components
        self.assertEqual([c.span_identity for c in component.network.compounds],[(1,11,12)])
        self.assertEqual(component.compound_decisions[0]['candidate_limit_m'],11.)
        self.assertTrue(seen)
        self.assertEqual(source,before)

    def test_authored_radii_override_category_defaults_without_shortening_them(self):
        source=self.short_fixture(10.);before=deepcopy(source)
        planner=m.OriginalSpanPlanner(('residential',),
            category_properties={'residential':{'radius':100.}},radius_for_end=lambda inventory,end:2.)
        component,=planner.prepare(source).components
        self.assertFalse(component.network.compounds)
        self.assertFalse(component.compound_decisions)
        self.assertEqual(source,before)

    def test_explicit_compound_is_not_duplicated_by_short_span_policy(self):
        planner=m.OriginalSpanPlanner(('residential',),category_properties={'residential':{'radius':2.}})
        result=planner.prepare(self.short_fixture(),compound_spans=((1,11,12),))
        self.assertEqual(len(result.components[0].network.compounds),1)
        self.assertFalse(result.components[0].compound_decisions)

    def test_unrelated_scope_issue_does_not_mark_accepted_component_partial(self):
        source=fixture();component=self.planner.prepare(source).components[0]
        outside=dict(kind='PHYSICAL',id=999,reason='Unsupported road elsewhere')
        self.assertFalse(self.planner.issues_for(source,component,(outside,)))
        source.junction_origins['components']=[dict(physical_junction_ids=(1001,21),logical_junction_ids=(1234,))]
        related=dict(kind='LOGICAL',id=1234,reason='Conflicting physical ownership')
        self.assertEqual(self.planner.issues_for(source,component,(outside,related)),(related,))

    def free_source(self):
        source=fixture()
        street=deepcopy(source.streets[7]);street['id']=90
        street['endpoints']=(None,None);street['bundle_id']=None
        street['items']=[dict(kind='Section',id=900,category='residential',width=4.,offset=0.,
                              points=((120.,0.),(130.,2.),(140.,0.)))]
        source.streets[90]=street
        return source

    def test_unattached_street_keeps_its_original_ids_without_a_junction(self):
        source=self.free_source();before=deepcopy(source)
        result=self.planner.prepare(source,standalone=True)
        component=next(c for c in result.components if c.ports.source_streets==(90,))
        self.assertEqual(source,before)
        self.assertFalse(component.ports.sites)
        self.assertEqual(component.network.routes[0].kind,'STANDALONE')
        span=component.ports.spans[0]
        self.assertEqual(span.section_ids,(900,))
        self.assertEqual((span.start.owner_kind,span.end.owner_kind),('FREE','FREE'))
        self.assertFalse(component.network.unresolved)

    def test_standalone_capability_and_explicit_site_selection_are_respected(self):
        source=self.free_source()
        for options in ({},{'standalone':True,'physical_ids':(1001,)}):
            self.assertNotIn(90,{sid for c in self.planner.prepare(source,**options).components for sid in c.ports.source_streets})

    def test_actual_physical_owner_is_not_replaced_by_a_free_end(self):
        source=self.free_source()
        site=deepcopy(source.junction_origins['physical'][0]);site['junction_id']=1999
        endpoint=site['ring'][0]['endpoints'][0]
        endpoint.update(street_id=90,end=0,terminal_kind='Section',terminal_id=900)
        source.junction_origins['physical'].append(site)
        result=self.planner.prepare(source,standalone=True)
        self.assertFalse(any(c.ports.source_streets==(90,) for c in result.components))
        self.assertTrue(any(row.get('id')==90 for row in result.issues))

    def test_unattached_bundle_or_dimension_change_is_explicitly_pending(self):
        for kind in ('bundle','dimensions'):
            source=self.free_source()
            if kind=='bundle':source.streets[90]['bundle_id']=123
            else:
                section=deepcopy(source.streets[90]['items'][0]);section.update(id=901,width=5.,points=((140.,0.),(160.,0.)))
                source.streets[90]['items'].append(section)
            result=self.planner.prepare(source,standalone=True)
            self.assertFalse(any(c.ports.source_streets==(90,) for c in result.components))
            self.assertTrue(any(row.get('id')==90 for row in result.issues))

    def test_dimensional_run_requires_capability_and_explicit_seam(self):
        source=self.free_source();first=source.streets[90]['items'][0]
        section=deepcopy(first);section.update(id=901,width=5.,points=(first['points'][-1],(160.,0.)))
        item=dict(kind='AuthoredTransition',id=12,points=[],surface_treatment='WIDTH_TRANSITION',
            transition=dict(before_section=900,after_section=901,source_kind='AuthoredTransition',source_id=12,
                            minimum_half_length=1.5,maximum_edge_slope=.3))
        source.streets[90]['items'].extend((item,section));before=deepcopy(source)
        ordinary=self.planner.prepare(source,standalone=True)
        self.assertFalse(any(c.ports.source_streets==(90,) for c in ordinary.components))
        planner=m.OriginalSpanPlanner(('residential',),dimension_runs=True)
        admitted=planner.prepare(source,standalone=True)
        component,=[c for c in admitted.components if c.ports.source_streets==(90,)]
        self.assertEqual(component.network.routes[0].kind,'STANDALONE')
        self.assertEqual(component.ports.spans[0].section_ids,(900,901));self.assertEqual(source,before)
        source.streets[90]['items'].remove(item)
        rejected=planner.prepare(source,standalone=True)
        self.assertFalse(any(c.ports.source_streets==(90,) for c in rejected.components))
        self.assertTrue(any('undeclared' in row['reason'] for row in rejected.issues))


if __name__=='__main__':unittest.main()
