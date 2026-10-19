"""A selected middle road run keeps exact source/end identity and scope."""
from copy import deepcopy
from types import SimpleNamespace
import unittest
from test_original_span_plan import m


def fixture(categories=('footway', 'residential', 'footway')):
    items=[dict(kind='Section',id=900+i,category=category,width=6. if category=='residential' else 2.,
                offset=0.,points=((20.*i,0.),(20.*i+10.,1.),(20.*(i+1),0.)))
           for i,category in enumerate(categories)]
    return SimpleNamespace(streets={90:dict(id=90,status='UNSUPPORTED',
        reason='Reserved component/category: footway',endpoints=(None,None),bundle_id=None,items=items)},
        bundles={},junction_origins={'physical':[]})


def cross_fixture(own_end,peer_end,y=0.):
    source=fixture(('residential',));road=source.streets[90]
    road.update(status='SUPPORTED',reason=None)
    road['items'][0]['points']=((0.,y),(10.,y+1.),(20.,y))
    location=road['items'][0]['points'][0 if own_end==0 else -1]
    points=(location,(location[0]+10.,y+10.),(location[0]+20.,y+20.))
    if peer_end:points=points[::-1]
    source.streets[91]=dict(id=91,status='UNSUPPORTED',reason='Reserved component/category: footway',
        endpoints=(None,None),bundle_id=None,items=[dict(kind='Section',id=910,category='footway',
                                                      width=2.,offset=0.,points=points)])
    source.junction_origins['physical']=[dict(junction_id=2000,location=location,ring=[])]
    source.junction_origins['seams']=[dict(junction_id=2000,location=location,endpoints=[
        dict(street_id=90,terminal_kind='Section',terminal_id=900,end=own_end),
        dict(street_id=91,terminal_kind='Section',terminal_id=910,end=peer_end)])]
    return source


class BoundaryRunTests(unittest.TestCase):
    def planner(self,**kwargs):
        return m.OriginalSpanPlanner(('residential',),section_scopes=True,boundary_runs=True,**kwargs)

    def test_two_category_boundaries_keep_complete_source_and_ids(self):
        source=fixture();before=deepcopy(source)
        result=self.planner().prepare(source)
        self.assertEqual(source,before);self.assertFalse(result.issues)
        component,=result.components;run,=component.network.boundary_runs
        self.assertEqual(run.span_identity,(90,901,901))
        self.assertEqual((run.roles,run.owners,run.next_sections),((4,4),(90,90),(900,902)))
        self.assertEqual(component.ports.source_streets,(90,))
        self.assertEqual(len(component.ports.spans),3)
        self.assertFalse(component.ports.sites);self.assertFalse(component.network.unresolved)

    def test_both_authoring_directions_of_free_category_run(self):
        for categories,expected in ((('residential','footway'),(1,4)),(('footway','residential'),(4,1))):
            result=self.planner().prepare(fixture(categories))
            run,=result.components[0].network.boundary_runs
            self.assertEqual(run.roles,expected)

    def test_recorded_free_end_on_excluded_first_section_keeps_category_scope(self):
        for categories in (('footway','residential'),('residential','footway')):
            source=fixture(categories)
            source.junction_origins['free_endpoints']=[dict(endpoint=dict(street_id=90,end=end,
                terminal_kind='Section',terminal_id=source.streets[90]['items'][0 if end==0 else -1]['id']),
                location=source.streets[90]['items'][0 if end==0 else -1]['points'][0 if end==0 else -1])
                for end in (0,1)]
            before=deepcopy(source);result=self.planner().prepare(source)
            self.assertFalse(result.issues)
            self.assertEqual(len(result.components[0].network.boundary_runs),1)
            self.assertEqual(source,before)

    def test_two_disconnected_runs_keep_one_editable_street(self):
        source=fixture(('residential','footway','residential'))
        result=self.planner().prepare(source);component,=result.components
        self.assertEqual([r.span_identity for r in component.network.boundary_runs],[(90,900,900),(90,902,902)])
        self.assertFalse(component.network.unresolved)

    def test_legacy_or_explicit_site_selection_does_not_add_unrequested_sources(self):
        source=fixture()
        self.assertFalse(m.OriginalSpanPlanner(('residential',),section_scopes=True).prepare(source).components)
        self.assertFalse(self.planner().prepare(source,physical_ids=()).components)
        with self.assertRaisesRegex(ValueError,'category scopes'):
            m.OriginalSpanPlanner(('residential',),boundary_runs=True)

    def test_unresolved_owner_is_never_reclassified_as_free(self):
        source=fixture(('residential','footway'));source.streets[90]['endpoints']=(123,None)
        result=self.planner().prepare(source)
        self.assertFalse(result.components)

    def test_bundle_and_dimensional_runs_remain_explicitly_pending(self):
        for kind in ('bundle','dimensions'):
            source=fixture(('footway','residential','residential','footway'))
            if kind=='bundle':source.streets[90]['bundle_id']=5
            else:source.streets[90]['items'][2]['width']=8.
            result=self.planner().prepare(source);component,=result.components
            self.assertFalse(component.network.boundary_runs)
            self.assertTrue(component.network.unresolved)
            self.assertTrue(any(i['kind']=='BOUNDARY_RUN' for i in result.issues))

    def test_all_cross_street_endpoint_directions_retain_the_original_node(self):
        for own_end in (0,1):
            for peer_end in (0,1):
                source=cross_fixture(own_end,peer_end);before=deepcopy(source)
                result=self.planner(category_connections=True,seam_connections=True).prepare(source)
                self.assertFalse(result.issues);component,=result.components
                run,=component.network.boundary_runs
                self.assertEqual(run.roles[own_end],5)
                self.assertEqual((run.owners[own_end],run.next_sections[own_end]),(2000,910))
                self.assertEqual(component.ports.source_streets,(90,))
                self.assertEqual(source,before)
                self.assertFalse(component.ports.sites)


if __name__=='__main__':unittest.main()
