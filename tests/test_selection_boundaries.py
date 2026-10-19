"""Selection trims retain original ownership; they never invent free ends."""
from copy import deepcopy
from dataclasses import replace
import unittest

from test_original_span_plan import m
from test_physical_ports import fixture
from test_original_boundary_runs import fixture as run_fixture


def declare(source,street,section,end,trim=2.5):
    source.junction_origins.setdefault('selection_boundaries',[]).append(dict(
        endpoint=dict(street_id=street,section_id=section,end=end),trim_m=trim))


class SelectionBoundaryTests(unittest.TestCase):
    def planner(self):
        return m.OriginalSpanPlanner(('residential',),section_scopes=True,boundary_runs=True)

    def test_approaches_in_both_source_directions_keep_original_identity(self):
        source=fixture()
        declare(source,2,20,0);declare(source,3,30,1,4.)
        before=deepcopy(source);result=self.planner().prepare(source)
        self.assertFalse(result.issues);self.assertEqual(source,before)
        component,=result.components
        routes={r.span.street_id:r for r in component.network.routes if r.span.street_id in (2,3)}
        for sid,end in ((2,0),(3,1)):
            route=routes[sid]
            self.assertEqual(route.kind,'SELECTION_APPROACH')
            bound=(route.span.start,route.span.end)[end]
            self.assertEqual(bound.owner_kind,'SELECTION_BOUNDARY')
            self.assertIsNone(bound.junction_id)
            self.assertEqual((bound.endpoint.street_id,bound.endpoint.section_id,bound.endpoint.end),(sid,sid*10,end))
            self.assertEqual(component.ports.selection_boundary(bound.endpoint).trim_m,2.5 if sid==2 else 4.)
        self.assertFalse(component.network.unresolved)

    def test_unattached_trimmed_run_keeps_all_source_controls_and_end_roles(self):
        for ends in ((0,),(1,),(0,1)):
            source=run_fixture(('residential',))
            source.streets[90].update(status='SUPPORTED',reason=None)
            for end in ends:declare(source,90,900,end)
            before=deepcopy(source);result=self.planner().prepare(source,standalone=True)
            self.assertFalse(result.issues);self.assertEqual(source,before)
            component,=result.components;run,=component.network.boundary_runs
            self.assertEqual(run.span_identity,(90,900,900))
            self.assertEqual(run.roles,tuple(6 if e in ends else 1 for e in (0,1)))
            self.assertEqual(run.owners,tuple(900 if e in ends else -1 for e in (0,1)))
            self.assertEqual(run.next_sections,(-1,-1))
            self.assertFalse(component.network.unresolved)

    def test_invalid_trim_values_are_not_accepted(self):
        for value in (0.,-1.,float('inf'),float('nan'),True,'2'):
            source=fixture();declare(source,2,20,0,value)
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'trim'):
                self.planner().reader.prepare(source,street_ids=(2,))

    def test_duplicate_or_existing_ownership_is_rejected(self):
        for kind in ('duplicate','junction','remote','free'):
            source=fixture();declare(source,2,20,0)
            if kind=='duplicate':declare(source,2,20,0)
            if kind=='junction':declare(source,2,20,1)
            if kind=='remote':source.streets[2]['endpoints']=(999,21)
            if kind=='free':source.junction_origins['free_endpoints']=[dict(
                endpoint=dict(street_id=2,terminal_id=20,terminal_kind='Section',end=0),location=(0.,30.))]
            with self.subTest(kind=kind),self.assertRaisesRegex(ValueError,'conflict'):
                self.planner().reader.prepare(source,street_ids=(2,))

    def test_interior_section_cut_requires_explicit_split(self):
        source=fixture();declare(source,1,11,1)
        with self.assertRaises(ValueError):
            self.planner().reader.prepare(source,street_ids=(1,))

    def test_paired_interior_cuts_keep_one_complete_street_and_all_sections(self):
        source=run_fixture(('residential',)*3)
        source.streets[90].update(status='SUPPORTED',reason=None)
        for section,end,trim in ((900,1,2.5),(901,0,3.),(901,1,1.),(902,0,4.)):
            declare(source,90,section,end,trim)
        before=deepcopy(source);result=self.planner().prepare(source,standalone=True)
        self.assertFalse(result.issues);self.assertEqual(source,before)
        component,=result.components;plan=component.ports
        self.assertEqual(plan.source_streets,(90,))
        self.assertEqual([s.section_ids for s in plan.spans],[(900,),(901,),(902,)])
        self.assertEqual([(s.first_control,s.last_control,s.source_control_count) for s in plan.spans],
                         [(0,2,7),(2,4,7),(4,6,7)])
        self.assertEqual([r.roles for r in component.network.boundary_runs],[(1,6),(6,6),(6,1)])
        for boundary in plan.selection_boundaries:
            peer=plan.selection_peer(boundary.endpoint)
            self.assertEqual(plan.selection_peer(peer),boundary.endpoint)
            self.assertEqual(peer.street_id,90)
        self.assertFalse(component.network.unresolved)
        dimensional=m.OriginalSpanPlanner(('residential',),section_scopes=True,
            boundary_runs=True,dimension_runs=True).prepare(source,standalone=True)
        self.assertFalse(dimensional.issues)
        self.assertEqual(dimensional.components[0].ports.spans,plan.spans)

    def test_pair_cannot_skip_an_intermediate_section_or_replace_a_seam(self):
        for defect in ('unpaired','nonconsecutive','seam','category'):
            source=run_fixture(('residential',)*3)
            source.streets[90].update(status='SUPPORTED',reason=None)
            declare(source,90,900,1)
            if defect!='unpaired':declare(source,90,902 if defect=='nonconsecutive' else 901,0)
            if defect=='seam':
                source.streets[90]['items'].insert(1,dict(kind='SymLane',id=4,
                    surface_treatment='COMPATIBLE_SEAM',points=[]))
            if defect=='category':source.streets[90]['items'][1]['category']='footway'
            with self.subTest(defect=defect),self.assertRaisesRegex(ValueError,'selection cut|paired declarations'):
                self.planner().reader.prepare(source,street_ids=(90,))

    def test_terminal_selection_has_no_interior_peer(self):
        source=run_fixture(('residential',));source.streets[90].update(status='SUPPORTED',reason=None)
        declare(source,90,900,0)
        plan=self.planner().reader.prepare(source,street_ids=(90,))
        self.assertIsNone(plan.selection_peer(plan.selection_boundaries[0].endpoint))

    def test_missing_or_duplicate_declaration_cannot_silently_become_free(self):
        source=fixture();declare(source,2,20,0)
        plan=self.planner().reader.prepare(source,street_ids=(2,))
        endpoint=plan.selection_boundaries[0].endpoint
        for declarations in ((),plan.selection_boundaries*2):
            with self.assertRaisesRegex(ValueError,'one exact original endpoint'):
                replace(plan,selection_boundaries=declarations).selection_boundary(endpoint)


if __name__=='__main__':unittest.main()
