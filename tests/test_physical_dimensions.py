"""Prepared dimensional seam ownership, independent of Blender geometry."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

spec = importlib.util.spec_from_file_location('_physical_dimension_tests', Path(__file__).resolve().parents[1]/'way/street_network/physical_dimensions.py')
module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module; spec.loader.exec_module(module)


def snapshot():
    a = dict(kind='Section', id=10, points=[(0,0,0), (20,0,0)], width=8., offset=0.)
    b = dict(kind='Section', id=20, points=[(20,0,0), (60,0,0)], width=10., offset=-1.)
    item = dict(kind='SideLane', id=2, points=[], surface_treatment='WIDTH_TRANSITION',
        transition=dict(before_section=10, after_section=20, source_kind='SideLane', source_id=2,
                        minimum_half_length=1.5, maximum_edge_slope=.3))
    return NS(streets={7:dict(id=7, status='SUPPORTED', items=[a,item,b])})


class PhysicalDimensionsTest(unittest.TestCase):
    def test_authored_transition_keeps_its_own_provenance(self):
        source=snapshot();item=source.streets[7]['items'][1]
        item['kind']='AuthoredTransition';item['transition']['source_kind']='AuthoredTransition'
        before=deepcopy(source)
        row,=module.PreparedPhysicalDimensions().read(source,[7])
        self.assertEqual((row.source_kind,row.source_id),('AuthoredTransition',2))
        self.assertEqual(source,before)
        item['transition']['source_kind']='SymLane'
        with self.assertRaisesRegex(ValueError,'identity'):
            module.PreparedPhysicalDimensions().read(source,[7])

    def test_category_seam_retains_identity_without_rendering_the_excluded_taper(self):
        from test_physical_ports import scoped_fixture,m as ports
        source=scoped_fixture();before=deepcopy(source)
        plan=ports.PreparedPhysicalPorts(('residential',),section_scopes=True).prepare(source,physical_ids=(1001,))
        reader=module.PreparedPhysicalDimensions()
        self.assertEqual(reader.read_spans(source,plan.spans,source_plan=plan),())
        self.assertEqual(source,before)
        spans=tuple(replace(s,end=replace(s.end,owner_kind='FREE'))
                    if s.street_id==4 and s.in_scope else s for s in plan.spans)
        with self.assertRaises(ValueError):reader.read_spans(source,spans,source_plan=plan)

    def test_original_ownership_and_input_preservation(self):
        source = snapshot(); before = deepcopy(source)
        row, = module.PreparedPhysicalDimensions().read(source, [7])
        self.assertEqual((row.street_id,row.before_section,row.after_section,row.source_kind,row.source_id), (7,10,20,'SideLane',2))
        self.assertEqual((row.minimum_half_length,row.maximum_edge_slope), (1.5,.3))
        self.assertEqual(source, before)

    def test_absent_transition_is_not_invented(self):
        source = snapshot(); source.streets[7]['items'].pop(1)
        with self.assertRaisesRegex(ValueError, 'undeclared'):
            module.PreparedPhysicalDimensions().read(source, [7])

    def test_dormant_declaration_survives_equal_dimensions(self):
        source = snapshot(); source.streets[7]['items'][-1].update(width=8.,offset=0.)
        self.assertEqual(len(module.PreparedPhysicalDimensions().read(source, [7])),1)

    def test_source_identity_and_adjacency(self):
        for key, value in (('source_id',99),('before_section',20),('after_section',99)):
            with self.subTest(key=key):
                source = snapshot(); source.streets[7]['items'][1]['transition'][key] = value
                with self.assertRaises(ValueError):module.PreparedPhysicalDimensions().read(source, [7])

    def test_separate_extent_and_disconnected_source_reject(self):
        source = snapshot(); source.streets[7]['items'][1]['points'] = [(20,0,0),(22,0,0)]
        with self.assertRaisesRegex(ValueError, 'extent'):module.PreparedPhysicalDimensions().read(source,[7])
        source = snapshot(); source.streets[7]['items'][-1]['points'][0] = (21,0,0)
        with self.assertRaisesRegex(ValueError, 'disconnected'):module.PreparedPhysicalDimensions().read(source,[7])

    def test_invalid_numeric_contract(self):
        for key, value in (('minimum_half_length',0),('minimum_half_length',float('nan')),
                           ('maximum_edge_slope',1.1),('maximum_edge_slope',float('inf'))):
            with self.subTest(key=key,value=value):
                source = snapshot(); source.streets[7]['items'][1]['transition'][key] = value
                with self.assertRaisesRegex(ValueError,'settings'):module.PreparedPhysicalDimensions().read(source,[7])

    def test_selection_owners(self):
        for ids in ([7,7],[8]):
            with self.assertRaises(ValueError):module.PreparedPhysicalDimensions().read(snapshot(),ids)


class SpanDimensionsTest(unittest.TestCase):
    def setUp(self):
        from test_physical_ports import fixture, m as ports
        self.source = fixture()
        items = self.source.streets[1]['items']
        for at, width in ((0, 4.), (2, 6.), (4, 8.), (6, 10.)):
            items[at]['width'] = width
        items[3] = dict(kind='SymLane', id=99, points=[], surface_treatment='WIDTH_TRANSITION',
            transition=dict(before_section=11, after_section=12, source_kind='SymLane', source_id=99,
                            minimum_half_length=1.5, maximum_edge_slope=.3))
        self.plan = ports.PreparedPhysicalPorts(('residential',)).prepare(self.source, interior_ids=(21, 22))
        self.reader = module.PreparedPhysicalDimensions()
        self.middle = next(s for s in self.plan.spans if s.section_ids == (11, 12))

    def test_junction_width_changes_are_not_longitudinal_tapers(self):
        before = deepcopy(self.source)
        row, = self.reader.read_spans(self.source, self.plan.spans)
        self.assertEqual((row.street_id, row.before_section, row.after_section, row.source_id), (1, 11, 12, 99))
        self.assertEqual(self.source, before)
        with self.assertRaisesRegex(ValueError, 'undeclared'):
            self.reader.read(self.source, (1,))

    def test_unselected_span_declarations_do_not_enter_selected_span(self):
        outer = next(s for s in self.plan.spans if s.section_ids == (10,))
        self.assertEqual(self.reader.read_spans(self.source, (outer,)), ())
        self.assertEqual(len(self.reader.read_spans(self.source, (self.middle,))), 1)

    def test_deferred_attachment_does_not_hide_missing_taper(self):
        self.source.streets[1]['items'][3] = dict(kind='Intersection', id=23, surface_treatment='DEFERRED_ATTACHMENT')
        with self.assertRaisesRegex(ValueError, 'undeclared'):
            self.reader.read_spans(self.source, (self.middle,))

    def test_forged_cut_cannot_hide_a_width_change(self):
        shortened = replace(self.middle, section_ids=(11,), last_control=3)
        with self.assertRaisesRegex(ValueError, 'original owned Section range'):
            self.reader.read_spans(self.source, (shortened,))

    def test_stale_control_and_endpoint_ownership_reject(self):
        for changed in (replace(self.middle, first_control=1),
                        replace(self.middle, source_control_count=8),
                        replace(self.middle, start=replace(self.middle.start, junction_id=999)),
                        replace(self.middle, end=self.middle.start)):
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                self.reader.read_spans(self.source, (changed,))

    def test_duplicate_spans_and_wrong_declarations_reject(self):
        with self.assertRaisesRegex(ValueError, 'repeats'):
            self.reader.read_spans(self.source, (self.middle, self.middle))
        self.source.streets[1]['items'][3]['transition']['before_section'] = 10
        with self.assertRaisesRegex(ValueError, 'consecutive'):
            self.reader.read_spans(self.source, (self.middle,))


if __name__ == '__main__':unittest.main()
