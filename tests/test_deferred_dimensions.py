"""Deferred crossing ownership and through-carriageway dimensional declarations."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest


ROOT = Path(__file__).resolve().parents[1]/'way/street_network'


def load(name):
    spec = importlib.util.spec_from_file_location('_deferred_'+name, ROOT/(name+'.py'))
    result = importlib.util.module_from_spec(spec);sys.modules[spec.name] = result
    spec.loader.exec_module(result)
    return result


transitions, dimensions, bundles = (load(name) for name in ('transitions', 'physical_dimensions', 'bundles'))


def fixture():
    before = dict(kind='Section', id=10, points=[(-30.,0.,0.),(0.,0.,0.)], width=6., offset=0.)
    after = dict(kind='Section', id=20, points=[(0.,0.,0.),(30.,0.,0.)], width=9., offset=-.4)
    attachment = dict(kind='Intersection', id=7, points=[], surface_treatment='DEFERRED_ATTACHMENT',
        location=(0.,0.), through=((10,1),(20,0)), attachments=[dict(street_id=99, category='footway',
            section_id=50, connector_id=3, end=0, side='LEFT')])
    return [before, attachment, after]


class DeferredDimensionsTest(unittest.TestCase):
    def test_crossing_identity_and_sections_are_preserved(self):
        items = fixture();before = deepcopy(items)
        transitions.PreparedTransitionAdapter().read(items)
        item = items[1]
        self.assertEqual(items[0],before[0]);self.assertEqual(items[2],before[2])
        self.assertEqual({k:item[k] for k in before[1]},before[1])
        self.assertEqual(item['transition']['source_kind'],'Intersection')
        self.assertEqual(item['transition']['source_id'],7)
        source = NS(streets={1:dict(id=1, items=items, status='SUPPORTED')})
        row, = dimensions.PreparedPhysicalDimensions().read(source,(1,))
        self.assertEqual((row.before_section,row.after_section,row.source_id),(10,20,7))
        self.assertEqual(row.source_kind,'Intersection')
        self.assertIsNone(bundles.PreparedBundleReview().member_reason(source.streets[1]))

    def test_native_carriageway_junction_keeps_its_dimensions(self):
        items = fixture();items[1]['surface_treatment'] = 'NATIVE_JUNCTION'
        before = deepcopy(items);transitions.PreparedTransitionAdapter().read(items)
        self.assertEqual(items,before)

    def test_equal_attachment_has_no_new_taper_but_existing_one_can_be_dormant(self):
        items = fixture();items[2].update(width=6.,offset=0.)
        transitions.PreparedTransitionAdapter().read(items)
        self.assertNotIn('transition',items[1])
        items[2]['width']=8.;transitions.PreparedTransitionAdapter().read(items)
        record=deepcopy(items[1]['transition'])
        items[2]['width']=6.;transitions.PreparedTransitionAdapter().read(items)
        self.assertEqual(record,items[1]['transition'])

    def test_bad_ownership_is_not_repaired_by_proximity(self):
        for key,value in (('through',((10,0),(20,0))),('through',((10,1),(20,0),(30,0))),
                          ('location',(1.,0.)),('location',(float('nan'),0.)),('attachments',[]),
                          ('points',[(0.,0.,0.),(1.,0.,0.)])):
            with self.subTest(key=key,value=value):
                items=fixture();items[1][key]=value
                transitions.PreparedTransitionAdapter().read(items)
                self.assertIn('transition_reason',items[1]);self.assertNotIn('transition',items[1])

    def test_physical_reader_rejects_stale_attachment_bindings(self):
        for key,value in (('through',((10,0),(20,0))),('location',(2.,0.)),('attachments',[])):
            with self.subTest(key=key):
                items=fixture();transitions.PreparedTransitionAdapter().read(items)
                items[1][key]=value
                source=NS(streets={1:dict(id=1,status='SUPPORTED',items=items)})
                with self.assertRaisesRegex(ValueError,'deferred attachment'):
                    dimensions.PreparedPhysicalDimensions().read(source,(1,))

    def test_policies_and_repeated_adaptation(self):
        items=fixture();policy=transitions.PreparedTransitionAdapter(minimum_half_length=2.,maximum_edge_slope=.2)
        policy.read(items);before=deepcopy(items);policy.read(items);self.assertEqual(items,before)
        self.assertEqual(items[1]['transition']['minimum_half_length'],2.)
        for args in ({'minimum_half_length':0},{'maximum_edge_slope':2},{'maximum_edge_slope':float('nan')}):
            with self.assertRaises(ValueError):transitions.PreparedTransitionAdapter(**args)

    def test_missing_declaration_still_rejects_in_reader(self):
        source=NS(streets={1:dict(id=1,status='SUPPORTED',items=fixture())})
        with self.assertRaisesRegex(ValueError,'undeclared'):
            dimensions.PreparedPhysicalDimensions().read(source,(1,))

    def test_disconnected_or_nonadjacent_sections_reject(self):
        items=fixture();items[2]['points'][0]=(1.,0.,0.)
        transitions.PreparedTransitionAdapter().read(items);self.assertIn('disconnected',items[1]['transition_reason'])
        items=fixture();items.insert(1,dict(kind='Feature',id=8))
        transitions.PreparedTransitionAdapter().read(items);self.assertIn('ownership',items[2]['transition_reason'])


if __name__=='__main__':unittest.main()
