"""Runtime union admission uses complete visible original source fibres."""
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.bundle_angular_curves import AngularCoordinates, AngularCrown
from tested_python_streets.gapped_runtime_admission import visible_overlap, _ranges
from tested_python_streets.gapped_source_plan import GappedSourcePlan
from tested_python_streets.gapped_source_inventory import GappedSourceUnavailable
from tested_python_streets.model import (Section, Network, Settings, Key, Traversal,
    Bundle, BundleMember, Intersection, End)
from tested_python_streets.paths import Path
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import Crown
from tested_python_streets.updates import BuildContext


def fixture(first=None, second=None, *, reverse=False, angle=0., shift=(0., 0.)):
    first = ((0., 3.), (10., 3.)) if first is None else first
    second = ((0., -3.), (10., -3.)) if second is None else second
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    frame = AngularCoordinates.__new__(AngularCoordinates)
    frame.center = np.array((0., -1000.)) @ rotation.T + shift
    frame.scale, frame.angle, frame.direction = 1000., np.pi*.5+angle, -1.
    products, sections, spans, curves = {}, [], [], []
    for role, controls in enumerate((first, second)):
        xy = np.asarray(controls) @ rotation.T + shift
        if reverse:xy = xy[::-1]
        section = Section(role+1, role+1, tuple((x, y, 0.) for x, y in xy), width=4., interpolation=0)
        path = Path(section)
        peak = (1 if role == 0 else -1)*(-1 if reverse else 1)
        profile = Crown(.08, 0., peak)
        products[Key('path', section.id)] = path
        products[Key('profile', section.id)] = profile
        visit = Traversal(section.id, not reverse)
        spans.append((SimpleNamespace(role=role, visits=(visit,), source_intervals=(
            dict(section=section.id, forward=visit.forward, first=0., last=path.length),)),))
        samples = path.samples(0., path.length, .5)
        own = [path.point(samples, .5*section.width*side) for side in (1., -1.)]
        if reverse:own = [points[::-1] for points in own[::-1]]
        # Authored rightward travel has Aouter/Ainner then Binner/Bouter.
        curves.extend((own[1], own[0]))
        sections.append(section)
    network = Network('gapped-admission', tuple(sections), settings=Settings(curve_tolerance=.02))
    context = BuildContext(network, products, PatternLibrary())
    inventory = SimpleNamespace(frame=frame, field=AngularCrown(frame, curves, datum=0., rise=.08),
                                spans_by_role=tuple(spans), owned_spans=())
    return context, inventory


def restrict(inventory, role, intervals):
    original = inventory.spans_by_role[role][0]
    row = original.source_intervals[0]
    spans = tuple(SimpleNamespace(role=role, visits=original.visits,
        source_intervals=(dict(row, first=first, last=last),)) for first, last in intervals)
    rows = list(inventory.spans_by_role);rows[role] = spans
    return SimpleNamespace(**dict(inventory.__dict__, spans_by_role=tuple(rows)))


class GappedRuntimeAdmissionTests(unittest.TestCase):
    def test_pure_separated_sources_and_approximation_scale_contacts_stay_unselected(self):
        for half_separation in (3., 2., 1.985):
            context, inventory = fixture(((0., half_separation), (10., half_separation)),
                                         ((0., -half_separation), (10., -half_separation)))
            self.assertIsNone(visible_overlap(context, inventory))

    def test_clear_overlap_is_reinverted_on_original_curbs_under_reversal_and_rigid_transforms(self):
        for reverse, angle, shift in ((False, 0., (0., 0.)), (True, .713, (1400., -1900.))):
            context, inventory = fixture(((0., 1.5), (10., 1.5)), ((0., -1.5), (10., -1.5)),
                reverse=reverse, angle=angle, shift=shift)
            source = context.network.to_dict()
            report = visible_overlap(context, inventory)
            self.assertIsNotNone(report)
            self.assertGreater(report['overlap_depth'], report['combined_source_allowance'])
            self.assertEqual(report['combined_source_allowance'], .04)
            self.assertTrue(report['analytical_curbs_inverted'] and report['complete_transverse_fibres'])
            for row in report['original_analytic_curbs']:
                path = context[Key('path', row['section'])]
                first, last = row['visible_interval']
                for side in ('inner', 'outer'):
                    actual = row[side]
                    self.assertLess(first, actual['source_station'])
                    self.assertLess(actual['source_station'], last)
                    np.testing.assert_array_equal(path.point(actual['source_station'], actual['offset']), actual['point'])
                    self.assertLessEqual(actual['inversion_residual'], report['numerical_margin'])
            self.assertEqual(context.network.to_dict(), source)

    def test_a_sampled_field_overlap_cannot_replace_original_analytic_separation(self):
        context, inventory = fixture()
        _, misleading = fixture(((0., 1.5), (10., 1.5)), ((0., -1.5), (10., -1.5)))
        inventory.field = misleading.field
        self.assertIsNone(visible_overlap(context, inventory))

    def test_owned_prefix_overlap_and_owned_gaps_never_supply_a_visible_witness(self):
        first = ((0., 1.5), (2., 1.5), (4., 3.), (10., 3.))
        second = tuple((x, -y) for x, y in first)
        context, inventory = fixture(first, second)
        self.assertIsNotNone(visible_overlap(context, inventory))
        for role in (0, 1):
            path = context[Key('path', role+1)]
            inventory = restrict(inventory, role, ((float(path.knots[2]), path.length),))
        self.assertIsNone(visible_overlap(context, inventory))
        first = ((0., 3.), (5., 3.), (10., 1.5), (20., 1.5), (25., 3.), (30., 3.))
        context, inventory = fixture(first, tuple((x, -y) for x, y in first))
        self.assertIsNotNone(visible_overlap(context, inventory))
        path = context[Key('path', 1)]
        inventory = restrict(inventory, 0, ((0., float(path.knots[1])),
                                           (float(path.knots[4]), path.length)))
        self.assertIsNone(visible_overlap(context, inventory))

    def test_whole_transverse_profile_planes_bound_the_eligible_fibre_range(self):
        context, inventory = fixture(((0., 1.5), (10., 1.5)), ((0., -1.5), (10., -1.5)))
        for role in (0, 1):inventory = restrict(inventory, role, ((2., 8.),))
        margin = 1e-7
        ranges = _ranges(context, inventory, margin)
        for role in (0, 1):
            source = ranges[role][0];path = source.path
            head = inventory.frame.local(path.point(2., np.array((source.inner, source.outer))))[:, 0]
            tail = inventory.frame.local(path.point(8., np.array((source.inner, source.outer))))[:, 0]
            self.assertAlmostEqual(source.lower, float(head.max())+margin)
            self.assertAlmostEqual(source.upper, float(tail.min())-margin)
            self.assertGreater(source.lower, float(head.min()))
            self.assertLess(source.upper, float(tail.max()))

    def test_runtime_requires_the_witness_before_meshing_and_explicit_api_stays_opt_in(self):
        context, inventory = fixture()
        source = context.network.section_index[1]
        continuation = replace(source, id=3, points=((-10., 3., 0.), (0., 3., 0.)))
        stem = replace(source, id=4, street_id=4, points=((0., 3., 0.), (0., 12., 0.)))
        bundle = Bundle(1, (BundleMember((Traversal(3), Traversal(1))), BundleMember((Traversal(2),))))
        site = Intersection(1, (End(3, 1), End(1, 0), End(4, 0)))
        context.network = replace(context.network, sections=(*context.network.sections, continuation, stem),
            intersections=(site,), bundles=(bundle,))
        with patch('tested_python_streets.gapped_source_plan.physical_context', return_value=context), \
             patch('tested_python_streets.gapped_source_plan.GappedSourceInventory.build', return_value=inventory), \
             patch('tested_python_streets.gapped_source_plan._coverage', side_effect=AssertionError('meshing reached')):
            with self.assertRaisesRegex(GappedSourceUnavailable, 'certified visible member overlap'):
                GappedSourcePlan.build(context, bundle, require_visible_overlap=True)
            with self.assertRaisesRegex(AssertionError, 'meshing reached'):
                GappedSourcePlan.build(context, bundle)

    def test_normal_previously_successful_separated_fixture_keeps_all_output_bytes(self):
        from test_python_street_merged_entrances import network
        from test_python_street_merged_run_inventory import side_junction
        from tested_python_streets.backend import PythonMeshBackend
        with patch('test_python_street_merged_run_inventory.network',
                   side_effect=lambda mode:network(mode, state='separated')):
            source = side_junction()
        with patch('tested_python_streets.gapped_source_plan.GappedSourcePlan.build',
                   side_effect=GappedSourceUnavailable('established separated provider comparison')):
            old = PythonMeshBackend(source.settings).prepare(source)
        current = PythonMeshBackend(source.settings).prepare(source)
        key = Key('bundle_plan', source.bundles[0].id)
        self.assertEqual(old.products[key].mode, 'separated')
        self.assertEqual(current.products[key].mode, 'separated')
        outputs = lambda update:{key:product.mesh for key, product in update.products.items()
            if getattr(product, 'mesh', None) is not None}
        expected, actual = outputs(old), outputs(current)
        self.assertEqual(expected.keys(), actual.keys())
        for key, mesh in expected.items():
            np.testing.assert_array_equal(actual[key].xyz, mesh.xyz)
            np.testing.assert_array_equal(actual[key].quads, mesh.quads)
            self.assertEqual(actual[key].material, mesh.material)
            self.assertEqual(actual[key].uv_layers.keys(), mesh.uv_layers.keys())
            for name, values in mesh.uv_layers.items():
                np.testing.assert_array_equal(actual[key].uv_layers[name], values)

    def test_only_typed_pre_admission_unavailability_declines_runtime_discovery(self):
        from tested_python_streets.bundle_axis import AxialChartUnavailable
        from tested_python_streets.ownership import CurbChainUnavailable
        context, inventory = fixture()
        source = context.network.section_index[1]
        continuation = replace(source, id=3, points=((-10., 3., 0.), (0., 3., 0.)))
        stem = replace(source, id=4, street_id=4, points=((0., 3., 0.), (0., 12., 0.)))
        bundle = Bundle(1, (BundleMember((Traversal(3), Traversal(1))), BundleMember((Traversal(2),))))
        site = Intersection(1, (End(3, 1), End(1, 0), End(4, 0)))
        context.network = replace(context.network, sections=(*context.network.sections, continuation, stem),
            intersections=(site,), bundles=(bundle,))
        for failure in (AxialChartUnavailable('finite frame unavailable'),
                        CurbChainUnavailable('another actual mouth interrupts the curb')):
            with patch('tested_python_streets.gapped_source_plan.physical_context', side_effect=failure):
                with self.assertRaisesRegex(GappedSourceUnavailable, 'runtime source chart is unavailable'):
                    GappedSourcePlan.build(context, bundle, require_visible_overlap=True)
                with self.assertRaises(type(failure)) as raised:GappedSourcePlan.build(context, bundle)
                self.assertIs(raised.exception, failure)
        failure = ValueError('invalid original source certificate')
        with patch('tested_python_streets.gapped_source_plan.physical_context', side_effect=failure):
            with self.assertRaises(ValueError) as raised:
                GappedSourcePlan.build(context, bundle, require_visible_overlap=True)
            self.assertIs(raised.exception, failure)
        with patch('tested_python_streets.gapped_source_plan.physical_context', return_value=context), \
             patch('tested_python_streets.gapped_source_plan.GappedSourceInventory.build', return_value=inventory), \
             patch('tested_python_streets.gapped_source_plan.visible_overlap', return_value={'proved':True}), \
             patch('tested_python_streets.gapped_source_plan._coverage', side_effect=AxialChartUnavailable('admitted coverage error')):
            with self.assertRaisesRegex(AxialChartUnavailable, 'admitted coverage error') as raised:
                GappedSourcePlan.build(context, bundle, require_visible_overlap=True)
            self.assertNotIsInstance(raised.exception, GappedSourceUnavailable)


if __name__ == '__main__':unittest.main()
