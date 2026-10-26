"""Original clipping provenance, finite support and persistent endpoint semantics."""
from copy import deepcopy
from dataclasses import replace
import json
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import clipped_fixture_source, fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.paths import Path
from tested_python_streets.source_support import ClippedApproachSupport
from tested_python_streets.routes import RouteTopology, route_key


class PythonSourceSupportTests(unittest.TestCase):
    def test_original_same_way_intervals_extend_either_authored_end_and_keep_all_old_controls(self):
        for reverse in (False, True):
            for mode in (0, 3):
                with self.subTest(reverse=reverse, interpolation=mode):
                    network, snapshot, context = clipped_fixture_source(model, reverse=reverse, interpolation=mode)
                    before = network.to_dict()
                    prepared = ClippedApproachSupport(snapshot, context).prepare(network)
                    self.assertFalse(prepared.error)
                    self.assertEqual(prepared.report['status'], 'CONTACTS_READY')
                    self.assertEqual(network.to_dict(), before)
                    current = prepared.network.section_index[2]
                    old = network.section_index[2]
                    self.assertTrue(all(p in current.points for p in old.points))
                    self.assertEqual(current.id, old.id)
                    self.assertEqual(current.street_id, old.street_id)
                    self.assertEqual((current.width, current.offset, current.rise), (old.width, old.offset, old.rise))
                    self.assertEqual(current.selection_ends, (0 if reverse else 1,))
                    extension, = prepared.report['extensions']
                    self.assertEqual(len(extension['extra_controls']), 2)
                    self.assertLessEqual(extension['added_polyline_length_m'], 40.)
                    self.assertEqual({r[0] for r in extension['source_intervals']}, {'way'})
                    np.testing.assert_allclose(current.points[0 if reverse else -1], (0., 8.4, 0.))
                    # Anchor station follows the same original centerline point
                    # even when prepending changes the directed metric origin.
                    old_point = Path(old).frame(network.furniture[0].station)[0]
                    anchor = prepared.network.furniture[0]
                    new_point = Path(current).frame(anchor.station)[0]
                    np.testing.assert_allclose(new_point, old_point, atol=1e-9)
                    update = PythonMeshBackend(prepared.network.settings).prepare(prepared.network)
                    part = update.products[model.Key('junction_mesh', 1)]
                    self.assertEqual(part.metadata['fitted_corner_radii'], [])
                    self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
                    self.assertLessEqual(part.mesh.xyz[:, 2].max(), .8 + 1e-8)

    def test_selection_cut_is_not_a_genuine_free_end_after_roundtrip(self):
        network, snapshot, context = clipped_fixture_source(model)
        prepared = ClippedApproachSupport(snapshot, context).prepare(network).network
        restored = model.Network.from_dict(json.loads(json.dumps(prepared.to_dict())))
        self.assertEqual(prepared, restored)
        candidate = RouteTopology(restored).trace(model.End(2, 0))
        self.assertTrue(candidate.terminal_selection)
        update = PythonMeshBackend(restored.settings).prepare(restored)
        self.assertFalse(update.products[route_key(model.End(2, 0))].terminal_free)
        with self.assertRaisesRegex(ValueError, 'selection cut'):
            replace(restored, sections=tuple(replace(s, selection_ends=(0,)) if s.id == 2 else s
                                            for s in restored.sections))

    def test_curved_outside_support_keeps_the_original_inner_controls_and_anchor_location(self):
        for reverse in (False, True):
            network, snapshot, context = clipped_fixture_source(model, reverse=reverse, curved=True)
            prepared = ClippedApproachSupport(snapshot, context).prepare(network)
            self.assertFalse(prepared.error)
            current = prepared.network.section_index[2]
            added = prepared.report['extensions'][0]['extra_controls']
            self.assertGreater(added[-1][0], 0.)
            self.assertLess(added[-1][0], 2.)
            before = Path(network.section_index[2]).frame(network.furniture[0].station)[0]
            after = Path(current).frame(prepared.network.furniture[0].station)[0]
            np.testing.assert_allclose(before, after, atol=1e-9)
            update = PythonMeshBackend(prepared.network.settings).prepare(prepared.network)
            self.assertLessEqual(update.products[model.Key('junction_mesh', 1)].mesh.validate()['maximum_edge_ratio'], 8.)

    def test_no_source_support_is_invented_across_a_way_or_tag_change(self):
        for field, value in (('source_way_id', 'another-way'), ('tags', {'highway': 'residential'})):
            network, snapshot, context = clipped_fixture_source(model)
            network = replace(network, intersections=(replace(network.intersections[0],
                radius_modes=('STRICT',) * 3),))
            context['nodes'][0]['legs'][1][field] = value
            prepared = ClippedApproachSupport(snapshot, context).prepare(network)
            self.assertTrue(prepared.error)
            extension, = prepared.report['extensions']
            self.assertEqual(len(extension['source_intervals']), 1)
            self.assertLessEqual(Path(prepared.network.section_index[2]).length, 4. + 1e-8)

    def test_wrong_cut_or_duplicate_ownership_is_diagnostic_and_keeps_original_geometry(self):
        for kind in ('position', 'interval', 'duplicate'):
            network, snapshot, context = clipped_fixture_source(model)
            if kind == 'position':
                context['nodes'][0]['incidences'][0]['clipped_point'] = (1., 2.5)
            elif kind == 'interval':
                context['nodes'][0]['open_ends'][0]['parameter'] = .5
            else:
                context['nodes'].append(deepcopy(context['nodes'][0]))
            prepared = ClippedApproachSupport(snapshot, context).prepare(network)
            self.assertTrue(prepared.error)
            self.assertEqual(network, prepared.network)
            self.assertEqual(prepared.report['extensions'], [])

    def test_extension_budget_is_enforced_and_an_infeasible_footprint_remains_rejected(self):
        network, snapshot, context = clipped_fixture_source(model)
        prepared = ClippedApproachSupport(snapshot, context, maximum_extension=1.).prepare(network)
        self.assertTrue(prepared.error)
        extension, = prepared.report['extensions']
        self.assertLessEqual(extension['added_polyline_length_m'], 1.)
        self.assertLessEqual(Path(prepared.network.section_index[2]).length, 3.5 + 1e-8)

    def test_original_outside_junction_is_not_replaced_by_an_inferred_free_road(self):
        network, snapshot, context = clipped_fixture_source(model)
        context['nodes'][0]['topology_kind'] = 'JUNCTION'
        prepared = ClippedApproachSupport(snapshot, context).prepare(network)
        self.assertIn('complete original incidence', prepared.error)
        self.assertEqual(prepared.report['extensions'], [])

    def test_proven_original_dead_end_can_be_restored_without_a_selection_alias(self):
        network, snapshot, context = clipped_fixture_source(model)
        node = context['nodes'][0]
        node['position'] = (0., 6.)
        node['legs'] = node['legs'][:1]
        node['topology_kind'] = 'ENDPOINT'
        node['incidences'][0]['parameter'] = 2.5 / 6.
        node['open_ends'][0]['parameter'] = 2.5 / 6.
        snapshot.streets[2]['items'][0]['source_segments'] = [('way', 'root', 'outside', 0., 2.5 / 6.)]
        prepared = ClippedApproachSupport(snapshot, context).prepare(network)
        self.assertFalse(prepared.error)
        self.assertEqual(prepared.network.section_index[2].selection_ends, ())
        self.assertEqual(prepared.report['extensions'][0]['output_end_kind'], 'SOURCE_ENDPOINT')

    def test_unchanged_feasible_source_without_context_keeps_the_exact_mesh(self):
        network = fixture_networks(model)[0]
        from types import SimpleNamespace
        snapshot = SimpleNamespace(streets={s.id: dict(items=[dict(kind='Section', id=s.id)]) for s in network.sections})
        prepared = ClippedApproachSupport(snapshot, {}).prepare(network)
        self.assertEqual(prepared.network, network)
        self.assertEqual(prepared.report['extensions'], [])

    def test_explicit_recovered_support_trim_has_a_finite_output_and_persistent_provenance(self):
        network, snapshot, context = clipped_fixture_source(model)
        network = replace(network, sections=tuple(replace(s, points=((0., 0., 0.), (0., 40., 0.)))
                                                 if s.id == 2 else s for s in network.sections))
        snapshot.streets[2]['items'][0]['clip_end'] = 25.
        prepared = ClippedApproachSupport(snapshot, {}).prepare(network)
        self.assertFalse(prepared.error)
        self.assertEqual(prepared.network.section_index[2].points[-1], (0., 15., 0.))
        self.assertEqual(prepared.network.section_index[2].selection_ends, (1,))
        self.assertEqual(prepared.report['context_trims'][0]['original_controls'][-1], (0., 40., 0.))
        update = PythonMeshBackend(prepared.network.settings).prepare(prepared.network)
        road = update.products[model.Key('road', 2)].mesh
        self.assertLessEqual(road.xyz[:, 1].max(), 15. + 1e-8)

    def test_original_node_positions_and_source_direction_must_agree_before_recovery(self):
        for kind in ('position', 'direction'):
            network, snapshot, context = clipped_fixture_source(model)
            if kind == 'position':
                context['nodes'][0]['position'] = (1., 4.)
            else:
                context['nodes'][0]['legs'][0]['source_end'] = 0
            prepared = ClippedApproachSupport(snapshot, context).prepare(network)
            self.assertTrue(prepared.error)
            self.assertEqual(prepared.report['extensions'], [])

    def test_failed_multi_section_context_trim_cannot_publish_partial_provenance(self):
        network, snapshot, context = clipped_fixture_source(model)
        snapshot.streets[1]['items'][0]['clip_end'] = 1.
        snapshot.streets[2]['items'][0]['clip_end'] = 3.
        prepared = ClippedApproachSupport(snapshot, {}).prepare(network)
        self.assertTrue(prepared.error)
        self.assertEqual(prepared.network, network)
        self.assertEqual(prepared.report['context_trims'], [])

    def test_an_explicit_radius_is_not_reduced_when_original_support_can_fit_it(self):
        network, snapshot, context = clipped_fixture_source(model)
        site = replace(network.intersections[0], radius_modes=('STRICT',) * 3)
        prepared = ClippedApproachSupport(snapshot, context).prepare(replace(network, intersections=(site,)))
        self.assertFalse(prepared.error)
        update = PythonMeshBackend(prepared.network.settings).prepare(prepared.network)
        layout = update.products[model.Key('mouths', 1)]
        self.assertTrue(all(c.radius == 2. for c in layout.corners))
        self.assertEqual(update.products[model.Key('junction_mesh', 1)].metadata['fitted_corner_radii'], [])

    def test_style_support_does_not_extend_a_selection_tail_to_chase_its_preferred_maximum(self):
        network, snapshot, context = clipped_fixture_source(model)
        neighbor = replace(network.sections[3], points=((4.5, 18., 0.), (4.5, 40., 0.)), width=4.)
        site = replace(network.intersections[0], radius=10.)
        network = replace(network, sections=(*network.sections[:3], neighbor), intersections=(site,))
        prepared = ClippedApproachSupport(snapshot, context).prepare(network)
        self.assertFalse(prepared.error)
        north = prepared.network.section_index[2]
        self.assertLess(north.points[-1][1], 18.)
        self.assertEqual(prepared.network.intersections[0].radius, 10.)
        self.assertEqual(prepared.report['style_support'], 'FIT_WITHIN_ALLOWED_RADIUS_RANGE')
        update = PythonMeshBackend(network.settings).prepare(prepared.network)
        fits = update.products[model.Key('junction_mesh', 1)].metadata['fitted_corner_radii']
        self.assertTrue(fits)
        self.assertTrue(all(2.5 <= fit['actual_radius'] < 10. for fit in fits))
        # The same explicitly requested 10 m radius retains its reach. It may
        # not silently become a smaller turn to hide a real unrelated overlap.
        strict = replace(network, intersections=(replace(site, radius_modes=('STRICT',) * 3),))
        strict_prepared = ClippedApproachSupport(snapshot, context).prepare(strict)
        self.assertFalse(strict_prepared.error)
        self.assertGreater(strict_prepared.network.section_index[2].points[-1][1], 18.)
        with self.assertRaisesRegex(ValueError, 'Pavement interiors overlap'):
            PythonMeshBackend(network.settings).prepare(strict_prepared.network)

    def test_path_metric_parameter_inverse_handles_both_modes_and_boundaries(self):
        for mode in (0, 3):
            section = model.Section(1, 1, ((0., 0., 0.), (20., 3., 0.), (40., 0., 0.)), interpolation=mode)
            path = Path(section)
            stations = np.linspace(0., path.length, 67)
            np.testing.assert_allclose(path.station(path.parameter(stations)), stations, atol=1e-7)
            parameters = np.r_[path.knots, np.linspace(0., path.knots[-1], 39)]
            np.testing.assert_allclose(path.parameter(path.station(parameters)), parameters, atol=1e-7)
            with self.assertRaises(ValueError):
                path.station(path.knots[-1] + 1.)


if __name__ == '__main__':
    unittest.main()
