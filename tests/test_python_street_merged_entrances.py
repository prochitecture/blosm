"""Normal early preparation, source aliases and physical junction ownership."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import internal_angular_state_network
from tested_python_streets import model
from tested_python_streets.bundle_merged_entrances import PairedAngularEntrance, preparation_key
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.ownership import PLAN, surface_key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import PathProvider, CarriagewayProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.routes import RouteTopology
from tested_python_streets.separator_ports import SeparatorPortIndex
from tested_python_streets.updates import BuildContext, UpdateCoordinator


def network(mode, state='overlap'):
    source = internal_angular_state_network(model, mode, state=state)
    return replace(source, sections=tuple(replace(section,
        points=tuple(map(tuple, np.asarray(section.points, dtype=np.float32).astype(float))),
        rise=float(np.float32(section.rise))) for section in source.sections),
        settings=replace(source.settings, exaggeration=1.))


def updater(source):
    return UpdateCoordinator((PathProvider(), ProfileProvider(), IntersectionProvider()),
                             PatternLibrary(source.settings.pattern_cache_bytes))


def context_for(source, update, patterns):
    return BuildContext(source, update.products, patterns, update.task_fingerprints,
                        update.source_fingerprints)


class MergedEntranceTests(unittest.TestCase):
    def qualify(self, source, update, patterns):
        ctx = context_for(source, update, patterns)
        spec = update.products[preparation_key(source.bundles[0])]
        self.assertIsNotNone(spec)
        plan = update.products[PLAN]
        self.assertEqual(plan.owner(model.Key('intersection', 500)), model.Key('intersection', 500))
        self.assertEqual(plan.owner(model.Key('intersection', 501)), model.Key('intersection', 500))
        product = update.products[surface_key(model.Key('intersection', 500))]
        self.assertIsNone(update.products[surface_key(model.Key('intersection', 501))])
        self.assertEqual(product.pavement.metadata['consumed_sections'], [301])
        self.assertLessEqual(product.pavement.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertEqual(product.group.members, tuple(sorted(spec.crossing.sites)))
        self.assertEqual(set(product.peaks), {(17, 1), (117, 0), (29, 1), (129, 0), (300, 1), (302, 0)})
        self.assertEqual(product.peaks[17, 1][:2], product.peaks[29, 1][:2])
        self.assertEqual(product.peaks[117, 0][:2], product.peaks[129, 0][:2])
        data = BoundaryConditions().build(product.group, product.group.loops[0],
                                         {p.region.owner: p.layout for p in product.group.footprints})
        paired_ports = [port for port in data.ports if len(port.aliases) == 2]
        self.assertEqual(len(paired_ports), 2)
        for port in paired_ports:
            with self.assertRaisesRegex(ValueError, 'separate original source parameter'):
                _ = port.source_u
            approach = product.group.physical_layout.approaches[port.approach_index]
            self.assertIsInstance(approach, PairedAngularEntrance)
            self.assertFalse(hasattr(approach.path, 'trims'))
            for alias in port.aliases:
                path = update.products[model.Key('path', alias.endpoint.section_id)]
                original_xy = path.frame(alias.station)[0]
                expected_station = approach.chart_station(product.group.physical_layout.outward_stations[port.approach_index])
                self.assertAlmostEqual(float(spec.frame.local(original_xy)[0]), float(expected_station), places=8)
                profile = update.products[model.Key('profile', alias.endpoint.section_id)].oriented(alias.endpoint.end)
                selected = np.asarray(port.outward_u)[alias.contains(port.outward_u)]
                original_q = 2.*alias.outward_parameter(selected)-1.
                physical_q = 2.*selected-1.
                np.testing.assert_allclose(profile.value(original_q), approach.profile.value(physical_q), atol=1e-12, rtol=0.)
                for outer in (alias.physical_outer_u, alias.physical_peak_u):
                    original_u = float(alias.outward_parameter(outer))
                    self.assertIn((alias.endpoint.section_id, alias.endpoint.end, round(original_u, 12)), product.port_points)
        ports = SeparatorPortIndex.build(ctx)
        for sid, end in ((17, 1), (29, 1), (117, 0), (129, 0)):
            visible = ports.resolve(model.End(sid, end))
            expected_site = 500 if sid in (17, 117) else 501
            self.assertEqual(visible.site.id, expected_site)
        self.assertIsNone(ports.resolve(model.End(301, 0)))
        self.assertIsNone(ports.resolve(model.End(301, 1)))
        consumed = CarriagewayProvider.span_metadata(ctx, source.section_index[301],
                     RouteTopology(source).attachments(source)[301])
        self.assertEqual(consumed['start'], update.products[model.Key('path', 301)].length)
        self.assertEqual(consumed['end'], 0.)
        self.assertEqual({span['owner'] for span in consumed['compound_spans']}, {'intersection:500'})
        return spec, product

    def test_both_modes_keep_two_original_endpoints_per_physical_mouth(self):
        for mode in (0, 3):
            source = network(mode)
            saved = source.to_dict()
            engine = updater(source)
            result = engine.prepare(source)
            with self.subTest(mode=mode):
                self.qualify(source, result, engine.patterns)
                self.assertEqual(source.to_dict(), saved)
                key = preparation_key(source.bundles[0])
                self.assertTrue(all(dep.kind not in ('junction_surface', 'bundle_plan', 'road', 'separator')
                                    for dep in result.dependencies[key]))

    def test_source_edit_rebuilds_preparation_and_both_logical_layouts(self):
        source = network(3)
        engine = updater(source)
        initial = engine.prepare(source)
        engine.commit(initial)
        edited = replace(source, sections=tuple(replace(section, width=section.width+.1)
                         if section.id == 17 else section for section in source.sections))
        local = engine.prepare(edited)
        full = updater(edited).prepare(edited)
        self.qualify(edited, local, engine.patterns)
        for key in (preparation_key(source.bundles[0]), model.Key('mouths', 500), model.Key('mouths', 501)):
            self.assertIn(key, local.rebuilt)
        product_key = surface_key(model.Key('intersection', 500))
        np.testing.assert_array_equal(local.products[product_key].pavement.mesh.xyz,
                                      full.products[product_key].pavement.mesh.xyz)
        np.testing.assert_array_equal(local.products[product_key].pavement.mesh.quads,
                                      full.products[product_key].pavement.mesh.quads)
        engine.commit(local)
        restored = engine.prepare(source)
        np.testing.assert_array_equal(restored.products[product_key].pavement.mesh.xyz,
                                      initial.products[product_key].pavement.mesh.xyz)
        self.assertIs(local.products[surface_key(model.Key('intersection', 210))],
                      initial.products[surface_key(model.Key('intersection', 210))])

    def test_separated_internal_branch_keeps_its_original_layout(self):
        source = network(3, 'multiple')
        engine = updater(source)
        result = engine.prepare(source)
        self.assertIsNone(result.products[preparation_key(source.bundles[0])])
        self.assertIsNot(result.products[model.Key('mouths', 500)], result.products[model.Key('mouths', 501)])

    def test_input_record_order_cannot_change_port_owners_or_signatures(self):
        source = network(3)
        engine = updater(source)
        first = engine.prepare(source)
        original = SeparatorPortIndex.build(context_for(source, first, engine.patterns))
        reordered = replace(source, intersections=source.intersections[::-1])
        engine.commit(first)
        result = engine.prepare(reordered)
        current = SeparatorPortIndex.build(context_for(reordered, result, engine.patterns))
        self.assertEqual(original.signature(source.bundles[0].sections), current.signature(source.bundles[0].sections))
        for endpoint in original.ports:
            self.assertEqual(original.resolve(endpoint).site.key, current.resolve(endpoint).site.key)
        self.assertEqual(result.rebuilt, ())


if __name__ == '__main__':
    unittest.main()
