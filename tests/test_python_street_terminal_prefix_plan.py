"""A paired original mouth replaces only its finite overlap prefix."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_entrances import network, updater, context_for
from test_python_street_dimensional_bundles import reverse_sections
from tested_python_streets.model import BundleMember, Key
from tested_python_streets.bundle_axis import SeparatedPlan
from tested_python_streets.bundle_merged_entrances import PairedAngularEntrance, preparation_key
from tested_python_streets.bundle_terminal_unions import TerminalMergedBundlePlan, _terminal_support
from tested_python_streets.providers import CarriagewayProvider, PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets import model
from tested_python_streets.routes import RouteTopology
from tested_python_streets.separator_ports import PORTS, SeparatorPortIndex
from tested_python_streets.updates import BuildContext
from tested_python_streets.uv import PlanarMapping, apply_mapping


def prepared_run(mode=0, *, reverse=False, angle=0., exaggeration=1.):
    source = network(mode)
    if reverse: source = reverse_sections(source, {s.id for s in source.sections})
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    source = replace(source, settings=replace(source.settings, exaggeration=exaggeration),
        sections=tuple(replace(s, points=tuple((*point, z) for point, z in
            zip(np.asarray(s.points)[:, :2]@rotation.T+(1200., -1800.) if angle else np.asarray(s.points)[:, :2],
                np.asarray(s.points)[:, 2]))) for s in source.sections))
    engine = updater(source)
    update = engine.prepare(source)
    context = context_for(source, update, engine.patterns)
    specification = context[preparation_key(source.bundles[0])]
    entrance = next(a for a in specification.layout.approaches
                    if isinstance(a, PairedAngularEntrance) and a.end_index == 0)
    # A consumer fixture selects one already prepared ORIGINAL source run.
    # All actual root/remote junctions and crossing-owned sources remain.
    bundle = replace(source.bundles[0], members=tuple(BundleMember(m.visits) for m in entrance.run.members),
                     interior_section_ids=())
    source = replace(source, bundles=(bundle,))
    context = BuildContext(source, dict(context.products), context.patterns,
                           context.task_digests, context.source_digests)
    context.products[PORTS] = SeparatorPortIndex.build(context)
    return context, bundle


class TerminalPrefixPlanTests(unittest.TestCase):
    def test_support_uses_only_exact_original_two_ended_c1_seams_and_stops_before_real_owners(self):
        sections = (model.Section(1, 1, ((0., 5., 0.), (20., 5., 0.)), width=10.),
                    model.Section(2, 1, ((20., 5., 0.), (40., 5., 0.)), width=10.),
                    model.Section(3, 1, ((40., 5., 0.), (60., 5., 0.)), width=10.),
                    model.Section(4, 4, ((40., 5., 0.), (40., 30., 0.)), width=3.),
                    model.Section(5, 5, ((0., -8., 0.), (60., -8., 0.)), width=10.))
        sections = tuple(replace(s, interpolation=0) for s in sections)
        bundle = model.Bundle(10, (model.BundleMember(tuple(model.Traversal(i) for i in (1, 2, 3))),
                                   model.BundleMember((model.Traversal(5),))))
        network = model.Network('original-support', sections,
            (model.Intersection(20, (model.End(2, 1), model.End(3, 0), model.End(4, 0))),),
            connections=(model.Connection(10, (model.End(1, 1), model.End(2, 0))),), bundles=(bundle,))
        def context(source):
            patterns = PatternLibrary()
            update = UpdateCoordinator((PathProvider(), ProfileProvider()), patterns).prepare(source)
            return BuildContext(source, update.products, patterns)
        self.assertEqual(tuple(v.section_id for v in _terminal_support(context(network), bundle, 0, 0)), (1, 2))
        self.assertEqual(tuple(v.section_id for v in _terminal_support(context(network), bundle, 0, 1)), (3,))
        dimensional = replace(network, sections=tuple(replace(s, width=7.) if s.id == 2 else s for s in sections))
        self.assertEqual(tuple(v.section_id for v in _terminal_support(context(dimensional), bundle, 0, 0)), (1,))
        bent = replace(network, sections=tuple(replace(s, points=((20., 5., 0.), (40., 9., 0.))) if s.id == 2 else s for s in sections))
        self.assertEqual(tuple(v.section_id for v in _terminal_support(context(bent), bundle, 0, 0)), (1,))

    def test_reversed_sources_rigid_transform_and_flat_display_preserve_physical_contract(self):
        context, bundle = prepared_run(0, reverse=True, angle=.713, exaggeration=0.)
        before = context.network.to_dict()
        plan = TerminalMergedBundlePlan.build(context, bundle)
        metadata = plan.pavement(context).metadata
        for check in metadata['stored_checks']:
            self.assertLessEqual(check['stored_quality']['maximum_edge_ratio'], 8.)
            self.assertEqual(check['drainage']['maximum_fill'], 0.)
            if check['material'] == 'pavement':
                self.assertGreaterEqual(check['stored_quality']['minimum_corner_sine'], .03)
                self.assertLessEqual(check['longitudinal_bound'], .03+1e-8)
        self.assertLessEqual(metadata['source_containment']['maximum_combined_error'], .002+1e-7)
        for part in (plan.pavement(context), plan.separator(context)):
            np.testing.assert_array_equal(part.mesh.xyz[:, 2], np.full(len(part.mesh.xyz), plan.datum))
        self.assertEqual(context.network.to_dict(), before)

    def test_original_complete_run_keeps_profiles_uv_query_and_cached_source_products(self):
        context, bundle = prepared_run()
        before = context.network.to_dict()
        identities = {key: id(value) for key, value in context.products.items()}
        plan = TerminalMergedBundlePlan.build(context, bundle)
        self.assertEqual(plan.mode, 'union')
        self.assertEqual(set(plan.source_ids), set(bundle.sections))
        for part in (plan.pavement(context), plan.separator(context)):
            self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
            mapped = apply_mapping(part.mesh, PlanarMapping())
            self.assertTrue(np.isfinite(mapped.uv_layers['UVMap']).all())
        metadata = plan.pavement(context).metadata
        self.assertTrue(metadata['computational_interfaces_only'])
        self.assertFalse(metadata['physical_aliases_created'])
        self.assertFalse(metadata['junction_geometry_republished'])
        samples = np.concatenate([row[0] for row in plan.evaluation_profiles.values()])
        np.testing.assert_allclose(plan.height(samples[:, :2]), samples[:, 2], atol=1e-8)
        query = plan.height(samples[:, :2])
        plan.release_geometry()
        self.assertIsNone(plan.pavement(context).mesh)
        np.testing.assert_array_equal(plan.height(samples[:, :2]), query)
        self.assertEqual(context.network.to_dict(), before)
        self.assertEqual({key: id(value) for key, value in context.products.items()}, identities)

    def test_computational_road_subset_preserves_original_absolute_station_uv_and_full_profile(self):
        context, bundle = prepared_run()
        sid = bundle.members[0].path[0].section_id
        products = dict(context.products)
        products[Key('bundle_plan', bundle.id)] = SeparatedPlan()
        context = BuildContext(context.network, products, context.patterns)
        section = context.network.section_index[sid]
        attached = RouteTopology(context.network).attachments(context.network)[sid]
        original = CarriagewayProvider.build(context, section, attached)
        first, last = original.metadata['start'], original.metadata['end']
        interval = (first+.25*(last-first), last-.25*(last-first))
        clipped = CarriagewayProvider.build(context, section, attached, source_interval=interval)
        self.assertFalse(clipped.metadata['computational_source_interval']['physical_aliases_created'])
        uv = clipped.mesh.uv_layers['UVMap']
        self.assertAlmostEqual(float(uv[..., 0].min()), interval[0])
        self.assertAlmostEqual(float(uv[..., 0].max()), interval[1])
        path = context[Key('path', sid)]
        profile = context[Key('profile', sid)]
        q = np.linspace(-1., 1., clipped.metadata['transverse_intervals']+1)
        for end, station in enumerate(interval):
            expected = np.c_[path.point(station, section.offset+.5*section.width*q),
                             profile.datum+profile.value(q)]
            np.testing.assert_array_equal(clipped.mesh.xyz[clipped.mesh.seams[end]], expected)
        for invalid in ((first-1., last), (first, last+1.), (first, first), (np.nan, last)):
            with self.assertRaisesRegex(ValueError, 'positive subset'):
                CarriagewayProvider.build(context, section, attached, source_interval=invalid)
        reference = CarriagewayProvider.build(context, section, attached)
        np.testing.assert_array_equal(reference.mesh.xyz, original.mesh.xyz)
        np.testing.assert_array_equal(reference.mesh.quads, original.mesh.quads)
        self.assertEqual(reference.metadata, original.metadata)


if __name__ == '__main__': unittest.main()
