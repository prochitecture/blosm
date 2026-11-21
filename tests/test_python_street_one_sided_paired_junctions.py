"""One-sided physical preparation preserves original source incidence."""
from dataclasses import fields, replace
from copy import copy
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.intersections import Approach, IntersectionProvider, JunctionLayout
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.one_sided_paired_junctions import OneSidedPairedSpecification, certify
from tested_python_streets.ownership import CompoundProduct
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.separator_ports import SeparatorPortIndex


def fixture(mode=0, reverse=False, angle=0., translation=(0., 0.), separation=8., middle_dip=.1, near_bend=.1):
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    shapes = ((1, ((-30., 4.), (-15., 4.+near_bend), (0., 4.)), 10.),
              (2, ((0., 4.), (15., 4.+near_bend), (30., 4.)), 10.),
              (3, ((-30., 4.-separation), (0., 4.-separation-middle_dip), (30., 4.-separation)), 10.),
              (4, ((0., 4.), (0., 30.)), 2.))
    sections = []
    for sid, points, width in shapes:
        moved = np.asarray(points)@rotation.T+translation
        xyz = tuple((float(x), float(y), 0.) for x, y in moved)
        sections.append(model.Section(sid, sid, xyz[::-1] if reverse else xyz,
                                      width=width, interpolation=mode))
    ends = (model.End(1, int(not reverse)), model.End(2, int(reverse)), model.End(4, int(reverse)))
    site = model.Intersection(100, ends, radius=2.)
    before, after, opposite = (model.Traversal(sid, not reverse) for sid in (1, 2, 3))
    bundle = model.Bundle(10, (model.BundleMember((before, after)), model.BundleMember((opposite,))))
    network = model.Network('one-sided-paired', tuple(sections), (site,), bundles=(bundle,))
    context = footprint_context(network)
    layout = IntersectionProvider.footprint(context, site)
    return context, site, layout, (before, after), opposite


class OneSidedPairedJunctionTests(unittest.TestCase):
    def test_retained_outer_curbs_join_the_actual_merged_and_ordinary_profiles(self):
        from tested_python_streets.one_sided_curb_heights import exterior_profile
        from tested_python_streets.terminal_curb_ancestry import OriginalTerminalBank
        from tested_python_streets.source_curb_conditions import SourceCurbConditionUnavailable
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.surface_fields import BoundaryJetField

        context, site, original, visits, opposite = fixture(mode=3, reverse=True, angle=.4)
        spec = OneSidedPairedSpecification.build(context, site, original, visits, opposite, 0)
        layout = spec.group.physical_layout
        through_curbs = []
        for i, curve in enumerate(layout.boundaries):
            self.assertEqual(len(curve.one_sided_edges), len(curve.points)-1)
            z, first, last, phase = exterior_profile(context, layout, i)
            self.assertEqual(z, 0.)
            for j, side, slope in ((i, -1., first), ((i+1)%3, 1., last)):
                actual = layout.approaches[j].physical_gradient(layout.outward_stations[j], side)
                self.assertAlmostEqual(slope, np.linalg.norm(actual))
            if all(isinstance(edge, OriginalTerminalBank) for edge in curve.one_sided_edges):
                self.assertEqual(phase, (0., 1.))
                through_curbs.append(i)
            else:
                self.assertLess(phase[1]-phase[0], 1.)
        self.assertEqual(len(through_curbs), 1)
        opposite_curb = through_curbs[0]
        sampled = AnalyticalBoundarySampler().build(spec.group)
        builder = CompoundSurfaceBuilder(context.patterns, source_heights=True)
        conditions, _, _ = builder.boundary_conditions(sampled, context)
        BoundaryJetField(conditions)
        bad = copy(layout)
        curve = copy(layout.boundaries[opposite_curb])
        first, *rest = curve.one_sided_edges
        curve.one_sided_edges = (replace(first, side=-first.side), *rest)
        bad.boundaries = tuple(curve if i == opposite_curb else c for i, c in enumerate(layout.boundaries))
        with self.assertRaisesRegex(SourceCurbConditionUnavailable, 'inner carriageway bank'):
            exterior_profile(context, bad, opposite_curb)

    def test_separator_keeps_the_complete_inner_plateau_and_rejects_an_outer_bank(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.one_sided_curb_heights import separator_level
        from tested_python_streets.provenance_rim_sampling import resample, strict_sampling
        from tested_python_streets.surface_fields import BoundaryJetField
        from tested_python_streets.terminal_curb_ancestry import OriginalTerminalBank
        from tested_python_streets.source_curb_conditions import SourceCurbConditionUnavailable

        context, site, original, visits, opposite = fixture(middle_dip=2.2, near_bend=0.)
        spec = OneSidedPairedSpecification.build(context, site, original, visits, opposite, 0)
        layout = spec.group.physical_layout
        pair, = layout.one_sided_separator_pairs
        height = separator_level(context, layout, pair)
        self.assertAlmostEqual(height, .08)
        builder = CompoundSurfaceBuilder(context.patterns, source_heights=True)
        strict = strict_sampling(spec.group)
        for group in (spec.group, strict, resample(strict, context)):
            conditions, _, _ = builder.boundary_conditions(group, context)
            hole, = (boundary for boundary in conditions if boundary.hole)
            np.testing.assert_allclose(hole.heights, height, rtol=0., atol=1e-12)
            np.testing.assert_array_equal(hole.gradients, np.zeros_like(hole.xy))
            BoundaryJetField(conditions)

        index, i = next((index, i) for index in pair
                        for i, edge in enumerate(layout.boundaries[index].one_sided_edges)
                        if isinstance(edge, OriginalTerminalBank))
        curb = copy(layout.boundaries[index])
        edges = list(curb.one_sided_edges)
        edges[i] = replace(edges[i], side=-edges[i].side)
        curb.one_sided_edges = tuple(edges)
        altered = copy(layout)
        altered.boundaries = tuple(curb if i == index else c for i, c in enumerate(layout.boundaries))
        with self.assertRaisesRegex(SourceCurbConditionUnavailable, 'original inner bank'):
            separator_level(context, altered, pair)

    def test_physical_layout_implements_the_normal_junction_and_road_signature_protocol(self):
        from tested_python_streets.ownership import PLAN, OwnershipPlan
        from tested_python_streets.providers import CarriagewayProvider
        from tested_python_streets.routes import RouteTopology
        context, site, original, visits, opposite = fixture()
        spec = OneSidedPairedSpecification.build(context, site, original, visits, opposite, 0)
        layout = spec.group.physical_layout
        self.assertTrue(all(hasattr(layout, field.name) for field in fields(JunctionLayout)))
        self.assertEqual(layout.mesh_strategy, 'boundary_chart')
        self.assertIsNone(layout.transition)
        self.assertEqual(layout.absorbed_transitions, ())
        self.assertIsNone(layout.absorbed_by)
        self.assertEqual(layout.finite_scope_ports, ())
        self.assertEqual(layout.finite_free_pads, ())
        np.testing.assert_array_equal(layout.center, original.center)
        self.assertIsNot(layout.center, original.center)
        context.products[site.layout_key] = layout
        context.products[PLAN] = OwnershipPlan((spec.group,), {site.key: spec.group})
        attachments = RouteTopology(context.network).attachments(context.network)[visits[0].section_id]
        _, mouths = CarriagewayProvider.signature(context, (), attachments, visits[0].section_id)
        self.assertTrue(mouths)
        self.assertTrue(all(row[-2] is False for row in mouths))

    def test_complete_original_regions_profiles_and_nonincidental_aliases_survive(self):
        for mode, reverse, angle, offset in ((0, False, 0., (0., 0.)),
                                           (3, True, .713, (1200., -1800.))):
            with self.subTest(mode=mode, reverse=reverse):
                context, site, layout, visits, opposite = fixture(mode, reverse, angle, offset)
                before = context.network.to_dict()
                spec = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
                sampled = AnalyticalBoundarySampler().build(spec.group)
                self.assertEqual(certify(sampled)['clipped_original_full_width_strips'], 3)
                self.assertEqual(spec.group.members, (site.key,))
                self.assertFalse(OneSidedPairedSpecification.candidates(context.network))
                self.assertTrue(spec.proof['original_near_mouths_internal'])
                self.assertTrue(spec.proof['original_root_caps_internal'])
                for entrance in spec.entrances:
                    sid = entrance.source.id
                    for station in entrance.protected_source_parameters[sid]:
                        self.assertTrue(np.any(entrance.source_curve_stations[entrance.near_role] == station))
                self.assertEqual({interval.section_id for interval in spec.source_intervals}, {1, 2, 3})
                data = BoundaryConditions().build(sampled, sampled.loops[0],
                    {f.region.owner: f.layout for f in sampled.footprints})
                self.assertEqual(len(data.ports), 3)
                paired = [p for p in data.ports if len(p.aliases) == 2]
                self.assertEqual(len(paired), 2)
                for port in paired:
                    self.assertEqual(sum(a.original_incidence for a in port.aliases), 1)
                    derived = next(a for a in port.aliases if not a.original_incidence)
                    self.assertEqual(derived.endpoint.section_id, 3)
                    self.assertGreater(derived.station, 0.)
                    self.assertLess(derived.station, context[model.Key('path', 3)].length)
                    approach = spec.group.physical_layout.approaches[port.approach_index]
                    for alias in port.aliases:
                        u = np.asarray(port.outward_u)
                        u = u[alias.contains(u)]
                        original_profile = context[model.Key('profile', alias.endpoint.section_id)].oriented(alias.endpoint.end)
                        np.testing.assert_allclose(original_profile.value(2.*alias.outward_parameter(u)-1.),
                            approach.profile.value(2.*u-1.), atol=1e-12, rtol=0.)
                context.products[site.layout_key] = spec.group.physical_layout
                context.task_digests = {site.layout_key: 'explicit-one-sided-test-layout'}
                index = SeparatorPortIndex.build(context)
                self.assertIsNone(index.resolve(model.End(3, 0)))
                self.assertIsNone(index.resolve(model.End(3, 1)))
                for key in spec.opposite_cut_keys:
                    self.assertIsNotNone(index.resolve_derived(key.owner, key.section_id, key.end))
                for endpoint in site.ends:
                    self.assertIsNotNone(index.resolve(endpoint))
                self.assertEqual(context.network.to_dict(), before)

    def test_no_authored_opposite_endpoint_or_impossible_merged_fibre_is_invented(self):
        context, site, layout, visits, opposite = fixture(separation=14.)
        with self.assertRaisesRegex(ValueError, 'no merged finite outward fibre'):
            OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
        with self.assertRaisesRegex(ValueError, 'no authored incidence'):
            OneSidedPairedSpecification.build(context, site, layout, visits, visits[0], 0)

    def test_normal_product_keeps_derived_source_samples_in_a_separate_registry(self):
        context, site, layout, visits, opposite = fixture()
        spec = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
        context.products[site.layout_key] = spec.group.physical_layout
        product = CompoundProduct.build(context, site.key, spec.group)
        self.assertEqual(set(product.peaks), {(e.section_id, e.end) for e in site.ends})
        self.assertTrue(all(sid != 3 for sid, end, u in product.port_points))
        self.assertEqual(len(product.derived_peaks), 2)
        self.assertGreater(len(product.derived_port_points), 4)
        for (owner, sid, end, u), (li, vertex, station) in product.derived_port_points.items():
            self.assertEqual((owner, sid), (site.key, 3))
            self.assertEqual(product.derived_peaks[owner, sid, end][2], station)
            self.assertTrue(np.isfinite(product.rims[li][vertex]).all())
        self.assertLessEqual(product.pavement.mesh.validate()['maximum_edge_ratio'], 8.)
        quality = product.pavement.metadata['one_sided_source_quality']
        self.assertLessEqual(quality['stored']['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(quality['stored']['minimum_corner_sine'], .03)
        self.assertIsNotNone(product.pavement.metadata['source_cycle_reconciliation'])

    def test_older_transport_rechecks_original_hidden_mouth_labels(self):
        context, site, layout, visits, opposite = fixture()
        spec = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
        physical = copy(spec.group.physical_layout)
        physical.one_sided_original_seams_hidden = False
        group = replace(spec.group, physical_layout=physical)
        self.assertTrue(certify(group)['original_local_footprint'])
        physical.one_sided_original_seams_hidden = False
        physical.one_sided_original_regions = physical.one_sided_original_regions[:1]
        with self.assertRaisesRegex(ValueError, 'internal original near-side mouth'):
            certify(group)

    def test_opposite_interior_alias_cannot_shadow_its_actual_remote_source_end(self):
        context, site, layout, visits, opposite = fixture()
        spec = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
        path = context[model.Key('path', 3)]
        stem = model.Section(5, 5, (tuple((*path.frame(0.)[0], 0.)), (-30., -30., 0.)), width=2.)
        remote = model.Intersection(200, (model.End(3, 0), model.End(5, 0)))
        context.network = replace(context.network, sections=(*context.network.sections, stem),
            intersections=(*context.network.intersections, remote))
        context.products[site.layout_key] = spec.group.physical_layout
        context.products[remote.layout_key] = SimpleNamespace(approaches=(
            Approach(path, 0, context[model.Key('profile', 3)]),
            Approach(Path(stem), 0, Crown(.08, 0.))), outward_stations=(1., 1.))
        context.task_digests = {site.layout_key: 'paired-prepared', remote.layout_key: 'original-remote'}
        index = SeparatorPortIndex.build(context)
        original = index.resolve(model.End(3, 0))
        self.assertEqual(original.site, remote)
        derived = index.resolve_derived(site.key, 3, 0)
        self.assertEqual(derived.owner, site.key)
        self.assertFalse(derived.alias.original_incidence)
        self.assertGreater(derived.alias.station, 1.)

    def test_real_separator_hole_survives_bounded_sampling_and_complete_source_proofs(self):
        from tested_python_streets.one_sided_rim_sampling import normalize_holes
        context, site, layout, visits, opposite = fixture(middle_dip=2.2, near_bend=0.)
        spec = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
        self.assertEqual(sum(loop.hole for loop in spec.group.loops), 1)
        originals = [loop.xy.copy() for loop in spec.group.loops]
        sampled = normalize_holes(spec.group)
        self.assertEqual(certify(sampled)['holes'], 1)
        np.testing.assert_array_equal(sampled.loops[0].xy, originals[0])
        for loop, original in zip(spec.group.loops, originals):
            np.testing.assert_array_equal(loop.xy, original)
        for row in sampled.physical_layout.source_hole_normalization:
            self.assertLessEqual(row['maximum_projection_error'], 1e-6)
        with self.assertRaises(ValueError):
            normalize_holes(spec.group, maximum_projection=1e-5)

    def test_source_cuts_preserve_whole_profiles_and_genuine_hole_topology(self):
        from tested_python_streets.one_sided_paired_cuts import OneSidedPairedCuts
        from tested_python_streets.partition import CompoundPartitioner
        for dip, holes in ((.1, 0), (2.2, 1)):
            context, site, layout, visits, opposite = fixture(middle_dip=dip, near_bend=0.)
            spec = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
            group = AnalyticalBoundarySampler().build(spec.group)
            partition = CompoundPartitioner(policy=OneSidedPairedCuts(spec.source_intervals[:2]),
                require_visibility=False).build(group, context)
            self.assertLess(partition.area_error, 1e-7)
            self.assertEqual(sum(loop.hole for loop in group.loops), holes)
            for li, loop in enumerate(group.loops):
                boundary = partition.xy[list(partition.boundary_loops[li])]
                for point in loop.xy:
                    self.assertTrue(np.any(np.all(boundary == point, axis=1)))
            with self.assertRaisesRegex(ValueError, 'owned original near-source'):
                tuple(OneSidedPairedCuts(spec.source_intervals[-1:]).cuts(group, context))


if __name__ == '__main__':
    unittest.main()
