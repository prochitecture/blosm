"""Authored Section splits and internal profile changes in close junctions."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_seam_fixture
from test_python_street_dimensional_bundles import meshes, reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.compound_coverage import CompoundSourceCoverage
from tested_python_streets.compounds import rim_label
from tested_python_streets.ownership import PLAN, footprint_key
from tested_python_streets.partition import BoundaryOwnershipCuts, CompoundPartitioner, CutPolicy
from tested_python_streets.paths import Path
from tested_python_streets.mesh import Mesh
from tested_python_streets.updates import BuildContext


class PythonCompoundRouteTests(unittest.TestCase):
    def qualify(self, network):
        saved = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        group = result.products[PLAN].by_site[model.Key('intersection', 1)]
        self.assertEqual(group.members, (model.Key('intersection', 1), model.Key('intersection', 2)))
        for sid in (1, 6):
            self.assertIsNone(getattr(result.products[model.Key('road', sid)], 'mesh', None))
        self.assertIsNone(result.products[model.Key('junction_mesh', 2)])
        mesh = result.products[model.Key('junction_mesh', 1)].mesh
        datum = network.section_index[2].datum
        self.assertGreaterEqual(float(mesh.xyz[:, 2].min()), datum - 1e-8)
        self.assertLessEqual(float(mesh.xyz[:, 2].max()), datum + .08 * network.settings.exaggeration + 1e-8)
        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertEqual(saved, network.to_dict())
        return result, group

    def test_splits_widths_and_internal_crowns_preserve_the_physical_profiles(self):
        for mode in (0, 3):
            for width in (6., 6.05, 6.4):
                for rise in (.08, .18, .26, .5):
                    with self.subTest(mode=mode, width=width, rise=rise):
                        self.qualify(compound_seam_fixture(model, mode=mode, width=width, rise=rise))

    def test_single_owner_prefix_is_not_a_false_exterior_port(self):
        for mode in (0, 3):
            for width in (6., 6.4):
                for rise in (.18, .5):
                    with self.subTest(mode=mode, width=width, rise=rise):
                        n = compound_seam_fixture(model, mode=mode, width=width, rise=rise, length=10., split=2.)
                        result, group = self.qualify(n)
                        self.assertEqual({claim.section_id for claim in group.claims}, {6})
                        layouts = {p.region.owner: p for p in group.footprints}
                        exposed = {layouts[span.owner].layout.approaches[label.index].source.id
                                   for loop in group.loops for span in loop.spans
                                   for label in (rim_label(span, layouts),) if label.kind == 'mouth'}
                        self.assertEqual(exposed, {2, 3, 4, 5})

    def test_authoring_directions_and_rigid_transforms_retain_source_ownership(self):
        for mode in (0, 3):
            for ids in ({1}, {6}, {1, 6}, {1, 2, 3, 4, 5, 6}):
                with self.subTest(mode=mode, reversed=ids):
                    self.qualify(reverse_sections(compound_seam_fixture(model, mode=mode, rise=.5), ids))
            n = compound_seam_fixture(model, mode=mode, rise=.5)
            angle = .731
            rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
            sections = []
            for s in n.sections:
                points = np.array(s.points)
                points[:, :2] = points[:, :2] @ rotation.T + (1234., -5678.)
                sections.append(replace(s, points=tuple(map(tuple, points)), datum=37.))
            self.qualify(replace(n, sections=tuple(sections)))

    def test_coincident_computational_cut_retains_both_original_section_ids(self):
        n = compound_seam_fixture(model, width=6., rise=.08)
        result, group = self.qualify(n)
        self.assertIsNone(group.physical_layout)
        self.assertEqual({c.section_id for c in group.claims}, {1, 6})
        context = BuildContext(n, result.products, None)
        cuts = tuple(BoundaryOwnershipCuts().cuts(group, context))
        self.assertEqual(len(cuts), 1)
        self.assertEqual(cuts[0].source_sections, (1, 6))
        self.assertEqual(cuts[0].section_id, 1)
        partition = CompoundPartitioner().build(group, context)
        self.assertEqual(len(partition.blocks), 2)
        self.assertEqual(partition.cuts[0].source_sections, (1, 6))

        class RepeatedAlias(CutPolicy):
            def cuts(self, group, context):
                return (cuts[0], replace(cuts[0], section_id=6, section_ids=()))
        with self.assertRaisesRegex(ValueError, 'same source Section'):
            CompoundPartitioner(policy=RepeatedAlias()).build(group, context)

    def test_coverage_needs_all_original_intervals_and_another_member_terminal(self):
        n = compound_seam_fixture(model)
        paths = {model.Key('path', s.id): Path(s) for s in n.sections}
        context = BuildContext(n, paths, None)
        first, last = model.Key('intersection', 1), model.Key('intersection', 2)
        group = SimpleNamespace(members=(first, last))
        mouth = lambda station: SimpleNamespace(station=station)
        # Both internal Sections are fully covered by only the first member.
        # A fully covered free branch must still keep its exterior profile.
        layouts = {first: SimpleNamespace(mouths={(sid, 0): mouth(paths[model.Key('path', sid)].length)
                                                   for sid in (1, 6, 2)}),
                   last: SimpleNamespace(mouths={})}
        coverage = CompoundSourceCoverage(context, group, layouts)
        self.assertEqual(coverage.fully_owned, {1, 6, 2})
        self.assertEqual(coverage.internal_sections(), {1, 6})
        length = paths[model.Key('path', 1)].length
        layouts[first].mouths[1, 0] = mouth(length * .25)
        layouts[last].mouths[1, 1] = mouth(length * .75)
        coverage = CompoundSourceCoverage(context, group, layouts)
        self.assertNotIn(1, coverage.fully_owned)
        self.assertEqual(coverage.internal_sections(), set())

    def test_no_deferral_at_a_free_terminal_or_without_real_compound_ownership(self):
        for mode in (0, 3):
            n = compound_seam_fixture(model, mode=mode, rise=.5)
            with self.subTest(mode=mode, terminal='free'), self.assertRaises(ValueError):
                PythonMeshBackend(n.settings).prepare(replace(n, intersections=n.intersections[:1]))
            n = compound_seam_fixture(model, mode=mode, rise=.5, length=12., split=2.)
            with self.subTest(mode=mode, terminal='separate'):
                with self.assertRaisesRegex(ValueError, 'actual opposite-junction compound ownership'):
                    PythonMeshBackend(n.settings).prepare(n)

    def test_incompatible_seams_cannot_hide_their_failure_or_advance_the_manifest(self):
        n = compound_seam_fixture(model, rise=.5)
        backend = PythonMeshBackend(n.settings)
        backend.commit(backend.prepare(n))
        saved = backend.manifest()
        for fields in (dict(datum=1.), dict(street_id=88), dict(interpolation=0)):
            changed = replace(n, sections=tuple(replace(s, **fields) if s.id == 6 else s for s in n.sections))
            with self.subTest(fields=fields), self.assertRaises(ValueError):
                backend.prepare(changed)
            self.assertEqual(backend.manifest(), saved)

    def test_explicit_local_full_and_restore_after_releasing_published_buffers(self):
        n = compound_seam_fixture(model, rise=.5)
        backend = PythonMeshBackend(n.settings)
        first = backend.prepare(n)
        saved = meshes(first)
        remote = first.products[model.Key('road', 99)]
        backend.commit(first)
        backend.release_published_meshes()
        changed = replace(n, sections=tuple(replace(s, rise=.26) if s.id == 6 else s for s in n.sections))
        local = backend.prepare(changed)
        full = PythonMeshBackend(n.settings).prepare(changed)
        for key in local.rebuilt:
            a, b = local.products[key], full.products[key]
            if getattr(a, 'mesh', None) is not None:
                np.testing.assert_array_equal(a.mesh.xyz, b.mesh.xyz)
                np.testing.assert_array_equal(a.mesh.quads, b.mesh.quads)
        self.assertIs(local.products[model.Key('road', 99)], remote)
        backend.commit(local)
        backend.release_published_meshes()
        restored = backend.prepare(n)
        for key, arrays in saved.items():
            part = restored.products[key]
            if part.mesh is not None:
                np.testing.assert_array_equal(part.mesh.xyz, arrays[0])
                np.testing.assert_array_equal(part.mesh.quads, arrays[1])
        self.assertIs(restored.products[model.Key('road', 99)], remote)

    def test_slope_bound_survives_actual_float32_publication_without_moving_the_rim(self):
        for exaggeration in (0., 1., 10.):
            for datum in (0., 37.):
                with self.subTest(exaggeration=exaggeration, datum=datum):
                    n = compound_seam_fixture(model, rise=.5, length=10., split=2.)
                    n = replace(n, settings=replace(n.settings, exaggeration=exaggeration),
                                sections=tuple(replace(s, datum=datum) for s in n.sections))
                    update, group = self.qualify(n)
                    part = update.products[model.Key('junction_mesh', 1)]
                    rounded = Mesh(part.mesh.xyz.astype(np.float32).astype(float), part.mesh.quads)
                    self.assertLessEqual(float(rounded.quad_slopes(exaggeration).max()), .1 + 1e-7)
                    correction = part.metadata['publication_roundoff']
                    if correction is not None:
                        self.assertTrue(correction['exact_boundary'])
                        self.assertTrue(correction['exact_xy'])
                        self.assertTrue(correction['exact_quads'])
                        self.assertLess(correction['maximum_height_change'], 1e-4)


if __name__ == '__main__':
    unittest.main()
