from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.providers import RoadSpan
from _python_bundle_fixtures import fixture_records, make_network


def mesh_parts(update):
    return {k: p for k, p in update.products.items() if hasattr(p, 'mesh')}


def area(mesh):
    p = mesh.xyz[mesh.quads, :2]
    return .5 * float(np.sum(p[:, :, 0] * np.roll(p[:, :, 1], -1, axis=1)
                             - p[:, :, 1] * np.roll(p[:, :, 0], -1, axis=1)))


class PythonBundleUnionTests(unittest.TestCase):
    def test_reference_families_have_dynamic_islands_and_no_duplicate_road_meshes(self):
        for fixture in fixture_records():
            with self.subTest(case=fixture['name']):
                n = make_network(model, fixture)
                if not fixture['valid']:
                    with self.assertRaisesRegex(ValueError, 'centerlines cross or touch'):
                        PythonMeshBackend(n.settings).prepare(n)
                    continue
                update = PythonMeshBackend(n.settings).prepare(n)
                plan = update.products[model.Key('bundle_plan', 7)]
                if fixture['name'] == 'curved_separated':
                    self.assertEqual(plan.mode, 'separated')
                    continue
                self.assertEqual(plan.mode, 'union')
                self.assertIsInstance(update.products[model.Key('road', 1)], RoadSpan)
                self.assertIsInstance(update.products[model.Key('road', 2)], RoadSpan)
                pavement = update.products[model.Key('bundle_pavement', 7)]
                if fixture['islands'] is not None:
                    self.assertEqual(pavement.metadata['islands'], fixture['islands'])
                for part in mesh_parts(update).values():
                    self.assertTrue(np.all(part.mesh.xyz[:, 2] >= 0.))
                    self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), 1.2 + 1e-10)
                    self.assertEqual(part.mesh.quads.shape[1], 4)

    def test_straight_union_area_profile_and_containment_are_analytically_correct(self):
        fixtures = {f['name']: f for f in fixture_records()}
        for name, width in (('exact_line_contact', 12.), ('contained_member', 14.)):
            n = make_network(model, fixtures[name])
            update = PythonMeshBackend(n.settings).prepare(n)
            p = update.products[model.Key('bundle_pavement', 7)]
            self.assertAlmostEqual(area(p.mesh), 96. * width, places=8)
            self.assertIsNone(update.products[model.Key('separator', 7)])
            plan = update.products[model.Key('bundle_plan', 7)]
            a, b = plan.peaks(np.asarray(48.))
            self.assertEqual(a, b)
            xy = plan.world(np.full(3, 48.), np.array((float(a) - width / 2., float(a), float(a) + width / 2.)))
            np.testing.assert_allclose(plan.height(xy), (0., 1.2, 0.), atol=1e-12)

    def test_switching_separate_union_and_back_matches_full_rebuild(self):
        f = fixture_records()[0]
        n = make_network(model, f)
        engine = PythonMeshBackend(n.settings)
        original = engine.prepare(n)
        engine.commit(original)
        sections = tuple(replace(s, points=tuple((x, y + (-2. if s.id == 1 else 2.), z) for x, y, z in s.points)) for s in n.sections)
        edited = replace(n, sections=sections)
        update = engine.prepare(edited)
        full = PythonMeshBackend(n.settings).prepare(edited)
        self.assertEqual(update.products[model.Key('bundle_plan', 7)].mode, 'union')
        for key, part in mesh_parts(update).items():
            np.testing.assert_array_equal(part.mesh.xyz, full.products[key].mesh.xyz)
            np.testing.assert_array_equal(part.mesh.quads, full.products[key].mesh.quads)
        engine.commit(update)
        restored = engine.prepare(n)
        self.assertIsNone(restored.products[model.Key('bundle_pavement', 7)])
        for key, part in mesh_parts(restored).items():
            np.testing.assert_array_equal(part.mesh.xyz, original.products[key].mesh.xyz)
            np.testing.assert_array_equal(part.mesh.quads, original.products[key].mesh.quads)

    def test_source_relative_furniture_uses_the_union_surface(self):
        f = fixture_records()[1]
        n = make_network(model, f)
        n = replace(n, furniture=(model.Placement(81, 1, 30., side='RIGHT', offset=-.2),))
        result = PythonMeshBackend(n.settings).prepare(n)
        anchor = result.products[model.Key('placement', 81)]
        plan = result.products[model.Key('bundle_plan', 7)]
        self.assertAlmostEqual(anchor.position[2], float(plan.height(np.asarray(anchor.position[:2]))))

    def test_reversing_lateral_roles_without_crossing_is_a_diagnostic(self):
        n = make_network(model, fixture_records()[1])
        b = n.bundles[0]
        n = replace(n, bundles=(replace(b, members=b.members[::-1]),))
        with self.assertRaisesRegex(ValueError, 'left-to-right order'):
            PythonMeshBackend(n.settings).prepare(n)

    def test_offset_can_reorder_physical_carriageways_without_filling_the_gap(self):
        # Centerlines remain ordered and disjoint, but the allowed carriageway
        # offset places A below B. The actual two pavement intervals are still
        # disjoint: 12 m pavement plus an 8 m separator, not 20 m of asphalt.
        n = make_network(model, fixture_records()[7])
        n = replace(n, sections=(replace(n.sections[0], offset=20.), n.sections[1]))
        update = PythonMeshBackend(n.settings).prepare(n)
        pavement = update.products[model.Key('bundle_pavement', 7)]
        separator = update.products[model.Key('separator', 7)]
        self.assertIsNotNone(separator)
        self.assertAlmostEqual(area(pavement.mesh), 96. * 12., places=8)
        self.assertAlmostEqual(area(separator.mesh), 96. * 8., places=8)


if __name__ == '__main__':
    unittest.main()
