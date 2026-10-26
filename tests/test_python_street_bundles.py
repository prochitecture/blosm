from dataclasses import replace
import itertools
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundles import certify_disjoint
from tested_python_streets.paths import Path
from tested_python_streets.patterns import PatternLibrary
from _python_street_fixtures import divided_network


def free_bundle(points_a, points_b, *, width_a=6., width_b=6.):
    sections = (model.Section(1, 10, points_a, width=width_a), model.Section(2, 20, points_b, width=width_b))
    bundle = model.Bundle(7, (model.BundleMember((model.Traversal(1),)), model.BundleMember((model.Traversal(2),))))
    return model.Network('free_bundle', sections, bundles=(bundle,))


def boundary_edges(mesh):
    edges, counts = np.unique(np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1),
                              return_counts=True, axis=0)
    assert np.all(counts <= 2)
    return mesh.xyz[edges[counts == 1]]


class PythonStreetBundleTests(unittest.TestCase):
    def test_general_quad_atlas_retains_rim_counts_and_has_no_t_junctions(self):
        library = PatternLibrary()
        for counts in itertools.product((4, 8, 12), repeat=4):
            pattern = library.quad(*counts)
            self.assertEqual(tuple(len(library.boundary(pattern, i)) - 1 for i in range(4)), counts)
            xyz = np.column_stack((pattern.uv, np.zeros(len(pattern.uv))))
            from tested_python_streets.mesh import Mesh
            mesh = Mesh(xyz, pattern.quads)
            mesh.validate()
            edges = boundary_edges(mesh)
            self.assertEqual(len(edges), sum(counts))
            # A connected disk of quads with precisely the requested boundary.
            self.assertEqual(len(xyz) - (2 * len(pattern.quads) + sum(counts) // 2) + len(pattern.quads), 1)
            area = .5 * np.sum(np.cross(pattern.uv[pattern.quads], np.roll(pattern.uv[pattern.quads], -1, axis=1)))
            self.assertAlmostEqual(area, 1.)

    def test_rounded_divided_t_uses_shared_boundary_and_bounded_plateau(self):
        network = divided_network(model)
        result = PythonMeshBackend(network.settings).prepare(network)
        separator = result.products[model.Key('separator', 7)]
        mesh = separator.mesh
        self.assertEqual(mesh.material, 'grass')
        self.assertEqual(separator.metadata['attached_caps'], [True, False])
        self.assertLess(separator.metadata['quality']['maximum_edge_ratio'], 4.)
        self.assertLess(len(mesh.quads), 3000)
        self.assertGreaterEqual(mesh.xyz[:, 2].min(), 0.)
        self.assertLessEqual(mesh.xyz[:, 2].max(), .8 + 1e-12)
        self.assertTrue(np.all(abs(mesh.xyz[mesh.xyz[:, 1] > 15., 2] - .8) < 1e-12))
        self.assertLess(mesh.slope_statistics(network.settings.exaggeration)['maximum'], .08)
        pavement_edges = np.concatenate([boundary_edges(part.mesh) for part in result.products.values()
                                        if hasattr(part, 'mesh') and part is not separator])
        shared = boundary_edges(mesh)
        # The one unconnected grass end has no neighboring pavement. Every
        # remaining grass boundary segment matches a whole pavement segment.
        shared = shared[~np.all(abs(shared[:, :, 1] - 45.) < 1e-8, axis=1)]
        for edge in shared:
            delta = np.max(abs(pavement_edges - edge), axis=(1, 2))
            reversed_delta = np.max(abs(pavement_edges[:, ::-1] - edge), axis=(1, 2))
            self.assertLess(float(np.minimum(delta, reversed_delta).min()), 1e-8)

    def test_curved_varying_gap_and_unequal_widths_generate_no_pavement_overlap(self):
        network = free_bundle(((-12., 0., 0.), (-14., 25., 0.), (-10., 50., 0.)),
                              ((8., 0., 0.), (14., 25., 0.), (10., 50., 0.)), width_b=8.)
        update = PythonMeshBackend(network.settings).prepare(network)
        part = update.products[model.Key('separator', 7)]
        self.assertEqual(part.metadata['attached_caps'], [False, False])
        self.assertGreater(part.metadata['quality']['minimum_corner'], 0.)
        self.assertTrue(np.all(abs(part.mesh.xyz[:, 2] - .8) < 1e-12))
        self.assertEqual(update.products[model.Key('profile', 1)].half, 1)
        self.assertEqual(update.products[model.Key('profile', 2)].half, -1)

    def test_member_orientation_is_not_inferred_from_polyline_order(self):
        original = free_bundle(((-10., 0., 0.), (-10., 30., 0.)), ((10., 0., 0.), (10., 30., 0.)))
        sections = (original.sections[0], replace(original.sections[1], points=original.sections[1].points[::-1]))
        b = original.bundles[0]
        reversed_network = replace(original, sections=sections, bundles=(replace(b, members=(b.members[0],
            model.BundleMember((model.Traversal(2, False),)))),))
        first = PythonMeshBackend(original.settings).prepare(original)
        second = PythonMeshBackend(original.settings).prepare(reversed_network)
        np.testing.assert_allclose(first.products[model.Key('separator', 7)].mesh.xyz,
                                   second.products[model.Key('separator', 7)].mesh.xyz, atol=1e-8)
        self.assertEqual(second.products[model.Key('profile', 2)].half, 1)

    def test_bundle_dependency_rebuild_matches_full_build_and_preserves_old_failure(self):
        network = divided_network(model)
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        moved = replace(network.sections[0], points=((-10., 0., 0.), (-12., 22., 0.), (-10., 45., 0.)))
        edited = replace(network, sections=(moved, *network.sections[1:]))
        incremental = engine.prepare(edited)
        full = PythonMeshBackend(network.settings).prepare(edited)
        self.assertIn(model.Key('separator', 7), incremental.rebuilt)
        self.assertNotIn(model.Key('road', 5), incremental.rebuilt)
        for key, part in incremental.products.items():
            if hasattr(part, 'mesh'):
                np.testing.assert_array_equal(part.mesh.xyz, full.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, full.products[key].mesh.quads)
        manifest = engine.manifest()
        bad = replace(network, sections=(replace(network.sections[0], points=((10., 0., 0.), (10., 45., 0.))), *network.sections[1:]))
        with self.assertRaises(ValueError):
            engine.prepare(bad)
        self.assertEqual(manifest, engine.manifest())

    def test_contact_certificate_rejects_crossing_touching_and_between_sample_contacts(self):
        a = Path(model.Section(1, 1, ((0., 0., 0.), (100., 0., 0.)), width=.01))
        for points in (((37.123, -10., 0.), (37.123, 10., 0.)),
                       ((37.123, 0., 0.), (37.123, 10., 0.)),
                       ((0., 1., 0.), (50., -.2, 0.), (100., 1., 0.))):
            b = Path(model.Section(2, 2, points, width=.01))
            with self.assertRaisesRegex(ValueError, 'cross or touch'):
                certify_disjoint((a,), (b,))
        b = Path(model.Section(2, 2, ((0., 1., 0.), (50., .2, 0.), (100., 1., 0.)), width=.01))
        self.assertGreater(certify_disjoint((a,), (b,)), 0)

    def test_multiple_sections_of_member_keep_explicit_section_seams(self):
        network = free_bundle(((-10., 0., 0.), (-10., 40., 0.)), ((10., 0., 0.), (10., 40., 0.)))
        sections = (replace(network.sections[0], points=((-10., 0., 0.), (-10., 20., 0.))),
                    network.sections[1], model.Section(3, 10, ((-10., 20., 0.), (-10., 40., 0.))))
        bundle = replace(network.bundles[0], members=(model.BundleMember((model.Traversal(1), model.Traversal(3))),
                                                       network.bundles[0].members[1]))
        network = replace(network, sections=sections, bundles=(bundle,),
            connections=(model.Connection(10, (model.End(1, 1), model.End(3, 0))),))
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertGreater(len(result.products[model.Key('separator', 7)].mesh.quads), 0)

    def test_shared_rise_policy_and_profile_derivatives_survive_roundtrip(self):
        network = divided_network(model)
        bundle = replace(network.bundles[0], profile_rise=.1)
        network = replace(network, bundles=(bundle,))
        import json
        reopened = model.Network.from_dict(json.loads(json.dumps(network.to_dict())))
        self.assertEqual(reopened, network)
        update = PythonMeshBackend(network.settings).prepare(reopened)
        for sid, half in ((1, 1), (2, -1)):
            profile = update.products[model.Key('profile', sid)]
            self.assertEqual(profile.half, half)
            self.assertEqual(float(profile.value(-half)), 0.)
            self.assertEqual(float(profile.value(half)), .1)
            self.assertEqual(float(profile.derivative(half)), 0.)
        mesh = update.products[model.Key('separator', 7)].mesh
        self.assertLessEqual(float(mesh.xyz[:, 2].max()), 1. + 1e-12)

    def test_pattern_size_budget_is_checked_before_allocating(self):
        patterns = PatternLibrary()
        with self.assertRaisesRegex(ValueError, 'budget'):
            patterns.quad(4_000_000, 4, 4, 4)
        with self.assertRaisesRegex(ValueError, 'budget'):
            patterns.strip(4, 4_000_000, 4)


if __name__ == '__main__':
    unittest.main()
