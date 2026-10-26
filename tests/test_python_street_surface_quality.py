"""Dynamic profile, seam, guide, slope and dependency proofs after XY grading."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_bundle_fixtures import bundle_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_axis import AxialPlan
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import UpdateCoordinator, BuildContext
from tested_python_streets.surface_query import QuadHeightField
from tested_python_streets.surface_quality import ProfiledSlopeCorrection
from tested_python_streets.mesh import Mesh
from tested_python_streets.slope_projection import corner_matrices


def context_for(network):
    patterns = PatternLibrary()
    update = UpdateCoordinator((PathProvider(), ProfileProvider()), patterns).prepare(network)
    return BuildContext(network, update.products, patterns)


def rim(mesh):
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return np.unique(unique[counts == 1])


class PythonProfiledQuadCorrectionTests(unittest.TestCase):
    @staticmethod
    def longitudinal(mesh):
        gradients = np.einsum('fcij,fj->fci', corner_matrices(mesh), mesh.xyz[mesh.quads, 2])
        return float(abs(gradients[:, :, 0]).max())

    def test_profile_interpolation_amplification_repairs_only_interior_shape(self):
        pattern = PatternLibrary().grid(4, 4)
        xyz = np.c_[pattern.uv, np.zeros(len(pattern.uv))]
        middle = np.flatnonzero(np.all(pattern.uv == (.5, .5), axis=1))[0]
        xyz[middle, :2] = (.31, .45)
        field = lambda xy: .02*(1.-(2.*xy[:, 1]-1.)**2)
        xyz[:, 2] = field(xyz[:, :2])
        original = Mesh(xyz, pattern.quads)
        saved = original.xyz.copy()
        self.assertLessEqual(original.validate()['maximum_edge_ratio'], 8.)
        self.assertGreater(float(original.quad_slopes().max()), .1)
        corrected, report = ProfiledSlopeCorrection().improve(original, field,
            datum=0., ceiling=.02, maximum_slope=.1, longitudinal_bound=self.longitudinal,
            maximum_longitudinal_slope=.03)
        self.assertTrue(report['applied'] and report['exact_fixed_xyz'])
        self.assertLessEqual(float(corrected.quad_slopes().max()), .1)
        self.assertLessEqual(self.longitudinal(corrected), .03)
        self.assertLessEqual(corrected.validate()['maximum_edge_ratio'], 8.)
        np.testing.assert_array_equal(corrected.quads, original.quads)
        np.testing.assert_array_equal(corrected.xyz[rim(original)], saved[rim(original)])
        movable = np.setdiff1d(np.arange(len(xyz)), rim(original))
        np.testing.assert_array_equal(corrected.xyz[movable, 2], field(corrected.xyz[movable, :2]))
        np.testing.assert_array_equal(original.xyz, saved)

    def test_qualified_profile_mesh_is_returned_exactly(self):
        pattern = PatternLibrary().grid(3, 3)
        xyz = np.c_[pattern.uv, .01*pattern.uv[:, 1]]
        original = Mesh(xyz, pattern.quads)
        corrected, report = ProfiledSlopeCorrection().improve(original, lambda xy: .01*xy[:, 1],
            datum=0., ceiling=.01, maximum_slope=.1, longitudinal_bound=self.longitudinal,
            maximum_longitudinal_slope=.03)
        self.assertIs(corrected, original)
        self.assertFalse(report['applied'])

    def test_fixed_profile_exceeding_physical_limit_is_still_rejected(self):
        pattern = PatternLibrary().grid(1, 1)
        xyz = np.c_[pattern.uv, .2*pattern.uv[:, 0]]
        original = Mesh(xyz, pattern.quads)
        with self.assertRaisesRegex(ValueError, 'unchanged physical limits'):
            ProfiledSlopeCorrection().improve(original, lambda xy: .2*xy[:, 0],
                datum=0., ceiling=.2, maximum_slope=.1, longitudinal_bound=self.longitudinal,
                maximum_longitudinal_slope=.03)
        np.testing.assert_array_equal(original.xyz, xyz)

    def test_four_bundle_regressions_keep_profiles_and_boundaries_without_extra_faces(self):
        names = ('reversed_member', 'smooth_islands', 'three_islands', 'unequal_widths')
        for network in bundle_networks(model, interpolation=3):
            if network.uuid.split(' / ')[-1] not in names:
                continue
            with self.subTest(network=network.uuid):
                context = context_for(network)
                plan = AxialPlan(context, network.bundles[0])
                raw, guides = plan.pavement_chart(context)
                self.assertGreater(raw.validate()['maximum_edge_ratio'], 8.)
                result = PythonMeshBackend(network.settings).prepare(network)
                pavement = result.products[model.Key('bundle_pavement', 7)]
                mesh = pavement.mesh
                self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
                np.testing.assert_array_equal(mesh.quads, raw.quads)
                np.testing.assert_array_equal(mesh.xyz[rim(raw)], raw.xyz[rim(raw)])
                np.testing.assert_array_equal(mesh.xyz[guides], raw.xyz[guides])
                self.assertEqual(mesh.material, raw.material)
                self.assertTrue(pavement.metadata['quad_correction']['exact_fixed_xyz'])
                self.assertLessEqual(mesh.slope_statistics(network.settings.exaggeration)['maximum'], .10)
                movable = np.ones(len(mesh.xyz), bool)
                movable[np.r_[rim(raw), guides]] = False
                np.testing.assert_array_equal(mesh.xyz[movable, 2], plan.height(mesh.xyz[movable, :2]))

    def test_common_datum_and_zero_exaggeration_keep_exact_guides_and_the_elevation_envelope(self):
        original = next(n for n in bundle_networks(model, interpolation=3) if n.uuid.endswith('unequal_widths'))
        for exaggeration in (0., 1., 10.):
            network = replace(original, sections=tuple(replace(s, datum=37.) for s in original.sections),
                              settings=replace(original.settings, exaggeration=exaggeration))
            result = PythonMeshBackend(network.settings).prepare(network)
            part = result.products[model.Key('bundle_pavement', 7)]
            self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
            self.assertGreaterEqual(part.mesh.xyz[:, 2].min(), 37.)
            self.assertLessEqual(part.mesh.xyz[:, 2].max(), 37. + .12 * exaggeration)
            self.assertTrue(part.metadata['quad_correction']['exact_fixed_xyz'])

    def test_profiled_correction_edit_matches_full_build_restores_and_moves_its_anchor(self):
        source = next(n for n in bundle_networks(model, interpolation=3) if n.uuid.endswith('three_islands'))
        network = replace(source, sections=(*source.sections, model.Section(3, 3, ((120., 0., 0.), (140., 0., 0.)))),
            furniture=(model.Placement(81, 1, 30., side='RIGHT', offset=-1.),))
        backend = PythonMeshBackend(network.settings)
        before = backend.prepare(network);backend.commit(before)
        points = list(network.sections[0].points)
        points[4] = (points[4][0], points[4][1] + .02, points[4][2])
        edited = replace(network, sections=(replace(network.sections[0], points=tuple(points)), *network.sections[1:]))
        local = backend.prepare(edited)
        full = PythonMeshBackend(edited.settings).prepare(edited)
        self.assertIs(local.products[model.Key('road', 3)], before.products[model.Key('road', 3)])
        for key, part in local.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, full.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, full.products[key].mesh.quads)
        anchor = local.products[model.Key('placement', 81)]
        self.assertNotEqual(anchor.position, before.products[model.Key('placement', 81)].position)
        mesh = local.products[model.Key('bundle_pavement', 7)].mesh
        physical = type(mesh)(mesh.xyz.copy(), mesh.quads)
        physical.xyz[:, 2] /= edited.settings.exaggeration
        field = QuadHeightField(physical, datum=0., ceiling=.12)
        surface_height = float(field(np.array((anchor.position[:2],)))[0]) * edited.settings.exaggeration
        self.assertLess(abs(surface_height - anchor.position[2]), .005 * edited.settings.exaggeration)
        backend.commit(local)
        restored = backend.prepare(network)
        for key, part in before.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, restored.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, restored.products[key].mesh.quads)
        self.assertEqual(restored.products[model.Key('placement', 81)], before.products[model.Key('placement', 81)])


if __name__ == '__main__':
    unittest.main()
