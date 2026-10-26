"""Preferred-radius fitting must not weaken authored geometry contracts."""
from dataclasses import replace
import json
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from _python_street_fixtures import fixture_networks, radius_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.corner_radii import CornerContactError, CornerRadiusPolicy, fitted_radius_report
from tested_python_streets.intersections import Approach, _corner
from tested_python_streets.paths import Path
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import UpdateCoordinator


def curved_pair():
    # Two 5.4 m arms diverge at the node and then become nearly parallel.
    # Their finite inner-curb separation cannot accommodate the 4 m fillet.
    sections = (model.Section(1, 1, ((0., 0., 0.), (-13., -16., 0.), (-32., -56., 0.)), width=5.4),
                model.Section(2, 2, ((0., 0., 0.), (-7., -24., 0.), (-21., -56., 0.)), width=5.4))
    network = model.Network('bent-fork', sections)
    products = UpdateCoordinator((PathProvider(), ProfileProvider()), PatternLibrary()).prepare(network).products
    return tuple(Approach(products[model.Key('path', i)], 0, products[model.Key('profile', i)]) for i in (1, 2))


class PythonCornerRadiusTests(unittest.TestCase):
    def test_exact_miter_curbs_use_unit_radius_normals_and_tangent_contacts(self):
        for reverse in (False, True):
            for angle in (0., .63):
                for second_angle in (1.2, 1.9):
                    with self.subTest(reverse=reverse, angle=angle, second_angle=second_angle):
                        def rotated(points, theta):
                            c, s = np.cos(theta), np.sin(theta)
                            return tuple((c*x-s*y, s*x+c*y, 0.) for x,y in points)
                        first = rotated(((0., 0.), (12., .2), (24., 0.)), angle)
                        second = rotated(((0., 0.), (12., -.3), (24., 0.)), angle+second_angle)
                        sections = (model.Section(1, 1, first[::-1] if reverse else first,
                                                  width=6., offset=.3 if reverse else -.3, interpolation=0),
                                    model.Section(2, 2, second, width=5., offset=.4, interpolation=0))
                        # Reversing the authored Section must not change the
                        # physical corner approached from its connected end.
                        from tested_python_streets.profiles import Crown
                        approaches = (Approach(Path(sections[0]), int(reverse), Crown(.12, 0.)),
                                      Approach(Path(sections[1]), 0, Crown(.12, 0.)))
                        corner = _corner(*approaches, 2.)
                        distances = (approaches[0].offset-.5*approaches[0].width,
                                     approaches[1].offset+.5*approaches[1].width)
                        for approach, station, distance in zip(approaches, corner.stations, distances):
                            contact, tangent = approach.jet(station, distance)
                            radial = contact-corner.center
                            self.assertAlmostEqual(float(np.linalg.norm(radial)), 2., delta=2e-8)
                            self.assertAlmostEqual(float(radial@tangent)/np.linalg.norm(tangent), 0., delta=2e-8)
                        self.assertEqual(tuple(s.points for s in sections),
                                         (first[::-1] if reverse else first, second))

    def test_old_records_remain_strict_and_explicit_override_changes_only_one_incidence(self):
        source = fixture_networks(model)[0]
        data = source.to_dict()
        for site in data['intersections']:
            site.pop('radius_modes')
        reopened = model.Network.from_dict(json.loads(json.dumps(data)))
        site = reopened.intersections[0]
        self.assertTrue(all(site.corner_radius(e) == (site.radius, 'STRICT') for e in site.ends))
        preferred = replace(site, radius_modes=('STYLE',) * len(site.ends))
        override = preferred.with_corner_radius(site.ends[1], 4.)
        self.assertEqual(override.corner_radius(site.ends[1]), (4., 'STRICT'))
        for end in (site.ends[0], site.ends[2]):
            self.assertEqual(override.corner_radius(end), (site.radius, 'STYLE'))
        network = replace(source, intersections=(override,))
        self.assertEqual(model.Network.from_dict(json.loads(json.dumps(network.to_dict()))), network)
        with self.assertRaises(ValueError):
            preferred.with_corner_radius(model.End(999, 0), 2.)

    def test_invalid_mode_dimensions_or_values_are_rejected(self):
        site = fixture_networks(model)[0].intersections[0]
        for modes in (('STYLE',), ('STYLE', 'unknown', 'STRICT')):
            with self.assertRaisesRegex(ValueError, 'radius modes'):
                replace(site, radius_modes=modes)
        with self.assertRaisesRegex(ValueError, 'positive'):
            replace(site, radius=0.)
        with self.assertRaises(ValueError):
            site.with_corner_radius(site.ends[0], -1.)

    def test_fitting_preserves_the_exact_circle_contact_and_correct_sector(self):
        a, b = curved_pair()
        controls = tuple(approach.source.points for approach in (a, b))
        with self.assertRaises(CornerContactError):
            _corner(a, b, 4.)
        corner = CornerRadiusPolicy.solve(a, b, 4., 'STYLE')
        self.assertLess(corner.radius, 4.)
        self.assertGreaterEqual(corner.radius, 1.)
        self.assertEqual(corner.requested_radius, 4.)
        self.assertLess(corner.sweep, 0.)
        self.assertGreaterEqual(corner.sweep, -np.pi)
        contacts = np.array((a.point(corner.stations[0], -.5 * a.width),
                             b.point(corner.stations[1], .5 * b.width)))
        np.testing.assert_allclose(np.linalg.norm(contacts - corner.center, axis=1), corner.radius, atol=1e-7)
        report, = fitted_radius_report(model.Key('intersection', 8), SimpleNamespace(corners=(corner,)))
        self.assertEqual(report['requested_radius'], 4.)
        self.assertEqual(report['actual_radius'], corner.radius)
        self.assertEqual((report['first_section'], report['second_section']), (1, 2))
        self.assertEqual(tuple(approach.source.points for approach in (a, b)), controls)

    def test_explicit_infeasible_override_is_never_reduced(self):
        a, b = curved_pair()
        with self.assertRaisesRegex(CornerContactError, 'Requested-radius'):
            CornerRadiusPolicy.solve(a, b, 4., 'STRICT')

    def test_style_fitting_is_bounded_and_does_not_erase_a_turn_with_a_tiny_fillet(self):
        a, b = curved_pair()
        with self.assertRaisesRegex(CornerContactError, 'Preferred style radius could not be fitted'):
            CornerRadiusPolicy.solve(a, b, 20., 'STYLE')
        candidates = tuple(CornerRadiusPolicy.candidates(20., 'STYLE'))
        self.assertLessEqual(len(candidates), 9)
        self.assertEqual(min(candidates), 5.)
        self.assertTrue(all(x > y for x, y in zip(candidates, candidates[1:])))

    def test_a_feasible_style_radius_preserves_all_existing_mesh_coordinates(self):
        source = fixture_networks(model)[0]
        preferred = replace(source, intersections=tuple(replace(site, radius_modes=('STYLE',) * len(site.ends))
                                                        for site in source.intersections))
        strict = PythonMeshBackend(source.settings).prepare(source)
        fit = PythonMeshBackend(preferred.settings).prepare(preferred)
        for key, part in strict.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fit.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fit.products[key].mesh.quads)
        self.assertEqual(fit.products[source.intersections[0].mesh_key].metadata['fitted_corner_radii'], [])

    def test_fitting_does_not_relax_indistinguishable_approach_directions(self):
        a, _ = curved_pair()
        with self.assertRaisesRegex(ValueError, 'indistinguishable'):
            CornerRadiusPolicy.solve(a, a, 4., 'STYLE')

    def test_style_fit_includes_profile_entrance_reach_after_a_valid_circle_contact(self):
        source, = radius_fixture_networks(model)
        north = replace(source.sections[1], points=((0., 0., 0.), (0., 5.1, 0.)))
        site = replace(source.intersections[0], radius=2.)
        network = replace(source, sections=(source.sections[0], north, *source.sections[2:]),
                          intersections=(site,))
        strict = replace(network, intersections=(replace(site, radius_modes=('STRICT',) * 3),))
        with self.assertRaisesRegex(ValueError, 'complete compatible approach'):
            PythonMeshBackend(network.settings).prepare(strict)
        update = PythonMeshBackend(network.settings).prepare(network)
        layout = update.products[model.Key('mouths', site.id)]
        fitted = fitted_radius_report(site.key, layout)
        self.assertEqual(len(fitted), 2)
        self.assertTrue(all(fit['actual_radius'] == 1.6 for fit in fitted))
        self.assertTrue(all(fit['attempts'] == [2., 1.6] for fit in fitted))
        self.assertTrue(all(s < a.length for a, s in zip(layout.approaches, layout.outward_stations)))
        self.assertEqual(network.intersections[0].radius, 2.)
        self.assertLessEqual(update.products[site.mesh_key].mesh.validate()['maximum_edge_ratio'], 8.)

    def test_complete_fitted_mesh_preserves_profiles_and_matches_incremental_full_restoration(self):
        network, = radius_fixture_networks(model)
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        part = initial.products[model.Key('junction_mesh', 1)]
        self.assertEqual(len(part.metadata['fitted_corner_radii']), 2)
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(part.mesh.xyz[:, 2].min(), 0.)
        self.assertLessEqual(part.mesh.xyz[:, 2].max(), .8 + 1e-8)
        layout = initial.products[model.Key('mouths', 1)]
        for a, station in zip(layout.approaches, layout.outward_stations):
            q = np.linspace(-1., 1., 2 * layout.radial + 1)
            cap = np.column_stack((a.point(station, a.offset + .5 * a.width * q),
                                   a.profile.value(q) * network.settings.exaggeration))
            self.assertLess(np.linalg.norm(part.mesh.xyz[:, None] - cap[None], axis=-1).min(axis=0).max(), 1e-8)
        for boundary in layout.boundaries:
            xy = boundary(np.linspace(0., 1., layout.rows + 1))
            distances = np.linalg.norm(part.mesh.xyz[:, None, :2] - xy[None], axis=-1)
            nearest = distances.argmin(axis=0)
            self.assertLess(distances.min(axis=0).max(), 1e-8)
            np.testing.assert_allclose(part.mesh.xyz[nearest, 2], 0., atol=1e-9)
        strict = replace(network, intersections=(network.intersections[0].with_corner_radius(model.End(1, 0), 10.),))
        with self.assertRaisesRegex(ValueError, 'Requested-radius'):
            backend.prepare(strict)
        self.assertIs(backend.updates.products[model.Key('junction_mesh', 1)], part)
        longer = replace(network, sections=(network.sections[0],
            replace(network.sections[1], points=((0., 0., 0.), (0., 12.5, 0.))), *network.sections[2:]))
        incremental = backend.prepare(longer)
        full = PythonMeshBackend(longer.settings).prepare(longer)
        self.assertIs(incremental.products[model.Key('road', 4)], initial.products[model.Key('road', 4)])
        for key, product in incremental.products.items():
            if getattr(product, 'mesh', None) is not None:
                np.testing.assert_array_equal(product.mesh.xyz, full.products[key].mesh.xyz)
                np.testing.assert_array_equal(product.mesh.quads, full.products[key].mesh.quads)
        backend.commit(incremental)
        restored = backend.prepare(network)
        np.testing.assert_array_equal(restored.products[model.Key('junction_mesh', 1)].mesh.xyz, part.mesh.xyz)


if __name__ == '__main__':
    unittest.main()
