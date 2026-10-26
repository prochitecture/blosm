from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.drainage import BoundaryDrainage
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.slope_projection import SlopeProjection
from tested_python_streets.providers import Part


def basin():
    pattern = PatternLibrary().grid(6, 6)
    xy = 6. * pattern.uv
    fixed = np.any((xy == 0.) | (xy == 6.), axis=1)
    heights = np.full(len(xy), .4)
    heights[fixed] = 0.
    inside = np.all((xy >= 2.) & (xy <= 4.), axis=1)
    heights[inside] = .1
    return Mesh(np.column_stack((xy, heights)), pattern.quads), fixed, inside


def refined_cross():
    origin = (-321.45941162109375, 117.5704345703125, 0.)
    rise = .07999999821186066
    sections = (
        model.Section(5, 159, ((-324.4237060546875, 174.3323516845703, 0.), origin),
                      width=10.800000190734863, rise=rise, category='primary'),
        model.Section(9, 159, (origin, (-319.3850402832031, 59.46158218383789, 0.)),
                      width=10.800000190734863, rise=rise, category='primary'),
        model.Section(15, 164, ((-312.9646911621094, 117.79276275634766, 0.), origin),
                      width=2., rise=rise, category='service'),
        model.Section(10, 10, (origin, (-409.9042053222656, 113.45535278320312, 0.)),
                      width=2., rise=rise, category='service', selection_ends=(1,)),
        model.Section(75, 164, ((-198.56924438476562, 121.52974700927734, 0.),
                               (-312.9646911621094, 117.79276275634766, 0.)),
                      width=2., rise=rise, category='service'))
    site = model.Intersection(29, tuple(model.End(i, e) for i, e in ((5, 1), (9, 0), (15, 1), (10, 0))),
        radius=4., radii=(4., 4., 1., 1.), radius_modes=('STYLE',) * 4)
    return model.Network('refined_narrow_cross', sections, (site,),
        connections=(model.Connection(3, (model.End(75, 1), model.End(15, 0))),),
        settings=model.Settings(exaggeration=1.))


class PythonDrainageTests(unittest.TestCase):
    def test_a_flat_basin_has_a_nonascending_path_to_the_boundary_after_minimal_fill(self):
        mesh, fixed, inside = basin()
        original = mesh.xyz.copy()
        heights, details = BoundaryDrainage().fill(mesh, fixed)
        np.testing.assert_array_equal(heights[inside], np.full(np.sum(inside), .4))
        np.testing.assert_array_equal(heights[fixed], mesh.xyz[fixed, 2])
        np.testing.assert_array_equal(mesh.xyz, original)
        self.assertEqual(details['raised_vertices'], 9)
        self.assertAlmostEqual(details['maximum_fill'], .3)
        # Independent reachability check, following only nondecreasing heights
        # outward from the fixed outlets (the reverse of drainage direction).
        reached = set(np.flatnonzero(fixed))
        while True:
            prior = len(reached)
            for face in mesh.quads:
                for a, b in zip(face, np.roll(face, -1)):
                    if a in reached and heights[b] >= heights[a] - 1e-12: reached.add(int(b))
                    if b in reached and heights[a] >= heights[b] - 1e-12: reached.add(int(a))
            if len(reached) == prior: break
        self.assertEqual(len(reached), len(heights))

    def test_face_order_and_common_datum_do_not_change_the_spill_solution(self):
        mesh, fixed, _ = basin()
        first, _ = BoundaryDrainage().fill(mesh, fixed)
        changed = Mesh(mesh.xyz + (0., 0., 37.), mesh.quads[::-1])
        other, _ = BoundaryDrainage().fill(changed, fixed)
        np.testing.assert_allclose(other - 37., first, atol=1e-14, rtol=0.)
        changed.xyz[:, 2] = other
        self.assertEqual(BoundaryDrainage().fill(changed, fixed)[1]['raised_vertices'], 0)

    def test_missing_outlets_and_budgets_are_explicit_errors(self):
        mesh, fixed, _ = basin()
        for solver in (BoundaryDrainage(max_cells=1), BoundaryDrainage(memory_bytes=4096)):
            with self.assertRaisesRegex(ValueError, 'budget'): solver.fill(mesh, fixed)
        with self.assertRaisesRegex(ValueError, 'Boolean'): BoundaryDrainage().fill(mesh, np.flatnonzero(fixed))
        with self.assertRaisesRegex(ValueError, 'outlets'): BoundaryDrainage().fill(mesh, np.zeros(len(fixed), bool))
        disconnected = Mesh(np.vstack((mesh.xyz, mesh.xyz + (10., 0., 0.))),
            np.vstack((mesh.quads, mesh.quads + len(mesh.xyz))))
        with self.assertRaisesRegex(ValueError, 'no fixed drainage outlet'):
            BoundaryDrainage().fill(disconnected, np.r_[fixed, np.zeros(len(fixed), bool)])

    def test_floor_projection_keeps_the_boundary_and_rejects_infeasible_floors(self):
        pattern = PatternLibrary().grid(2, 2)
        mesh = Mesh(np.column_stack((2. * pattern.uv, np.zeros(len(pattern.uv)))), pattern.quads)
        fixed = np.ones(len(mesh.xyz), bool); fixed[4] = False
        floor = np.zeros(len(mesh.xyz)); floor[4] = .03
        candidate, _ = SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08,
            maximum_slope=.05, minimum_heights=floor)
        self.assertAlmostEqual(candidate.xyz[4, 2], .03)
        np.testing.assert_array_equal(candidate.xyz[fixed], mesh.xyz[fixed])
        self.assertLessEqual(candidate.quad_slopes().max(), .05 + 1e-9)
        floor[4] = .08
        with self.assertRaisesRegex(ValueError, 'iteration budget'):
            SlopeProjection(max_iterations=20).solve(mesh, fixed, datum=0., ceiling=.08,
                maximum_slope=.05, minimum_heights=floor)
        for bad in (np.full(len(mesh.xyz), -.1), np.full(len(mesh.xyz), np.nan), np.zeros(2)):
            with self.assertRaisesRegex(ValueError, 'minimum heights'):
                SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08,
                    maximum_slope=.05, minimum_heights=bad)

    def test_real_unequal_width_cross_uses_bounded_regrading_and_preserves_profiles(self):
        network = refined_cross()
        result = PythonMeshBackend(network.settings).prepare(network)
        part = result.products[model.Key('junction_mesh', 29)]
        self.assertIsNotNone(part.metadata.get('refined_grading'))
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertLessEqual(part.metadata['physical_slopes']['maximum'], part.metadata['measured_slope_limit'] + 1e-8)
        details = part.metadata['refined_grading']
        self.assertTrue(details['exact_boundary'])
        self.assertEqual(details['joins']['isolated_interior_minima'], 0)
        self.assertLessEqual(details['joins']['maximum_entrance_longitudinal_slope'], .03)
        self.assertLessEqual(details['rounds'][-1]['drainage']['maximum_fill'], 1e-6)

    def test_regraded_field_queries_and_local_edits_survive_publication_release(self):
        network = refined_cross()
        network = replace(network, furniture=(model.Placement(81, 5, 55.),))
        backend = PythonMeshBackend(network.settings)
        first = backend.prepare(network)
        surface = first.products[model.Key('junction_surface', 29)]
        mesh = first.products[model.Key('junction_mesh', 29)].mesh
        centers = mesh.xyz[mesh.quads].mean(axis=1)
        np.testing.assert_allclose(surface.field(centers[:, :2]), centers[:, 2], atol=1e-10, rtol=0.)

        def buffers(update, previous=None):
            # Publication owns unchanged buffers after the numerical cache is
            # released. Apply the proposed changes to that saved scene state.
            result = dict(previous or {})
            for key in update.removed:
                result.pop(key, None)
            result.update({k: (p.mesh.xyz.copy(), p.mesh.quads.copy()) for k, p in update.products.items()
                           if isinstance(p, Part) and p.mesh is not None})
            return result

        original = buffers(first)
        anchor = first.products[model.Key('placement', 81)]
        backend.commit(first)
        backend.release_published_meshes()
        np.testing.assert_allclose(surface.field(centers[:, :2]), centers[:, 2], atol=1e-10, rtol=0.)
        point = network.sections[0].points[0]
        changed_section = replace(network.sections[0], points=((point[0] + .05, *point[1:]),
                                                              network.sections[0].points[1]))
        changed = replace(network, sections=(changed_section, *network.sections[1:]))
        local = backend.prepare(changed)
        full = PythonMeshBackend(changed.settings).prepare(changed)
        a, b = buffers(local, original), buffers(full)
        self.assertEqual(a.keys(), b.keys())
        for key in a:
            np.testing.assert_array_equal(a[key][0], b[key][0])
            np.testing.assert_array_equal(a[key][1], b[key][1])
        self.assertNotEqual(local.products[model.Key('placement', 81)], anchor)
        self.assertEqual(local.products[model.Key('placement', 81)], full.products[model.Key('placement', 81)])
        backend.commit(local)
        restored = buffers(backend.prepare(network), a)
        for key in original:
            np.testing.assert_array_equal(original[key][0], restored[key][0])
            np.testing.assert_array_equal(original[key][1], restored[key][1])


if __name__ == '__main__':
    unittest.main(verbosity=2)
