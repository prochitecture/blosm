"""Published-curb free closures retain source boundaries and physical heights."""
import unittest
from dataclasses import replace
import numpy as np
import _python_street_package
from tested_python_streets.finite_free_closure import FiniteFreeCap, inside_vertex
from tested_python_streets.finite_free_heights import reconcile, boundary_mask
from tested_python_streets.mesh import Mesh
from tested_python_streets.drainage import BoundaryDrainage
from tested_python_streets.model import Section, Intersection, End, Key
from tested_python_streets.backend import PythonMeshBackend
from test_python_street_partial_bundle_caps import bergamo_ground_ramps


def physical_stopping_curb():
    """Original airport controls: the actual187 curb crosses190's free plane."""
    network = bergamo_ground_ramps()
    sections = (
        Section(91, 91, ((-1640.246337890625, 817.523681640625, 0.),
            (-1649.06787109375, 812.3607788085938, 0.),
            (-1668.99853515625, 797.7609252929688, 0.),
            (-1681.733642578125, 787.68994140625, 0.),
            (-1687.280517578125, 783.3388671875, 0.),
            (-1696.84228515625, 773.0778198242188, 0.)), width=14., rise=.08),
        Section(94, 94, ((-1729.1622314453125, 847.32568359375, 0.),
            (-1727.39794921875, 842.7833251953125, 0.),
            (-1715.3770751953125, 829.3214721679688, 0.),
            (-1706.3104248046875, 821.05908203125, 0.),
            (-1700.0042724609375, 814.7901000976562, 0.),
            (-1696.529052734375, 809.6795654296875, 0.),
            (-1692.9923095703125, 801.5411376953125, 0.),
            (-1692.5430908203125, 794.505615234375, 0.),
            (-1692.7008056640625, 786.5463256835938, 0.),
            (-1693.152587890625, 784.3089599609375, 0.),
            (-1693.8768310546875, 780.7135009765625, 0.),
            (-1696.84228515625, 773.0778198242188, 0.)), width=5.400000095367432, rise=.08))
    site = Intersection(187, (End(91, 1), End(94, 1), End(93, 0)), 10.,
        (10., 2., 10.), ('STYLE', 'STYLE', 'STYLE'))
    return replace(network, sections=(*network.sections, *sections),
                   intersections=(*network.intersections, site))


class FiniteFreeClosureTests(unittest.TestCase):
    def test_unique_plane_crossing_uses_existing_inside_vertex_under_rigid_transforms(self):
        points = np.array(((3., 0., 0.), (1., .5, 0.), (.2, 1., 0.), (-.3, 1.5, 0.), (-1., 2., 0.)))
        for angle, translation in ((0., (0., 0.)), (.71, (2100., -1700.))):
            rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
            xyz = points.copy(); xyz[:, :2] = points[:, :2] @ rotation.T + translation
            node, proof = inside_vertex(xyz, np.asarray(translation), rotation[:, 0])
            self.assertEqual(node, 2)
            self.assertLess(proof['plane_offset'], proof['source_chord_length'])
            self.assertFalse(proof['physical_boundary_samples_added'])
            np.testing.assert_array_equal(xyz[node], xyz[2])

    def test_plane_requires_a_unique_transverse_crossing(self):
        for values in ((1., .5, .1), (1., 0., 1.), (1., -1., 1., -1.), (1., 0., 0., -1.)):
            points = np.c_[values, np.arange(len(values)), np.zeros(len(values))]
            with self.assertRaises(ValueError):
                inside_vertex(points, np.zeros(2), np.array((1., 0.)))
        for middle in (0., -1e-10):
            points = np.array(((1., 0., 0.), (middle, 1., 0.), (-1., 2., 0.)))
            node, proof = inside_vertex(points, np.zeros(2), np.array((1., 0.)))
            self.assertEqual(node, 1 if middle == 0. else 0)
            self.assertLessEqual(proof['source_chord_length'], np.sqrt(2.)+1e-8)

    def test_only_computational_connector_gets_one_collinear_parity_vertex(self):
        curve = np.array(((0., 0., .08), (1., 0., 0.)))
        curb = np.array(((1., 0., 0.), (2., 0., 0.), (3., 1., 0.)))
        free = np.array((3., 3., .08))
        connector = np.linspace(curb[-1], free, 5)
        for role in (0, 1):
            cap = FiniteFreeCap(curve, curb, connector, free, role)
            cap.continuation = dict(finite_free_closure={})
            refined = cap.refine_exterior_connector()
            self.assertEqual(len(refined.rim), len(cap.rim)+1)
            for point in cap.rim:
                self.assertTrue(np.any(np.all(refined.rim == point, axis=1)))
            np.testing.assert_array_equal(refined._curve, curve)
            np.testing.assert_array_equal(refined._curb, curb)
            delta = refined._connector-free
            direction = (free-curb[-1])[:2]
            self.assertLess(abs(delta[:, 0]*direction[1]-delta[:, 1]*direction[0]).max(), 1e-12)
            with self.assertRaisesRegex(ValueError, 'bounded to one'):
                refined.refine_exterior_connector()

    @staticmethod
    def height_seed():
        x, y = np.meshgrid(np.linspace(0., 4., 9), np.linspace(0., 3., 7))
        xyz = np.c_[x.ravel()+1700., y.ravel()-1300., np.zeros(x.size)]
        faces = []
        for row in range(6):
            for col in range(8):
                a = row*9+col;faces.append((a, a+1, a+10, a+9))
        mesh = Mesh(xyz, np.asarray(faces), 'grass')
        fixed = boundary_mask(mesh)
        # Low fixed source curbs surround a deliberately excessive seed crest.
        mesh.xyz[:, 2] = .01 * (mesh.xyz[:, 0]-1700.)
        mesh.xyz[~fixed, 2] = .08
        return mesh, fixed

    def test_physical_height_reconciliation_preserves_complete_rim_and_query(self):
        mesh, fixed = self.height_seed()
        self.assertGreater(mesh.quad_slopes().max(), .1)
        output, query, physical, report = reconcile(mesh, datum=0., ceiling=.08)
        np.testing.assert_array_equal(output.xyz[fixed], mesh.xyz[fixed])
        np.testing.assert_array_equal(output.xyz[:, :2], mesh.xyz[:, :2])
        np.testing.assert_array_equal(output.quads, mesh.quads)
        self.assertLessEqual(report['maximum_stored_physical_slope'], .1+1e-7)
        self.assertLessEqual(BoundaryDrainage().fill(physical, fixed)[1]['maximum_fill'], 1e-6)
        samples = output.xyz[output.quads].mean(axis=1)
        np.testing.assert_allclose(query.evaluate(samples[:, :2]), samples[:, 2], atol=1e-10)

    def test_zero_and_exaggerated_display_require_explicit_physical_buffer(self):
        mesh, fixed = self.height_seed()
        for exaggeration in (0., 10.):
            display = mesh.with_positions(mesh.xyz * np.array((1., 1., exaggeration)))
            with self.assertRaisesRegex(ValueError, 'physical height buffer'):
                reconcile(display, datum=0., ceiling=.08, exaggeration=exaggeration)
            output, query, physical, report = reconcile(display, datum=0., ceiling=.08,
                exaggeration=exaggeration, physical_mesh=mesh)
            np.testing.assert_array_equal(output.xyz[fixed], display.xyz[fixed])
            np.testing.assert_array_equal(output.xyz[:, :2], mesh.xyz[:, :2])
            self.assertTrue(report['explicit_physical_reconstruction'])
            self.assertLessEqual(report['maximum_stored_physical_slope'], .1+1e-7)
            if exaggeration == 0.:
                self.assertTrue(np.all(output.xyz[:, 2] == 0.))
                self.assertGreater(physical.xyz[:, 2].max(), 0.)

    def test_infeasible_original_boundary_is_rejected(self):
        mesh, fixed = self.height_seed()
        mesh.xyz[fixed, 2] = .5 * (mesh.xyz[fixed, 0]-1700.)
        with self.assertRaisesRegex(ValueError, 'fixed boundary'):
            reconcile(mesh, datum=0., ceiling=2.)

    def test_actual_stopping_junction_keeps_next_road_and_all_published_bank_edges(self):
        network = physical_stopping_curb()
        update = PythonMeshBackend(network.settings).prepare(network)
        part = update.products[Key('separator', 4)]
        records = [r for r in part.metadata['cap_continuations'] if 'finite_free_closure' in r]
        self.assertEqual(len(records), 1)
        proof = records[0]['finite_free_closure']
        self.assertEqual(proof['extra_physical_curb']['next_mouth'], [94, 1])
        self.assertFalse(proof['extra_physical_curb']['next_street_traversed'])
        self.assertFalse(proof['new_source_end'])
        self.assertLessEqual(proof['plane_offset'], proof['source_chord_length'])
        self.assertLessEqual(part.metadata['finite_free_reconciliation']['maximum_stored_physical_slope'], .1+1e-7)
        for sid in (91, 93, 94):
            self.assertGreater(len(update.products[Key('road', sid)].mesh.quads), 0)
            self.assertEqual(update.products[Key('profile', sid)].half, 0)
        # All grass edges on actual pavement are whole, existing source edges.
        def edges(mesh):
            ids = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)),
                                   axis=-1).reshape(-1, 2), axis=1)
            unique, counts = np.unique(ids, axis=0, return_counts=True)
            return {tuple(sorted((tuple(mesh.xyz[a]), tuple(mesh.xyz[b]))))
                    for a, b in unique[counts == 1]}
        grass = edges(part.mesh)
        road = edges(update.products[Key('road', 93)].mesh)
        self.assertGreater(len(grass & road), 5)
        from tested_python_streets.ownership import surface_key
        junction = update.products[surface_key(Key('intersection', 187))]
        self.assertGreater(len(grass & edges(junction.pavement.mesh)), 5)
        # The extra physical curb is represented by normal task inputs and by
        # the relevant port's surface digest; no undeclared cached road read.
        dependencies = update.dependencies[Key('separator', 4)]
        self.assertIn(Key('mouths', 187), dependencies)
        for sid in (91, 93, 94):
            self.assertIn(Key('path', sid), dependencies)
            self.assertIn(Key('road', sid), dependencies)
        from tested_python_streets.updates import BuildContext
        from tested_python_streets.patterns import PatternLibrary
        from tested_python_streets.separator_ports import PORTS, SeparatorPortIndex
        from tested_python_streets.bundles import BundleProvider
        context = BuildContext(network, update.products, PatternLibrary(),
            update.task_fingerprints, update.source_fingerprints)
        before = BundleProvider.signature(context, network.bundles[0], dependencies)
        digests = dict(update.task_fingerprints)
        digests[surface_key(Key('intersection', 187))] = 'changed published stopping curb'
        changed = BuildContext(network, dict(update.products), context.patterns,
            digests, update.source_fingerprints)
        changed.products[PORTS] = SeparatorPortIndex.build(changed)
        self.assertNotEqual(before, BundleProvider.signature(changed, network.bundles[0], dependencies))

    def test_actual_zero_and_exaggerated_caps_keep_the_original_physical_buffer(self):
        source = physical_stopping_curb()
        for exaggeration in (0., 10.):
            network = replace(source, settings=replace(source.settings, exaggeration=exaggeration))
            update = PythonMeshBackend(network.settings).prepare(network)
            part = update.products[Key('separator', 4)]
            plan = update.products[Key('bundle_plan', 4)]
            report = part.metadata['finite_free_reconciliation']
            self.assertTrue(report['explicit_physical_reconstruction'])
            self.assertLessEqual(report['maximum_stored_physical_slope'], .1+1e-7)
            self.assertIsNotNone(plan._finite_free_physical)
            if exaggeration == 0.:
                self.assertTrue(np.all(part.mesh.xyz[:, 2] == 0.))
                self.assertGreater(plan._finite_free_physical.xyz[:, 2].max(), 0.)
            points = part.mesh.xyz[part.mesh.quads].mean(axis=1)
            np.testing.assert_allclose(plan.height(points[:, :2]), points[:, 2], atol=1e-9)


if __name__ == '__main__':
    unittest.main()
