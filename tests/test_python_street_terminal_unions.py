"""Original nonparallel terminal cuts join a merged Bundle through direct quads."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import nonparallel_terminal_fixture, one_sided_free_cap_fixture
from test_python_street_dimensional_bundles import meshes, reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_terminal_caps import SkewMergedCapPlan, bounded_terminal_front, graded_segment


class PythonTerminalUnionTests(unittest.TestCase):
    def equal(self, first, second):
        self.assertEqual(first.keys(), second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0], second[key][0])
            np.testing.assert_array_equal(first[key][1], second[key][1])

    def qualify(self, network):
        authored = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(), authored)
        plan = result.products[model.Key('bundle_plan', 7)]
        self.assertIsInstance(plan, SkewMergedCapPlan)
        part = result.products[model.Key('bundle_pavement', 7)]
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
        if network.settings.resolution == 1:
            self.assertLess(len(part.mesh.quads), 12000)
        self.assertGreaterEqual(float(part.mesh.xyz[:, 2].min()), plan.datum-1e-8)
        self.assertLessEqual(float(part.mesh.xyz[:, 2].max()),
                             plan.datum+plan.rise*network.settings.exaggeration+1e-8)
        edges = np.sort(np.stack((part.mesh.quads, np.roll(part.mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        self.assertLessEqual(int(counts.max()), 2)
        self.assertEqual(len(part.mesh.xyz)-len(unique)+len(part.mesh.quads), 1)
        boundary = unique[counts == 1]
        for end, profile in plan.union_profiles.items():
            plane = (plan.first, plan.last)[end]
            x = (part.mesh.xyz[:, :2]-plan.origin)@plan.axis
            self.assertFalse(np.any(np.all(abs(x[boundary]-plane) < 1e-8, axis=1)))
            cap = plan.cap_meshes[-1, end]
            for vertex in profile:
                self.assertLess(float(np.linalg.norm(cap.xyz-vertex, axis=1).min()), 1e-8)
            self.assertLessEqual(plan.union_cap_metadata[end]['maximum_body_join_longitudinal_slope'], .03+1e-8)
        self.assertEqual(result.products[model.Key('separator', 7)].metadata['islands'], 1)
        return result

    def test_captured_move_and_its_opposite_keep_complete_boundary_and_profile(self):
        for mode in (0, 3):
            for move in (.05, -.05):
                with self.subTest(mode=mode, move=move):
                    network = nonparallel_terminal_fixture(model, move=move)
                    network = replace(network, sections=tuple(replace(s, interpolation=mode) for s in network.sections))
                    result = self.qualify(network)
                    plan = result.products[model.Key('bundle_plan', 7)]
                    cap = plan.cap_meshes[-1, 1]
                    # Both original outer terminal endpoints remain on the
                    # exterior of the merged cap, despite their different cuts.
                    for sid, side in ((2, -1.), (4, 1.)):
                        path = result.products[model.Key('path', sid)]
                        point = path.point(path.length, path.section.offset+side*.5*path.section.width)
                        at = int(np.argmin(np.linalg.norm(cap.xyz[:, :2]-point, axis=1)))
                        self.assertLess(float(np.linalg.norm(cap.xyz[at, :2]-point)), 1e-8)
                        self.assertAlmostEqual(float(cap.xyz[at, 2]), plan.datum, places=12)

    def test_both_authoring_directions_and_common_traversal_reversal(self):
        network = nonparallel_terminal_fixture(model)
        self.qualify(reverse_sections(network, {1, 4}))
        bundle = network.bundles[0]
        members = tuple(replace(m, path=tuple(replace(t, forward=not t.forward) for t in m.path[::-1]))
                        for m in bundle.members[::-1])
        self.qualify(replace(network, bundles=(replace(bundle, members=members),)))

    def test_small_and_intermediate_moves_preserve_quality_and_whole_body_edges(self):
        for mode in (0, 3):
            for move in (.001, -.001, -.01, .02, -.02):
                with self.subTest(mode=mode, move=move):
                    network = nonparallel_terminal_fixture(model, move=move)
                    network = replace(network, sections=tuple(replace(s, interpolation=mode) for s in network.sections))
                    self.qualify(network)

    def test_precision_front_preserves_corners_and_bounds_continuous_sampling_error(self):
        def samples(points):
            return np.concatenate([np.linspace(a,b,301) for a,b in zip(points,points[1:])])

        def distance(points, chain):
            a, edge = chain[:-1], np.diff(chain,axis=0)
            delta = points[:,None]-a
            t = np.clip(np.sum(delta*edge,axis=-1)/np.sum(edge*edge,axis=-1),0.,1.)
            return np.linalg.norm(delta-t[...,None]*edge,axis=-1).min(axis=1).max()

        angle=.73
        rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        for step in (.0001,.001,.0039,.015):
            with self.subTest(step=step):
                original=np.array(((0.,-9.),(0.,1.),(-step,1.),(step,9.)))@rotation.T+(123.,-456.)
                actual,meta=bounded_terminal_front(original,.005)
                np.testing.assert_array_equal(actual[[0,-1]],original[[0,-1]])
                error=max(distance(samples(original),actual),distance(samples(actual),original))
                self.assertLessEqual(float(error),meta['maximum_front_chord_error_bound']+1e-12)
                self.assertLessEqual(meta['maximum_front_chord_error_bound']+meta['source_curb_error_budget'],.005)
                if step==.015:
                    np.testing.assert_array_equal(actual,original)

    def test_local_boundary_grading_is_metric_bounded_and_allocation_limited(self):
        for length, first, last in ((10.,.08,.7),(8.,.7,.08),(.005,.7,.7),(10.,.01,.01)):
            t=graded_segment(length,.7,first,last)
            cells=length*np.diff(t)
            self.assertEqual(t[0],0.);self.assertEqual(t[-1],1.)
            self.assertTrue(np.all(cells>0.))
            self.assertLessEqual(float(cells.max()),.7+1e-12)
            self.assertLessEqual(float(cells[0]),first+1e-12)
            self.assertLessEqual(float(cells[-1]),last+1e-12)
            self.assertLessEqual(float(np.maximum(cells[1:]/cells[:-1],cells[:-1]/cells[1:]).max(initial=1.)),2.+1e-10)
        with self.assertRaises(ValueError):
            graded_segment(10000.,.7,.7,.7)
        with self.assertRaises(ValueError):
            graded_segment(10.,.7,.08,.7,extra=1000000)

    def test_finer_meshes_and_flat_display_keep_physical_join_checks(self):
        network = nonparallel_terminal_fixture(model)
        for resolution in (2, 3, 4):
            with self.subTest(resolution=resolution):
                self.qualify(replace(network, sections=tuple(replace(s, datum=37.) for s in network.sections),
                    settings=replace(network.settings, resolution=resolution,
                                     exaggeration=0. if resolution==2 else 10.)))

    def test_updates_restore_exactly_and_release_retains_cap_queries(self):
        initial = one_sided_free_cap_fixture(model)
        changed = nonparallel_terminal_fixture(model)
        changed = replace(changed, uuid=initial.uuid, furniture=initial.furniture)
        backend = PythonMeshBackend(initial.settings)
        base = backend.prepare(initial);backend.commit(base)
        saved = meshes(base)
        remote = base.products[model.Key('road', 99)]
        local = backend.prepare(changed)
        self.equal(meshes(local), meshes(PythonMeshBackend(changed.settings).prepare(changed)))
        self.assertIs(local.products[model.Key('road', 99)], remote)
        backend.commit(local)
        plan = local.products[model.Key('bundle_plan', 7)]
        query = plan.cap_queries[-1]
        position = query.xy[0]+query.origin
        height = float(plan.height(position))
        backend.release_published_meshes()
        self.assertFalse(plan.cap_meshes)
        self.assertAlmostEqual(float(plan.height(position)), height, places=12)
        self.equal(meshes(backend.prepare(changed, force=True)), meshes(PythonMeshBackend(changed.settings).prepare(changed)))
        restored = backend.prepare(initial)
        self.assertNotIn(model.Key('road', 99), restored.rebuilt)
        # Published unchanged buffers live in Blender, rather than a second
        # numerical cache. Combine retained published data with replacements.
        self.equal({**saved, **meshes(restored)}, saved)

    def test_rigid_motion_preserves_the_same_construction(self):
        network = nonparallel_terminal_fixture(model)
        angle = .713
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        sources = []
        for source in network.sections:
            p = np.asarray(source.points).copy()
            p[:, :2] = p[:, :2]@rotation.T + (1234., -5678.)
            sources.append(replace(source, points=tuple(map(tuple, p))))
        self.qualify(replace(network, sections=tuple(sources)))


if __name__ == '__main__':
    unittest.main(verbosity=2)
