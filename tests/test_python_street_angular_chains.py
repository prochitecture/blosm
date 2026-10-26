"""Real angular Section routes: geometry, identities and Update dependencies."""
from dataclasses import replace
from functools import partial
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundles import _beziers
from tested_python_streets.bundle_angular_curves import polyline_distance
from _python_angular_fixtures import angular_chain_network as make_network
from test_python_street_angular_bundles import buffers


network = partial(make_network, model)


def island_boundary_error(plan,mesh):
    """Inspect actual mesh edges against the original physical island loop."""
    edges = np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
    unique,counts = np.unique(edges,axis=0,return_counts=True)
    adjacent = {}
    for a,b in unique[counts==1]:
        adjacent.setdefault(int(a),[]).append(int(b))
        adjacent.setdefault(int(b),[]).append(int(a))
    assert all(len(neighbors)==2 for neighbors in adjacent.values())
    first = min(adjacent)
    chain,previous,current = [first],None,first
    while True:
        following = next(v for v in adjacent[current] if v!=previous)
        if following==first:
            break
        chain.append(following)
        previous,current = current,following
        assert len(chain)<=len(adjacent)
    assert len(chain)==len(adjacent)
    actual,intended = mesh.xyz[chain,:2],plan.islands[0].xy
    return max(float(polyline_distance(points,np.vstack((curve,curve[:1]))).max())
               for points,curve in ((actual,intended),(intended,actual)))


class PythonAngularChainTests(unittest.TestCase):
    def qualify(self, source, backend=None):
        saved = source.to_dict()
        backend = backend or PythonMeshBackend(source.settings)
        update = backend.prepare(source)
        self.assertEqual(saved, source.to_dict())
        plan = update.products[model.Key('bundle_plan', 13)]
        self.assertEqual(type(plan).__name__, 'AngularBundlePlan')
        expected = tuple(t.section_id for m in source.bundles[0].members for t in m.path)
        self.assertEqual(plan.source_ids, expected)
        part = update.products[model.Key('bundle_pavement', 13)]
        self.assertEqual(part.metadata['source_sections'], list(expected))
        self.assertIn(part.metadata['islands'], (0, 1))
        self.assertLessEqual(part.metadata['quality']['maximum_edge_ratio'], 8.)
        self.assertLessEqual(part.metadata['angular_longitudinal_bound'], .03+1e-8)
        self.assertLessEqual(part.metadata['physical_slopes']['maximum'], plan.maximum_slope+1e-8)
        self.assertTrue(part.metadata['physical_boundary_contract']['degree_two'])
        self.assertLessEqual(len(part.mesh.quads), 100000)
        self.assertGreaterEqual(float(part.mesh.xyz[:,2].min()), plan.datum-1e-8)
        self.assertLessEqual(float(part.mesh.xyz[:,2].max()),
                             plan.datum+plan.rise*source.settings.exaggeration+1e-8)
        for role, route in enumerate(plan.paths):
            self.assertEqual([p.section.id for p in route.pieces],
                             [t.section_id for t in source.bundles[0].members[role].path])
            self.assertEqual(part.metadata['member_routes'][role], route.source_spans())
            for piece in route.pieces:
                self.assertIs(piece.path, update.products[model.Key('path', piece.section.id)])
                # Curve certification refers to original polynomials, not a
                # replacement cubic fitted through concatenated controls.
                self.assertTrue(_beziers(piece.path))
                owned = update.products[model.Key('road', piece.section.id)]
                self.assertEqual(owned.metadata['composite_owner'], 'bundle:13')
                self.assertFalse(hasattr(owned, 'mesh'))
        for item in source.furniture:
            if item.id == 999:
                continue
            placement = update.products[model.Key('placement', item.id)]
            self.assertAlmostEqual(placement.position[2],
                float(plan.height(np.asarray(placement.position[:2]))), places=10)
        return backend, update, plan

    def equal(self, first, second):
        self.assertEqual(first.keys(), second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0], second[key][0])
            np.testing.assert_array_equal(first[key][1], second[key][1])

    def test_unequal_section_partitions_both_modes_all_sweeps_and_islands(self):
        for sweep in (180., 225., 270.):
            for variable in (False, True):
                for mode in (0, 3):
                    with self.subTest(sweep=sweep, variable=variable, mode=mode):
                        source = network(sweep, variable, mode)
                        _, update, _ = self.qualify(source)
                        self.assertEqual(update.products[model.Key('bundle_pavement',13)].metadata['islands'],
                                         int(variable))
                        grass = update.products[model.Key('separator', 13)]
                        self.assertEqual(grass is not None, variable)
                        if grass is not None:
                            self.assertLessEqual(grass.mesh.validate()['maximum_edge_ratio'], 8.)

    def test_reversed_middle_sections_offsets_datum_and_rigid_transform(self):
        self.qualify(network(270., True, reversed_sections=(17, 117, 229),
            widths=(5., 7.), offsets=(.25, -.15), rises=(.04, .10), datum=37.,
            rotation=.713, translation=(1234., -5678.)))

    def test_float32_source_storage_preserves_collinear_exact_seams(self):
        for sweep in (180.,225.,270.):
            for variable in (False,True):
                with self.subTest(sweep=sweep,variable=variable):
                    source = network(sweep,variable,0,reversed_sections=(117,229))
                    source = replace(source,sections=tuple(replace(s,
                        points=tuple(map(tuple,np.asarray(s.points,dtype=np.float32).astype(float))))
                        for s in source.sections))
                    _,update,plan = self.qualify(source)
                    part = update.products[model.Key('bundle_pavement',13)]
                    errors = part.metadata['member_route_curb_discrepancy']
                    limits = part.metadata['member_route_curb_tolerance']
                    self.assertTrue(all(error<=limit<=1e-5 for error,limit in zip(errors,limits)))
                    self.assertTrue(any(error>1e-8 for error in errors))

    def test_blender_middle_control_edit_retains_island_boundary_accuracy(self):
        for sweep in (180.,225.,270.):
            for mode in (0,3):
                with self.subTest(sweep=sweep,mode=mode):
                    source = network(sweep,True,mode,reversed_sections=(117,229))
                    sections = []
                    for s in source.sections:
                        points = np.asarray(s.points,dtype=np.float32)
                        if s.id==229:
                            # Same interior control selected by the ordinary
                            # Blender gallery helper, including float32 storage.
                            at = 1+(len(points)-2)//2
                            points[at,0] = np.float32(float(points[at,0])+.1)
                        sections.append(replace(s,points=tuple(map(tuple,points.astype(float)))))
                    _,update,plan = self.qualify(replace(source,sections=tuple(sections)))
                    separator = update.products[model.Key('separator',13)]
                    chart = separator.metadata['chart']
                    self.assertLessEqual(chart['boundary_chord_error']+chart.get('source_boundary_error',0.),
                                         source.settings.curve_tolerance+1e-8)
                    self.assertLessEqual(island_boundary_error(plan,separator.mesh),
                                         source.settings.curve_tolerance+1e-8)

    def test_blender_width_edit_uses_valid_fixed_rim_quad_atlas(self):
        for sweep in (180.,225.,270.):
            for mode in (0,3):
                with self.subTest(sweep=sweep,mode=mode):
                    source = network(sweep,True,mode,reversed_sections=(117,229))
                    source = replace(source,sections=tuple(replace(s,
                        points=tuple(map(tuple,np.asarray(s.points,dtype=np.float32).astype(float))),
                        width=float(np.float32(6.1)) if s.street_id==41 else s.width)
                        for s in source.sections))
                    _,update,plan = self.qualify(source)
                    part = update.products[model.Key('separator',13)]
                    mesh = part.mesh
                    self.assertLessEqual(mesh.validate()['maximum_edge_ratio'],8.)
                    self.assertLessEqual(island_boundary_error(plan,mesh),source.settings.curve_tolerance+1e-8)
                    self.assertLessEqual(part.metadata['chart'].get('complete_curve_error_bound',0.),
                                         source.settings.curve_tolerance)

    def test_seam_and_middle_control_edits_equal_full_and_keep_remote_outputs(self):
        source = network()
        backend, initial, _ = self.qualify(source)
        backend.commit(initial)
        saved, remote = buffers(initial), initial.products[model.Key('road', 99)]
        middle = source.section_index[117]
        control = list(middle.points)
        control[1] = (control[1][0]+.1, control[1][1], 0.)
        interior = replace(source, sections=tuple(replace(s, points=tuple(control)) if s.id==117 else s
                                                  for s in source.sections))
        seam = source.connections[0]
        endpoints = {e.section_id:e.end for e in seam.ends}
        moved = []
        for section in source.sections:
            if section.id in endpoints:
                points = list(section.points)
                at = 0 if endpoints[section.id]==0 else -1
                x,y,z = points[at]
                points[at] = (x+.1,y,z)
                section = replace(section, points=tuple(points))
            moved.append(section)
        changes = (interior, replace(source, sections=tuple(moved)),
            replace(source, sections=tuple(replace(s, width=6.1) if s.street_id==41 else s
                                          for s in source.sections)))
        for changed in changes:
            with self.subTest(changed=changed.section_index[117]):
                _, local, _ = self.qualify(changed, backend)
                full = PythonMeshBackend(changed.settings).prepare(changed)
                self.equal(buffers(local), buffers(full))
                self.assertIn(model.Key('bundle_plan', 13), local.rebuilt)
                self.assertIn(model.Key('placement', 910), local.rebuilt)
                self.assertIs(local.products[model.Key('road', 99)], remote)
                self.assertNotIn(model.Key('placement', 999), local.rebuilt)
                backend.commit(local)
                restored = backend.prepare(source)
                self.equal(buffers(restored), saved)
                backend.commit(restored)

    def test_release_keeps_middle_section_queries_and_failed_edit_is_transactional(self):
        source = network(reversed_sections=(117, 229))
        backend, initial, plan = self.qualify(source)
        backend.commit(initial)
        backend.release_published_meshes()
        self.assertIsNone(initial.products[model.Key('bundle_pavement', 13)].mesh)
        changed = replace(source, furniture=tuple(replace(p, station=5.1) if p.id==910 else p
                                                 for p in source.furniture))
        local = backend.prepare(changed)
        full = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(local.products[model.Key('placement', 910)], full.products[model.Key('placement', 910)])
        self.assertNotIn(model.Key('bundle_plan', 13), local.rebuilt)
        old = backend.manifest()
        broken = replace(source, connections=source.connections[1:])
        with self.assertRaisesRegex(ValueError, 'explicit two-ended Section seams'):
            backend.prepare(broken)
        self.assertEqual(backend.manifest(), old)
        wrong = replace(source, sections=tuple(replace(s, datum=2.) if s.id==117 else s
                                               for s in source.sections))
        with self.assertRaises(ValueError):
            backend.prepare(wrong)
        self.assertEqual(backend.manifest(), old)

    def test_finer_resolution_and_zero_exaggeration_keep_source_contracts(self):
        for resolution, exaggeration in ((2,0.), (4,10.)):
            with self.subTest(resolution=resolution):
                self.qualify(network(180., True,
                    settings=model.Settings(resolution=resolution, exaggeration=exaggeration)))


if __name__ == '__main__':
    unittest.main(verbosity=2)
