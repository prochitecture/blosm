"""Curved changing-width source routes, ownership and attachment dependencies.

These tests cover complete-overlap surfaces. The separate twelve-case diagnostic
retains the currently failing island cases; they are not counted as supported.
"""
from dataclasses import replace
from functools import partial
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import angular_dimensional_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_boundaries import dimensional_sites
from tested_python_streets.bundle_angular_curves import polyline_distance
from tested_python_streets.transitions import output_key
from test_python_street_angular_bundles import buffers

network = partial(angular_dimensional_network,model)


class PythonAngularDimensionTests(unittest.TestCase):
    def qualify(self,source,backend=None):
        before = source.to_dict()
        backend = backend or PythonMeshBackend(source.settings)
        update = backend.prepare(source)
        self.assertEqual(before,source.to_dict())
        plan = update.products[model.Key('bundle_plan',13)]
        self.assertEqual(type(plan).__name__,'AngularBundlePlan')
        sites = dimensional_sites(source,source.bundles[0])
        self.assertEqual(set(plan.transition_owners),{site.key for site in sites})
        self.assertTrue(sites)
        for site in sites:
            self.assertIsNone(update.products[output_key(site.key)])
            self.assertIn(site.layout_key,update.dependencies[model.Key('bundle_plan',13)])
        part = update.products[model.Key('bundle_pavement',13)]
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'],8.)
        self.assertLessEqual(part.metadata['physical_slopes']['maximum'],plan.maximum_slope+1e-8)
        self.assertLessEqual(part.metadata['angular_longitudinal_bound'],.03+1e-8)
        self.assertTrue(all(item['lower']>1e-8 for item in part.metadata['angular_transition_certificates']))
        self.assertEqual(set(part.metadata['source_sections']),set(source.bundles[0].sections))
        for route in plan.paths:
            for piece in route.pieces:
                self.assertIs(piece.path,update.products[model.Key('path',piece.section.id)])
                self.assertFalse(hasattr(update.products[model.Key('road',piece.section.id)],'mesh'))
        return backend,update,plan

    def equal(self,first,second):
        self.assertEqual(first.keys(),second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0],second[key][0])
            np.testing.assert_array_equal(first[key][1],second[key][1])

    def test_constant_overlap_both_modes_all_three_curved_sweeps(self):
        for sweep in (180.,225.,270.):
            for mode in (0,3):
                with self.subTest(sweep=sweep,mode=mode):
                    _,update,plan = self.qualify(network(sweep,False,mode))
                    self.assertIsNone(update.products[model.Key('separator',13)])
                    self.assertEqual(len(plan.transition_owners),4)
                    self.assertLess(len(update.products[model.Key('bundle_pavement',13)].mesh.quads),10000)

    def test_analytical_quintic_curbs_and_traversal_have_independent_dense_checks(self):
        source = network(225.,False,3)
        _,update,plan = self.qualify(source)
        for role,member in enumerate(source.bundles[0].members):
            for previous,following in zip(member.path,member.path[1:]):
                site = next((j for j in source.connections if
                    {e.section_id for e in j.ends}=={previous.section_id,following.section_id}),None)
                transition = update.products[site.layout_key].transition
                if transition is None:
                    continue
                direction = 1 if transition.source_index(previous.section_id)==0 else -1
                for curb in transition.curbs:
                    points = curb(np.linspace(0.,1.,2001))
                    angles = np.unwrap(np.arctan2(points[:,1]-plan.frame.center[1],
                                                  points[:,0]-plan.frame.center[0]))
                    self.assertTrue(np.all(np.diff(angles)*direction*plan.frame.direction>0.))
                    error = min(float(polyline_distance(points,curve.points).max())
                                for curve in plan.field.curves[2*role:2*role+2])
                    self.assertLessEqual(error,.11*source.settings.curve_tolerance)

    def test_float32_source_coordinates_keep_dimensional_seams_connected(self):
        for sweep in (180.,225.,270.):
            for mode in (0,3):
                with self.subTest(sweep=sweep,mode=mode):
                    source = network(sweep,False,mode)
                    source = replace(source,sections=tuple(replace(s,
                        points=tuple(map(tuple,np.asarray(s.points,dtype=np.float32).astype(float))))
                        for s in source.sections))
                    _,update,plan = self.qualify(source)
                    metadata = update.products[model.Key('bundle_pavement',13)].metadata
                    self.assertTrue(all(error<=limit<=1e-5 for error,limit in zip(
                        metadata['member_route_curb_discrepancy'],metadata['member_route_curb_tolerance'])))

    def test_local_width_and_vertex_changes_equal_full_restore_and_keep_remote(self):
        source = network()
        backend,initial,_ = self.qualify(source)
        backend.commit(initial)
        saved = buffers(initial)
        remote = initial.products[model.Key('road',99)]
        wider = replace(source,sections=tuple(replace(s,width=7.1) if s.id==117 else s for s in source.sections))
        middle = source.section_index[229]
        points = list(middle.points)
        x,y,z = points[1]
        points[1] = (x+.05,y,z)
        moved = replace(source,sections=tuple(replace(s,points=tuple(points)) if s.id==229 else s
                                             for s in source.sections))
        uniform = replace(source,sections=tuple(replace(s,width=6.) if s.id!=99 else s for s in source.sections))
        for changed in (wider,moved,uniform):
            local = backend.prepare(changed)
            full = PythonMeshBackend(changed.settings).prepare(changed)
            self.equal(buffers(local),buffers(full))
            self.assertIs(local.products[model.Key('road',99)],remote)
            self.assertNotIn(model.Key('placement',999),local.rebuilt)
            for item in changed.furniture:
                self.assertEqual(local.products[model.Key('placement',item.id)],
                                 full.products[model.Key('placement',item.id)])
            backend.commit(local)
            restored = backend.prepare(source)
            self.equal(buffers(restored),saved)
            backend.commit(restored)

    def test_transition_anchor_queries_survive_released_buffers(self):
        source = network()
        source = replace(source,furniture=(*source.furniture,model.Placement(930,117,.5),model.Placement(931,229,.5)))
        backend,initial,plan = self.qualify(source)
        backend.commit(initial)
        for placement in source.furniture:
            if placement.id==999:
                continue
            result = initial.products[model.Key('placement',placement.id)]
            self.assertAlmostEqual(result.position[2],float(plan.height(np.asarray(result.position[:2]))),places=10)
        backend.release_published_meshes()
        changed = replace(source,furniture=tuple(replace(p,station=.6) if p.id==930 else p for p in source.furniture))
        local = backend.prepare(changed)
        full = PythonMeshBackend(source.settings).prepare(changed)
        self.assertNotIn(model.Key('bundle_plan',13),local.rebuilt)
        self.assertEqual(local.products[model.Key('placement',930)],full.products[model.Key('placement',930)])

    def test_invalid_incidence_preserves_the_committed_build(self):
        source = network()
        backend,initial,_ = self.qualify(source)
        backend.commit(initial)
        manifest = backend.manifest()
        broken = replace(source,connections=source.connections[1:])
        with self.assertRaises(ValueError):
            backend.prepare(broken)
        self.assertEqual(backend.manifest(),manifest)


if __name__=='__main__':
    unittest.main()
