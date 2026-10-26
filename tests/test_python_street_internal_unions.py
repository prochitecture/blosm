"""Original compound mouths connect overlapping curved Bundle runs."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import internal_angular_state_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from tested_python_streets.paths import Path
from tested_python_streets.regions import Region, RegionUnion, signed_area
from tested_python_streets.slope_projection import corner_matrices
from test_python_street_angular_bundles import buffers
from test_python_street_angular_island_chain import edges, components


def network(mode, *, datum=0., exaggeration=10.):
    source = internal_angular_state_network(model, mode, state='multiple')
    return replace(source, sections=tuple(replace(section,
        points=tuple(map(tuple, np.asarray(section.points, dtype=np.float32).astype(float))),
        rise=float(np.float32(section.rise)), datum=datum) for section in source.sections),
        furniture=(*source.furniture, model.Placement(904, 301, 5.)),
        settings=replace(source.settings, exaggeration=exaggeration))


def distance(query, segments):
    first = segments[:, 0]
    delta = segments[:, 1]-first
    length2 = np.sum(delta*delta, axis=1)
    result = []
    for start in range(0, len(query), 64):
        d = query[start:start+64, None]-first
        t = np.clip(np.sum(d*delta, axis=-1)/length2, 0., 1.)
        result.extend(np.sqrt(np.min(np.sum((d-t[..., None]*delta)**2, axis=-1), axis=1)))
    return np.asarray(result)


class InternalAngularUnionTests(unittest.TestCase):
    def qualify(self, source, result):
        plan = result.products[model.Key('bundle_plan', 13)]
        self.assertEqual(plan.mode, 'union')
        self.assertEqual(len(plan.children), 2)
        self.assertEqual(set(plan.source_ids), {17, 29, 117, 129})
        self.assertEqual({s.id for s in source.sections if s.street_id == 800}, {300, 301, 302})
        self.assertEqual(source.bundles[0].interior_section_ids, (301,))
        for sid in plan.source_ids:
            road = result.products[model.Key('road', sid)]
            self.assertIsNone(getattr(road, 'mesh', None))
            self.assertEqual(road.metadata['composite_owner'], 'bundle:13')
        self.assertIn(model.Key('placement', 904), result.products)
        self.assertEqual(len(result.products[model.Key('separator', 13)].metadata['runs']), 2)
        grass = result.products[model.Key('separator', 13)].mesh
        self.assertEqual(components(edges(grass)[0]), 4)
        paving, parts = set(), []
        for key, part in result.products.items():
            if getattr(part, 'mesh', None) is None:
                continue
            mesh = Mesh(part.mesh.xyz.astype(np.float32), part.mesh.quads)
            self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
            datum = source.sections[0].datum
            self.assertGreaterEqual(float(mesh.xyz[:, 2].min()), datum-1e-6)
            self.assertLessEqual(float(mesh.xyz[:, 2].max()),
                datum+source.settings.exaggeration*source.sections[0].rise+1e-6)
            if source.settings.exaggeration:
                physical = mesh.xyz.copy()
                physical[:, 2] = datum+(physical[:, 2]-datum)/source.settings.exaggeration
                gradient = np.einsum('fcij,fj->fci', corner_matrices(Mesh(physical, mesh.quads)),
                                     physical[mesh.quads, 2]-datum)
                self.assertLessEqual(float(np.linalg.norm(gradient, axis=2).max()), .1+1e-6)
            if key.kind != 'separator':
                paving.update(tuple(sorted((tuple(mesh.xyz[a]), tuple(mesh.xyz[b]))))
                              for a, b in edges(mesh)[1])
            if key != model.Key('road', 99):
                parts.append(mesh)
        xyz = grass.xyz.astype(np.float32)
        self.assertTrue(all(tuple(sorted((tuple(xyz[a]), tuple(xyz[b])))) in paving for a, b in edges(grass)[1]))
        vertex = np.vstack([part.xyz for part in parts])
        offset = np.cumsum([0, *[len(part.xyz) for part in parts]])
        vertex, inverse = np.unique(vertex, axis=0, return_inverse=True)
        faces = np.vstack([inverse[offset[i]:offset[i+1]][part.quads] for i, part in enumerate(parts)])
        whole = Mesh(vertex, faces)
        unique, boundary = edges(whole)
        self.assertEqual(components(unique), 1)
        self.assertEqual(len(vertex)-len(unique)+len(faces), 1)
        degree = np.bincount(boundary.ravel(), minlength=len(vertex))
        self.assertTrue(np.all(degree[np.unique(boundary)] == 2))
        replacement = result.products[model.Key('bundle_pavement', 13)].mesh
        self.assertLess(len(grass.quads)+len(replacement.quads), 8000)
        return plan

    def footprint(self, source, plan):
        for run, child in plan.children:
            regions = []
            for sid in child.source_ids:
                section = source.section_index[sid]
                path = Path(section, source.settings.curve_tolerance)
                first, last = child.ranges[sid]
                samples = path.samples(first, last, source.settings.spacing,
                                       tolerance=.1*source.settings.curve_tolerance)
                curves = [path.point(samples, section.offset+q*.5*section.width) for q in (-1., 1.)]
                regions.append(Region(model.Key('section', sid), np.vstack((curves[1], curves[0][::-1]))))
            loops = RegionUnion().build(regions)
            reference = np.concatenate([np.stack((loop.xy, np.roll(loop.xy, -1, axis=0)), axis=1) for loop in loops])
            published = child._pavement.mesh
            xy = published.xyz.astype(np.float32)[:, :2].astype(float)
            actual = xy[edges(published)[1]]
            source_queries = np.concatenate([(1.-t)*reference[:, 0]+t*reference[:, 1] for t in np.linspace(0., 1., 9)])
            actual_queries = np.concatenate([(1.-t)*actual[:, 0]+t*actual[:, 1] for t in np.linspace(0., 1., 9)])
            exact = all(source.section_index[sid].interpolation == 0 for sid in child.source_ids)
            tolerance = 2e-5+(0. if exact else source.settings.curve_tolerance)
            self.assertLessEqual(float(distance(actual_queries, reference).max()), tolerance)
            self.assertLessEqual(float(distance(source_queries, actual).max()), tolerance)
            expected = float(sum(signed_area(loop.xy) for loop in loops))
            cells = xy[published.quads]
            area = float(.5*np.sum(cells[:, :, 0]*np.roll(cells[:, :, 1], -1, axis=1)
                                  -cells[:, :, 1]*np.roll(cells[:, :, 0], -1, axis=1)))
            perimeter = float(np.linalg.norm(reference[:, 1]-reference[:, 0], axis=1).sum())
            self.assertLessEqual(abs(area-expected), perimeter*tolerance+len(loops)*np.pi*tolerance*tolerance)

    def test_exact_original_ribbon_union_and_published_interfaces(self):
        source = network(0)
        saved = source.to_dict()
        result = PythonMeshBackend(source.settings).prepare(source)
        self.footprint(source, self.qualify(source, result))
        self.assertEqual(source.to_dict(), saved)

    def test_smooth_original_ribbon_union_and_published_interfaces(self):
        source = network(3)
        saved = source.to_dict()
        result = PythonMeshBackend(source.settings).prepare(source)
        self.footprint(source, self.qualify(source, result))
        self.assertEqual(source.to_dict(), saved)

    def test_combined_edit_equals_full_restores_and_keeps_remote_and_all_anchors(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = network(mode)
                backend = PythonMeshBackend(source.settings)
                initial = backend.prepare(source)
                original = buffers(initial)
                backend.commit(initial)
                sections = []
                for section in source.sections:
                    if section.id == 29:
                        points = np.asarray(section.points).copy()
                        points[len(points)//2, 0] += .05
                        section = replace(section, points=tuple(map(tuple, points)))
                    if section.id == 17:
                        section = replace(section, width=section.width+.05)
                    sections.append(section)
                changed = replace(source, sections=tuple(sections))
                actual = backend.prepare(changed)
                self.qualify(changed, actual)
                expected = PythonMeshBackend(changed.settings).prepare(changed)
                for key, (xyz, quads) in buffers(actual).items():
                    np.testing.assert_array_equal(xyz, buffers(expected)[key][0])
                    np.testing.assert_array_equal(quads, buffers(expected)[key][1])
                self.assertIs(actual.products[model.Key('road', 99)], initial.products[model.Key('road', 99)])
                self.assertNotIn(model.Key('placement', 999), actual.rebuilt)
                for key in actual.products:
                    if key.kind == 'placement':
                        self.assertEqual(actual.products[key], expected.products[key])
                backend.commit(actual)
                restored = buffers(backend.prepare(source))
                for key, (xyz, quads) in original.items():
                    np.testing.assert_array_equal(xyz, restored[key][0])
                    np.testing.assert_array_equal(quads, restored[key][1])

    def test_zero_exaggeration_and_translated_datum_keep_complete_interfaces(self):
        source = network(3, datum=1.7, exaggeration=0.)
        result = PythonMeshBackend(source.settings).prepare(source)
        self.qualify(source, result)
        self.assertTrue(all(np.all(part.mesh.xyz[:, 2] == 1.7) for part in result.products.values()
                            if getattr(part, 'mesh', None) is not None))

    def test_reversed_source_storage_keeps_original_incidence_and_footprint(self):
        source = network(3)
        reversed_ids = {29, 117}
        sections = tuple(replace(section, points=section.points[::-1], offset=-section.offset)
                         if section.id in reversed_ids else section for section in source.sections)
        sites = tuple(replace(site, ends=tuple(model.End(end.section_id, 1-end.end)
            if end.section_id in reversed_ids else end for end in site.ends)) for site in source.intersections)
        bundle = replace(source.bundles[0], members=tuple(replace(member,
            path=tuple(model.Traversal(visit.section_id, not visit.forward)
                       if visit.section_id in reversed_ids else visit for visit in member.path))
            for member in source.bundles[0].members))
        source = replace(source, sections=sections, intersections=sites, bundles=(bundle,))
        result = PythonMeshBackend(source.settings).prepare(source)
        self.footprint(source, self.qualify(source, result))


if __name__ == '__main__':
    unittest.main()
