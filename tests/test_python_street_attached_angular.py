"""Actual rounded ports, physical profiles and transactional curved Bundle edits."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import attached_angular_network, attached_angular_state_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from tested_python_streets.bundle_angular_curves import polyline_distance
from test_python_street_angular_bundles import buffers
from test_python_street_angular_island_chain import edges, components


def source_network(mode=3, *, resolution=1, exaggeration=10., datum=0., state=None, radius=40.):
    source = (attached_angular_network(model, mode, radius=radius) if state is None
              else attached_angular_state_network(model, mode, state=state))
    return replace(source, sections=tuple(replace(section,
        points=tuple(map(tuple, np.asarray(section.points, dtype=np.float32).astype(float))),
        rise=float(np.float32(section.rise)), datum=datum) for section in source.sections),
        settings=replace(source.settings, resolution=resolution, exaggeration=exaggeration))


class AttachedAngularTests(unittest.TestCase):
    def qualify(self, source, update, *, islands=2):
        plan = update.products[model.Key('bundle_plan', 13)]
        self.assertEqual(type(plan).__name__, 'AttachedAngularPlan')
        self.assertGreater(plan.proof_checks, 0)
        pavement = update.products[model.Key('bundle_pavement', 13)]
        grass = update.products[model.Key('separator', 13)]
        self.assertLessEqual(pavement.metadata['physical_slopes']['maximum'], .1+1e-8)
        self.assertLessEqual(pavement.metadata['angular_longitudinal_bound'], .03+1e-8)
        self.assertEqual(grass.metadata['islands'], islands)
        self.assertEqual(grass.metadata['attached_caps'], [True, True])
        self.assertEqual(components(edges(pavement.mesh)[0]), 2 if islands == 1 else 1)
        self.assertEqual(components(edges(grass.mesh)[0]), islands)
        self.assertEqual({key.id for key in update.products if key.kind == 'junction_mesh'}, {200, 201, 210, 211})
        for sid in (17, 29):
            road = update.products[model.Key('road', sid)]
            self.assertFalse(hasattr(road, 'mesh'))
            self.assertEqual(road.metadata['composite_owner'], 'bundle:13')
        # Independent published-precision whole-edge audit across all real parts.
        available = set()
        for key, part in update.products.items():
            if getattr(part, 'mesh', None) is None:
                continue
            mesh = Mesh(part.mesh.xyz.astype(np.float32), part.mesh.quads)
            self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
            self.assertEqual(len(np.unique(mesh.quads)), len(mesh.xyz))
            self.assertGreaterEqual(mesh.xyz[:, 2].min(), plan.datum-1e-6)
            self.assertLessEqual(mesh.xyz[:, 2].max(), plan.datum+source.settings.exaggeration*plan.rise+1e-6)
            if key.kind == 'separator':
                continue
            for a, b in edges(mesh)[1]:
                available.add(tuple(sorted((tuple(mesh.xyz[a]), tuple(mesh.xyz[b])))))
        published_grass = grass.mesh.xyz.astype(np.float32)
        for a, b in edges(grass.mesh)[1]:
            self.assertIn(tuple(sorted((tuple(published_grass[a]), tuple(published_grass[b])))), available)
        # Each actual full mouth remains a complete shared edge partition.
        bundle = pavement.mesh.xyz.astype(np.float32)
        bundle_edges = {tuple(sorted((tuple(bundle[a]), tuple(bundle[b])))) for a, b in edges(pavement.mesh)[1]}
        for junction in source.intersections:
            mouth_end = next(end for end in junction.ends if end.section_id in (17, 29))
            layout = update.products[junction.layout_key]
            mouth = layout.mouths[mouth_end.section_id, mouth_end.end]
            path = update.products[model.Key('path', mouth_end.section_id)]
            profile = update.products[model.Key('profile', mouth_end.section_id)]
            q = np.linspace(-1., 1., mouth.transverse_intervals+1)
            xy = path.point(mouth.station, path.section.offset+.5*path.section.width*q)
            xyz = np.c_[xy, profile.datum+source.settings.exaggeration*profile.value(q)].astype(np.float32)
            for a, b in zip(xyz, xyz[1:]):
                self.assertIn(tuple(sorted((tuple(a), tuple(b)))), bundle_edges)
        if source.settings.resolution == 1:
            self.assertLess(len(pavement.mesh.quads)+len(grass.mesh.quads), 10000)
        identities = {item.id for item in source.furniture}
        if 901 in identities:
            self.assertIn(model.Key('separator', 13), update.dependencies[model.Key('placement', 901)])
        for identity in identities.intersection((901, 902)):
            placement = update.products[model.Key('placement', identity)]
            self.assertAlmostEqual(placement.position[2], float(plan.height(np.asarray(placement.position[:2]))), places=10)
        return plan

    def test_exact_original_rounded_ports_and_float32_interfaces(self):
        source = source_network(0)
        original = source.to_dict()
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))
        self.assertEqual(original, source.to_dict())

    def test_smooth_original_rounded_ports_and_float32_interfaces(self):
        source = source_network(3)
        original = source.to_dict()
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))
        self.assertEqual(original, source.to_dict())

    def test_vertex_and_width_updates_equal_full_restore_remote_and_anchors(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = source_network(mode)
                backend = PythonMeshBackend(source.settings)
                initial = backend.prepare(source)
                saved = buffers(initial)
                backend.commit(initial)
                for action in ('vertex', 'width'):
                    sections = []
                    for section in source.sections:
                        if action == 'vertex' and section.id == 29:
                            points = list(section.points)
                            x, y, z = points[8]
                            points[8] = (float(np.float32(x+.05)), y, z)
                            section = replace(section, points=tuple(points))
                        elif action == 'width' and section.id == 17:
                            section = replace(section, width=float(np.float32(section.width+.1)))
                        sections.append(section)
                    changed = replace(source, sections=tuple(sections))
                    local = backend.prepare(changed)
                    self.qualify(changed, local)
                    full = PythonMeshBackend(changed.settings).prepare(changed)
                    for key, (xyz, quads) in buffers(local).items():
                        np.testing.assert_array_equal(xyz, buffers(full)[key][0])
                        np.testing.assert_array_equal(quads, buffers(full)[key][1])
                    self.assertIs(local.products[model.Key('road', 99)], initial.products[model.Key('road', 99)])
                    self.assertNotIn(model.Key('placement', 999), local.rebuilt)
                    for key in local.products:
                        if key.kind == 'placement':
                            self.assertEqual(local.products[key], full.products[key])
                    backend.commit(local)
                    restored = backend.prepare(source)
                    for key, (xyz, quads) in saved.items():
                        np.testing.assert_array_equal(xyz, buffers(restored)[key][0])
                        np.testing.assert_array_equal(quads, buffers(restored)[key][1])
                    backend.commit(restored)

    def test_reversed_original_storage_keeps_directed_ports(self):
        source = source_network(3)
        reverse = {17, 29}
        changed = replace(source,
            sections=tuple(replace(s, points=s.points[::-1], offset=-s.offset) if s.id in reverse else s for s in source.sections),
            intersections=tuple(replace(j, ends=tuple(model.End(e.section_id, 1-e.end) if e.section_id in reverse else e for e in j.ends))
                                for j in source.intersections),
            bundles=(replace(source.bundles[0], members=tuple(model.BundleMember(tuple(model.Traversal(v.section_id, not v.forward)
                for v in member.path)) for member in source.bundles[0].members)),), furniture=())
        # Geometry incidence is tested independently; reversed stations would
        # need corresponding source-relative placement station conversion.
        result = PythonMeshBackend(changed.settings).prepare(changed)
        self.qualify(changed, result)

    def test_zero_exaggeration_retains_physical_caps_on_translated_datum(self):
        source = source_network(3, exaggeration=0., datum=1.7)
        result = PythonMeshBackend(source.settings).prepare(source)
        plan = self.qualify(source, result)
        self.assertGreater(plan.grass_query.heights.max(), plan.datum+.07)
        self.assertTrue(all(np.all(part.mesh.xyz[:, 2] == plan.datum) for part in result.products.values()
                            if getattr(part, 'mesh', None) is not None))

    def test_resolution_two_retains_complete_native_partitions(self):
        source = source_network(3, resolution=2)
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))

    def test_tighter_entrances_retain_whole_source_curbs_and_profiles(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = source_network(mode, radius=20.)
                update = PythonMeshBackend(source.settings).prepare(source)
                plan = self.qualify(source, update)
                mesh = update.products[model.Key('bundle_pavement', 13)].mesh
                points = mesh.xyz[edges(mesh)[1], :2]
                queries = np.concatenate([(1.-t)*points[:, 0]+t*points[:, 1] for t in np.linspace(0., 1., 9)])
                references = [curve.points for curve in plan.field.curves]
                for junction in source.intersections:
                    endpoint = next(end for end in junction.ends if end.section_id in (17, 29))
                    mouth = update.products[junction.layout_key].mouths[endpoint.section_id, endpoint.end]
                    path = update.products[model.Key('path', endpoint.section_id)]
                    references.append(path.point(mouth.station, path.section.offset+.5*path.section.width*np.array((-1., 1.))))
                error = np.min(np.stack([polyline_distance(queries, reference) for reference in references]), axis=0)
                self.assertLessEqual(float(error.max()), 1e-7 if mode == 0 else source.settings.curve_tolerance)

    def test_continuous_separator_in_both_modes_retains_complete_interfaces(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = source_network(mode, state='separated')
                self.qualify(source, PythonMeshBackend(source.settings).prepare(source), islands=1)

    def test_one_three_two_island_updates_equal_full_restore_remote_and_anchors(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = source_network(mode, state='separated')
                backend = PythonMeshBackend(source.settings)
                initial = backend.prepare(source)
                original = buffers(initial)
                remote = initial.products[model.Key('road', 99)]
                backend.commit(initial)
                for state, islands in (('multiple', 3), ('overlap', 2), ('separated', 1)):
                    changed = replace(source_network(mode, state=state), uuid=source.uuid)
                    local = backend.prepare(changed)
                    self.qualify(changed, local, islands=islands)
                    full = PythonMeshBackend(changed.settings).prepare(changed)
                    reference = buffers(full)
                    for key, (xyz, quads) in buffers(local).items():
                        np.testing.assert_array_equal(xyz, reference[key][0])
                        np.testing.assert_array_equal(quads, reference[key][1])
                    self.assertIs(local.products[model.Key('road', 99)], remote)
                    self.assertNotIn(model.Key('placement', 999), local.rebuilt)
                    for key in local.products:
                        if key.kind == 'placement':
                            self.assertEqual(local.products[key], full.products[key])
                    backend.commit(local)
                restored = buffers(local)
                self.assertEqual(original.keys(), restored.keys())
                for key, (xyz, quads) in original.items():
                    np.testing.assert_array_equal(xyz, restored[key][0])
                    np.testing.assert_array_equal(quads, restored[key][1])
