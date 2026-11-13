"""Real-extract and translated-terrain qualification for native building import."""
from pathlib import Path
import importlib
import json
import math
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
from release_support import arguments, enable, catalog, configure_osm, run, PACKAGE

args = arguments()
enable(args)
runtime = importlib.import_module(PACKAGE + '.app.osm_import')
gn_runtime = importlib.import_module(PACKAGE + '.building.padl_integration.gn_runtime')
get_gn_input = importlib.import_module(PACKAGE + '.util.geometry_nodes').getGnInput
alive = importlib.import_module(PACKAGE + '.progressive.ownership').alive
metrics = []


def evaluated_summary(results):
    """Bound memory by footprint count, retaining evaluated topology/shape checks.

    Error-only groups can legitimately produce empty meshes. Per-footprint face
    and corner counts, bounds, corner centroid and area still detect missing,
    duplicated or differently placed geometry without retaining every face.
    """
    bpy.context.view_layer.update()
    graph = bpy.context.evaluated_depsgraph_get()
    summary = {}
    for result in results:
        obj = result.output.evaluated_get(graph)
        mesh = obj.to_mesh()
        try:
            if not mesh.polygons:
                continue
            ids = mesh.attributes.get('footprint_id')
            assert ids is not None, 'Generated mesh lost footprint identities'
            for face in mesh.polygons:
                identity = ids.data[face.index].value
                row = summary.setdefault(identity, dict(faces=0, corners=0, area=0.,
                    low=[float('inf')] * 3, high=[-float('inf')] * 3, total=[0.] * 3))
                row['faces'] += 1
                row['corners'] += len(face.vertices)
                row['area'] += face.area
                for index in face.vertices:
                    point = obj.matrix_world @ mesh.vertices[index].co
                    assert all(math.isfinite(value) for value in point)
                    for axis in range(3):
                        row['low'][axis] = min(row['low'][axis], point[axis])
                        row['high'][axis] = max(row['high'][axis], point[axis])
                        row['total'][axis] += point[axis]
        finally:
            obj.to_mesh_clear()
    for row in summary.values():
        row['center'] = [value / row['corners'] for value in row.pop('total')]
    return summary


def terrain_state(terrain):
    return (terrain.data.as_pointer(), tuple(tuple(row) for row in terrain.matrix_world),
            tuple(tuple(vertex.co) for vertex in terrain.data.vertices),
            tuple(tuple(face.vertices) for face in terrain.data.polygons),
            tuple((modifier.name, modifier.type) for modifier in terrain.modifiers),
            tuple(material.as_pointer() for material in terrain.data.materials))


class PadlDatasetChecks(unittest.TestCase):
    def prepare(self, backend, fixture, progressive):
        entry = next(item for item in catalog() if item['id'] == fixture)
        settings = configure_osm(args, entry, '3Dsimple', False)
        settings.buildingBackend = backend
        settings.buildingSeed = 113
        settings.progressiveImport = progressive
        # Check that progressive grouping ignores the saved blocking preference.
        settings.singleObject = progressive
        return settings

    def session(self):
        session = runtime.OsmImport(bpy.context)

        def cleanup():
            if not session.finished:
                session.cancel()
                session.finish(silent=True)

        self.addCleanup(cleanup)
        return session

    def record(self, backend, fixture, mode):
        metrics.append(dict(backend=backend, fixture=fixture, mode=mode,
                            groups=len(runtime.app.padlBuildingResults),
                            **json.loads(json.dumps(runtime.last_result))))

    def assert_no_readers(self):
        from padl_status_support import assert_status_resources
        assert_status_resources(runtime)

    def assert_summary_equal(self, expected, actual):
        self.assertTrue(expected, 'Dataset generated no geometry')
        self.assertEqual(set(expected), set(actual), 'Generated footprint identities changed')
        for identity, before in expected.items():
            after = actual[identity]
            self.assertEqual(before['faces'], after['faces'], identity)
            self.assertEqual(before['corners'], after['corners'], identity)
            for field in ('low', 'high', 'center'):
                for left, right in zip(before[field], after[field]):
                    self.assertAlmostEqual(left, right, delta=.01, msg=(identity, field))
            self.assertAlmostEqual(before['area'], after['area'],
                                   delta=max(.001, abs(before['area']) * .0001), msg=identity)

    def report_snapshot(self):
        results = runtime.app.padlBuildingResults
        report = runtime.app.padlBuildingReport
        self.assertEqual(report['statuses'], [row for result in results for row in result.statuses])
        for result in results:
            text = result.owner.get('blosm_building_report')
            self.assertIsInstance(text, bpy.types.Text)
            self.assertEqual(json.loads(text.as_string()), report)
        return json.loads(json.dumps(report))

    def drain(self, session, on_commit=None):
        count = session.committed_units
        for _ in range(200000):
            if session.job.done:
                break
            session.job.step()
            if session.committed_units != count:
                self.assertEqual(session.committed_units, count + 1)
                count = session.committed_units
                self.assert_no_readers()
                if on_commit:
                    on_commit()
        self.assertTrue(session.job.done)
        self.assertIsNone(session.job.error, str(session.job.error))
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'success', result)
        self.assertFalse(result['errors'], result)
        self.assertEqual(result['committed'], len(runtime.app.padlBuildingResults))
        self.assert_no_readers()

    def terrain(self, settings):
        extent = 1000.
        corners = [(-extent, -extent), (extent, -extent),
                   (extent, extent), (-extent, extent)]
        mesh = bpy.data.meshes.new('Translated sloped user terrain')
        mesh.from_pydata([(x, y, 20 + .05*x + .03*y) for x, y in corners], [], [(0, 1, 2, 3)])
        terrain = bpy.data.objects.new('Translated sloped user terrain', mesh)
        bpy.context.scene.collection.objects.link(terrain)
        terrain.location = (30, -20, 15)
        bpy.context.view_layer.update()
        settings.terrainObject = terrain.name
        return terrain

    def assert_terrain_helper(self, backend, terrain):
        helpers = [obj for obj in bpy.data.objects if obj.type == 'MESH'
                   and len(obj.modifiers) == 1 and obj.modifiers[0].type == 'NODES'
                   and obj.modifiers[0].node_group
                   and obj.modifiers[0].node_group.get('prochitecture_terrain_z_limits_v1')]
        self.assertEqual(len(helpers), 1 if backend == 'geometry_nodes' else 0)
        if helpers:
            helper = helpers[0]
            self.assertNotIn('_blosm_pending_import', helper)
            modifier = helper.modifiers[0]
            self.assertEqual(get_gn_input(modifier, gn_runtime.socket_id(modifier.node_group, 'terrain')), terrain)
            for result in runtime.app.padlBuildingResults:
                modifier = result.owner.modifiers['PADL Buildings']
                self.assertEqual(get_gn_input(modifier, gn_runtime.socket_id(modifier.node_group, 'Terrain limits')), helper)
            return helper
        return None

    def assert_projected(self, summary):
        # The translated terrain is around z=34 here; unprojected inputs are z=0.
        self.assertEqual(len(summary), 3)
        self.assertTrue(all(row['low'][2] > 25 for row in summary.values()), summary)
        self.assertGreater(max(row['low'][2] for row in summary.values()) -
                           min(row['low'][2] for row in summary.values()), .1)

    def test_karl_marx_allee_full_native_geometry_and_status_parity(self):
        for backend in ('geometry_nodes', 'python'):
            with self.subTest(backend=backend):
                self.prepare(backend, 'karl_marx_allee', False)
                self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
                self.assertEqual(runtime.last_result['outcome'], 'success')
                self.record(backend, 'karl_marx_allee', 'blocking_separate')
                before_report = self.report_snapshot()
                before = evaluated_summary(runtime.app.padlBuildingResults)
                group_count = len(runtime.app.padlBuildingResults)
                self.assertGreater(group_count, 3)
                self.prepare(backend, 'karl_marx_allee', True)
                session = self.session()
                self.drain(session)
                self.record(backend, 'karl_marx_allee', 'progressive')
                self.assertEqual(len(runtime.app.padlBuildingResults), group_count)
                # Legitimate unsupported source/roof diagnostics are part of
                # the baseline contract, not a reason to skip the dataset.
                self.assertEqual(self.report_snapshot(), before_report)
                self.assert_summary_equal(before, evaluated_summary(runtime.app.padlBuildingResults))

    def test_translated_sloped_terrain_parity_and_shared_helpers(self):
        for backend in ('geometry_nodes', 'python'):
            with self.subTest(backend=backend):
                settings = self.prepare(backend, 'shared_edges', False)
                terrain = self.terrain(settings)
                source = terrain_state(terrain)
                self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
                self.record(backend, 'shared_edges_terrain', 'blocking_separate')
                self.assertEqual(terrain_state(terrain), source)
                self.assert_terrain_helper(backend, terrain)
                before_report = self.report_snapshot()
                before = evaluated_summary(runtime.app.padlBuildingResults)
                self.assert_projected(before)
                settings = self.prepare(backend, 'shared_edges', True)
                terrain = self.terrain(settings)
                source = terrain_state(terrain)
                helper_ids = []

                def verify_commit():
                    self.assertEqual(terrain_state(terrain), source)
                    helper = self.assert_terrain_helper(backend, terrain)
                    if helper:
                        helper_ids.append(helper.as_pointer())

                self.drain(self.session(), verify_commit)
                self.record(backend, 'shared_edges_terrain', 'progressive')
                self.assertEqual(terrain_state(terrain), source)
                if backend == 'geometry_nodes':
                    self.assertEqual(len(helper_ids), 3)
                    self.assertEqual(len(set(helper_ids)), 1, 'Terrain helper duplicated between groups')
                self.assertEqual(self.report_snapshot(), before_report)
                after = evaluated_summary(runtime.app.padlBuildingResults)
                self.assert_projected(after)
                self.assert_summary_equal(before, after)

    def test_cancel_preserves_terrain_and_accepted_helper(self):
        for backend in ('geometry_nodes', 'python'):
            with self.subTest(backend=backend):
                settings = self.prepare(backend, 'shared_edges', True)
                terrain = self.terrain(settings)
                source = terrain_state(terrain)
                session = self.session()
                for _ in range(20000):
                    if session.committed_units or session.job.done:
                        break
                    session.job.step()
                self.assertIsNone(session.job.error, str(session.job.error))
                self.assertEqual(session.committed_units, 1)
                self.assertFalse(session.job.done)
                helper = self.assert_terrain_helper(backend, terrain)
                before_report = self.report_snapshot()
                before = evaluated_summary(runtime.app.padlBuildingResults)
                session.cancel()
                result = session.finish(silent=True)
                self.record(backend, 'shared_edges_terrain', 'cancel_after_first')
                self.assertEqual(result['outcome'], 'cancelled')
                self.assertFalse(result['errors'])
                self.assertEqual(terrain_state(terrain), source)
                self.assertEqual(self.report_snapshot(), before_report)
                self.assert_summary_equal(before, evaluated_summary(runtime.app.padlBuildingResults))
                retained_helper = self.assert_terrain_helper(backend, terrain)
                if helper:
                    self.assertTrue(alive(helper))
                    self.assertEqual(retained_helper, helper)
                self.assertFalse(session.ownership.pending)
                self.assert_no_readers()


try:
    run(args, PadlDatasetChecks)
finally:
    (args.artifacts / 'padl_dataset_metrics.json').write_text(json.dumps(dict(
        description='Observed import timings; no speedup claim. Geometry inspection is outside import timing.',
        runs=metrics), indent=2), encoding='utf-8')
