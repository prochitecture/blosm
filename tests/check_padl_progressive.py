"""Native simple-3D buildings: atomic progressive groups on both backends.

Run through the release runner with Blender 5.2 and the packaged PADL assets.
All fixture lookups and output locations follow release_support arguments.
"""
from pathlib import Path
import importlib
import json
import math
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
from release_support import (arguments, enable, catalog, configure_osm,
                             assert_geometry, run, PACKAGE)

args = arguments()
addon = enable(args)
runtime = importlib.import_module(PACKAGE + '.app.osm_import')
ownership = importlib.import_module(PACKAGE + '.progressive.ownership')
backends = importlib.import_module(PACKAGE + '.building.padl_integration.backends')
read_footprints = importlib.import_module(PACKAGE + '.padl.persistence').read_footprints
FaceType = importlib.import_module(PACKAGE + '.padl.bindings').FaceType
BACKENDS = ('geometry_nodes', 'python')


def evaluated_geometry(results):
    """Compare final volumes, including GN evaluation and excluding Python inputs."""
    bpy.context.view_layer.update()
    graph = bpy.context.evaluated_depsgraph_get()
    faces = []
    for result in results:
        evaluated = result.output.evaluated_get(graph)
        mesh = evaluated.to_mesh()
        try:
            assert mesh.vertices and mesh.polygons, 'Generated building is empty'
            for vertex in mesh.vertices:
                assert all(math.isfinite(value) for value in vertex.co)
            for face in mesh.polygons:
                assert face.area > 1e-9, 'Generated building has a degenerate face'
                faces.append(sorted(tuple(evaluated.matrix_world @ mesh.vertices[i].co)
                                    for i in face.vertices))
        finally:
            evaluated.to_mesh_clear()
    return sorted(faces)


def data_ids():
    return {kind: {item.as_pointer() for item in getattr(bpy.data, kind) if not item.get(ownership.TRANSIENT)}
            for kind in ownership.KINDS}


class PadlProgressiveChecks(unittest.TestCase):
    def prepare(self, backend, *, fixture='shared_edges', progressive=True, single=True):
        entry = next(item for item in catalog() if item['id'] == fixture)
        settings = configure_osm(args, entry, '3Dsimple', single)
        settings.buildingBackend = backend
        settings.buildingSeed = 113
        settings.progressiveImport = progressive
        return settings

    def session(self):
        session = runtime.OsmImport(bpy.context)

        def cleanup():
            if not session.finished:
                session.cancel()
                session.finish(silent=True)

        self.addCleanup(cleanup)
        return session

    def assert_no_temporary_readers(self):
        from padl_status_support import assert_status_resources
        assert_status_resources(runtime)

    def assert_report(self, count, input_count):
        results = runtime.app.padlBuildingResults
        self.assertEqual(len(results), count)
        self.assertIs(runtime.app.padlBuildingResult, results[0] if count == 1 else None)
        report = runtime.app.padlBuildingReport
        self.assertEqual(report['input_count'], input_count)
        self.assertEqual(report['statuses'], [row for result in results for row in result.statuses])
        self.assertTrue(report['statuses'], 'Completed buildings have no generation statuses')
        self.assertTrue(all(not row.get('error') for row in report['statuses']), report)
        text = results[0].owner.get('blosm_building_report')
        self.assertIsInstance(text, bpy.types.Text)
        self.assertEqual(json.loads(text.as_string()), report)
        for result in results:
            self.assertEqual(result.owner.get('blosm_building_report'), text)
            self.assertTrue(ownership.alive(result.owner))
            self.assertTrue(ownership.alive(result.output))
            self.assertTrue(read_footprints(result.owner), 'Editable inputs were lost')
            self.assertNotIn(ownership.PENDING, result.owner)
            self.assertNotIn(ownership.PENDING, result.output)
            if bpy.context.scene.blosm.buildingBackend == 'python':
                self.assertNotEqual(result.owner, result.output)
                self.assertTrue(result.owner.hide_get())
                self.assertTrue(result.owner.hide_render)
            else:
                self.assertEqual(result.owner, result.output)
        self.assert_no_temporary_readers()
        # A snapshot prevents later mutations from making a comparison vacuous.
        return json.loads(json.dumps(report))

    def first_group(self, session, input_count=3):
        for _ in range(20000):
            if session.committed_units or session.job.done:
                break
            session.job.step()
        self.assertIsNone(session.job.error, str(session.job.error))
        self.assertEqual(session.committed_units, 1)
        self.assertFalse(session.job.done, 'First commit did not yield to the UI')
        self.assert_report(1, input_count)
        return runtime.app.padlBuildingResults[0]

    def drain(self, session, *, input_count=3, expected_groups=3):
        commits = []
        previous = session.committed_units
        for _ in range(20000):
            if session.job.done:
                break
            session.job.step()
            if session.committed_units != previous:
                self.assertEqual(session.committed_units, previous + 1,
                                 'Multiple building groups committed in one UI turn')
                previous = session.committed_units
                commits.append(previous)
                self.assert_report(previous, input_count)
        self.assertTrue(session.job.done, 'Progressive building import did not terminate')
        self.assertIsNone(session.job.error, str(session.job.error))
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'success', result)
        self.assertFalse(result['errors'], result)
        self.assertEqual(result['committed'], expected_groups)
        self.assert_report(expected_groups, input_count)
        self.assertFalse(session.ownership.pending)
        return commits

    def blocking(self, backend, *, single, fixture='shared_edges', groups=3, inputs=3):
        self.prepare(backend, fixture=fixture, progressive=False, single=single)
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertEqual(runtime.last_result['outcome'], 'success')
        self.assert_report(1 if single else groups, inputs)
        return evaluated_geometry(runtime.app.padlBuildingResults)

    def test_grouped_and_separate_blocking_match_progressive_outputs(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                grouped = self.blocking(backend, single=True)
                separate = self.blocking(backend, single=False)
                assert_geometry(self, grouped, separate)
                settings = self.prepare(backend)
                session = self.session()
                self.assertEqual(self.drain(session), [1, 2, 3])
                self.assertTrue(settings.singleObject, 'Saved blocking preference was changed')
                assert_geometry(self, separate, evaluated_geometry(runtime.app.padlBuildingResults))

    def test_cancel_keeps_complete_group_metadata_report_and_user_edits(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                self.prepare(backend)
                session = self.session()
                first = self.first_group(session)
                first.owner.name = 'Inputs renamed during import'
                first.output.name = 'Output renamed during import'
                first.output.location.z += 7
                sentinel = bpy.data.objects.new('Unrelated user object', None)
                bpy.context.scene.collection.objects.link(sentinel)
                note = bpy.data.texts.new('Unrelated user notes')
                note.write('Keep this edit')
                before = data_ids()
                report = self.assert_report(1, 3)
                geometry = evaluated_geometry((first,))
                session.cancel()
                outcome = session.finish(silent=True)
                self.assertEqual(outcome['outcome'], 'cancelled')
                self.assertEqual(outcome['committed'], 1)
                self.assertFalse(outcome['errors'])
                self.assertEqual(data_ids(), before, 'Cancellation removed committed or user data')
                self.assertEqual(self.assert_report(1, 3), report)
                self.assertEqual(note.as_string(), 'Keep this edit')
                assert_geometry(self, geometry, evaluated_geometry((first,)))
                self.assertFalse(session.ownership.pending)
                self.assertIsNone(session.finish(silent=True))

    def test_failed_next_group_rolls_back_all_pending_data_only(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                self.prepare(backend)
                session = self.session()
                first = self.first_group(session)
                first.owner.name = 'Previously committed editable inputs'
                sentinel_mesh = bpy.data.meshes.new('User mesh between turns')
                sentinel = bpy.data.objects.new('User object between turns', sentinel_mesh)
                bpy.context.scene.collection.objects.link(sentinel)
                sentinel['user_edit'] = 17
                note = bpy.data.texts.new('User text between turns')
                note.write('Preserve on failure')
                before = data_ids()
                report = self.assert_report(1, 3)
                geometry = evaluated_geometry((first,))
                generated = []
                original = backends.BuildingGenerationService.generate

                def fail_after_generation(service, *positional, **keywords):
                    # Exercise rollback after mesh/output and persistence Text
                    # allocations, not merely a failure before work starts.
                    result = original(service, *positional, **keywords)
                    generated.append((result.owner.name, result.output.name))
                    raise RuntimeError('Injected failure after native building generation')

                with patch.object(backends.BuildingGenerationService, 'generate', fail_after_generation):
                    for _ in range(20000):
                        if session.job.done:
                            break
                        session.job.step()
                self.assertTrue(session.job.done)
                self.assertEqual(len(generated), 1)
                self.assertIn('Injected failure', str(session.job.error))
                outcome = session.finish(silent=True)
                self.assertEqual(outcome['outcome'], 'partial', outcome)
                self.assertEqual(outcome['committed'], 1)
                self.assertTrue(any('Injected failure' in error for error in outcome['errors']))
                self.assertEqual(data_ids(), before, 'Failed group leaked data or removed user data')
                self.assertEqual(self.assert_report(1, 3), report)
                self.assertEqual(sentinel['user_edit'], 17)
                self.assertEqual(note.as_string(), 'Preserve on failure')
                assert_geometry(self, geometry, evaluated_geometry((first,)))
                self.assertFalse(session.ownership.pending)

    def test_first_gn_group_failure_reclaims_recursive_assets_and_texts(self):
        self.prepare('geometry_nodes')
        # The common reset retains fake-user assets from earlier test cases.
        # Start this isolated test with a cold GN library so rollback must deal
        # with the root, its recursive dependencies and their document Texts.
        bpy.data.batch_remove(ids=list(bpy.data.node_groups))
        bpy.data.orphans_purge(do_recursive=True)
        note = bpy.data.texts.new('Pre-existing user document')
        note.write('Keep the complete unrelated dependency graph')
        helper = bpy.data.node_groups.new('Pre-existing user helper', 'GeometryNodeTree')
        helper['document'] = note
        root = bpy.data.node_groups.new('Pre-existing user generator', 'GeometryNodeTree')
        root.nodes.new('GeometryNodeGroup').node_tree = helper
        root.use_fake_user = True
        sentinel = bpy.data.objects.new('Pre-existing user object', None)
        bpy.context.scene.collection.objects.link(sentinel)
        sentinel['generator'] = root
        before = data_ids()
        # user_map includes ID kinds beyond Ownership.KINDS (e.g. libraries).
        all_before = {item.as_pointer(): (item.bl_rna.identifier, item.name)
                      for item in bpy.data.user_map()}
        session = self.session()
        generated = []
        original = backends.BuildingGenerationService.generate

        def fail_after_generation(service, *positional, **keywords):
            result = original(service, *positional, **keywords)
            tree = result.owner.modifiers['PADL Buildings'].node_group
            graph, pending = set(), [tree]
            while pending:
                group = pending.pop()
                pointer = group.as_pointer()
                if pointer in graph:
                    continue
                graph.add(pointer)
                pending.extend(node.node_tree for node in group.nodes
                               if getattr(node, 'node_tree', None) is not None)
            current = data_ids()
            generated.append({
                'graph': graph,
                'new_groups': current['node_groups'] - before['node_groups'],
                'new_texts': current['texts'] - before['texts'],
            })
            raise RuntimeError('Injected failure before the first native building commit')

        with patch.object(backends.BuildingGenerationService, 'generate', fail_after_generation):
            for _ in range(20000):
                if session.job.done:
                    break
                session.job.step()
        self.assertTrue(session.job.done)
        self.assertEqual(len(generated), 1, 'Failure did not reach real GN generation')
        self.assertGreater(len(generated[0]['graph']), 1, 'No recursive GN dependencies exercised')
        self.assertTrue(generated[0]['graph'] <= generated[0]['new_groups'],
                        'The cold-library test reused a pre-existing PADL dependency')
        self.assertTrue(generated[0]['new_texts'], 'No persistence Text allocations exercised')
        self.assertEqual(session.committed_units, 0)
        self.assertIn('Injected failure', str(session.job.error))
        outcome = session.finish(silent=True)
        self.assertEqual(outcome['outcome'], 'failed', outcome)
        self.assertEqual(outcome['committed'], 0)
        self.assertTrue(any('Injected failure' in error for error in outcome['errors']))
        self.assert_no_temporary_readers()
        self.assertFalse(session.ownership.pending)
        self.assertFalse(runtime.app.padlBuildingResults)
        self.assertEqual(data_ids(), before, 'First failed group leaked assets or removed existing data')
        all_after = {item.as_pointer(): (item.bl_rna.identifier, item.name)
                     for item in bpy.data.user_map()}
        self.assertEqual(all_after, all_before, 'Datablock inventory changed after first-group failure')
        self.assertEqual(note.as_string(), 'Keep the complete unrelated dependency graph')
        self.assertEqual(sentinel['generator'], root)
        self.assertEqual(root.nodes[0].node_tree, helper)
        self.assertEqual(helper['document'], note)

    def test_building_outline_and_nine_parts_commit_together(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                baseline = self.blocking(backend, single=False, fixture='building_parts', groups=1, inputs=10)
                self.prepare(backend, fixture='building_parts')
                session = self.session()
                first = self.first_group(session, input_count=10)
                footprints = read_footprints(first.owner)
                self.assertEqual(len(footprints), 10)
                parents = [row for row in footprints if row.face_type == FaceType.BUILDING_WITH_PARTS]
                parts = [row for row in footprints if row.face_type == FaceType.PART]
                self.assertEqual(len(parents), 1)
                self.assertEqual(len(parts), 9)
                self.assertTrue(all(row.parent_id == parents[0].footprint_id for row in parts))
                self.assertEqual(len(first.statuses), 9, 'Outline became an extra volume')
                assert_geometry(self, baseline, evaluated_geometry((first,)))
                self.assertEqual(self.drain(session, input_count=10, expected_groups=1), [])
                self.assertEqual(len(runtime.app.padlBuildingResults), 1)

    def assert_gn_edit_stops_import(self, edit):
        self.prepare('geometry_nodes')
        handlers = tuple(bpy.app.handlers.depsgraph_update_post)
        session = self.session()
        first = self.first_group(session)
        print('GN asset edit:', edit, 'first group ready', flush=True)
        tree = first.owner.modifiers['PADL Buildings'].node_group
        groups, pending = {}, [tree]
        while pending:
            group = pending.pop()
            if group.as_pointer() in groups:
                continue
            groups[group.as_pointer()] = group
            pending.extend(node.node_tree for node in group.nodes
                           if getattr(node, 'node_tree', None) is not None)
        nested = next(group for group in groups.values() if group != tree)
        changed_node, changed_text, original_text = None, None, None
        changed_default = None
        try:
            if edit in ('root', 'nested'):
                target = tree if edit == 'root' else nested
                changed_node = (target, target.nodes.new('ShaderNodeMath'))
            elif edit == 'interface':
                socket = next(item for item in tree.interface.items_tree
                              if item.item_type == 'SOCKET' and item.in_out == 'INPUT'
                              and item.name == 'Seed')
                changed_default = socket, socket.default_value
                socket.default_value += 1
            elif edit == 'document':
                changed_text = tree['padl_document']
                original_text = changed_text.as_string()
                changed_text.write('\nUser changed the persisted program\n')
            else:
                bpy.data.node_groups.remove(tree if edit == 'delete_root' else nested, do_unlink=True)
            print('GN asset edit:', edit, 'mutation complete', flush=True)
            # Real UI edits are evaluated between timer callbacks. Do not mark
            # the guard dirty directly: verify the actual Blender notification.
            # Deletion is checked independently before the importer evaluates
            # a graph with missing helpers; forcing that invalid graph to
            # evaluate can stall Blender itself. Text edits need no notification.
            if edit in ('root', 'nested', 'interface'):
                bpy.context.view_layer.update()
            print('GN asset edit:', edit, 'evaluation complete', flush=True)
            for _ in range(20000):
                if session.job.done:
                    break
                session.job.step()
            self.assertTrue(session.job.done)
            self.assertIsNotNone(session.job.error, 'Modified shared GN asset was silently reused')
            outcome = session.finish(silent=True)
            print('GN asset edit:', edit, 'import finished', flush=True)
            self.assertEqual(outcome['outcome'], 'partial', outcome)
            self.assertEqual(outcome['committed'], 1)
            self.assertTrue(ownership.alive(first.owner), 'Completed user-editable owner was removed')
            self.assertEqual(len(runtime.app.padlBuildingResults), 1)
            self.assertFalse(session.ownership.pending)
            self.assert_no_temporary_readers()
            self.assertEqual(tuple(bpy.app.handlers.depsgraph_update_post), handlers,
                             'Progressive GN invalidation handler leaked after failure')
        finally:
            if not session.finished:
                session.cancel()
                session.finish(silent=True)
            if changed_node is not None:
                changed_node[0].nodes.remove(changed_node[1])
            if changed_text is not None:
                changed_text.clear()
                changed_text.write(original_text)
            if changed_default is not None:
                changed_default[0].default_value = changed_default[1]
            if edit.startswith('delete_'):
                # A deleted helper leaves a deliberately invalid published
                # graph. Remove this test's remaining closure before reuse.
                bpy.data.batch_remove(ids=[group for group in groups.values() if ownership.alive(group)])
            print('GN asset edit:', edit, 'cleanup complete', flush=True)

    def test_progressive_gn_rejects_edited_root_graph(self):
        self.assert_gn_edit_stops_import('root')

    def test_progressive_gn_rejects_changed_interface_default(self):
        self.assert_gn_edit_stops_import('interface')

    def test_progressive_gn_unedited_batch_reuses_asset_verification(self):
        self.prepare('geometry_nodes')
        handlers = tuple(bpy.app.handlers.depsgraph_update_post)
        checks = []
        original = backends.BuildingProgramLibrary.tree

        def verified_tree(library, role):
            checks.append(role)
            return original(library, role)

        with patch.object(backends.BuildingProgramLibrary, 'tree', verified_tree):
            self.drain(self.session())
        self.assertEqual(checks, ['simple_3d'], 'Own generation invalidated the cached GN verification')
        self.assertEqual(tuple(bpy.app.handlers.depsgraph_update_post), handlers,
                         'Progressive GN invalidation handler leaked after success')

    def test_progressive_gn_rejects_edited_nested_graph(self):
        self.assert_gn_edit_stops_import('nested')

    def test_progressive_gn_rejects_persisted_text_edits_without_graph_update(self):
        self.assert_gn_edit_stops_import('document')

    def test_progressive_gn_rejects_deleted_root_or_dependency(self):
        for edit in ('delete_root', 'delete_nested'):
            with self.subTest(edit=edit):
                self.assert_gn_edit_stops_import(edit)


    def test_finish_unit_failure_restores_committed_report(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                self.prepare(backend)
                session = self.session()
                first = self.first_group(session)
                report = self.assert_report(1, 3)
                before = data_ids()
                with patch.object(runtime.Renderer, 'finishUnit',
                                  side_effect=RuntimeError('Injected unit finalization failure')):
                    for _ in range(20000):
                        if session.job.done:
                            break
                        session.job.step()
                self.assertTrue(session.job.done)
                self.assertIn('Injected unit finalization failure', str(session.job.error))
                # The failed unit has already written into the accepted Text;
                # its metadata must be restored before pending IDs are removed.
                self.assertEqual(self.assert_report(1, 3), report)
                outcome = session.finish(silent=True)
                self.assertEqual(outcome['outcome'], 'partial', outcome)
                self.assertEqual(outcome['committed'], 1)
                self.assertEqual(data_ids(), before)
                self.assertEqual(self.assert_report(1, 3), report)
                self.assertIs(runtime.app.padlBuildingResults[0], first)

    def test_blocking_later_group_failure_does_not_publish_pending_results(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                self.prepare(backend, progressive=False, single=False)
                before = data_ids()
                session = self.session()
                generated = []
                original = backends.BuildingGenerationService.generate

                def fail_second_group(service, *positional, **keywords):
                    result = original(service, *positional, **keywords)
                    generated.append(result)
                    if len(generated) == 2:
                        raise RuntimeError('Injected later blocking group failure')
                    return result

                with patch.object(backends.BuildingGenerationService, 'generate', fail_second_group):
                    for _ in range(20000):
                        if session.job.done:
                            break
                        session.job.step()
                self.assertTrue(session.job.done)
                self.assertEqual(len(generated), 2)
                self.assertIn('Injected later blocking group failure', str(session.job.error))
                outcome = session.finish(silent=True)
                self.assertEqual(outcome['outcome'], 'failed', outcome)
                self.assertEqual(outcome['committed'], 0)
                self.assertEqual(runtime.app.padlBuildingResults, ())
                self.assertIsNone(runtime.app.padlBuildingResult)
                self.assertEqual(runtime.app.padlBuildingReport['statuses'], [])
                self.assertEqual(data_ids(), before)
                self.assertTrue(all(not ownership.alive(result.owner) for result in generated))
                self.assert_no_temporary_readers()

    def test_deleted_report_after_final_group_preserves_warning(self):
        for backend in BACKENDS:
            with self.subTest(backend=backend):
                self.prepare(backend)
                session = self.session()
                first = self.first_group(session)
                manager = next(manager for manager in runtime.app.managers
                               if manager.id == 'padl_buildings')
                manager.report.issue('test:source', 'test_notice', 'Source mapping notice')
                for _ in range(20000):
                    if session.committed_units == 3 or session.job.done:
                        break
                    session.job.step()
                self.assertEqual(session.committed_units, 3)
                self.assertFalse(session.job.done, 'Final group did not expose a UI boundary')
                text = first.owner['blosm_building_report']
                name = text.name
                bpy.data.texts.remove(text)
                for _ in range(20000):
                    if session.job.done:
                        break
                    session.job.step()
                self.assertTrue(session.job.done)
                self.assertIsNone(session.job.error, str(session.job.error))
                outcome = session.finish(silent=True)
                self.assertEqual(outcome['outcome'], 'success', outcome)
                self.assertEqual(outcome['committed'], 3)
                self.assertEqual(len(outcome['warnings']), 1)
                self.assertIn('1 mapping/source notices.', outcome['warnings'][0])
                self.assertNotIn('See ', outcome['warnings'][0])
                self.assertIsNone(bpy.data.texts.get(name), 'Cleanup recreated user-deleted data')
                self.assertEqual(len(runtime.app.padlBuildingResults), 3)
                self.assert_no_temporary_readers()


run(args, PadlProgressiveChecks)
