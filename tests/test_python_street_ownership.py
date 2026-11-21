from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, merging_bundle_network, short_junction_pair
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import PLAN, StandaloneProduct, surface_key
from tested_python_streets.providers import Part, RoadSpan


class PythonJunctionOwnershipTests(unittest.TestCase):
    def test_same_owner_corrected_mouth_epoch_rebuilds_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_SAMPLING_POLICY')

    def test_exact_owned_ranges_epoch_rebuilds_plan_and_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_RANGE_POLICY')

    def test_native_bank_parameter_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_CURB_COORDINATE_POLICY')

    def test_source_bank_rounding_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_BANK_POLICY')

    def test_retained_arc_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_ARC_POLICY')

    def test_owned_connector_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_OWNED_CONNECTOR_POLICY')

    def test_original_parent_chain_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_PARENT_CHAIN_POLICY')

    def test_residual_gap_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_RESIDUAL_GAP_POLICY')

    def test_owned_bundle_normal_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_OWNED_BUNDLE_NORMAL_POLICY')

    def test_owned_prefix_corner_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_OWNED_PREFIX_CORNER_POLICY')

    def test_native_curvature_root_epoch_rebuilds_all_projected_consumers(self):
        self._assert_projected_policy_epoch('SOURCE_ROOT_POLICY')

    def _assert_projected_policy_epoch(self,policy):
        from types import SimpleNamespace
        from unittest.mock import patch
        from tested_python_streets import ownership
        from tested_python_streets.compounds import CompoundGroup
        from tested_python_streets.intersections import Mouth
        from tested_python_streets.junction_transitions import JunctionTransitions
        from tested_python_streets.patterns import PatternLibrary
        from tested_python_streets.updates import Task,UpdateCoordinator
        network=short_junction_pair(model);before=network.to_dict()
        members=tuple(j.key for j in network.intersections)
        raw={ownership.footprint_key(j.key):object() for j in network.intersections}
        seam=model.Key('test_seam',0);dependent=model.Key('test_dependent',0)
        class ProjectedConsumers:
            station=5.
            def tasks(self,network):
                for junction in network.intersections:
                    key=ownership.footprint_key(junction.key)
                    yield Task(key,(junction.key,),lambda context,key=key:raw[key])
                def plan(context):
                    layout=SimpleNamespace(approaches=(),mouths={(3,0):Mouth(self.station,8)})
                    group=CompoundGroup(members=members,claims=(),footprints=(),loops=(),physical_layout=layout)
                    return ownership.OwnershipPlan((group,),{key:group for key in members})
                tasks=tuple(ownership.JunctionOwnership.tasks(network))
                yield replace(next(t for t in tasks if t.key==PLAN),build=plan)
                for junction in network.intersections:
                    yield Task(junction.layout_key,(ownership.footprint_key(junction.key),PLAN),
                        lambda context,j=junction:ownership.JunctionOwnership.layout(context,j),
                        signature=lambda context,j=junction:ownership.JunctionOwnership.layout_signature(context,j))
                    task=next(t for t in tasks if t.key==surface_key(junction.key))
                    yield replace(task,build=lambda context,j=junction:context[j.layout_key].mouths[3,0].station)
                yield Task(seam,(*raw,PLAN),lambda context:self.station,
                    signature=lambda context:JunctionTransitions.seam_signature(context,tuple(raw),network.intersections))
                yield Task(dependent,(surface_key(members[0]),network.intersections[0].layout_key,seam),
                    lambda context:(context[surface_key(members[0])],context[seam]))
        provider=ProjectedConsumers();coordinator=UpdateCoordinator((provider,),PatternLibrary())
        with patch.object(ownership,policy,'legacy-source-ownership-policy'):
            legacy=coordinator.prepare(network);coordinator.commit(legacy)
        provider.station=7.
        current=coordinator.prepare(network)
        for key in (PLAN,*(j.layout_key for j in network.intersections),
                    *(surface_key(j.key) for j in network.intersections),seam,dependent):
            self.assertIn(key,current.rebuilt)
            self.assertNotEqual(legacy.task_fingerprints[key],current.task_fingerprints[key])
        self.assertEqual(legacy.products[PLAN].groups[0].members,current.products[PLAN].groups[0].members)
        self.assertEqual(legacy.products[PLAN].groups[0].physical_layout.mouths[3,0].station,5.)
        self.assertEqual(current.products[PLAN].groups[0].physical_layout.mouths[3,0].station,7.)
        self.assertEqual(current.products[dependent],(7.,7.))
        for key in raw:
            self.assertNotIn(key,current.rebuilt)
            self.assertIs(legacy.products[key],current.products[key])
            self.assertEqual(legacy.task_fingerprints[key],current.task_fingerprints[key])
        self.assertEqual(network.to_dict(),before)
        coordinator.commit(current)
        self.assertEqual(coordinator.prepare(network).rebuilt,())

    def test_complete_source_sampling_epoch_rebuilds_plan_and_reuses_raw_inputs(self):
        network=short_junction_pair(model)
        backend=PythonMeshBackend(network.settings)
        current=backend.updates.providers
        class LegacyPlanProvider:
            def __init__(self,provider):
                self.provider=provider
                self.source_kinds=provider.source_kinds
            def tasks(self,network):
                for task in self.provider.tasks(network):
                    yield replace(task,version='1') if task.key==PLAN else task
        backend.updates.providers=tuple(LegacyPlanProvider(p) for p in current)
        legacy=backend.prepare(network)
        backend.commit(legacy)
        backend.updates.providers=current
        update=backend.prepare(network)
        self.assertIn(PLAN,update.rebuilt)
        self.assertNotEqual(legacy.task_fingerprints[PLAN],update.task_fingerprints[PLAN])
        raw={key for key in legacy.products if key.kind in ('path','profile','physical_intersection')}
        self.assertTrue(raw)
        for key in raw:
            self.assertNotIn(key,update.rebuilt)
            self.assertIs(legacy.products[key],update.products[key])
            self.assertEqual(legacy.task_fingerprints[key],update.task_fingerprints[key])
        backend.commit(update)
        self.assertEqual(backend.prepare(network).rebuilt,())

    def test_divided_cross_keeps_four_external_separator_runs_and_one_central_island(self):
        network = compound_cross_network(model)
        update = PythonMeshBackend(network.settings).prepare(network)
        meshes = {key: p for key, p in update.products.items() if isinstance(p, Part)}
        self.assertEqual(sum(key.kind == 'junction_mesh' for key in meshes), 1)
        self.assertEqual(sum(key.kind == 'junction_grass' for key in meshes), 1)
        self.assertEqual(sum(isinstance(p, RoadSpan) for p in update.products.values()), 4)
        for bundle in network.bundles:
            separator = meshes[model.Key('separator', bundle.id)]
            self.assertEqual(len(separator.metadata['runs']), 2)
            self.assertLess(len(separator.mesh.quads), 3000)
            np.testing.assert_allclose(separator.mesh.xyz[:, 2], .8, atol=1e-12)
        root = meshes[model.Key('junction_mesh', 4)]
        self.assertLess(root.mesh.slope_statistics(network.settings.exaggeration)['maximum'], .1)
        query = update.products[surface_key(model.Key('intersection', 4))]
        np.testing.assert_allclose(query.height(np.array(((0., 0.), (100., 100.)))), (.8, 0.), atol=1e-12)

    def test_merging_bundle_uses_its_common_junction_as_the_rounded_separator_cap(self):
        network = merging_bundle_network(model)
        update = PythonMeshBackend(network.settings).prepare(network)
        separator = update.products[model.Key('separator', 1)]
        self.assertEqual(separator.metadata['attached_caps'], [True, False])
        self.assertGreater(separator.metadata['centerline_proof_checks'], 0)
        self.assertLess(len(separator.mesh.quads), 2000)
        np.testing.assert_allclose(separator.mesh.xyz[:, 2], .8, atol=1e-12)

    def test_default_backend_publishes_one_compound_and_consumes_the_short_link(self):
        network = short_junction_pair(model)
        prepared = PythonMeshBackend(network.settings).prepare(network)
        owner = model.Key('intersection', 1)
        self.assertEqual(prepared.products[PLAN].owner(owner), owner)
        part = prepared.products[model.Key('junction_mesh', 1)]
        self.assertIsInstance(part, Part)
        self.assertEqual(part.owner, model.Key('compound', 1))
        self.assertIsNone(prepared.products[model.Key('junction_mesh', 2)])
        self.assertIsInstance(prepared.products[model.Key('road', 1)], RoadSpan)
        self.assertEqual(part.metadata['member_junctions'], ['intersection:1', 'intersection:2'])
        self.assertTrue(part.metadata['default_provider_integration'])

    def test_edit_is_local_matches_clean_rebuild_and_moves_a_consumed_link_anchor(self):
        network = short_junction_pair(model)
        network = replace(network, sections=(*network.sections,
            model.Section(6, 6, ((60., 50., 0.), (90., 50., 0.)))),
            furniture=(model.Placement(81, 1, 3., offset=-2.),))
        engine = PythonMeshBackend(network.settings)
        original = engine.prepare(network)
        engine.commit(original)
        changed = replace(network, sections=tuple(replace(s, width=6.4) if s.id == 3 else s for s in network.sections))
        update = engine.prepare(changed)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertIn(surface_key(model.Key('intersection', 1)), update.rebuilt)
        self.assertIn(model.Key('placement', 81), update.rebuilt)
        self.assertNotIn(model.Key('road', 6), update.rebuilt)
        self.assertIs(update.products[model.Key('road', 6)], original.products[model.Key('road', 6)])
        for key, product in update.products.items():
            if isinstance(product, Part):
                np.testing.assert_array_equal(product.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(product.mesh.quads, fresh.products[key].mesh.quads)
        anchor = update.products[model.Key('placement', 81)]
        compound = update.products[surface_key(model.Key('intersection', 1))]
        self.assertAlmostEqual(anchor.position[2], compound.height(np.array((anchor.position[:2],)))[0])
        self.assertNotEqual(anchor.position[2], original.products[model.Key('placement', 81)].position[2])

    def test_released_compound_buffers_do_not_survive_through_query_aliases(self):
        network = short_junction_pair(model)
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        product = engine.updates.products[surface_key(model.Key('intersection', 1))]
        height = product.height(np.array(((3., 0.),)))
        engine.release_published_meshes()
        self.assertIsNone(product.pavement.mesh)
        self.assertIn('bounds', product.pavement.metadata)
        np.testing.assert_array_equal(height, product.height(np.array(((3., 0.),))))
        self.assertEqual(engine.prepare(network).rebuilt, ())

    def test_unrelated_compound_surfaces_retain_identity_after_global_plan_refresh(self):
        first, remote = short_junction_pair(model), short_junction_pair(model)
        remote_sections = tuple(replace(s, id=s.id + 100, street_id=s.street_id + 100,
            points=tuple((x + 100., y + 100., z) for x, y, z in s.points)) for s in remote.sections)
        remote_junctions = tuple(replace(j, id=j.id + 100,
            ends=tuple(model.End(e.section_id + 100, e.end) for e in j.ends)) for j in remote.intersections)
        network = replace(first, sections=(*first.sections, *remote_sections),
                          intersections=(*first.intersections, *remote_junctions))
        engine = PythonMeshBackend(network.settings)
        initial = engine.prepare(network)
        engine.commit(initial)
        changed = replace(network, sections=tuple(replace(s, width=6.2) if s.id == 3 else s for s in network.sections))
        update = engine.prepare(changed)
        self.assertIn(PLAN, update.rebuilt)
        self.assertNotIn(surface_key(model.Key('intersection', 101)), update.rebuilt)
        self.assertNotIn(model.Key('junction_mesh', 101), update.rebuilt)
        self.assertIs(initial.products[model.Key('junction_mesh', 101)], update.products[model.Key('junction_mesh', 101)])

    def test_compound_can_split_back_into_independent_junctions_without_stale_pavement(self):
        network = short_junction_pair(model)
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        separated = short_junction_pair(model, length=20.)
        update = engine.prepare(separated)
        self.assertFalse(update.products[PLAN].groups)
        for jid in (1, 2):
            self.assertEqual(update.products[model.Key('junction_mesh', jid)].owner, model.Key('intersection', jid))
            self.assertIsInstance(update.products[surface_key(model.Key('intersection', jid))], StandaloneProduct)
        self.assertIsInstance(update.products[model.Key('road', 1)], Part)
        engine.commit(update)
        restored = engine.prepare(network)
        fresh = PythonMeshBackend(network.settings).prepare(network)
        self.assertIsNone(restored.products[model.Key('junction_mesh', 2)])
        np.testing.assert_array_equal(restored.products[model.Key('junction_mesh', 1)].mesh.xyz,
                                      fresh.products[model.Key('junction_mesh', 1)].mesh.xyz)


if __name__ == '__main__':
    unittest.main()
