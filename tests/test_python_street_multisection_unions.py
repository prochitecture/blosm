"""Bundle unions retain Section chains and their junction-consumed prefixes."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import split_attached_fixture, attached_overlap_fixture_networks
from test_python_street_attached_unions import meshes
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.paths import Path


def reverse_middle(network):
    length = Path(network.section_index[21]).length
    swap = lambda e: model.End(21, 1-e.end) if e.section_id == 21 else e
    return replace(network,
        sections=tuple(replace(s, points=s.points[::-1]) if s.id == 21 else s for s in network.sections),
        connections=tuple(replace(c, ends=tuple(map(swap, c.ends))) for c in network.connections),
        bundles=(replace(network.bundles[0], members=tuple(replace(m,
            path=tuple(replace(t, forward=not t.forward) if t.section_id == 21 else t for t in m.path))
            for m in network.bundles[0].members)),),
        furniture=tuple(replace(p, station=length-p.station,
            side='RIGHT' if p.side == 'LEFT' else 'LEFT') if p.section_id == 21 else p for p in network.furniture))


class PythonMultiSectionUnionTests(unittest.TestCase):
    def equal(self, first, second):
        self.assertEqual(first.keys(), second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0], second[key][0])
            np.testing.assert_array_equal(first[key][1], second[key][1])

    def test_attached_chains_merge_pavement_and_keep_every_source_section(self):
        original = attached_overlap_fixture_networks(model)[0]
        reference = PythonMeshBackend(original.settings).prepare(original)
        cap_slope = reference.products[model.Key('separator', 7)].mesh.slope_statistics(10.)['maximum']
        for short in (False, True):
            with self.subTest(consumed=short):
                network = split_attached_fixture(model, short_ends=short)
                source = network.to_dict()
                update = PythonMeshBackend(network.settings).prepare(network)
                plan = update.products[model.Key('bundle_plan', 7)]
                self.assertEqual(plan.mode, 'union')
                self.assertEqual(set(plan.ranges), {1, 11, 12, 2, 21, 22})
                self.assertEqual(network.to_dict(), source)
                published = meshes(update)
                self.assertFalse(any(model.Key('road', sid) in published for sid in plan.ranges))
                grass = update.products[model.Key('separator', 7)]
                self.assertEqual(grass.metadata['islands'], 2)
                self.assertEqual(grass.metadata['welded_cap_body_vertices'], 34)
                for kind in ('bundle_pavement', 'separator'):
                    p = update.products[model.Key(kind, 7)]
                    self.assertLessEqual(p.mesh.validate()['maximum_edge_ratio'], 8.)
                    self.assertLessEqual(p.mesh.slope_statistics(10.)['maximum'], .1 if kind == 'bundle_pavement' else cap_slope + 1e-8)
                self.assertTrue(np.all(np.isfinite(update.products[model.Key('placement', 82)].position)))
                if short:
                    for sid in (1, 2, 12, 22):
                        self.assertLessEqual(plan.ranges[sid][1] - plan.ranges[sid][0], 1e-7)

    def test_both_interpolation_modes_and_reversed_middle_sections_keep_visible_ports(self):
        for mode in (0, 3):
            for short in (False, True):
                for reverse in (False, True):
                    with self.subTest(mode=mode, consumed=short, reverse=reverse):
                        network = split_attached_fixture(model, short_ends=short)
                        network = replace(network, sections=tuple(replace(s, interpolation=mode) for s in network.sections))
                        if reverse:
                            network = reverse_middle(network)
                        result = PythonMeshBackend(network.settings).prepare(network)
                        plan = result.products[model.Key('bundle_plan', 7)]
                        self.assertEqual(plan.mode, 'union')
                        self.assertEqual([[t.section_id for t in m.path] for m in plan.visible_members],
                                         [[11], [21]] if short else [[1, 11, 12], [2, 21, 22]])
                        self.assertEqual(result.products[model.Key('separator', 7)].metadata['islands'], 2)
                        self.assertLessEqual(result.products[model.Key('bundle_pavement', 7)].mesh.validate()['maximum_edge_ratio'], 8.)

    def test_finer_resolutions_reversed_consumed_members_and_translated_datum(self):
        for mode in (0, 3):
            for resolution in (2, 3, 4):
                with self.subTest(mode=mode, resolution=resolution):
                    network = split_attached_fixture(model, short_ends=True)
                    network = replace(network, sections=tuple(replace(s, interpolation=mode, datum=37.) for s in network.sections),
                        settings=replace(network.settings, resolution=resolution, exaggeration=0. if resolution == 2 else 10.))
                    network = reverse_middle(network)
                    update = PythonMeshBackend(network.settings).prepare(network)
                    for kind in ('bundle_pavement', 'separator'):
                        part = update.products[model.Key(kind, 7)]
                        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
                        self.assertGreaterEqual(float(part.mesh.xyz[:, 2].min()), 37.-1e-8)
                        self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), 37.+.12*network.settings.exaggeration+1e-8)

    def test_consumed_prefix_edits_match_full_build_restore_and_release_caches(self):
        network = split_attached_fixture(model, short_ends=True)
        engine = PythonMeshBackend(network.settings)
        initial = engine.prepare(network); engine.commit(initial)
        saved = meshes(initial)
        remote = initial.products[model.Key('road', 99)]
        middle = network.section_index[11]
        points = tuple((x, y+.2 if abs(x) < 1e-8 else y, z) for x, y, z in middle.points)
        edits = (
            replace(network, sections=tuple(replace(s, points=points) if s.id == 11 else s for s in network.sections)),
            replace(network, sections=tuple(replace(s, width=6.2) if s.id in {1, 11, 12, 2, 21, 22} else s
                                           for s in network.sections)),
            replace(network, intersections=tuple(replace(j, radius=3.) if j.id in (1, 2) else j
                                                for j in network.intersections)))
        for changed in edits:
            local = engine.prepare(changed)
            self.equal(meshes(local), meshes(PythonMeshBackend(changed.settings).prepare(changed)))
            self.assertIs(local.products[model.Key('road', 99)], remote)
            self.assertNotIn(model.Key('placement', 99), local.rebuilt)
            engine.commit(local)
            restored = engine.prepare(network)
            self.equal(meshes(restored), saved); engine.commit(restored)
        engine.release_published_meshes()
        self.assertIsNone(engine.updates.products[model.Key('bundle_plan', 7)]._pavement_part)
        self.equal(meshes(engine.prepare(network, force=True)), saved)

    def test_coincident_controls_without_explicit_seams_and_out_of_range_anchors_are_rejected(self):
        network = split_attached_fixture(model)
        broken = replace(network, connections=tuple(c for c in network.connections if c.id != 40))
        with self.assertRaisesRegex(ValueError, 'explicit two-ended Section seams'):
            PythonMeshBackend(broken.settings).prepare(broken)
        network = split_attached_fixture(model, short_ends=True)
        broken = replace(network, furniture=tuple(replace(p, station=2.0001) if p.id == 82 else p for p in network.furniture))
        with self.assertRaisesRegex(ValueError, 'station lies outside'):
            PythonMeshBackend(broken.settings).prepare(broken)


if __name__ == '__main__':
    unittest.main(verbosity=2)
