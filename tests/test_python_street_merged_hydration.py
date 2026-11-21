"""Only source-relevant released junctions are privately rebuilt for a Bundle."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import short_junction_pair
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_merged_plan import physical_context
from tested_python_streets.ownership import surface_key
from tested_python_streets.updates import BuildContext


def pair_context(exaggeration):
    source = short_junction_pair(model)
    remote = short_junction_pair(model)
    sections = tuple(replace(s, id=s.id+100, street_id=s.street_id+100,
        points=tuple((x+200., y+200., z) for x, y, z in s.points)) for s in remote.sections)
    sites = tuple(replace(j, id=j.id+100,
        ends=tuple(model.End(e.section_id+100,e.end) for e in j.ends)) for j in remote.intersections)
    network = replace(source, sections=(*source.sections,*sections),
        intersections=(*source.intersections,*sites),settings=replace(source.settings,exaggeration=exaggeration))
    engine = PythonMeshBackend(network.settings)
    prepared = engine.prepare(network)
    return BuildContext(network,prepared.products,engine.patterns,prepared.task_fingerprints,prepared.source_fingerprints)


class MergedPhysicalHydrationTests(unittest.TestCase):
    def test_released_bundle_update_matches_fresh_and_restores_without_touching_unrelated_owner(self):
        from test_python_street_merged_bundles import source_for,snapshot
        source = source_for(3)
        engine = PythonMeshBackend(source.settings)
        initial = engine.prepare(source)
        original = snapshot(initial)
        engine.commit(initial)
        unrelated_key = surface_key(model.Key('intersection',210))
        unrelated = initial.products[unrelated_key]
        engine.release_published_meshes()
        changed = replace(source,sections=tuple(replace(s,width=s.width+.1) if s.id==17 else s for s in source.sections))
        local = engine.prepare(changed)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertTrue(snapshot(local))
        for key,(xyz,quads) in snapshot(local).items():
            np.testing.assert_array_equal(xyz,fresh.products[key].mesh.xyz)
            np.testing.assert_array_equal(quads,fresh.products[key].mesh.quads)
        for item in changed.furniture:
            self.assertEqual(local.products[model.Key('placement',item.id)],fresh.products[model.Key('placement',item.id)])
        self.assertIs(local.products[unrelated_key],unrelated)
        self.assertIsNone(unrelated.pavement.mesh)
        engine.commit(local)
        engine.release_published_meshes()
        restored = engine.prepare(source)
        for key,(xyz,quads) in snapshot(restored).items():
            np.testing.assert_array_equal(xyz,original[key][0])
            np.testing.assert_array_equal(quads,original[key][1])

    def test_unrelated_release_is_skipped_and_required_release_rebuilds_privately(self):
        context = pair_context(1.)
        local,remote = model.Key('intersection',1),model.Key('intersection',101)
        current = context[surface_key(local)]
        saved = current.physical_pavement
        unrelated = context[surface_key(remote)]
        unrelated.release_geometry()
        view = physical_context(context,required_owners=(local,))
        self.assertIs(view[surface_key(local)].pavement.mesh,saved)
        self.assertIs(view[surface_key(remote)],unrelated)
        self.assertIsNone(unrelated.pavement.mesh)
        self.assertEqual(view.hydrated_physical_owners,())
        current.release_geometry()
        rebuilt = physical_context(context,required_owners=(local,))
        self.assertEqual(rebuilt.hydrated_physical_owners,(local,))
        self.assertIs(context[surface_key(local)],current)
        self.assertIsNone(current.pavement.mesh)
        self.assertIsNone(current.physical_pavement)
        self.assertIsNone(unrelated.pavement.mesh)
        np.testing.assert_array_equal(rebuilt[surface_key(local)].pavement.mesh.xyz,saved.xyz)
        np.testing.assert_array_equal(rebuilt[surface_key(local)].pavement.mesh.quads,saved.quads)

    def test_zero_exaggeration_release_reconstructs_original_physical_heights(self):
        context = pair_context(0.)
        owner = model.Key('intersection',1)
        product = context[surface_key(owner)]
        reference = product.physical_pavement
        self.assertEqual(float(product.pavement.mesh.xyz[:,2].max()),0.)
        self.assertGreater(float(reference.xyz[:,2].max()),.07)
        product.release_geometry()
        view = physical_context(context,required_owners=(owner,))
        physical = view[surface_key(owner)]
        self.assertEqual(physical.exaggeration,1.)
        self.assertEqual(context.network.settings.exaggeration,0.)
        self.assertIsNone(product.physical_pavement)
        np.testing.assert_array_equal(physical.pavement.mesh.xyz,reference.xyz)
        np.testing.assert_array_equal(physical.pavement.mesh.quads,reference.quads)
        with self.assertRaisesRegex(ValueError,'original junction incidence'):
            physical_context(context,required_owners=(model.Key('intersection',999),))


if __name__=='__main__':unittest.main()
