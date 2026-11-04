"""Numerical corridor fit: topology, complete profiles, clearance and contacts."""
import unittest
from dataclasses import replace
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.terrain import PlacementSettings, CorridorComponent, fit
from tested_python_streets.terrain.corridor import prepare
from tested_python_streets.terrain.alignment import quad_samples
from tested_python_streets.providers import Part
from tested_python_streets.mesh import Mesh
from _python_street_fixtures import transition_fixture_networks, short_junction_pair


def terrain(xy):
    return 3. + .015*xy[:,0] + .025*xy[:,1] + .18*np.exp(-np.sum((xy-(10.,2.))**2,axis=1)/70.)


class TerrainTests(unittest.TestCase):
    def network(self):
        return model.Network('terrain_test', (
            model.Section(1,1,((-50.,0.,0.),(0.,0.,0.),(50.,0.,0.))),))

    def check_network(self, network, extra=None):
        backend = PythonMeshBackend(network.settings)
        update = backend.prepare(network)
        settings = PlacementSettings()
        refs, comps, _ = prepare(network, update.products, settings)
        if extra:
            extra(refs,comps)
        result = fit(refs, comps, terrain, settings)
        self.assertLessEqual(result.report['maximum_reference_grade'], settings.maximum_grade+1e-10)
        contacts = {}
        for c in comps:
            placed = result.coordinates[c.identity]
            np.testing.assert_array_equal(placed[:,:2],c.xyz[:,:2])
            for s in np.unique(c.stations):
                delta = (placed[:,2]-c.xyz[:,2])[c.stations == s]
                self.assertLess(float(np.ptp(delta)),1e-12)
            for s,contact in c.contacts:
                height = float(result.curves[c.reference](s))
                self.assertAlmostEqual(height, contacts.setdefault(contact,height), places=10)
            for ids, sample in quad_samples(placed,c.quads,settings.sample_spacing):
                clearance = sample[:,:,2]-terrain(sample[:,:,:2].reshape(-1,2)).reshape(sample.shape[:2])
                self.assertGreaterEqual(float(clearance.min()),settings.clearance-1e-8)
            curve = result.curves[c.reference]
            self.assertEqual(curve.derivatives[0],0.)
            self.assertEqual(curve.derivatives[-1],0.)
        return result

    def test_full_crown_and_side_slope(self):
        self.check_network(self.network())

    def test_junctions_and_shared_minimum(self):
        self.check_network(short_junction_pair(model))

    def test_divided_profiles_and_separator(self):
        self.check_network(transition_fixture_networks(model)[1])

    def test_future_raised_strip_constrains_same_reference(self):
        def side_strip(refs,comps):
            road = comps[0]
            xyz = road.xyz.copy()
            xyz[:,1] += 12.
            xyz[:,2] += .15
            comps.append(CorridorComponent('future sidewalk',road.reference,xyz,
                road.quads,road.stations,road.datum,clearance=.01))
        base = self.check_network(self.network())
        raised = self.check_network(self.network(),side_strip)
        self.assertGreater(float(raised.coordinates['road:1'][:,2].max()),float(base.coordinates['road:1'][:,2].max()))

    def test_missing_terrain_and_infeasible_lift(self):
        network=self.network(); update=PythonMeshBackend(network.settings).prepare(network)
        refs,comps,_=prepare(network,update.products,PlacementSettings())
        with self.assertRaisesRegex(ValueError,'cover'):
            fit(refs,comps,lambda xy:np.full(len(xy),np.nan))
        with self.assertRaisesRegex(ValueError,'lifts'):
            fit(refs,comps,terrain,PlacementSettings(maximum_lift=.01))

    def test_feature_metadata_binds_a_side_strip(self):
        network=self.network();update=PythonMeshBackend(network.settings).prepare(network)
        road=update.products[model.Key('road',1)]
        xyz=road.mesh.xyz.copy();xyz[:,1]+=12.;xyz[:,2]+=.15
        stations=np.repeat(road.metadata['terrain_stations'],road.metadata['transverse_intervals']+1)
        update.products[model.Key('sidewalk',1)]=Part(network.sections[0].key,Mesh(xyz,road.mesh.quads),
            dict(terrain_corridor=dict(section_id=1,stations=stations.tolist(),clearance=.01)))
        refs,comps,_=prepare(network,update.products,PlacementSettings())
        result=fit(refs,comps,terrain)
        a=result.coordinates['road:1'][:,2]-road.mesh.xyz[:,2]
        b=result.coordinates['sidewalk:1'][:,2]-xyz[:,2]
        np.testing.assert_allclose(a,b,atol=1e-12)

    def test_flat_terrain_preserves_relief_exactly(self):
        network=self.network(); update=PythonMeshBackend(network.settings).prepare(network)
        refs,comps,_=prepare(network,update.products,PlacementSettings())
        result=fit(refs,comps,lambda xy:np.full(len(xy),7.))
        np.testing.assert_allclose(result.coordinates[comps[0].identity][:,2],comps[0].xyz[:,2]+7.02,atol=1e-12)

    def test_reference_finer_than_mesh_has_no_between_row_dips(self):
        network=self.network();update=PythonMeshBackend(network.settings).prepare(network)
        settings=PlacementSettings(control_spacing=.1)
        refs,comps,_=prepare(network,update.products,settings)
        result=fit(refs,comps,lambda xy:np.full(len(xy),3.),settings)
        for curve in result.curves.values():
            samples=np.linspace(curve.stations[0],curve.stations[-1],3001)
            np.testing.assert_allclose(curve(samples),3.02,atol=1e-12)


if __name__ == '__main__':
    unittest.main()
