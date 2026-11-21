"""Finite shared scope mouths retain source incidence and physical guards."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.finite_scope_ports import _contains_segment
from tested_python_streets.mesh import Mesh
from tested_python_streets.ownership import PLAN, surface_key
from tested_python_streets.regions import _classify


def fixture():
    # Authored synthetic finite Y: two retained arms have distinct known scope
    # incidences; the third road really has a free remote end in this fixture.
    # Its finite widths/lengths reproduce the physical overlapping-tail class.
    ends=((11.828857421875,4.1549072265625,0.),
          (8.9735107421875,8.5625,0.),(-22.5750732421875,-14.74359130859375,0.))
    sections=tuple(model.Section(i,i,((0.,0.,0.),point),width=5.400000095367432,
        rise=.07999999821186066,scope_ends=(1,) if i<3 else ()) for i,point in enumerate(ends,1))
    junction=model.Intersection(1,tuple(model.End(i,0) for i in (1,2,3)),
        radius_modes=('STYLE',)*3)
    return model.Network('finite-ground-scope-pool',sections,(junction,),(),(),model.Settings())


def reverse(network,ids):
    swap=lambda e:model.End(e.section_id,1-e.end) if e.section_id in ids else e
    return replace(network,sections=tuple(replace(s,points=s.points[::-1],
        scope_ends=tuple(1-e for e in s.scope_ends)) if s.id in ids else s for s in network.sections),
        intersections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in network.intersections))


class FiniteScopePortTests(unittest.TestCase):
    def qualify(self,network):
        source=network.to_dict()
        prepared=PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(),source)
        owner=model.Key('intersection',1)
        group=prepared.products[PLAN].by_site[owner]
        layout=group.physical_layout
        self.assertEqual(len(layout.finite_scope_ports),1)
        port=layout.finite_scope_ports[0]
        self.assertEqual({a.scope_end for a in port.aliases},
            {model.End(s.id,s.scope_ends[0]) for s in network.sections if s.scope_ends})
        for alias in port.aliases:
            path=prepared.products[model.Key('path',alias.junction_end.section_id)]
            self.assertAlmostEqual(alias.length,path.length)
            self.assertEqual(alias.scope_station,0. if alias.scope_end.end==0 else path.length)
            self.assertIsNone(getattr(prepared.products[model.Key('road',alias.junction_end.section_id)],'mesh',None))
        self.assertEqual([c.radius for c in layout.corners],[2.,2.])
        product=prepared.products[surface_key(owner)]
        self.assertEqual(product.pavement.metadata['consumed_sections'],[1,2])
        checks=product.pavement.metadata['finite_port_reconciliation']
        self.assertTrue(checks['exact_boundary'])
        self.assertEqual(checks['joins']['isolated_interior_minima'],0)
        self.assertLessEqual(checks['joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(checks['rounds'][-1]['drainage']['maximum_fill'],0.)
        stored=Mesh(product.pavement.mesh.xyz.astype(np.float32).astype(float),product.pavement.mesh.quads)
        stored.xyz[:,2]/=network.settings.exaggeration
        fixed=np.zeros(len(stored.xyz),bool)
        for nodes in product.pavement.mesh.seams.values():fixed[nodes]=True
        source_layouts={p.region.owner:p.layout for p in product.group.footprints}
        boundaries=[replace(BoundaryConditions().build(product.group,loop,source_layouts),
            xy=loop.xy.astype(np.float32).astype(float)) for loop in product.group.loops]
        stored_joins=CompoundSurfaceBuilder.projection_checks(stored,fixed,boundaries)
        self.assertEqual(stored_joins['isolated_interior_minima'],0)
        self.assertLessEqual(stored_joins['maximum_entrance_longitudinal_slope'],.03)
        for p in prepared.products.values():
            if getattr(p,'mesh',None) is None:continue
            self.assertLessEqual(p.mesh.validate()['maximum_edge_ratio'],8.)
            self.assertIn('UVMap',p.mesh.uv_layers)
            self.assertTrue(np.isfinite(p.mesh.uv_layers['UVMap'].astype(np.float32)).all())
            stored=Mesh(p.mesh.xyz.astype(np.float32).astype(float),p.mesh.quads)
            self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
            self.assertLessEqual(float(stored.quad_slopes(network.settings.exaggeration).max()),.1+1e-7)
        # The finite outline keeps every original scope corner or covers it
        # inside the union, including corners not on a surviving terminal cut.
        for a in port.original_approaches:
            p=a.point(a.length,np.array((-.5,.5))*a.width)
            inside,boundary=_classify(p,group.loops[0].xy,1e-7)
            self.assertTrue(np.all(inside | (boundary>=0)))
        physical=product.pavement.mesh
        for corner in port.points:
            self.assertLess(float(np.linalg.norm(physical.xyz[:,:2]-corner,axis=1).min()),1e-8)
        q=np.linspace(-1.,1.,2001)
        self.assertAlmostEqual(float(port.profile.value(0.)),port.profile.rise)
        self.assertGreaterEqual(float(port.profile.value(q).min()),0.)
        self.assertLessEqual(float(port.profile.value(q).max()),port.profile.rise)
        np.testing.assert_allclose(port.profile.derivative(2*port.parameters[1:-1]-1.),0.,atol=1e-12)
        for bank, side in ((port.original_approaches[1],-1.),(port.original_approaches[0],1.)):
            np.testing.assert_allclose(port.physical_gradient(0.,side),
                                       bank.physical_gradient(bank.length,side), atol=1e-10)
        with self.assertRaisesRegex(ValueError,'extrapolated'):port.point(.01,0.)
        return prepared

    def test_full_mesh_profiles_quality_uv_and_stored_coordinates(self):
        for mode in (0,3):
            with self.subTest(interpolation=mode):
                n=fixture()
                self.qualify(replace(n,sections=tuple(replace(s,interpolation=mode) for s in n.sections)))

    def test_directions_resolution_and_rotated_translated_source(self):
        for ids,resolution,angle in (({1},1,0.),({2,3},1,.731),({1,2,3},2,2.1)):
            with self.subTest(reverse=ids,resolution=resolution):
                n=reverse(fixture(),ids)
                rot=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
                sources=[]
                for s in n.sections:
                    points=np.asarray(s.points).copy()
                    points[:,:2]=points[:,:2]@rot.T+(1366.,1029.)
                    sources.append(replace(s,points=tuple(map(tuple,points))))
                self.qualify(replace(n,sections=tuple(sources),settings=replace(n.settings,resolution=resolution)))

    def test_unknown_selection_or_strict_incidence_cannot_be_merged(self):
        n=fixture()
        variants=[replace(n,sections=tuple(replace(s,scope_ends=()) if s.id==1 else s for s in n.sections)),
            replace(n,sections=tuple(replace(s,scope_ends=(),selection_ends=(1,)) if s.id==1 else s for s in n.sections)),
            replace(n,intersections=(replace(n.intersections[0],radius_modes=('STRICT','STYLE','STYLE')),))]
        for network in variants:
            with self.assertRaises(ValueError):PythonMeshBackend(network.settings).prepare(network)

    def test_continuous_source_proof_detects_a_subsample_exit(self):
        # A 1 micrometre notch would evade a fixed-density centreline scan.
        polygon=np.array(((-1.,-1.),(11.,-1.),(11.,1.),(5.000001,1.),
            (5.000001,-.1),(5.,-.1),(5.,1.),(-1.,1.)))
        self.assertFalse(_contains_segment(np.array((0.,0.)),np.array((10.,0.)),polygon))


if __name__=='__main__':unittest.main()
