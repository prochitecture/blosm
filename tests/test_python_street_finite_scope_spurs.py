"""Open finite ground scope cuts beside a supported continuing road."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.mesh import Mesh
from tested_python_streets.ownership import PLAN,surface_key


def fixture(*,split=True,mode=3):
    # Explicit synthetic source topology: a finite scope spur beside a long
    # continuing ground road. Remote ends here are authored true free ends;
    # actual-import qualification separately keeps the complete source graph.
    middle=(-11.83251953125,3.1627197265625,0.)
    far=(-204.9942626953125,54.94976806640625,0.)
    sections=[model.Section(1,1,((0.,0.,0.),middle if split else far),width=7.,interpolation=mode),
        model.Section(2,2,((0.,0.,0.),(-14.2369384765625,-1.18963623046875,0.)),
            width=2.,scope_ends=(1,),interpolation=mode),
        model.Section(3,1,((0.,0.,0.),(18.6395263671875,-5.00018310546875,0.)),width=7.,interpolation=mode)]
    seams=()
    if split:
        sections.append(model.Section(4,1,(middle,far),width=7.,interpolation=mode))
        seams=(model.Connection(2,(model.End(1,1),model.End(4,0))),)
    junction=model.Intersection(1,tuple(model.End(i,0) for i in (1,2,3)),
        radii=(2.,1.,2.),radius_modes=('STYLE',)*3)
    return model.Network('finite-open-scope-spur',tuple(sections),(junction,),
        settings=model.Settings(exaggeration=1.),connections=seams)


def reverse(n,ids):
    swap=lambda e:model.End(e.section_id,1-e.end) if e.section_id in ids else e
    return replace(n,sections=tuple(replace(s,points=s.points[::-1],
        scope_ends=tuple(1-e for e in s.scope_ends)) if s.id in ids else s for s in n.sections),
        intersections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in n.intersections),
        connections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in n.connections))


class FiniteScopeSpurTests(unittest.TestCase):
    def qualify(self,network):
        before=network.to_dict()
        update=PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(before,network.to_dict())
        owner=model.Key('intersection',1)
        group=update.products[PLAN].by_site[owner]
        layout=group.physical_layout
        port,=layout.finite_scope_ports
        alias,=port.aliases
        self.assertEqual(alias.scope_end,model.End(2,network.section_index[2].scope_ends[0]))
        self.assertEqual(port.report()['boundary_kind'],'OPEN_SCOPE_BOUNDARY')
        self.assertTrue(port.report()['reconciled_inner_half'])
        self.assertFalse(port.report()['endpoint_extrapolation'])
        self.assertEqual(set(s.id for s in network.sections),
            {key.id for key in update.source_fingerprints if key.kind=='section'})
        self.assertIsNone(getattr(update.products[model.Key('road',2)],'mesh',None))
        spur=port.spur
        physical_peak=port.point(0.,.5*port.width*port.profile.peak)
        np.testing.assert_allclose(physical_peak,spur.frame(spur.length)[0],atol=1e-10)
        q=np.linspace(-1.,0.,101) if port.outer_side<0 else np.linspace(0.,1.,101)
        original=spur.point(spur.length,.5*spur.width*q)
        physical_q=2.*((original-port.origin)@port.normal)/port.width
        np.testing.assert_allclose(port.profile.value(physical_q),spur.profile.value(q),atol=1e-12)
        self.assertAlmostEqual(float(port.profile.derivative(-port.outer_side)),0.,places=12)
        self.assertLessEqual(float(port.profile.value(np.linspace(-1.,1.,1001)).max()),spur.profile.rise+1e-12)
        product=update.products[surface_key(owner)]
        part=update.products[model.Key('junction_mesh',1)]
        limit=part.metadata['measured_slope_limit']
        self.assertLessEqual(limit,.176000001)
        checks=part.metadata['finite_port_reconciliation']
        self.assertEqual(checks['joins']['isolated_interior_minima'],0)
        self.assertLessEqual(checks['joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(checks['rounds'][-1]['drainage']['maximum_fill'],0.)
        fixed=np.zeros(len(part.mesh.xyz),bool)
        for nodes in part.mesh.seams.values():fixed[nodes]=True
        stored=Mesh(part.mesh.xyz.astype(np.float32).astype(float),part.mesh.quads)
        data={f.region.owner:f.layout for f in product.group.footprints}
        conditions=[replace(BoundaryConditions().build(product.group,loop,data),
            xy=loop.xy.astype(np.float32).astype(float)) for loop in product.group.loops]
        joins=CompoundSurfaceBuilder.projection_checks(stored,fixed,conditions)
        self.assertEqual(joins['isolated_interior_minima'],0)
        self.assertLessEqual(joins['maximum_entrance_longitudinal_slope'],.03)
        for p in update.products.values():
            if getattr(p,'mesh',None) is None:continue
            self.assertIn('UVMap',p.mesh.uv_layers)
            self.assertTrue(np.isfinite(p.mesh.uv_layers['UVMap'].astype(np.float32)).all())
            rounded=Mesh(p.mesh.xyz.astype(np.float32).astype(float),p.mesh.quads)
            self.assertLessEqual(rounded.validate()['maximum_edge_ratio'],8.)
            self.assertLessEqual(float(rounded.quad_slopes().max()),limit+1e-7)
        # The continuing wide road still joins on its complete original crown.
        wide=port.wide
        endpoint,station=list(wide.path.trims(port.wide_station).items())[-1]
        road=update.products[model.Key('road',endpoint[0])]
        row=road.mesh.xyz[road.mesh.seams[endpoint[1]]]
        np.testing.assert_allclose(row[:,2],wide.profile.value(np.linspace(-1.,1.,len(row))),atol=1e-12)
        with self.assertRaisesRegex(ValueError,'extrapolated'):port.point(.1,0.)
        return update

    def test_original_topology_profiles_uvs_and_float32_mesh(self):
        self.qualify(fixture())

    def test_source_reversal_exact_and_reflected_ground_geometry(self):
        exact=fixture(split=False,mode=0)
        axis=np.array((-11.83251953125,3.1627197265625));axis/=np.linalg.norm(axis)
        exact=replace(exact,sections=tuple(replace(s,points=((0.,0.,0.),
            tuple(axis*212. if s.id==1 else -axis*19.3)+(0.,))) if s.id in (1,3) else s for s in exact.sections))
        variants=[reverse(fixture(),{1,2,3,4}),exact]
        n=fixture();rot=np.array(((np.cos(.713),-np.sin(.713)),(np.sin(.713),np.cos(.713))))
        sections=[]
        for s in n.sections:
            p=np.asarray(s.points).copy();p[:,1]*=-1.
            p[:,:2]=p[:,:2]@rot.T+(1300.,-600.)
            sections.append(replace(s,points=tuple(map(tuple,p))))
        variants.append(replace(n,sections=tuple(sections)))
        for n in variants:
            with self.subTest(source=n.to_dict()):self.qualify(n)

    def test_unknown_free_selection_and_strict_are_not_scope_handoffs(self):
        n=fixture()
        variants=[replace(n,sections=tuple(replace(s,scope_ends=()) if s.id==2 else s for s in n.sections)),
            replace(n,sections=tuple(replace(s,scope_ends=(),selection_ends=(1,)) if s.id==2 else s for s in n.sections)),
            replace(n,intersections=(replace(n.intersections[0],radius_modes=('STRICT','STYLE','STYLE')),))]
        for n in variants:
            with self.assertRaises(ValueError):PythonMeshBackend(n.settings).prepare(n)

    def test_local_scope_edit_after_release_matches_fresh_and_failure_is_transactional(self):
        n=fixture();engine=PythonMeshBackend(n.settings)
        engine.commit(engine.prepare(n));engine.release_published_meshes()
        changed=replace(n,sections=tuple(replace(s,points=(s.points[0],(-14.20,-1.18,0.)))
            if s.id==2 else s for s in n.sections))
        local=engine.prepare(changed)
        fresh=PythonMeshBackend(changed.settings).prepare(changed)
        for key,p in local.products.items():
            if getattr(p,'mesh',None) is None:continue
            np.testing.assert_array_equal(p.mesh.xyz,fresh.products[key].mesh.xyz)
            np.testing.assert_array_equal(p.mesh.quads,fresh.products[key].mesh.quads)
        engine.commit(local);manifest=engine.manifest()
        unsupported=replace(changed,sections=tuple(replace(s,scope_ends=())
            if s.id==2 else s for s in changed.sections))
        with self.assertRaises(ValueError):engine.prepare(unsupported)
        self.assertEqual(engine.manifest(),manifest)


if __name__=='__main__':unittest.main()
