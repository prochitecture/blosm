"""A finite wide free cap belongs to one profiled source-backed junction."""
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


def fixture():
    # Authored synthetic finite wide pad with two external ground approaches.
    # Full original-import qualification retains the remote source junctions.
    a=((15.0621337890625,150.10589599609375,0.),(10.0977783203125,75.73248291015625,0.),
       (8.335205078125,58.34405517578125,0.),(0.,0.,0.))
    b=((0.,0.,0.),(7.91845703125,-2.77044677734375,0.),(17.1048583984375,-5.6185302734375,0.),
       (21.943115234375,-8.05560302734375,0.),(28.3685302734375,-11.7724609375,0.),
       (30.8814697265625,-15.71270751953125,0.),(31.6448974609375,-21.93536376953125,0.),
       (30.5269775390625,-33.37921142578125,0.),(28.3443603515625,-50.23333740234375,0.))
    return model.Network('finite-wide-free-cap',(
        model.Section(1,1,a,width=7.),
        model.Section(2,1,((0.,0.,0.),(-.3262939453125,-2.293212890625,0.)),width=7.),
        model.Section(3,3,b,width=2.)),
        (model.Intersection(1,(model.End(1,1),model.End(2,0),model.End(3,0)),
            radii=(2.,2.,1.),radius_modes=('STYLE',)*3),),settings=model.Settings(exaggeration=1.))


def reverse(n,ids):
    swap=lambda e:model.End(e.section_id,1-e.end) if e.section_id in ids else e
    return replace(n,sections=tuple(replace(s,points=s.points[::-1]) if s.id in ids else s for s in n.sections),
        intersections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in n.intersections))


class FiniteFreeCapTests(unittest.TestCase):
    def qualify(self,n):
        original=n.to_dict();update=PythonMeshBackend(n.settings).prepare(n)
        self.assertEqual(n.to_dict(),original)
        owner=model.Key('intersection',1)
        group=update.products[PLAN].by_site[owner]
        claim,=group.physical_layout.finite_free_pads
        self.assertTrue(claim.source_union_cap)
        self.assertEqual(claim.pad.source.id,2)
        self.assertFalse(claim.pad.source.scope_ends)
        free_curb = group.physical_layout.boundaries[0]
        self.assertEqual(len(free_curb.free_cap_banks), len(free_curb.points)-1)
        self.assertIn(None, free_curb.free_cap_banks)
        for bank in filter(None, free_curb.free_cap_banks):
            role, side, first, last = bank
            self.assertIn(side, (-1,1))
            self.assertTrue(0. <= min(first,last) < max(first,last) <= claim.distances[role]+1e-8)
        self.assertIsNone(getattr(update.products[model.Key('road',2)],'mesh',None))
        end=claim.pad.frame(claim.pad.length)[0]
        np.testing.assert_allclose((claim.cap_points-end)@claim.pad.direction,0.,atol=1e-10)
        self.assertGreaterEqual(claim.radius,.5)
        for curve in claim.side_curves:self.assertLessEqual(curve.curvature_bound,(1.+1e-8)/claim.radius)
        product=update.products[surface_key(owner)]
        part=update.products[model.Key('junction_mesh',1)]
        closest=int(np.argmin(np.linalg.norm(part.mesh.xyz[:,:2]-end,axis=1)))
        self.assertLess(np.linalg.norm(part.mesh.xyz[closest,:2]-end),1e-8)
        self.assertEqual(part.mesh.xyz[closest,2],0.)
        check=part.metadata['finite_port_reconciliation']
        self.assertEqual(check['joins']['isolated_interior_minima'],0)
        self.assertLessEqual(check['joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(check['rounds'][-1]['drainage']['maximum_fill'],0.)
        fixed=np.zeros(len(part.mesh.xyz),bool)
        for ids in part.mesh.seams.values():fixed[ids]=True
        stored=Mesh(part.mesh.xyz.astype(np.float32).astype(float),part.mesh.quads)
        layouts={f.region.owner:f.layout for f in product.group.footprints}
        boundaries=[replace(BoundaryConditions().build(product.group,loop,layouts),
            xy=loop.xy.astype(np.float32).astype(float)) for loop in product.group.loops]
        joins=CompoundSurfaceBuilder.projection_checks(stored,fixed,boundaries)
        self.assertEqual(joins['isolated_interior_minima'],0)
        self.assertLessEqual(joins['maximum_entrance_longitudinal_slope'],.03)
        limit=part.metadata['measured_slope_limit']
        self.assertLessEqual(limit,.17600001)
        for p in update.products.values():
            if getattr(p,'mesh',None) is None:continue
            self.assertIn('UVMap',p.mesh.uv_layers)
            self.assertTrue(np.isfinite(p.mesh.uv_layers['UVMap'].astype(np.float32)).all())
            rounded=Mesh(p.mesh.xyz.astype(np.float32).astype(float),p.mesh.quads)
            self.assertLessEqual(rounded.validate()['maximum_edge_ratio'],8.)
            self.assertLessEqual(float(rounded.quad_slopes().max()),limit+1e-7)
        for approach,station in zip(group.physical_layout.approaches,group.physical_layout.outward_stations):
            endpoint,at=list(approach.path.trims(station).items())[-1]
            road=update.products[model.Key('road',endpoint[0])]
            row=road.mesh.xyz[road.mesh.seams[endpoint[1]]]
            np.testing.assert_allclose(row[:,2],approach.profile.value(np.linspace(-1.,1.,len(row))),atol=1e-12)
        return update

    def test_whole_finite_cap_profiles_and_stored_quad_quality(self):self.qualify(fixture())

    def test_reversed_sources_and_rigid_transform(self):
        n=reverse(fixture(),{1,2,3})
        angle=.713;rot=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        sections=[]
        for s in n.sections:
            p=np.asarray(s.points).copy();p[:,:2]=p[:,:2]@rot.T+(1111.,689.)
            sections.append(replace(s,points=tuple(map(tuple,p))))
        self.qualify(replace(n,sections=tuple(sections)))

    def test_scope_selection_and_strict_inputs_remain_failures(self):
        n=fixture()
        variants=[replace(n,sections=tuple(replace(s,scope_ends=(1,)) if s.id==2 else s for s in n.sections)),
            replace(n,sections=tuple(replace(s,selection_ends=(1,)) if s.id==2 else s for s in n.sections)),
            replace(n,intersections=(replace(n.intersections[0],radius_modes=('STYLE','STRICT','STYLE')),))]
        for n in variants:
            with self.assertRaises(ValueError):PythonMeshBackend(n.settings).prepare(n)


if __name__=='__main__':unittest.main()
