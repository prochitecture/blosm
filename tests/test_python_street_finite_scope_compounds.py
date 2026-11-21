"""Compose a finite scope with its actual adjacent deferred source junction."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.finite_scope_compounds import ExpandedScopePort,ScopePlaneCuts
from tested_python_streets.finite_containment import contains_segment
from tested_python_streets.ownership import PLAN,footprint_key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.updates import UpdateCoordinator
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from tested_python_streets.mesh import cross
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder


class OwnershipOnly(IntersectionProvider):
    def tasks(self,network):
        for task in super().tasks(network):
            if task.key.kind in ('approach','physical_intersection','compound_plan'):yield task


def fixture():
    # Authored synthetic complete ground network: two close sites, a recorded
    # scope, and a third real junction beyond the narrow continuing profile.
    # Remote ends of the other authored roads are explicitly true free ends.
    a=(0.,0.,0.);b=(18.6395263671875,-5.00018310546875,0.)
    far=(-25.50335693359375,-15.114501953125,0.)
    middle=(-11.83251953125,3.1627197265625,0.)
    specs=((1,7.,(a,middle)),(2,2.,(a,(-14.2369384765625,-1.18963623046875,0.))),
        (3,7.,(a,b)),(4,7.,(middle,(-204.9942626953125,54.94976806640625,0.))),
        (5,7.,(b,(71.9366455078125,-20.033447265625,0.))),
        (6,2.,(b,(4.908447265625,-4.74273681640625,0.),(-3.92156982421875,-5.777099609375,0.),
            (-15.2880859375,-9.85015869140625,0.),far)),
        (7,2.,(far,(far[0]-20.,far[1],0.))),(8,2.,(far,(far[0],far[1]-20.,0.))))
    sections=tuple(model.Section(i,1 if i in (1,3,4,5) else i,points,width=w,
        scope_ends=(1,) if i==2 else ()) for i,w,points in specs)
    sites=(model.Intersection(1,(model.End(1,0),model.End(2,0),model.End(3,0)),
        radii=(2.,1.,2.),radius_modes=('STYLE',)*3),
        model.Intersection(2,(model.End(3,1),model.End(5,0),model.End(6,0)),
            radii=(2.,2.,1.),radius_modes=('STYLE',)*3),
        model.Intersection(3,(model.End(6,1),model.End(7,0),model.End(8,0)),
            radii=(1.,1.,1.),radius_modes=('STYLE',)*3))
    return model.Network('composed-finite-ground-scope',sections,sites,
        settings=model.Settings(exaggeration=1.),connections=(model.Connection(9,(model.End(1,1),model.End(4,0))),))


def prepare(n):return UpdateCoordinator((PathProvider(),ProfileProvider(),OwnershipOnly()),PatternLibrary()).prepare(n)


class FiniteScopeCompoundTests(unittest.TestCase):
    def qualify(self,n):
        before=n.to_dict();result=prepare(n);self.assertEqual(before,n.to_dict())
        g=result.products[PLAN].by_site[model.Key('intersection',1)]
        self.assertIs(g,result.products[PLAN].by_site[model.Key('intersection',2)])
        self.assertEqual(g.members,(model.Key('intersection',1),model.Key('intersection',2)))
        physical=g.physical_layout;port,=physical.finite_scope_ports
        self.assertIsInstance(port,ExpandedScopePort)
        self.assertEqual({a.source.id for a in physical.approaches},{1,2,5,6})
        self.assertEqual({c.section_id for c in g.claims},{1,2,3})
        self.assertTrue(physical.finite_scope_source_proof)
        self.assertTrue(physical.finite_scope_curb_sampling)
        self.assertTrue(physical.finite_scope_rounded_contacts)
        self.assertTrue(any(c.fit_policy=='COMPOSED_SCOPE_ORIGINAL_CURB_CONTACT' for c in physical.corners))
        for curve in physical.boundaries:
            self.assertTrue(curve.source_curves)
            self.assertLessEqual(curve.maximum_projection_error,1e-6)
            np.testing.assert_array_equal(curve.points[[0,-1]],curve.original_polyline[[0,-1]])
            edges=np.diff(curve.points,axis=0);unit=edges/np.linalg.norm(edges,axis=1)[:,None]
            turn=np.arctan2(cross(unit[:-1],unit[1:]),np.sum(unit[:-1]*unit[1:],axis=1))
            self.assertLess(float(abs(turn).max()),.06)
        self.assertEqual(port.report()['boundary_kind'],'OPEN_SCOPE_BOUNDARY')
        self.assertTrue(port.report()['reconciled_both_halves'])
        self.assertFalse(port.report()['endpoint_extrapolation'])
        self.assertEqual(port.aliases[0].scope_end,model.End(2,n.section_index[2].scope_ends[0]))
        center=port.spur.frame(port.spur.length)[0]
        np.testing.assert_allclose(port.point(0.,.5*port.width*port.profile.peak),center,atol=1e-10)
        np.testing.assert_allclose(port.profile.derivative(np.array((-1.,port.profile.peak,1.))),0.,atol=1e-12)
        self.assertLessEqual(port.profile.value(np.linspace(-1.,1.,1001)).max(),port.spur.profile.rise+1e-12)
        far=result.products[footprint_key(model.Key('intersection',3))]
        end=port.other.end
        mouth=physical.mouths[6,end].station;opposite=far.mouths[6,1-end].station
        self.assertLess(mouth,opposite) if end==0 else self.assertGreater(mouth,opposite)
        cuts=ScopePlaneCuts().cuts(g,None)
        self.assertEqual(len(cuts),2)
        for cut in cuts:
            def point(ref):
                p=g.loops[ref.loop].xy;return p[ref.edge]+ref.fraction*(p[(ref.edge+1)%len(p)]-p[ref.edge])
            self.assertTrue(contains_segment(point(cut.first),point(cut.last),g.loops[0].xy,tolerance=1e-7))
        with self.assertRaisesRegex(ValueError,'extrapolated'):port.point(.01,0.)
        return result

    def test_original_incidence_width_grade_scope_and_real_far_profile(self):self.qualify(fixture())

    def test_complete_mesh_original_profiles_uv_stored_slopes_joins_and_drainage(self):
        n=fixture();before=n.to_dict();result=PythonMeshBackend(n.settings).prepare(n)
        self.assertEqual(before,n.to_dict())
        product=result.products[model.Key('junction_surface',1)]
        part=product.pavement;limit=part.metadata['measured_slope_limit']
        self.assertEqual(part.metadata['finite_port_reconciliation']['joins']['isolated_interior_minima'],0)
        self.assertEqual(part.metadata['finite_port_reconciliation']['rounds'][-1]['drainage']['maximum_fill'],0.)
        fixed=np.zeros(len(part.mesh.xyz),bool)
        for ids in part.mesh.seams.values():fixed[ids]=True
        stored=Mesh(part.mesh.xyz.astype(np.float32).astype(float),part.mesh.quads)
        sources={f.region.owner:f.layout for f in product.group.footprints}
        conditions=[replace(BoundaryConditions().build(product.group,loop,sources),
            xy=loop.xy.astype(np.float32).astype(float)) for loop in product.group.loops]
        joins=CompoundSurfaceBuilder.projection_checks(stored,fixed,conditions)
        self.assertEqual(joins['isolated_interior_minima'],0)
        self.assertLessEqual(joins['maximum_entrance_longitudinal_slope'],.03)
        for p in result.products.values():
            if getattr(p,'mesh',None) is None:continue
            self.assertIn('UVMap',p.mesh.uv_layers)
            self.assertTrue(np.isfinite(p.mesh.uv_layers['UVMap'].astype(np.float32)).all())
            rounded=Mesh(p.mesh.xyz.astype(np.float32).astype(float),p.mesh.quads)
            self.assertLessEqual(rounded.validate()['maximum_edge_ratio'],8.)
            self.assertLessEqual(float(rounded.quad_slopes().max()),limit+1e-7)
        physical=product.group.physical_layout
        for a,s in zip(physical.approaches,physical.outward_stations):
            if isinstance(a,ExpandedScopePort):continue
            endpoint,_=list(a.path.trims(s).items())[-1]
            road=result.products[model.Key('road',endpoint[0])]
            seam=road.mesh.xyz[road.mesh.seams[endpoint[1]]]
            np.testing.assert_allclose(seam[:,2],a.profile.value(np.linspace(-1.,1.,len(seam))),atol=1e-12)

    def test_reversed_storage_and_reflected_rotated_complete_network(self):
        n=fixture();swap=lambda e:model.End(e.section_id,1-e.end)
        reverse=replace(n,sections=tuple(replace(s,points=s.points[::-1],scope_ends=tuple(1-e for e in s.scope_ends)) for s in n.sections),
            intersections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in n.intersections),
            connections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in n.connections))
        self.qualify(reverse)
        rotation=np.array(((np.cos(.37),-np.sin(.37)),(np.sin(.37),np.cos(.37))))
        sections=[]
        for s in n.sections:
            p=np.asarray(s.points).copy();p[:,1]*=-1.;p[:,:2]=p[:,:2]@rotation.T+(1200.,-700.)
            sections.append(replace(s,points=tuple(map(tuple,p))))
        self.qualify(replace(n,sections=tuple(sections)))

    def test_unknown_selection_strict_and_missing_real_incidence_are_not_composed(self):
        n=fixture()
        variants=(replace(n,sections=tuple(replace(s,scope_ends=()) if s.id==2 else s for s in n.sections)),
            replace(n,sections=tuple(replace(s,scope_ends=(),selection_ends=(1,)) if s.id==2 else s for s in n.sections)),
            replace(n,intersections=tuple(replace(j,radius_modes=('STRICT',)*3) if j.id==2 else j for j in n.intersections)),
            replace(n,intersections=tuple(j for j in n.intersections if j.id!=1)))
        for variant in variants:
            with self.assertRaises(ValueError):prepare(variant)


if __name__=='__main__':unittest.main()
