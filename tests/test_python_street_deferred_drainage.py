"""An unchanged ordinary fork must reconcile slopes without creating a basin."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Network,Section,Intersection,End,Key,Settings
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.intersections import Approach,IntersectionProvider
from tested_python_streets.ownership import PLAN,OwnershipPlan,JunctionOwnership,footprint_key
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.compounds import CompoundGroup,sample_footprint
from tested_python_streets.regions import RegionUnion
from tested_python_streets.drainage import BoundaryDrainage
from tested_python_streets.four_side_charts import quality


# Resolved real source geometry, with local identities. No external extract or
# saved preparation is required to reproduce the rejected projection.
SHAPES = (((-1323.4530029296875, -236.1239471435547),
  (-1336.9072265625, -218.6884307861328),
  (-1352.7115478515625, -205.6829071044922),
  (-1396.557861328125, -185.19044494628906),
  (-1409.2049560546875, -177.6623992919922),
  (-1422.23291015625, -168.33087158203125),
  (-1439.2734375, -153.57711791992188),
  (-1448.4193115234375, -144.02378845214844),
  (-1456.063720703125, -133.3909912109375),
  (-1542.7215576171875, 25.88353157043457),
  (-1555.914306640625, 57.77984619140625)),
 ((-1382.116943359375, -311.95318603515625),
  (-1385.8868408203125, -298.482666015625),
  (-1399.789794921875, -270.20440673828125),
  (-1401.7093505859375, -261.944091796875),
  (-1401.8240966796875, -253.3390655517578),
  (-1398.93603515625, -244.3116912841797),
  (-1394.4140625, -234.3273468017578),
  (-1384.804931640625, -225.90260314941406),
  (-1372.8084716796875, -221.20758056640625),
  (-1361.302734375, -219.5069122314453),
  (-1348.29638671875, -221.41331481933594),
  (-1336.6290283203125, -227.24896240234375),
  (-1323.4530029296875, -236.1239471435547)),
 ((-1323.4530029296875, -236.1239471435547),
  (-1314.88330078125, -250.16314697265625),
  (-1310.024169921875, -263.6672058105469),
  (-1303.5869140625, -282.2366638183594)))


def fixture():
    rise=float(np.float32(.08))
    sections=tuple(Section(i+1,i+1,tuple((x,y,0.) for x,y in points),width=3.5,rise=rise)
                   for i,points in enumerate(SHAPES))
    site=Intersection(1,(End(1,0),End(2,1),End(3,0)),radius=10.,radii=(10.,)*3,radius_modes=('STYLE',)*3)
    network=Network('curved_three_arm_drainage',sections,(site,),settings=Settings(exaggeration=1.))
    products={PLAN:OwnershipPlan((),{})}
    for section in sections:
        products[Key('path',section.id)]=Path(section)
        products[Key('profile',section.id)]=Crown(rise,0.)
    context=BuildContext(network,products,PatternLibrary())
    approaches=[Approach(Route((Piece(products[Key('path',end.section_id)],end.end==0),)),
                         end.end,products[Key('profile',end.section_id)]) for end in site.ends]
    physical=IntersectionProvider.curved_boundary(context,site,approaches)
    products[footprint_key(site.key)]=physical
    try:layout=IntersectionProvider.layout(context,site,physical)
    except ValueError:layout=replace(physical,mesh_strategy='boundary_chart')
    products[site.layout_key]=layout
    return context,site


class DeferredDrainageTests(unittest.TestCase):
    def test_exhausted_projection_can_recover_on_the_original_chart(self):
        from tested_python_streets.slope_projection import SlopeProjection
        context,site=fixture()
        footprint=sample_footprint(site.key,context[footprint_key(site.key)])
        group=CompoundGroup((site.key,),(),(footprint,),RegionUnion().build((footprint.region,)))
        calls=[]
        def projector(**kwargs):
            calls.append(None)
            return SlopeProjection(max_iterations=1 if len(calls)==1 else 1000,**kwargs)
        with patch('tested_python_streets.compound_surfaces.SlopeProjection',side_effect=projector):
            surface=CompoundSurfaceBuilder(context.patterns,max_refinements=0).build(group,context)
        report=surface.metadata['deferred_projection_reconciliation']
        self.assertIn('iteration budget',report['initial_projection_failure'])
        self.assertEqual(report['discarded_refinement_attempts'],[])
        fixed=np.zeros(len(surface.mesh.xyz),bool)
        fixed[list(surface.boundary_refs)]=True
        self.assertLessEqual(BoundaryDrainage().fill(surface.mesh,fixed)[1]['maximum_fill'],1e-6)
        self.assertLessEqual(report['joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertLessEqual(float(surface.mesh.quad_slopes().max()),surface.metadata['measured_slope_limit']+1e-8)

    def test_failed_projection_and_refinement_reconcile_the_original_chart(self):
        context,site=fixture()
        network=context.network.to_dict()
        paths={s.id:context[Key('path',s.id)].coefficients.copy() for s in context.network.sections}
        footprint=sample_footprint(site.key,context[footprint_key(site.key)])
        group=CompoundGroup((site.key,),(),(footprint,),RegionUnion().build((footprint.region,)))
        surface=CompoundSurfaceBuilder(context.patterns).build(group,context)
        metadata=surface.metadata
        recovery=metadata['deferred_projection_reconciliation']
        self.assertIn('interior depression',recovery['initial_projection_failure'])
        self.assertEqual(len(recovery['discarded_refinement_attempts']),3)
        self.assertEqual(metadata['refinements'],[])
        mesh=surface.mesh
        fixed=np.zeros(len(mesh.xyz),bool);fixed[list(surface.boundary_refs)]=True
        group=surface.chart.partition.group
        layouts={f.region.owner:f.layout for f in group.footprints}
        conditions=[BoundaryConditions().build(group,loop,layouts) for loop in group.loops]
        for node,ref in surface.boundary_refs.items():
            self.assertAlmostEqual(mesh.xyz[node,2],conditions[ref.loop].edge_samples(ref.edge,ref.fraction)[1],places=13)
        checks=CompoundSurfaceBuilder.projection_checks(mesh,fixed,conditions)
        self.assertEqual(checks['isolated_interior_minima'],0)
        self.assertLessEqual(checks['maximum_entrance_longitudinal_slope'],.03)
        self.assertLessEqual(BoundaryDrainage().fill(mesh,fixed)[1]['maximum_fill'],1e-6)
        stored=mesh.with_positions(mesh.xyz.astype(np.float32).astype(float))
        for candidate in (mesh,stored):
            q=quality(candidate)
            self.assertLessEqual(q['maximum_edge_ratio'],8.)
            self.assertGreaterEqual(q['minimum_corner_sine'],.03)
            self.assertLessEqual(float(candidate.quad_slopes().max()),metadata['measured_slope_limit']+1e-7)
            self.assertGreaterEqual(float(candidate.xyz[:,2].min()),surface.field.datum-1e-10)
            self.assertLessEqual(float(candidate.xyz[:,2].max()),surface.field.ceiling+1e-8)
        self.assertEqual(context.network.to_dict(),network)
        for sid,coefficients in paths.items():
            np.testing.assert_array_equal(context[Key('path',sid)].coefficients,coefficients)


if __name__=='__main__':unittest.main()
