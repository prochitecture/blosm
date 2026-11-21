"""A moderate-aspect four-arm fork needs cuts behind both rounded noses."""
import unittest
from dataclasses import replace
import numpy as np
import _python_street_package
from tested_python_streets.model import Network,Section,Intersection,Connection,End,Key,Settings
from tested_python_streets.paths import Path,unit
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.intersections import Approach,IntersectionProvider
from tested_python_streets.compounds import CompoundGroup,sample_footprint
from tested_python_streets.regions import RegionUnion
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.longitudinal_cuts import LongitudinalCuts
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.four_side_charts import quality

# Original source controls, with local identities. Shared endpoint tangents
# follow the original explicit two-ended seam; no saved preparation is read.
SOURCES=((1,
  3.5,
  1,
  ((2232.87255859375, -511.6051330566406),
   (2211.25927734375, -512.4949951171875),
   (2130.62109375, -532.9673461914062))),
 (2,
  10.800000190734863,
  -1,
  ((2130.62109375, -532.9673461914062),
   (1861.2762451171875, -540.4628295898438),
   (1638.1763916015625, -551.0003051757812))),
 (3,
  10.800000190734863,
  -1,
  ((2232.529296875, -530.288818359375),
   (2194.79931640625, -531.12353515625),
   (2130.62109375, -532.9673461914062))),
 (4,
  5.400000095367432,
  -1,
  ((2180.05322265625, -509.8879699707031),
   (2148.98681640625, -518.8802490234375),
   (2141.822998046875, -523.4797973632812),
   (2130.62109375, -532.9673461914062))),
 (5,
  5.400000095367432,
  -1,
  ((2232.87255859375, -502.2189636230469),
   (2206.61083984375, -503.4349365234375),
   (2180.05322265625, -509.8879699707031))))


def fixture():
    rise=float(np.float32(.08))
    sections=tuple(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width,rise=rise)
                   for sid,width,half,points in SOURCES)
    site=Intersection(1,(End(3,1),End(1,1),End(4,1),End(2,0)),4.,(4.,10.,2.,4.),('STYLE',)*4)
    seam=Connection(1,(End(4,0),End(5,1)))
    n=Network('two_rounded_fork_noses',sections,(site,),connections=(seam,),settings=Settings(exaggeration=1.))
    a,b=np.asarray(n.section_index[4].points),np.asarray(n.section_index[5].points)
    tangent=unit(unit(a[1,:2]-a[0,:2])-unit(b[-2,:2]-b[-1,:2]))
    p={}
    for section,(_,_,half,_) in zip(sections,SOURCES):
        directions={0:tangent} if section.id==4 else {1:tangent} if section.id==5 else None
        p[Key('path',section.id)]=Path(section,endpoint_directions=directions)
        p[Key('profile',section.id)]=Crown(rise,0.,half)
    context=BuildContext(n,p,PatternLibrary());topology=RouteTopology(n)
    approaches=[Approach(Route.build(context,topology.trace(e)),e.end,p[Key('profile',e.section_id)]) for e in site.ends]
    layout=IntersectionProvider.curved_boundary(context,site,approaches)
    p[site.layout_key]=replace(layout,mesh_strategy='boundary_chart')
    footprint=sample_footprint(site.key,layout)
    group=CompoundGroup((site.key,),(),(footprint,),RegionUnion().build((footprint.region,)))
    return context,group


class MultiNoseLongitudinalTests(unittest.TestCase):
    def test_all_nose_cuts_qualify_without_changing_source_curbs_or_profiles(self):
        context,group=fixture();original=context.network.to_dict()
        self.assertFalse(LongitudinalCuts.eligible(group,context.network.settings.spacing))
        with self.assertRaisesRegex(ValueError,'not an elongated'):
            LongitudinalCuts(context.network.settings.spacing).cuts(group,context)
        one=LongitudinalCuts(context.network.settings.spacing,minimum_aspect=2.).cuts(group,context)
        all_noses=LongitudinalCuts(context.network.settings.spacing,minimum_aspect=2.,all_noses=True).cuts(group,context)
        self.assertGreater(len(all_noses),len(one))
        surface=CompoundSurfaceBuilder(context.patterns).build(group,context)
        policy=surface.metadata['planar_chart']['fixed_rim_compound_recovery']
        self.assertEqual(policy['policy'],'ALL_NOSES_LONGITUDINAL')
        self.assertFalse(policy['physical_boundary_changed'])
        self.assertLess(len(surface.mesh.quads),6000)
        fixed=np.zeros(len(surface.mesh.xyz),bool);fixed[list(surface.boundary_refs)]=True
        prepared=surface.chart.partition.group
        conditions=[BoundaryConditions().build(prepared,loop,{f.region.owner:f.layout for f in prepared.footprints})
                    for loop in prepared.loops]
        for node,ref in surface.boundary_refs.items():
            self.assertAlmostEqual(surface.mesh.xyz[node,2],conditions[ref.loop].edge_samples(ref.edge,ref.fraction)[1],places=13)
        for mesh in (surface.mesh,surface.mesh.with_positions(surface.mesh.xyz.astype(np.float32).astype(float))):
            q=quality(mesh)
            self.assertLessEqual(q['maximum_edge_ratio'],8.)
            self.assertGreaterEqual(q['minimum_corner_sine'],.03)
            self.assertLessEqual(float(mesh.quad_slopes().max()),surface.metadata['measured_slope_limit']+1e-7)
        checks=CompoundSurfaceBuilder.projection_checks(surface.mesh,fixed,conditions)
        self.assertEqual(checks['isolated_interior_minima'],0)
        self.assertLessEqual(checks['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(context.network.to_dict(),original)


if __name__=='__main__':unittest.main()
