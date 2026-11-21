"""Typed chart recovery preserves a source-complete original compound."""
import unittest
from unittest.mock import patch
from dataclasses import replace
import numpy as np
import _python_street_package
from _python_street_fixtures import short_junction_pair
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import PLAN
from tested_python_streets.model import Network,Section,Intersection,End,Key,Settings,Bundle,BundleMember,Traversal
from tested_python_streets.paths import Path
from tested_python_streets.routes import RouteTopology,Route,route_key
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.quad_charts import QuadChartCoreError
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.four_side_charts import quality
from tested_python_streets.cycle_dimensional_retry import certify_layout
from tested_python_streets.source_strip_containment import SourceStripContainmentError

# Original curved source controls and physical half-profiles. Local IDs avoid
# any dependence on an OSM extract, installed asset or saved preparation.
SOURCES=((1,
  14.0,
  0,
  ((-1707.9224853515625, 756.8281860351562),
   (-1711.3333740234375, 741.6785278320312),
   (-1710.27099609375, 729.8672485351562),
   (-1707.3258056640625, 720.3375244140625),
   (-1701.363037109375, 709.704833984375),
   (-1694.7850341796875, 701.8995361328125),
   (-1687.51416015625, 696.0198974609375),
   (-1682.023681640625, 692.9794311523438))),
 (2,
  14.0,
  0,
  ((-1643.685546875, 694.3049926757812),
   (-1632.8416748046875, 699.4451293945312),
   (-1625.8406982421875, 702.771728515625),
   (-1615.3231201171875, 709.4592895507812))),
 (3,
  10.800000190734863,
  0,
  ((-1682.023681640625, 692.9794311523438), (-1677.82421875, 690.662841796875))),
 (4,
  3.5,
  0,
  ((-1779.298095703125, 629.0199584960938),
   (-1748.3897705078125, 650.2955322265625),
   (-1725.1060791015625, 664.9052734375),
   (-1702.7252197265625, 677.5450439453125),
   (-1686.979248046875, 689.2849731445312),
   (-1682.023681640625, 692.9794311523438))),
 (5,
  14.0,
  0,
  ((-1677.82421875, 690.662841796875),
   (-1659.118408203125, 689.3776245117188),
   (-1643.685546875, 694.3049926757812))),
 (6,
  7.0,
  -1,
  ((-1677.82421875, 690.662841796875),
   (-1665.921875, 668.8076171875),
   (-1660.44384765625, 648.5682983398438),
   (-1652.8873291015625, 516.8753051757812))),
 (7,
  7.0,
  -1,
  ((-1642.260498046875, 523.6852416992188),
   (-1650.0631103515625, 666.6771850585938),
   (-1648.737548828125, 679.04443359375),
   (-1645.4677734375, 690.2423095703125),
   (-1643.685546875, 694.3049926757812))))
SITES=((1, ((1, 1), (4, 1), (3, 0)), 10.0, (10.0, 10.0, 4.0)),
 (2, ((2, 0), (5, 1), (7, 1)), 10.0, (10.0, 10.0, 10.0)),
 (3, ((3, 1), (6, 0), (5, 0)), 4.0, (4.0, 10.0, 10.0)))


def _legacy_fixture():
    rise=float(np.float32(.08))
    sections=tuple(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width,rise=rise)
                   for sid,width,half,points in SOURCES)
    sites=tuple(Intersection(sid,tuple(End(*e) for e in ends),radius,radii,('STYLE',)*len(ends))
                for sid,ends,radius,radii in SITES)
    network=Network('three_physical_junctions',sections,sites,
        bundles=(Bundle(1,(BundleMember((Traversal(6,False),)),BundleMember((Traversal(7,True),)))),),
        settings=Settings(exaggeration=1.))
    products={}
    for section,(_,_,half,_) in zip(sections,SOURCES):
        products[Key('path',section.id)]=Path(section)
        products[Key('profile',section.id)]=Crown(rise,0.,half)
    context=BuildContext(network,products,PatternLibrary())
    topology=RouteTopology(network)
    for site in sites:
        for end in site.ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    layouts={site.key:IntersectionProvider.footprint(context,site) for site in sites}
    group=DeferredCompoundPlanner(context,sites,layouts).exterior(tuple(site.key for site in sites))
    return context,group


def fixture():
    # The previous three-site fixture's declared 14 m internal strip escaped
    # its provisional exterior by 47 mm. Keep that original input below as a
    # rejection regression. This maintained two-owner graph retains all its
    # original dimensions and complete source prefixes. Only a typed primary
    # chart rejection is injected; fallback charts and height remain real.
    network=short_junction_pair(model,length=2.527)
    network=replace(network,settings=replace(network.settings,exaggeration=1.))
    engine=PythonMeshBackend(network.settings)
    update=engine.prepare(network)
    return BuildContext(network,update.products,engine.patterns),update.products[PLAN].groups[0]


class CompoundChartRecoveryTests(unittest.TestCase):
    def test_typed_core_failure_keeps_complete_compound_rims_profiles_and_sources(self):
        context,group=fixture()
        original=context.network.to_dict()
        consumed={claim.section_id for claim in group.claims}
        before=certify_layout(context,group,consumed)
        self.assertEqual({row['section'] for row in before['sections']},{1,2,3,4,5})
        build=CompoundQuadChart._build_steps;attempts=[]
        def primary_failure(chart,prepared,ctx,*,require_quality):
            if not chart.four_sided:
                attempts.append('typed_primary_core')
                raise QuadChartCoreError('No positive six-vertex direct quad core')
            attempts.append('real_fixed_rim_chart')
            return build(chart,prepared,ctx,require_quality=require_quality)
        with patch.object(CompoundQuadChart,'_build_steps',primary_failure):
            with patch.object(CompoundQuadChart,'_reflex',
                    side_effect=QuadChartCoreError('Local radial core unavailable')):
                surface=CompoundSurfaceBuilder(context.patterns).build(group,context)
        self.assertEqual(attempts,['typed_primary_core','real_fixed_rim_chart'])
        self.assertIn('fixed_rim_compound_recovery',surface.metadata['planar_chart'])
        self.assertLess(len(surface.mesh.quads),8000)
        for mesh in (surface.mesh,surface.mesh.with_positions(surface.mesh.xyz.astype(np.float32).astype(float))):
            result=quality(mesh)
            self.assertLessEqual(result['maximum_edge_ratio'],8.)
            self.assertGreaterEqual(result['minimum_corner_sine'],.03)
            self.assertLessEqual(float(mesh.quad_slopes().max()),surface.metadata['measured_slope_limit']+1e-7)
        prepared=surface.chart.partition.group
        self.assertEqual(certify_layout(context,prepared,consumed)['maximum_error'],before['maximum_error'])
        conditions=[BoundaryConditions().build(prepared,loop,{f.region.owner:f.layout for f in prepared.footprints})
                    for loop in prepared.loops]
        fixed=np.zeros(len(surface.mesh.xyz),bool);fixed[list(surface.boundary_refs)]=True
        for node,ref in surface.boundary_refs.items():
            xy,height,_=conditions[ref.loop].edge_samples(ref.edge,ref.fraction)
            np.testing.assert_allclose(surface.mesh.xyz[node,:2],xy,rtol=0.,atol=1e-12)
            self.assertAlmostEqual(surface.mesh.xyz[node,2],height,places=13)
        checks=CompoundSurfaceBuilder.projection_checks(surface.mesh,fixed,conditions)
        self.assertEqual(checks['isolated_interior_minima'],0)
        self.assertLessEqual(checks['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(context.network.to_dict(),original)

    def test_original_three_site_width_escape_cannot_be_admitted(self):
        # Stop before admission, retaining the exact provisional geometry;
        # no proof result is fabricated and no source guard is bypassed.
        class CandidateCaptured(RuntimeError):pass
        captured=[]
        def capture(context,group,consumed):
            captured.append((context,group,consumed))
            raise CandidateCaptured()
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',side_effect=capture):
            with self.assertRaises(CandidateCaptured):_legacy_fixture()
        self.assertEqual(len(captured),1)
        context,group,consumed=captured[0]
        self.assertEqual(context[Key('path',5)].section.width,14.)
        with self.assertRaisesRegex(SourceStripContainmentError,'complete source strip of Section 5'):
            certify_layout(context,group,consumed)


if __name__=='__main__':unittest.main()
