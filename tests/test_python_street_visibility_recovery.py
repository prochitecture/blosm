"""A source-owned two-junction compound can need cuts beyond its owner seam."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compounds import CompoundPlanner
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.four_side_charts import quality
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.paths import Path
from tested_python_streets.quad_charts import VisibilityKernelError
from tested_python_streets.routes import RouteTopology, Route, route_key
from tested_python_streets.visibility_recovery import source_visibility
from test_python_street_compounds import footprint_context


CONTROLS=(
 (1,3.5,((-1323.4530029296875,-236.1239471435547),(-1314.88330078125,-250.16314697265625),
         (-1310.024169921875,-263.6672058105469),(-1303.5869140625,-282.2366638183594))),
 (2,3.5,((-1303.5869140625,-282.2366638183594),(-1301.6533203125,-298.5787658691406))),
 (3,3.5,((-1301.6533203125,-298.5787658691406),(-1301.2039794921875,-307.6180114746094),
         (-1303.5755615234375,-339.45489501953125))),
 (4,2.,((-1301.6533203125,-298.5787658691406),(-1293.5546875,-296.6768798828125))),
 (5,2.,((-1293.9801025390625,-285.0439147949219),(-1296.4149169921875,-284.3309631347656),
        (-1303.5869140625,-282.2366638183594))),
 (6,2.,((-1293.5546875,-296.6768798828125),(-1288.3348388671875,-296.789306640625),
        (-1284.37548828125,-297.88104248046875),(-1273.835693359375,-303.527099609375))),
 (7,2.,((-1222.8424072265625,-301.67828369140625),(-1240.2440185546875,-298.9363708496094),
        (-1256.6964111328125,-295.1816101074219),(-1278.0025634765625,-289.70037841796875),
        (-1293.9801025390625,-285.0439147949219))))


def fixture(angle=0., offset=(0.,0.)):
    rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
    rise=float(np.float32(.08))
    sections=tuple(model.Section(sid, {6:4,7:5}.get(sid,sid),
        tuple((float(x),float(y),0.) for x,y in np.asarray(points)@rotation.T+offset),
        width=width,rise=rise,interpolation=3) for sid,width,points in CONTROLS)
    sites=(model.Intersection(1,(model.End(1,1),model.End(2,0),model.End(5,1)),10.,(10.,10.,1.),('STYLE',)*3),
           model.Intersection(2,(model.End(2,1),model.End(3,0),model.End(4,0)),10.,(10.,10.,1.),('STYLE',)*3))
    connections=(model.Connection(3,(model.End(4,1),model.End(6,0))),
                 model.Connection(4,(model.End(7,1),model.End(5,0))))
    network=model.Network('curved-two-junction-kernel',sections,sites,connections=connections,
                          settings=model.Settings(exaggeration=1.))
    context=footprint_context(network)
    products=context.products
    topology=RouteTopology(network)
    for site in sites:
        for end in site.ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    layouts={site.key:IntersectionProvider.footprint(context,site) for site in sites}
    group,=CompoundPlanner().build(layouts)
    return context,group


class VisibilityRecoveryTests(unittest.TestCase):
    def test_original_profiles_source_width_and_stored_surface_survive_missing_kernel(self):
        context,group=fixture()
        before=context.network.to_dict()
        chart=CompoundQuadChart(context.patterns)
        with self.assertRaises(VisibilityKernelError):chart._build(group,context,require_quality=True)
        original=BoundaryConditions().build(group,group.loops[0],{f.region.owner:f.layout for f in group.footprints})
        surface=CompoundSurfaceBuilder(context.patterns).build(group,context)
        metadata=surface.metadata['planar_chart']['visibility_kernel_recovery']
        self.assertEqual(metadata['final_width_proof']['maximum_error'],.002)
        self.assertFalse(metadata['physical_boundary_changed'])
        self.assertIsNotNone(surface.metadata['source_cycle_reconciliation'])
        prepared=surface.chart.partition.group
        result=BoundaryConditions().build(prepared,prepared.loops[0],{f.region.owner:f.layout for f in prepared.footprints})
        for port in original.ports:
            other=next(p for p in result.ports if (p.section_id,p.end)==(port.section_id,port.end))
            np.testing.assert_array_equal(original.xy[list(port.vertices)],result.xy[list(other.vertices)])
            np.testing.assert_array_equal(original.heights[list(port.vertices)],result.heights[list(other.vertices)])
        stored=surface.mesh.with_positions(publication_coordinates(surface.mesh.xyz,1.,physical=True))
        self.assertLessEqual(quality(stored)['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(quality(stored)['minimum_corner_sine'],.03)
        fixed=np.zeros(len(surface.mesh.xyz),bool);fixed[list(surface.boundary_refs)]=True
        checks=CompoundSurfaceBuilder.projection_checks(surface.mesh,fixed,[result])
        self.assertEqual(checks['isolated_interior_minima'],0)
        self.assertLessEqual(checks['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(context.network.to_dict(),before)

    def test_rotated_source_chart_retains_claims_rims_and_original_metric_proof(self):
        context,group=fixture(.47,(1500.,-2100.))
        result=source_visibility(CompoundQuadChart(context.patterns),group,context,[])
        self.assertEqual(result.partition.group.members,group.members)
        self.assertEqual(result.partition.group.claims,group.claims)
        proof=result.metadata['source_width_visibility']
        self.assertEqual([row['section'] for row in proof['sections']],[2])
        self.assertEqual(proof['maximum_error'],.002)
        self.assertLessEqual(result.metadata['fixed_rim_compound_recovery']['stored_quality']['maximum_edge_ratio'],8.)

    def test_clipped_source_and_absent_claims_cannot_use_visibility_recovery(self):
        context,group=fixture();chart=CompoundQuadChart(context.patterns)
        with self.assertRaisesRegex(ValueError,'multi-junction'):
            source_visibility(chart,replace(group,claims=()),context,[])
        original=context[model.Key('path',2)].section
        context.products[model.Key('path',2)]=Path(replace(original,
            points=tuple((x+100.,y,z) for x,y,z in original.points)))
        with self.assertRaisesRegex(ValueError,'complete source strip'):
            source_visibility(chart,group,context,[])


if __name__=='__main__':unittest.main()
