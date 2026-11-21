"""Failed shape correction can use bounded energy seeds without topology edits."""
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from test_python_street_interior_quad_refinement import distorted_grid
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup,sample_footprint
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.four_side_charts import quality
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.mesh import Mesh,publication_coordinates
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.refined_grading_recovery import recover
from tested_python_streets.regions import RegionUnion


SOURCES=(
 (1,5.,((-733.142333984375,-616.5968017578125),(-731.7499389648438,-612.4558715820312))),
 (2,2.,((-733.3095092773438,-667.4586791992188),(-746.48388671875,-634.0504150390625),
        (-754.7672119140625,-627.325927734375),(-765.1515502929688,-621.6253051757812),
        (-841.5867309570312,-624.3112182617188),(-925.8790893554688,-628.1085815429688),
        (-1003.6723022460938,-631.1488647460938),(-1019.9674072265625,-633.6402587890625),
        (-1034.162353515625,-641.5309448242188),(-1043.7110595703125,-647.3182373046875))),
 (3,5.,((-733.3095092773438,-667.4586791992188),(-730.4793090820312,-632.3710327148438),
        (-733.142333984375,-616.5968017578125))),
 (4,5.,((-748.3426513671875,-690.400146484375),(-733.3095092773438,-667.4586791992188))),
 (5,5.,((-703.6390991210938,-717.683837890625),(-715.4287719726562,-694.9785766601562),
        (-733.3095092773438,-667.4586791992188))),
 (6,2.,((-714.4686889648438,-717.683837890625),(-721.049072265625,-702.892822265625),
        (-733.3095092773438,-667.4586791992188))))


def fixture():
    rise=float(np.float32(.08))
    sections=tuple(model.Section(sid,3 if sid==1 else sid,
        tuple((x,y,0.) for x,y in points),width=width,rise=rise,interpolation=3)
        for sid,width,points in SOURCES)
    site=model.Intersection(1,(model.End(3,0),model.End(2,0),model.End(4,1),model.End(6,1),model.End(5,1)),
        3.,(3.,1.,3.,1.,3.),('STYLE',)*5)
    network=model.Network('five_arm_refined_grading',sections,(site,),
        connections=(model.Connection(2,(model.End(3,1),model.End(1,0))),),
        settings=model.Settings(cell_size=float(np.float32(.7)),exaggeration=1.))
    context=footprint_context(network)
    layout=IntersectionProvider.footprint(context,site)
    source=sample_footprint(site.key,layout)
    return context,CompoundGroup((site.key,),(),(source,),RegionUnion().build((source.region,)))


class RefinedGradingRecoveryTests(unittest.TestCase):
    def test_energy_repair_retains_fixed_rim_topology_and_stored_limits(self):
        mesh=distorted_grid();before=mesh.xyz.copy();faces=mesh.quads.copy()
        output,report=recover(mesh,mesh.seams[0],100,dict(qualified=False))
        self.assertTrue(report['qualified'])
        self.assertTrue(report['refined_energy_recovery']['fixed_topology'])
        np.testing.assert_array_equal(output.quads,faces)
        np.testing.assert_array_equal(output.xyz[mesh.seams[0]],before[mesh.seams[0]])
        np.testing.assert_array_equal(mesh.xyz,before)
        self.assertLessEqual(report['refined_energy_recovery']['stored_quality']['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(report['refined_energy_recovery']['stored_quality']['minimum_corner_sine'],.03)

    def test_actual_five_arm_controls_pass_complete_height_and_boundary_guards(self):
        context,group=fixture();before=context.network.to_dict()
        result=CompoundSurfaceBuilder(context.patterns).build(group,context)
        repair=result.metadata['refined_grading']['planar']['refined_energy_recovery']
        self.assertFalse(repair['previous_correction']['qualified'])
        self.assertEqual(repair['energy_iterations'],500)
        self.assertEqual(context.network.to_dict(),before)
        prepared=result.chart.partition.group
        conditions=[BoundaryConditions().build(prepared,loop,{f.region.owner:f.layout for f in prepared.footprints})
                    for loop in prepared.loops]
        fixed=np.zeros(len(result.mesh.xyz),bool);fixed[list(result.boundary_refs)]=True
        for node,ref in result.boundary_refs.items():
            self.assertEqual(result.mesh.xyz[node,2],conditions[ref.loop].edge_samples(ref.edge,ref.fraction)[1])
        checks=CompoundSurfaceBuilder.projection_checks(result.mesh,fixed,conditions)
        self.assertEqual(checks['isolated_interior_minima'],0)
        self.assertLessEqual(checks['maximum_entrance_longitudinal_slope'],.03)
        stored=result.mesh.with_positions(publication_coordinates(result.mesh.xyz,1.,physical=True))
        self.assertLessEqual(quality(stored)['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(quality(stored)['minimum_corner_sine'],.03)
        self.assertLessEqual(float(stored.quad_slopes().max()),result.metadata['measured_slope_limit']+1e-7)

    def test_already_qualified_refined_stage_does_not_call_the_new_recovery(self):
        pattern=PatternLibrary().grid(2,2)
        mesh=Mesh(np.c_[pattern.uv,np.full(len(pattern.uv),.04)],pattern.quads,
                  seams={0:np.array((0,1,2,5,8,7,6,3))})
        class Field:
            datum=0.;ceiling=.08;metadata={}
            def __call__(self,xy):return np.full(len(xy),.04)
        refs={int(i):None for i in mesh.seams[0]}
        with patch('tested_python_streets.refined_grading_recovery.recover') as recovery:
            result,field,report=CompoundSurfaceBuilder(PatternLibrary()).correct_refined_grading(
                mesh,Field(),refs,[],.1,100)
        recovery.assert_not_called()
        np.testing.assert_array_equal(result.xyz,mesh.xyz)
        np.testing.assert_array_equal(result.quads,mesh.quads)

    def test_immutable_bad_quad_nonplanar_and_allocation_limits_still_reject(self):
        skinny=Mesh(np.array(((0.,0.,0.),(20.,0.,0.),(20.,1.,0.),(0.,1.,0.))),np.array(((0,1,2,3),)))
        with self.assertRaises(ValueError):recover(skinny,np.arange(4),100,{})
        mesh=distorted_grid()
        with self.assertRaises(ValueError):recover(mesh,mesh.seams[0],len(mesh.quads)-1,{})
        xyz=mesh.xyz.copy();xyz[0,2]=.1
        with self.assertRaises(ValueError):recover(mesh.with_positions(xyz),mesh.seams[0],100,{})


if __name__=='__main__':unittest.main()
