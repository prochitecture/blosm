"""Source-normal recovery needs real whole strips and preserves shared holes."""
from dataclasses import replace
import math
import unittest
import numpy as np
import _python_street_package
from test_python_street_source_cycles import triangle,prepare
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.physical_hole_recovery import build
from tested_python_streets.conditions import BoundaryConditions


class SourceNormalRecoveryTests(unittest.TestCase):
    def test_complete_sources_profiles_and_hole_edges_survive_rigid_transforms(self):
        for angle,offset in ((0.,(0.,0.)),(.47,(1200.,-1800.))):
            network=triangle(2.)
            rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
            network=replace(network,sections=tuple(replace(section,
                points=tuple((float(x),float(y),0.) for x,y in np.asarray(section.points)[:,:2]@rotation.T+offset))
                for section in network.sections))
            original=network.to_dict();context,group=prepare(network)
            result=build(CompoundQuadChart(context.patterns),group,context,[dict(reason='explicit recovery regression')])
            self.assertEqual(network.to_dict(),original)
            proof=result.metadata['source_normal_recovery']
            self.assertFalse(proof['physical_boundary_changed'])
            self.assertEqual(proof['final_width_proof']['maximum_error'],.002)
            self.assertEqual({row['section'] for row in proof['final_width_proof']['sections']},{1,2,3})
            self.assertEqual(len(result.holes),1)
            for index,hole in result.holes.items():
                np.testing.assert_array_equal(result.mesh.xyz[result.mesh.seams[index]][::-1],hole.mesh.xyz[hole.mesh.seams[0]])
            prepared=result.partition.group
            before=[BoundaryConditions().build(group,loop,{f.region.owner:f.layout for f in group.footprints}) for loop in group.loops]
            after=[BoundaryConditions().build(prepared,loop,{f.region.owner:f.layout for f in prepared.footprints}) for loop in prepared.loops]
            for first,last in zip(before,after):
                for port in first.ports:
                    other=next(p for p in last.ports if (p.section_id,p.end)==(port.section_id,port.end))
                    np.testing.assert_array_equal(first.xy[list(port.vertices)],last.xy[list(other.vertices)])
                    np.testing.assert_array_equal(first.heights[list(port.vertices)],last.heights[list(other.vertices)])

    def test_absent_claims_bad_source_coverage_and_cell_budget_do_not_qualify(self):
        context,group=prepare(triangle(2.));chart=CompoundQuadChart(context.patterns)
        with self.assertRaisesRegex(ValueError,'actual consumed'):
            build(chart,replace(group,claims=()),context,[])
        moved=replace(group,loops=tuple(replace(loop,xy=loop.xy+(1000.,1000.)) for loop in group.loops))
        with self.assertRaises(ValueError):build(chart,moved,context,[])
        with self.assertRaises(ValueError):build(CompoundQuadChart(context.patterns,max_cells=4),group,context,[])


if __name__=='__main__':unittest.main()
