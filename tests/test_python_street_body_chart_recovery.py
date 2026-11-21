"""Uneven fixed banks must retain their physical edges in a direct-quad body."""
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from tested_python_streets.body_chart_recovery import FixedRimBodyRecovery
from tested_python_streets.four_side_charts import quality
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.patterns import PatternLibrary


class Frame:
    def __init__(self, angle=0., origin=(0.,0.)):
        c,s=np.cos(angle),np.sin(angle)
        self.rotation=np.array(((c,-s),(s,c)))
        self.origin=np.asarray(origin)
    def world(self,xy):return np.asarray(xy)@self.rotation.T+self.origin
    def local(self,xy):return (np.asarray(xy)-self.origin)@self.rotation
    def longitudinal_bound(self,mesh):
        local=mesh.xyz.copy();local[:,:2]=self.local(local[:,:2])
        # This upper bound includes the transverse gradient; stricter than
        # the longitudinal-only bound for this simple analytical fixture.
        return mesh.with_positions(local).slope_statistics()['maximum']


def inputs(angle=0.,origin=(0.,0.)):
    frame=Frame(angle,origin);length,width=28.,7.
    y=np.linspace(width,0.,9)
    curb=np.c_[np.linspace(0.,length,25),np.zeros(25)]
    physical=np.vstack((np.c_[np.zeros(9),y],curb[1:],np.c_[np.full(8,length),y[-2::-1]]))
    kinds=[('cut','head')]*8+[('outer',0)]*24+[('cut','tail')]*8
    keys=[('physical',i) for i in range(len(physical))]
    stations=np.r_[0.,.018,.055,.16,.45,1.1,np.linspace(2.,26.,19),26.9,27.55,27.84,27.945,27.982,28.]
    # Semantic parity agrees without changing either prescribed profile.
    if (len(stations)-1)%2:stations=np.insert(stations,-1,27.991)
    seam=frame.world(np.c_[stations,np.full(len(stations),width)])
    shared_keys=[('bank',i) for i in range(len(seam))]
    field=SimpleNamespace(datum=0.,ceiling=.08)
    field.values=lambda xy:.08*frame.local(xy)[:,1]/width
    builder=SimpleNamespace(patterns=PatternLibrary(),spacing=.7,max_cells=4000,
        maximum_slope=.1,maximum_longitudinal_slope=.03)
    return builder,frame,field,(frame.world(physical),kinds,keys),(seam,shared_keys)


class BodyChartRecoveryTests(unittest.TestCase):
    def test_uneven_bank_and_full_profiles_are_exact_at_storage_precision(self):
        for angle,origin in ((0.,(0.,0.)),(.47,(300.,-700.))):
            with self.subTest(angle=angle):
                builder,frame,field,outer,seam=inputs(angle,origin)
                helper=FixedRimBodyRecovery(builder,frame,field)
                rim,expected_keys,values,_=helper.boundary(outer,seam,field)
                mesh,keys,report=helper.build(outer,seam,'test body')
                np.testing.assert_array_equal(mesh.xyz[:len(rim)],np.c_[rim,values])
                self.assertEqual(keys[:len(rim)],expected_keys)
                self.assertEqual(len(keys),len(set(keys)))
                self.assertLessEqual(len(mesh.quads),builder.max_cells)
                for candidate in (mesh,mesh.with_positions(publication_coordinates(mesh.xyz,1.,physical=True))):
                    measured=quality(candidate)
                    self.assertLessEqual(measured['maximum_edge_ratio'],8.)
                    self.assertGreaterEqual(measured['minimum_corner_sine'],.03)
                self.assertEqual(report['drainage']['maximum_fill'],0.)
                self.assertTrue(report['all_direct_quads'])

    def test_cell_budget_and_uncoordinated_profile_parity_are_rejected(self):
        builder,frame,field,outer,seam=inputs()
        builder.max_cells=4
        with self.assertRaises(ValueError):FixedRimBodyRecovery(builder,frame,field).build(outer,seam,'bounded')
        builder.max_cells=4000
        wrong=seam[0][:-1],seam[1][:-1]
        with self.assertRaises(ValueError):FixedRimBodyRecovery(builder,frame,field).build(outer,wrong,'unmatched')
        interior=np.insert(seam[0],2,.5*(seam[0][1]+seam[0][2]),axis=0)
        keys=seam[1][:2]+[('extra',)]+seam[1][2:]
        with self.assertRaises(ValueError):FixedRimBodyRecovery(builder,frame,field).build(outer,(interior,keys),'parity')


if __name__=='__main__':unittest.main()
