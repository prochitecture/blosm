"""An interior diagonal reaches only its first unchanged physical curb."""
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.regions import BoundaryLoop
from tested_python_streets.partition import BoundaryRef
from tested_python_streets.computational_cut_extension import extend
from tested_python_streets.decomposition import VisibilitySplitter


def fixture():
    outer=BoundaryLoop(np.array(((0.,0.),(10.,0.),(10.,10.),(0.,10.))),())
    hole=BoundaryLoop(np.array(((4.,4.),(4.,6.),(6.,6.),(6.,4.))),())
    return SimpleNamespace(loops=(outer,hole)),{(li,e) for li in (0,1) for e in range(4)}


class ComputationalCutExtensionTests(unittest.TestCase):
    def test_semantic_internal_corners_unlock_a_bowed_strip_without_extra_budget(self):
        x=np.linspace(11.,.5,101)
        boundary=np.vstack((np.c_[x,3.+.046*(x-2.5)**2],(0.,0.),(12.,6.)))
        physical=frozenset(range(101));policy=VisibilitySplitter(max_candidates=3072,max_cuts=32,minimum_clearance_ratio=.01)
        with self.assertRaisesRegex(ValueError,'another compatible cut family'):
            policy.split(boundary,endpoints=physical,preferred=(25,))
        cuts=policy.split(boundary,semantic_endpoints=(101,102),preferred=(25,))
        self.assertTrue(cuts)
        self.assertLessEqual(len(cuts),32)
        self.assertTrue(any(a in (101,102) or b in (101,102) for a,b in cuts))
        with self.assertRaisesRegex(ValueError,'semantic endpoints'):
            policy.split(boundary,endpoints=physical,semantic_endpoints=(102,))

    def test_both_internal_ends_extend_to_first_curbs_without_moving_rim(self):
        g,edges=fixture();before=[loop.xy.copy() for loop in g.loops]
        refs,proof=extend(g,(2.,2.),(7.,2.),curb_edges=edges)
        self.assertEqual(refs,(BoundaryRef(0,3,.8),BoundaryRef(0,1,.2)))
        self.assertAlmostEqual(proof['extended_length'],10.)
        self.assertFalse(proof['physical_rim_changed'])
        for loop,old in zip(g.loops,before):np.testing.assert_array_equal(loop.xy,old)

    def test_existing_endpoint_is_exact_and_hole_is_first_actual_curb(self):
        g,edges=fixture();known=BoundaryRef(0,3,.5)
        refs,proof=extend(g,(0.,5.),(3.,5.),curb_edges=edges,first_ref=known)
        self.assertIs(refs[0],known)
        self.assertEqual(refs[1],BoundaryRef(1,0,.5))
        self.assertAlmostEqual(proof['last_extension_m'],1.)

    def test_reverse_and_world_transform_keep_same_first_physical_contacts(self):
        rotation=np.array(((math.cos(.713),-math.sin(.713)),(math.sin(.713),math.cos(.713))))
        def world(q):return np.asarray(q)@rotation.T+(1200.,-1800.)
        for reverse in (False,True):
            g,edges=fixture()
            g=SimpleNamespace(loops=tuple(BoundaryLoop(world(loop.xy),()) for loop in g.loops))
            first,last=world((2.,2.)),world((7.,2.))
            if reverse:first,last=last,first
            refs,proof=extend(g,first,last,curb_edges=edges)
            expected=(BoundaryRef(0,3,.8),BoundaryRef(0,1,.2))
            if reverse:expected=expected[::-1]
            self.assertEqual([(r.loop,r.edge) for r in refs],[(r.loop,r.edge) for r in expected])
            np.testing.assert_allclose([r.fraction for r in refs],[r.fraction for r in expected],atol=1e-12)
            self.assertAlmostEqual(proof['extended_length'],10.)

    def test_profile_hole_crossing_boundary_overlap_or_unbounded_input_rejects(self):
        g,edges=fixture()
        with self.assertRaisesRegex(ValueError,'profile'):
            extend(g,(2.,2.),(7.,2.),curb_edges=edges-{(0,1)})
        with self.assertRaisesRegex(ValueError,'crosses a physical hole'):
            extend(g,(2.,5.),(8.,5.),curb_edges=edges)
        with self.assertRaisesRegex(ValueError,'boundary interval'):
            extend(g,(2.,0.),(7.,0.),curb_edges=edges)
        with self.assertRaisesRegex(ValueError,'segment budget'):
            extend(g,(2.,2.),(7.,2.),curb_edges=edges,max_segments=4)
        with self.assertRaisesRegex(ValueError,'moved an existing'):
            extend(g,(1.,2.),(7.,2.),curb_edges=edges,first_ref=BoundaryRef(0,3,.8))


if __name__=='__main__':unittest.main()
