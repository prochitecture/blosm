"""Full-corner equivalence and allocation bounds for imported street inspection."""
import importlib.util
from pathlib import Path
import tracemalloc
import unittest

import numpy as np

spec=importlib.util.spec_from_file_location('surface_validation',
    Path(__file__).resolve().parents[1]/'way/street_network/surface_validation.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
check=module.minimum_positive_quad_corner


def quads(count):
    rng=np.random.default_rng(70419)
    origin=rng.uniform(-1000.,1000.,(count,1,3))
    width=rng.uniform(.01,2.,(count,1))
    height=rng.uniform(.01,2.,(count,1))
    points=np.broadcast_to(origin,(count,4,3)).copy()
    points[:,1:3,0]+=width
    points[:,2:,1]+=height
    return points.reshape(-1,3),np.arange(count*4,dtype=np.int32).reshape(-1,4)


class SurfaceValidationTests(unittest.TestCase):
    def test_complete_minimum_is_exact_across_partial_batches(self):
        xyz,faces=quads(657)
        edges=np.roll(xyz[faces,:2],-1,axis=1)-xyz[faces,:2]
        following=np.roll(edges,-1,axis=1)
        expected=float(np.min(edges[:,:,0]*following[:,:,1]-edges[:,:,1]*following[:,:,0]))
        original=xyz.copy(),faces.copy()
        for size in (1,17,128,16384):
            self.assertEqual(check(xyz,faces,batch_size=size),expected)
        np.testing.assert_array_equal(xyz,original[0])
        np.testing.assert_array_equal(faces,original[1])

    def test_folded_and_collapsed_faces_are_not_lost_at_batch_boundaries(self):
        xyz,faces=quads(35)
        for index in (0,15,16,17,32,34):
            for kind in ('reversed','collapsed','concave'):
                with self.subTest(face=index,kind=kind):
                    changed=xyz.copy();ids=faces[index]
                    if kind=='reversed':changed[ids]=changed[ids[::-1]]
                    elif kind=='collapsed':changed[ids[2]]=changed[ids[1]]
                    else:changed[ids[2]]=changed[ids[0]]-.25*(changed[ids[2]]-changed[ids[0]])
                    with self.assertRaisesRegex(ValueError,'folded or collapsed'):
                        check(changed,faces,batch_size=16)

    def test_nonfinite_coordinates_and_invalid_inputs_reject(self):
        xyz,faces=quads(3)
        for value in (float('nan'),float('inf'),-float('inf')):
            changed=xyz.copy();changed[-1,2]=value
            with self.assertRaisesRegex(ValueError,'non-finite'):
                check(changed,faces)
        for size in (0,-1,1.5,True):
            with self.assertRaises(ValueError):check(xyz,faces,batch_size=size)
        with self.assertRaises(ValueError):check(xyz,faces[:0])
        with self.assertRaises(ValueError):check(xyz,faces[:,:3])

    def test_large_inspection_bounds_temporary_allocation(self):
        xyz,faces=quads(230000)
        tracemalloc.start()
        try:
            self.assertGreater(check(xyz,faces),0.)
            _,peak=tracemalloc.get_traced_memory()
        finally:tracemalloc.stop()
        # The former whole-array edge/corner expressions require over 50 MiB
        # on this count. Input arrays are allocated before measurement.
        self.assertLess(peak,8*1024*1024)


if __name__=='__main__':unittest.main()
