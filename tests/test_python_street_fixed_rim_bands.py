"""Standalone separator bands never insert a unilateral pavement seam sample."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.fixed_rim_bands import BandedSeparatorChart, GradedRimSeparatorChart
from tested_python_streets.separator_grading import LocalSeparatorBands
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.four_side_charts import quality
from tested_python_streets.validation import validate_overlaps
from tested_python_streets.model import Key


class FixedRimBandTests(unittest.TestCase):
    def test_cluster_at_rim_origin_retains_every_whole_edge(self):
        # Eight source intervals over 20 mm, between ordinary half-metre edges.
        x = np.r_[np.arange(0.,4.,.5),
            4.+np.array((-.01,-.0073,-.0037,-.0027,0.,.0027,.0037,.0073,.01)),
            np.arange(4.5,8.1,.5)]
        xy = np.vstack((np.c_[x,np.zeros_like(x)], ((8.,.5),(8.,1.),(8.,1.5)),
            np.c_[np.arange(8.,-.1,-.5),np.full(17,2.)], ((0.,1.5),(0.,1.),(0.,.5))))
        xy = np.roll(xy,-12,axis=0) @ np.array(((.8,.6),(-.6,.8))) + (423.,275.)
        result = LocalSeparatorBands(PatternLibrary(), bands=8, spacing=.7, max_cells=4000).build(xy)
        self.assertEqual(result.metadata['short_intervals'],8)
        np.testing.assert_array_equal(result.mesh.xyz[:len(xy),:2],xy)
        edges = np.sort(np.stack((result.mesh.quads,np.roll(result.mesh.quads,-1,axis=1)),
                                 axis=-1).reshape(-1,2),axis=1)
        edges, counts = np.unique(edges,axis=0,return_counts=True)
        rim = {tuple(edge) for edge in edges[counts==1]}
        self.assertEqual(rim,{tuple(sorted((i,(i+1)%len(xy)))) for i in range(len(xy))})

    def test_two_short_bank_edges_remain_shared_when_the_neighbouring_edges_are_long(self):
        bank = [(x, 0.) for x in (0., .5, 1., 1.494, 1.5, 1.506, 2., 2.5, 3.)]
        boundary = np.array(bank + [(3., .5), (3., 1.), (2.5, 1.), (2., 1.), (1.5, 1.),
                                    (1., 1.), (.5, 1.), (0., 1.), (0., .5)])
        for angle, shift in ((0., (0., 0.)), (.39, (975., -611.))):
            c, s = np.cos(angle), np.sin(angle)
            xy = boundary @ np.array(((c,s),(-s,c))) + shift
            chart = GradedRimSeparatorChart(PatternLibrary(), spacing=.7, max_cells=200, bands=4)
            result = chart.build(xy)
            self.assertEqual(result.metadata['spacing_jumps'], 2)
            self.assertEqual(result.metadata['rim_certificate']['added_boundary_samples'], 0)
            self.assertTrue(np.array_equal(result.mesh.xyz[result.mesh.seams[0], :2], xy))
            self.assertLessEqual(len(result.mesh.quads), 200)
            stored = result.mesh.with_positions(publication_coordinates(result.mesh.xyz, 1., physical=True))
            for mesh in (result.mesh, stored):
                self.assertLessEqual(quality(mesh)['maximum_edge_ratio'], 8.)
                self.assertGreater(quality(mesh)['minimum_corner_sine'], 0.)
                validate_overlaps({Key('grass', 0): mesh})
        with self.assertRaisesRegex(ValueError, 'cell budget'):
            GradedRimSeparatorChart(PatternLibrary(), max_cells=len(boundary)).build(boundary)

    def test_curved_banks_keep_exact_whole_edges_and_a_positive_quad_disk(self):
        x = np.linspace(0., 72., 65)
        center = .003*(x-36.)**2
        first = np.c_[x, center-2.]
        second = np.c_[x[::-1], (center+2.)[::-1]]
        boundary = np.vstack((first, second))
        output = BandedSeparatorChart(PatternLibrary(), bands=8, spacing=.7, max_cells=8000).build(boundary)
        mesh = output.mesh
        np.testing.assert_array_equal(mesh.xyz[:len(boundary), :2], boundary)
        self.assertTrue(output.metadata['exact_whole_boundary_edges'])
        self.assertEqual(output.metadata['rim_certificate']['added_boundary_samples'], 0)
        self.assertLessEqual(len(mesh.quads), 8000)
        self.assertEqual(mesh.material, 'grass')
        for candidate in (mesh, mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))):
            report = quality(candidate)
            self.assertLessEqual(report['maximum_edge_ratio'], 8.)
            self.assertGreater(report['minimum_corner_sine'], 0.)
        validate_overlaps({Key('separator', 10): mesh}, cell_size=2.)

    def test_bands_and_allocations_remain_bounded(self):
        with self.assertRaisesRegex(ValueError, 'sixty-four'):
            BandedSeparatorChart(PatternLibrary(), bands=65)
        with self.assertRaisesRegex(ValueError, 'Invalid bounded'):
            BandedSeparatorChart(PatternLibrary(), max_cells=50001).build(np.array(((0., 0.), (2., 0.), (2., 1.), (0., 1.))))


if __name__ == '__main__':
    unittest.main()
