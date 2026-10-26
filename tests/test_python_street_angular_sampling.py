"""Source-phase locality, physical curb error and bounded angular sampling."""
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import angular_dimensional_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_angular_boundaries import AngularMemberBoundary
from tested_python_streets.bundle_ends import PathRange
from tested_python_streets.paths import Path


class AngularSourceSamplingTests(unittest.TestCase):
    def test_trims_keep_remote_phase_exact_controls_and_original_physical_curve(self):
        for mode in (0,3):
            source=angular_dimensional_network(model,180.,False,mode)
            update=PythonMeshBackend(source.settings).prepare(source)
            plan=update.products[model.Key('bundle_plan',13)]
            original=update.products[model.Key('path',117)]
            collector=AngularMemberBoundary.__new__(AngularMemberBoundary)
            collector.context=SimpleNamespace(network=source)
            collector.frame=plan.frame;collector.exact_corners=[]
            traversal=model.Traversal(117,False)
            traces=[]
            for fraction in (.07,.14):
                first,last=original.length*fraction,original.length*.91
                samples=collector._path_boundary(PathRange(original,first,last),traversal,1.)
                phase=original.samples(0.,original.length,source.settings.spacing,
                                       tolerance=.1*source.settings.curve_tolerance)
                phase=np.r_[first,phase[(phase>first+1e-10)&(phase<last-1e-10)],last]
                expected=original.point(phase,original.section.offset-.5*original.section.width)[::-1]
                np.testing.assert_array_equal(samples,expected)
                traces.append(samples)
                if mode==0:
                    knots=original.stations[(original.stations>first)&(original.stations<last)]
                    for p in original.point(knots,original.section.offset-.5*original.section.width):
                        self.assertTrue(any(np.array_equal(p,q) for q in samples))
            # The second restriction retains every unchanged old phase row.
            self.assertTrue(all(any(np.array_equal(p,q) for q in traces[0]) for p in traces[1][1:-1]))

    def test_smooth_dimensional_islands_generate_from_float32_controls_for_every_sweep(self):
        for sweep in (180.,225.,270.):
            with self.subTest(sweep=sweep):
                source=angular_dimensional_network(model,sweep,True,3)
                source=replace(source,sections=tuple(replace(s,points=tuple(map(tuple,
                    np.asarray(s.points,dtype=np.float32).astype(float)))) for s in source.sections))
                before=source.to_dict()
                update=PythonMeshBackend(source.settings).prepare(source)
                plan=update.products[model.Key('bundle_plan',13)]
                self.assertEqual(type(plan).__name__,'AngularBundlePlan')
                self.assertEqual(len(plan.islands),1)
                self.assertEqual(len(plan.transition_owners),4)
                self.assertEqual(source.to_dict(),before)
                for key in (model.Key('bundle_pavement',13),model.Key('separator',13)):
                    mesh=update.products[key].mesh
                    self.assertLessEqual(mesh.validate()['maximum_edge_ratio'],8.)
                    self.assertLess(len(mesh.quads),10000)
                self.assertLessEqual(update.products[model.Key('bundle_pavement',13)].metadata['angular_longitudinal_bound'],.03+1e-8)

    def test_sample_budget_rejects_before_any_sampling_array_allocation(self):
        for mode in (0,3):
            section=model.Section(1,1,((0.,0.,0.),(1.,0.,0.),(2.,.1,0.)),width=1.,interpolation=mode)
            path=Path(section)
            with patch('tested_python_streets.paths.np.linspace',side_effect=AssertionError('Unexpected allocation')):
                with self.assertRaisesRegex(ValueError,'allocation budget'):
                    path.samples(0.,path.length,1e-9,max_samples=100)
            for limit in (0,2,3.5,True):
                with self.assertRaisesRegex(ValueError,'sample budget'):
                    path.samples(0.,path.length,.2,max_samples=limit)


if __name__=='__main__': unittest.main()
