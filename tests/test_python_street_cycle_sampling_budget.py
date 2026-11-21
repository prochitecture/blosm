"""Complete future grids and rendering-only rows share a bounded source phase."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Section
from tested_python_streets.paths import Path
from tested_python_streets.compound_cycles import _canonical_samples


def previous(paths,cuts,spacing,tolerance):
    result={}
    for sid,path in paths.items():
        compulsory=[]
        for value in sorted({0.,path.length,*cuts.get(sid,())}):
            if not compulsory or value-compulsory[-1]>1e-9:compulsory.append(value)
            elif value==path.length:compulsory[-1]=value
        values=[]
        for first,last in zip(compulsory,compulsory[1:]):
            if last<=first+1e-10:continue
            values.extend(path.samples(first,last,min(.5,spacing),tolerance=tolerance)[:-1])
        values=np.r_[values,path.length]
        result[sid]=values[np.r_[True,np.diff(values)>1e-9]]
    return result


class CycleSamplingBudgetTests(unittest.TestCase):
    def test_short_bend_preserves_long_straight_grid_under_shared_budget(self):
        points = ((0,0,0), (2,.15,0), (4,.5,0), (6,.5,0), (8,.5,0), (160,.5,0))
        for controls in (points, points[::-1]):
            path = Path(Section(1, 1, controls, width=2., interpolation=3))
            cuts = {1: {path.length/3}}
            self.assertGreater(len(previous({1:path}, cuts, .5, .00025)[1]), 500)
            samples = _canonical_samples({1:path}, cuts, .5, .00025, max_samples=500)
            stations = samples[1]
            self.assertTrue(samples.complete_source_ranges)
            self.assertLessEqual(len(stations), 500)
            self.assertEqual((stations[0], stations[-1]), (0., path.length))
            self.assertIn(path.length/3, stations)
            self.assertTrue(np.all(np.diff(stations) > 0.))
            for offset in (-1., 1.):
                edge = path.point(stations, offset)
                chords = np.diff(edge, axis=0)
                for t in (.25, .5, .75):
                    actual = path.point(stations[:-1] + t*np.diff(stations), offset)
                    # Offset length varies with curvature along the centerline.
                    # Measure the spatial distance to the finite curb chord.
                    u = np.sum((actual-edge[:-1])*chords, axis=1) / np.sum(chords*chords, axis=1)
                    chord = edge[:-1] + np.clip(u, 0., 1.)[:, None]*chords
                    self.assertLessEqual(np.linalg.norm(actual-chord, axis=1).max(), .00025)

    def test_owned_head_tail_and_compulsory_contacts_are_byte_exact(self):
        for mode in (0,3):
            path=Path(Section(1,1,((0.,0.,0.),(10.,1.,0.),(30.,8.,0.),(45.,8.,0.)),width=2.,interpolation=mode))
            first,last=3.45,path.length-4.25
            cuts={1:{.8,first,path.length/2.,last,path.length-1.1}}
            reference=previous({1:path},cuts,.5,.00025)[1]
            actual=_canonical_samples({1:path},cuts,.5,.00025,intervals={1:((0.,first),(last,path.length))})[1]
            for low,high in ((0.,first),(last,path.length)):
                np.testing.assert_array_equal(actual[(actual>=low)&(actual<=high)],reference[(reference>=low)&(reference<=high)])
            self.assertTrue(cuts[1]<={float(v) for v in actual})
            full=_canonical_samples({1:path},cuts,.5,.00025)[1]
            np.testing.assert_array_equal(full,reference)

    def test_unused_ten_kilometre_tail_never_reaches_the_sampler(self):
        path=Path(Section(1,1,((0.,0.,0.),(10001.,0.,0.)),width=2.,interpolation=0))
        controls=path.points.copy()
        with patch.object(path,'samples',wraps=path.samples) as sampled:
            actual=_canonical_samples({1:path},{1:{1.}},.5,.00025,
                intervals={1:((0.,1.),)},max_samples=4)[1]
        np.testing.assert_array_equal(actual,(0.,.5,1.,path.length))
        self.assertEqual(sampled.call_count,1)
        self.assertEqual(sampled.call_args.args[:2],(0.,1.))
        self.assertEqual(sampled.call_args.kwargs['max_samples'],4)
        np.testing.assert_array_equal(path.points,controls)

    def test_sparse_rendering_grid_cannot_certify_future_mouths_or_source_cells(self):
        from test_python_street_source_mouth_clearance import fixture
        from tested_python_streets.model import Key
        from tested_python_streets.finite_profile_search import _stations,advance
        from tested_python_streets.source_mouth_clearance import advance_past_overlap
        for reverse,angle in ((False,0.),(True,.713)):
            with self.subTest(reverse=reverse,angle=angle):
                planner,arms,outward,ports,sources,_,members,_=fixture(reverse=reverse,angle=angle)
                paths={sid:planner.context[Key('path',sid)] for sid in (1,2)}
                cuts={};intervals={}
                for arm,station in zip(arms,outward):
                    for (sid,end),at in arm.path.trims(station).items():
                        cuts.setdefault(sid,set()).add(at)
                        intervals.setdefault(sid,[]).append((0.,at) if end==0 else (at,paths[sid].length))
                full=_canonical_samples(paths,cuts,.5,.00025)
                owned=_canonical_samples(paths,cuts,.5,.00025,intervals=intervals)
                self.assertTrue(full.complete_source_ranges)
                self.assertFalse(owned.complete_source_ranges)
                # Rendering is unchanged; only the future proof grid differs.
                for sid,ranges in intervals.items():
                    for first,last in ranges:
                        np.testing.assert_array_equal(
                            owned[sid][(owned[sid]>=first)&(owned[sid]<=last)],
                            full[sid][(full[sid]>=first)&(full[sid]<=last)])
                self.assertGreater(len(_stations(planner.context,arms[0],full)),100)
                for call in (
                    lambda:_stations(planner.context,arms[0],owned),
                    lambda:advance(planner,arms,outward,(0,),ports,sources,owned,members),
                    lambda:advance_past_overlap(planner,arms,outward,ports,sources,owned,members)):
                    with self.assertRaisesRegex(ValueError,'complete certified'):
                        call()
                result,records=advance_past_overlap(planner,arms,outward,ports,sources,full,members)
                self.assertGreater(result[0],80.5)
                self.assertLess(result[0],83.)
                self.assertEqual(result[1],outward[1])
                self.assertEqual(len(records),1)
                self.assertGreater(records[0]['certified_remainder_cells'],100)

    def test_full_path_and_aggregate_budget_reject_before_large_array_allocation(self):
        path=Path(Section(1,1,((0.,0.,0.),(10001.,0.,0.)),width=2.,interpolation=0))
        with patch('numpy.linspace',side_effect=AssertionError('large array was allocated')):
            with self.assertRaisesRegex(ValueError,'allocation budget'):
                _canonical_samples({1:path},{},.5,.00025)
        short=Path(Section(2,2,((0.,0.,0.),(1.,0.,0.)),width=2.,interpolation=0))
        with self.assertRaisesRegex(ValueError,'allocation budget'):
            _canonical_samples({1:short,2:short},{},.5,.00025,max_samples=5)
        with self.assertRaisesRegex(ValueError,'finite allocation'):
            _canonical_samples({1:short},{},.5,.00025,max_samples=0)


if __name__=='__main__':unittest.main()
