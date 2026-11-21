"""Original terminal bank identities survive the unchanged polygon operations."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Key
from tested_python_streets.compounds import CompoundGroup
from tested_python_streets.cycle_terminal_ports import footprint, minimum_station, _clip
from tested_python_streets.terminal_curb_ancestry import original_strip, clip, validate, ClipEdge
from tested_python_streets.source_curb_conditions import SourceCurbConditions, SourceCurbConditionUnavailable
from test_python_street_terminal_ports import fixture


class TerminalCurbAncestryTests(unittest.TestCase):
    def test_linear_ray_uses_the_incident_chord_on_each_side_of_a_knot(self):
        from tested_python_streets.finite_terminal_entrances import LinearFrame,LinearRay
        frame=LinearFrame(np.array((1700.,-500.)),np.array((.6,.8)))
        ray=LinearRay(frame,frame.world(np.array(((0.,0.),(2.,1.),(4.,.5)))))
        a,b,c=ray.stations
        point=np.array([b])
        y1,d1=ray.jet(point,interval=(a,b))
        y2,d2=ray.jet(point,interval=(c,b))
        np.testing.assert_array_equal(y1,y2)
        self.assertAlmostEqual(d1,.5)
        self.assertAlmostEqual(d2,-.25)
        with self.assertRaises(ValueError):
            ray.jet(point,interval=(a,c))

    def test_geometry_and_original_stations_survive_exact_smooth_storage_and_rigid_changes(self):
        for mode in (0, 3):
            for reverse, reflection, angle, translation in ((False, False, 0., (0., 0.)),
                    (True, True, .713, (1800., -1200.))):
                context, pair, _ = fixture(mode, reverse, reflection, angle, translation)
                saved = context.network.to_dict()
                approach = pair.build(context, {})
                station = minimum_station(approach)+6.
                center, normal = approach.cut_plane(station)
                for role in (0, 1):
                    points, tags = original_strip(approach, role)
                    expected = _clip(points, center, normal)
                    actual, traced = clip(points, tags, center, normal)
                    np.testing.assert_array_equal(actual, expected)
                    self.assertEqual(len(traced), len(actual))
                source = footprint(Key('terminal_ancestry', 1), approach, station, 4)
                ancestry = source.layout.terminal_curb_ancestry
                self.assertEqual(set(ancestry), {label.index for label in source.segments if label.kind == 'curb'})
                for index, bank in ancestry.items():
                    error = validate(approach, bank, source.layout.boundaries[index].points)
                    self.assertLessEqual(error, approach.source_curve_tolerance+1e-8)
                    self.assertEqual(approach.paths[bank.role].section.id, bank.section)
                self.assertEqual(context.network.to_dict(), saved)

    def test_missing_or_misassigned_original_ancestry_is_rejected_without_bank_search(self):
        context, pair, _ = fixture()
        approach = pair.build(context, {})
        source = footprint(Key('terminal_ancestry', 1), approach, minimum_station(approach)+6., 4)
        index, bank = next(iter(source.layout.terminal_curb_ancestry.items()))
        points = source.layout.boundaries[index].points
        for changed in (None, replace(bank, section=999), replace(bank, side=-bank.side)):
            with self.assertRaises(ValueError):
                validate(approach, changed, points)
        with self.assertRaises(ValueError):
            ClipEdge('curb', (bank, bank))
        group = CompoundGroup((), (), (source,), ())
        source.layout.terminal_curb_ancestry.pop(index)
        with self.assertRaisesRegex(SourceCurbConditionUnavailable, 'missing its original bank ancestry'):
            SourceCurbConditions(context, group).source_rule(source.region.owner, index, 0., 1.)
        saved = approach.source_curve_stations
        approach.source_curve_stations = (saved[0][:-1], saved[1])
        try:
            with self.assertRaisesRegex(ValueError, 'station arrays'):
                original_strip(approach, 0)
        finally:
            approach.source_curve_stations = saved

    def test_original_bank_jets_are_used_only_outside_the_proven_common_field_domain(self):
        context, pair, _ = fixture()
        approach = pair.build(context, {})
        source = footprint(Key('terminal_ancestry', 1), approach, minimum_station(approach)+6., 4)
        group = CompoundGroup((), (), (source,), ())
        resolver = SourceCurbConditions(context, group)
        low = max(float(ray.stations[0]) for ray in approach.field.curves)
        high = min(float(ray.stations[-1]) for ray in approach.field.curves)
        outside = inside = 0
        for index, bank in source.layout.terminal_curb_ancestry.items():
            curve = source.layout.boundaries[index]
            rule = resolver.source_rule(source.region.owner, index, 0., 1.)
            for t in (0., .5, 1.):
                xy = curve(np.array((t,)))
                coordinate = float(approach.frame_chart.local(xy)[0, 0])
                height, gradient = rule.sample(t)
                # The endpoint jet belongs to the interior of this physical
                # curb. Identify its four actual ray chords independently.
                adjacent=t+(-1e-8 if t==1. else 1e-8)
                inward=float(approach.frame_chart.local(curve(np.array((adjacent,))))[0,0])
                if low <= inward <= high:
                    inside += 1
                    chords=[]
                    for ray in approach.field.curves:
                        index=np.searchsorted(ray.stations,inward,side='right')-1
                        chords.append(ray.stations[index:index+2])
                    interval=max(c[0] for c in chords),min(c[1] for c in chords)
                    self.assertLessEqual(interval[0]-1e-10,coordinate)
                    self.assertLessEqual(coordinate,interval[1]+1e-10)
                    np.testing.assert_allclose(height,approach.field.values(xy)[0],rtol=0.,atol=1e-13)
                    np.testing.assert_allclose(gradient,approach.field.gradient(xy,chord_interval=interval)[0],
                                               rtol=0.,atol=1e-13)
                else:
                    outside += 1
                    original = approach.original_approaches[bank.role]
                    path = approach.paths[bank.role]
                    station = bank.first+t*(bank.last-bank.first)
                    distance = station if original.end == 0 else path.length-station
                    q = bank.side*(1 if original.end == 0 else -1)
                    profile = context[Key('profile', bank.section)]
                    self.assertEqual(float(height), float(profile.datum+profile.value(bank.side)))
                    np.testing.assert_array_equal(gradient, original.physical_gradient(distance, q))
        self.assertGreater(inside, 0)
        self.assertGreater(outside, 0)


if __name__ == '__main__':
    unittest.main()
