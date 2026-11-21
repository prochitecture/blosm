"""Native side ledgers, including a closed-strip curve containing both banks."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown, ProfileProvider, memberships
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_cycles import _footprint
from tested_python_streets.compounds import CompoundGroup
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.source_curb_coordinates import CurbCoordinates
from tested_python_streets.source_curb_conditions import (SourceCurbConditions,
    SourceCurbConditionUnavailable, _NativeCoordinateRule)


def fixture(reverse=False, mode=0):
    controls = (((0., 0., 0.), (4., 0., 0.)), ((0., 4., 0.), (4., 4., 0.)))
    sections = tuple(model.Section(i, i, points[::-1] if reverse else points,
                                   width=2., interpolation=mode) for i, points in enumerate(controls, 1))
    bundle = model.Bundle(1, (model.BundleMember((model.Traversal(2, not reverse),)),
                             model.BundleMember((model.Traversal(1, not reverse),))))
    network = model.Network('native-bank-condition', sections, (), bundles=(bundle,))
    ownership = memberships(network); products = {}
    for section in sections:
        products[model.Key('path', section.id)] = Path(section, .005)
        products[model.Key('profile', section.id)] = ProfileProvider.build(section, ownership[section.id], network.section_index)
    context = BuildContext(network, products, PatternLibrary())
    path = products[model.Key('path', 1)]
    # Unequal valid subdivisions keep exactly the same two analytical banks.
    # The convenience split then puts right/left/caps in the SAME curve1.
    stations = np.array((0., .25, .5, .75, 1., 1., 0.))*path.length
    sides = np.array((1, 1, 1, 1, 1, -1, -1))
    coordinates = CurbCoordinates(path, stations, sides.astype(float), sides)
    points = path.point(stations, sides.astype(float))
    source = _footprint(model.Key('cycle_strip', 1), points, internal_edges=(4, 6),
                        source_coordinates=coordinates)
    group = CompoundGroup((), (), (source,), ())
    return context, group, source


class NativeCurbConditionTests(unittest.TestCase):
    def test_exact_miter_uses_the_actual_one_sided_offset_jacobian(self):
        from tested_python_streets.compound_cycles import PolylineCurb
        for reverse in (False,True):
            controls=((0.,0.,0.),(8.,0.,0.),(8.,8.,0.))
            if reverse:controls=controls[::-1]
            section=model.Section(1,1,controls,width=4.,offset=.25,interpolation=0)
            path=Path(section,.005);profile=Crown(.08,0.,half=1)
            context={model.Key('path',1):path,model.Key('profile',1):profile}
            for side in (-1,1):
                stations=np.array((2.,6.,10.,14.))
                if reverse:stations=stations[::-1]
                offset=np.full(len(stations),section.offset+.5*section.width*side)
                curve=PolylineCurb(path.point(stations,offset))
                curve.source_coordinates=CurbCoordinates(path,stations,offset,np.full(len(stations),side))
                rule=_NativeCoordinateRule(context,curve,0.,1.)
                branches=(0.,*rule.polynomial_breaks(),1.)
                self.assertGreater(len(branches),2)
                self.assertGreater(np.max(np.linalg.norm(path.normals,axis=1)),1.)
                for first,last in zip(branches,branches[1:]):
                    t=np.linspace(0.,1.,7)
                    height,gradient=rule.sample_piece(t,first,last)
                    metric=(first+t*(last-first))*curve.length
                    edge=int(np.clip(np.searchsorted(curve.stations,.5*(first+last)*curve.length,side='right')-1,0,len(stations)-2))
                    a,b=curve.stations[edge:edge+2]
                    native=stations[edge]+(metric-a)/(b-a)*(stations[edge+1]-stations[edge])
                    # One-sided differences use only the declared original
                    # segment, including either side of the authored miter.
                    lo,hi=sorted((native[0],native[-1]));step=1e-5
                    before=np.maximum(lo,native-step);after=np.minimum(hi,native+step)
                    velocity=(path.point(after,offset[0])-path.point(before,offset[0]))/(after-before)[:,None]
                    transverse=.5*section.width*path.frame(native)[2]
                    np.testing.assert_allclose(np.sum(gradient*velocity,axis=1),0.,rtol=0.,atol=1e-10)
                    np.testing.assert_allclose(np.sum(gradient*transverse,axis=1),profile.derivative(side),rtol=0.,atol=1e-13)
                    np.testing.assert_array_equal(height,np.full(len(t),profile.datum+profile.value(side)))
                    single=_NativeCoordinateRule(context,curve,first,last)
                    self.assertEqual(single.polynomial_breaks(),())
                    direct_height,direct_gradient=single.sample(t)
                    np.testing.assert_array_equal(direct_height,height)
                    np.testing.assert_allclose(direct_gradient,gradient,rtol=0.,atol=1e-14)
                if len(branches)>2:
                    with self.assertRaisesRegex(SourceCurbConditionUnavailable,'every declared'):
                        rule.sample_piece(np.array((0.,1.)),0.,1.)
                for invalid in (np.array((-.01,.5)),np.array((.5,1.01)),np.array((np.nan,))):
                    with self.assertRaisesRegex(ValueError,'finite local'):
                        rule.sample_piece(invalid,branches[0],branches[1])

    def test_smooth_breaks_are_original_metric_knots_and_preserve_the_transverse_jet(self):
        from tested_python_streets.compound_cycles import PolylineCurb
        section=model.Section(1,1,((0.,0.,0.),(12.,0.,0.),(24.,2.,0.)),width=2.,interpolation=3)
        path=Path(section,.005);profile=Crown(.08,0.)
        context={model.Key('path',1):path,model.Key('profile',1):profile}
        native_knots=path.station(path.knots)
        knot=native_knots[1]
        stations=np.array((knot-.1,knot+.1))
        curve=PolylineCurb(path.point(stations,1.))
        curve.source_coordinates=CurbCoordinates(path,stations,np.ones(2),np.ones(2))
        rule=_NativeCoordinateRule(context,curve,0.,1.)
        np.testing.assert_allclose(rule.polynomial_breaks(),(.5,),rtol=0.,atol=1e-13)
        np.testing.assert_array_equal(curve.points,path.point(stations,1.))
        for first,last in ((0.,.5),(.5,1.)):
            # Use the exact returned phase, rather than introducing a second
            # independently rounded parameter for the same original knot.
            first=0. if first==0. else rule.polynomial_breaks()[0]
            last=1. if last==1. else rule.polynomial_breaks()[0]
            t=np.linspace(0.,1.,17);at=stations[0]+(first+t*(last-first))*np.diff(stations)[0]
            _,gradient=rule.sample_piece(t,first,last)
            tangent=path.frame(at)[1];normal=path.frame(at)[2]
            np.testing.assert_allclose(np.sum(gradient*tangent,axis=1),0.,rtol=0.,atol=1e-14)
            np.testing.assert_allclose(np.sum(gradient*normal,axis=1),profile.derivative(1.),rtol=0.,atol=1e-14)

    def test_each_span_uses_true_native_side_even_when_one_curve_contains_both_banks(self):
        for reverse in (False, True):
            for mode in (0, 3):
                context, group, source = fixture(reverse, mode)
                resolver = SourceCurbConditions(context, group)
                native = context[model.Key('path', 1)]
                profile = context[model.Key('profile', 1)]
                arm = Approach(native, 0, profile)
                curve = source.layout.boundaries[1]
                self.assertEqual(set(curve.source_coordinates.sides), {-1, 1})
                represented = set()
                for label in source.segments:
                    if label.kind != 'curb': continue
                    current = source.layout.boundaries[label.index]
                    middle = .5*(label.first+label.last)*current.length
                    edge = int(np.clip(np.searchsorted(current.stations, middle, side='right')-1,
                                       0, len(current.points)-2))
                    piece, = current.source_coordinates.edge_sources(edge)
                    side = piece.authored_side
                    rule = resolver.source_rule(source.region.owner, label.index, label.first, label.last)
                    height, gradient = rule.sample(np.array((.5,)))
                    self.assertEqual(float(height[0]), profile.datum+float(profile.value(side)))
                    station = .5*(piece.first_station+piece.last_station)
                    np.testing.assert_allclose(gradient[0], arm.physical_gradient(station, side), rtol=0., atol=1e-12)
                    if label.index == 1: represented.add(side)
                self.assertEqual(represented, {-1, 1})

    def test_a_source_cap_is_not_a_curb_and_invalid_parameter_ranges_do_not_fall_back(self):
        context, group, source = fixture()
        resolver = SourceCurbConditions(context, group)
        cap = next(label for label in source.segments if label.kind == 'internal')
        with self.assertRaisesRegex(SourceCurbConditionUnavailable, 'source cap'):
            resolver.source_rule(source.region.owner, cap.index, cap.first, cap.last)
        for first, last in ((-.01, .1), (.2, .1), (0., np.nan), (0., 1.01)):
            with self.assertRaises(ValueError):
                resolver.source_rule(source.region.owner, 0, first, last)


if __name__ == '__main__': unittest.main()
