"""Finite width steps preserve real roots and exact incoming curb endpoints."""
from dataclasses import replace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compound_dimensional_corners import prepare
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.intersections import Approach, ApproachExtentError
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.regions import Region, RegionUnion, signed_area
from tested_python_streets.routes import Piece, Route
from tested_python_streets.sampling import CertifiedCurbSampler


def fixture(width=3.5, length=40., angle=0., translation=(0., 0.), reverse=False, routed=False):
    rotation = np.array(((math.cos(angle), -math.sin(angle)), (math.sin(angle), math.cos(angle))))
    result = []
    for sid, w, end in ((1, width, (length, 0.)), (2, 5., (-length, 0.))):
        points = np.array(((0., 0.), end)) @ rotation.T + translation
        if reverse:
            points = points[::-1]
        source = model.Section(sid, sid, tuple((x, y, 0.) for x, y in points), width=w)
        path = Path(source)
        if routed:
            path = Route((Piece(path, not reverse),))
        result.append(Approach(path, int(reverse), Crown(.08, 0.)))
    return tuple(result)


def strip(a, sid, fixed):
    at = np.unique(np.r_[np.linspace(0., a.length, 33), fixed])
    xy = np.vstack((a.point(at, -.5*a.width), a.point(at[::-1], .5*a.width)))
    return Region(model.Key('strip', sid), xy if signed_area(xy) > 0. else xy[::-1])


class DimensionalCornerTests(unittest.TestCase):
    def test_fixed_original_contacts_are_certified_without_a_coarse_reach_search(self):
        a,b=fixture()
        corner=CornerRadiusPolicy.solve(a,b,2.,'STRICT')
        before=(a.source.points,b.source.points,a.profile,b.profile)
        record=prepare(model.Key('fixed',1),corner,model.Settings(),(8.,.25),fixed_stations=True)
        self.assertEqual(record.stations,(8.,.25))
        self.assertEqual(record.report['reach_steps'],0)
        self.assertTrue(record.report['fixed_original_contact_candidate'])
        self.assertLessEqual(record.report['curvature_bound'],.5*(1.+1e-8))
        self.assertEqual(before,(a.source.points,b.source.points,a.profile,b.profile))

    def test_uncertified_fixed_contacts_decline_without_moving_to_another_plane(self):
        from tested_python_streets.continuations import ContinuationReachError
        a,b=fixture();corner=CornerRadiusPolicy.solve(a,b,2.,'STRICT')
        with self.assertRaises(ContinuationReachError):
            prepare(model.Key('fixed',1),corner,model.Settings(),(.175,.25),fixed_stations=True)
        # The ordinary finite search remains available with its old behavior.
        record=prepare(model.Key('searched',1),corner,model.Settings(),(.175,.25))
        self.assertGreater(record.stations[0],.175)
        self.assertNotIn('fixed_original_contact_candidate',record.report)

    def test_original_inner_peaks_need_explicit_positive_zero_jet_opt_in(self):
        settings=model.Settings(curve_tolerance=.00025)
        for reverse,angle,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            a,b=fixture(reverse=reverse,routed=reverse,angle=angle,translation=translation)
            corner=CornerRadiusPolicy.solve(a,b,2.,'STYLE')
            a.profile=Crown(.08,0.,-1);b.profile=Crown(.08,0.,1)
            before=(a.source.points,b.source.points,a.profile,b.profile)
            with self.assertRaisesRegex(ValueError,'full-profile'):
                prepare(model.Key('inner',1),corner,settings)
            record=prepare(model.Key('inner',1),corner,settings,original_inner_half_profiles=True)
            self.assertEqual(record.report['facing_inner_heights'],[.08,.08])
            self.assertEqual(record.report['facing_inner_derivatives'],[0.,0.])
            self.assertEqual(record.report['original_profile_halves'],[-1,1])
            self.assertLessEqual(record.report['curvature_bound'],.5*(1.+1e-8))
            self.assertEqual(before,(a.source.points,b.source.points,a.profile,b.profile))
            with self.assertRaisesRegex(ValueError,'both original'):
                prepare(model.Key('inner',1),corner,settings,original_inner_half_profiles=True,original_outer_half_profiles=True)
            with self.assertRaisesRegex(ValueError,'positive facing original'):
                prepare(model.Key('inner',1),replace(corner,radius_mode='STRICT'),settings,original_inner_half_profiles=True)
            saved=b.profile;b.profile=Crown(.08,0.,-1)
            try:
                with self.assertRaisesRegex(ValueError,'positive facing original'):
                    prepare(model.Key('inner',1),corner,settings,original_inner_half_profiles=True)
            finally:b.profile=saved

    def test_original_bundle_outer_half_at_actual_1747_needs_explicit_typed_opt_in(self):
        controls={827:(((-605.824462890625,-706.3928833007812),(-589.4397583007812,-700.3386840820312),(-584.7872924804688,-698.4467163085938)),3.5,0),
                  874:(((-654.4186401367188,-724.1100463867188),(-605.824462890625,-706.3928833007812)),7.,1)}
        settings=model.Settings(curve_tolerance=.00025)
        for reverse,angle,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
            approaches=[]
            for sid,(controls_xy,width,end) in controls.items():
                xy=np.array(controls_xy)@rotation.T+translation
                if reverse:xy=xy[::-1]
                source=model.Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=width)
                actual_end=1-end if reverse else end
                profile_half=(-1 if actual_end==0 else 1) if sid==874 else 0
                approaches.append(Approach(Path(source),actual_end,Crown(.08,0.,profile_half)))
            a,b=approaches;corner=CornerRadiusPolicy.solve(a,b,2.,'STYLE')
            before=[(x.source.points,x.profile) for x in approaches]
            self.assertIsNone(corner.center)
            with self.assertRaisesRegex(ValueError,'full-profile'):
                prepare(model.Key('patch',1747),corner,settings)
            record=prepare(model.Key('patch',1747),corner,settings,original_outer_half_profiles=True)
            self.assertEqual(record.report['minimum_radius'],2.)
            self.assertLessEqual(record.report['curvature_bound'],.5*(1.+1e-8))
            self.assertEqual(record.report['original_profile_halves'],[0,-1])
            self.assertEqual(record.report['facing_outer_heights'],[0.,0.])
            self.assertEqual(float(b.profile.value(-1.)),.08)
            self.assertEqual(before,[(x.source.points,x.profile) for x in approaches])
            np.testing.assert_allclose(record.curve(0.),a.point(record.stations[0],-.5*a.width),atol=1e-10)
            np.testing.assert_allclose(record.curve(1.),b.point(record.stations[1],.5*b.width),atol=1e-10)
            # Explicit STRICT and the opposite plateau side stay outside this
            # original-root STYLE extension even when a caller opts in.
            for invalid in (replace(corner,radius_mode='STRICT'),replace(corner,first=b,second=a)):
                with self.assertRaisesRegex(ValueError,'zero-height outer halves'):
                    prepare(model.Key('patch',1747),invalid,settings,original_outer_half_profiles=True)

    def test_actual_curved_337_and_424_use_finite_original_widths_and_requested_radius(self):
        controls = {
            168: (5., ((-1225.73388671875,-645.4879760742188),(-1222.2176513671875,-646.0897827148438),
                       (-1210.805419921875,-649.3203125),(-1163.52294921875,-662.2200927734375))),
            169: (2.5, ((-1240.212646484375,-649.1141357421875),(-1235.1865234375,-646.977783203125),
                        (-1231.125244140625,-645.7651977539062),(-1227.5699462890625,-645.5432739257812),
                        (-1225.73388671875,-645.4879760742188))),
            212: (2.5, ((-1163.52294921875,-662.2200927734375),(-1146.198486328125,-666.6427001953125),
                        (-1140.566162109375,-667.9906616210938)))}
        settings = model.Settings(curve_tolerance=.00025)
        for site, ends in ((337, ((168,0),(169,1))), (424, ((212,0),(168,1)))):
            for angle, translation in ((0.,(0.,0.)),(.713,(1200.,-1800.))):
                with self.subTest(site=site, angle=angle):
                    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
                    approaches=[]
                    for sid,end in ends:
                        width,xy=controls[sid]
                        points=np.asarray(xy)@rotation.T+translation
                        section=model.Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width)
                        approaches.append(Approach(Path(section),end,Crown(.08,0.)))
                    original=CornerRadiusPolicy.solve(*approaches,3.,'STYLE',
                        math.radians(settings.continuation_angle_degrees))
                    if site==337:
                        self.assertIsNotNone(original.center)
                        # Only the caller's explicit root-support failure may
                        # select this existing near-through STYLE alternative.
                        with self.assertRaisesRegex(ValueError,'request'):
                            prepare(model.Key('patch',site),original,settings)
                    corner=replace(original,center=None,continuation=True,stations=(0.,0.))
                    record=prepare(model.Key('patch',site),corner,settings)
                    self.assertEqual(record.report['minimum_radius'],3.)
                    self.assertLessEqual(record.report['curvature_bound'],(1.+1e-8)/3.)
                    self.assertLessEqual(record.report['reach_steps'],16)
                    for a,s in zip(approaches,record.stations):
                        self.assertLess(s,a.length)
                    for i,label in enumerate(record.source.segments):
                        curve=record.source.layout.boundaries[label.index]
                        np.testing.assert_allclose(curve(label.first),record.source.region.xy[i],atol=1e-10)

    def test_paired_width_steps_keep_full_original_roots_and_all_preparation_seams_internal(self):
        settings = model.Settings(curve_tolerance=.00025)
        for width in (2.5, 3.5):
            for reverse, routed, angle in ((False, False, 0.), (True, True, .713)):
                with self.subTest(width=width, reverse=reverse, routed=routed):
                    a, b = fixture(width=width, reverse=reverse, routed=routed, angle=angle,
                                   translation=(1200., -1800.))
                    before = (a.source.points, b.source.points)
                    records = []
                    for first, last in ((a, b), (b, a)):
                        corner = CornerRadiusPolicy.solve(first, last, 10., 'STRICT')
                        record = prepare(model.Key('patch', len(records)), corner, settings,
                                         minimum_stations=(.4, .6))
                        self.assertEqual(corner.radius, 10.)
                        self.assertLessEqual(record.report['curvature_bound'], (1.+1e-8)/10.)
                        np.testing.assert_allclose(record.curve(0.), first.point(record.stations[0], -.5*first.width), atol=1e-10)
                        np.testing.assert_allclose(record.curve(1.), last.point(record.stations[1], .5*last.width), atol=1e-10)
                        self.assertGreaterEqual(record.stations[0], .4)
                        self.assertGreaterEqual(record.stations[1], .6)
                        # The actual wider root is a mandatory exact piece join.
                        anchor = (last.point(0., .5*last.width) if last.width > first.width
                                  else first.point(0., -.5*first.width))
                        joins = np.array([p(0.) for p in record.curve.pieces]+[record.curve.pieces[-1](1.)])
                        self.assertLess(np.linalg.norm(joins-anchor, axis=1).min(), 1e-10)
                        sampler = CertifiedCurbSampler(settings.curve_tolerance, settings.spacing, 8192)
                        self.assertGreater(len(sampler.sample(record.curve, 0., 1.)), 4)
                        coarse = CertifiedCurbSampler(.02, settings.spacing, 8192)
                        np.testing.assert_array_equal(coarse.sample(record.curve, 0., 1.),
                                                      sampler.sample(record.curve, 0., 1.))
                        for i, label in enumerate(record.source.segments):
                            curve = record.source.layout.boundaries[label.index]
                            np.testing.assert_allclose(curve(label.first), record.source.region.xy[i], atol=1e-10)
                            np.testing.assert_allclose(curve(label.last), record.source.region.xy[(i+1)%len(record.source.region.xy)], atol=1e-10)
                        records.append(record)
                    loops = RegionUnion().build(tuple(r.source.region for r in records)+
                        (strip(a, 1, (records[0].stations[0], records[1].stations[1])),
                         strip(b, 2, (records[0].stations[1], records[1].stations[0]))))
                    self.assertEqual(len(loops), 1)
                    patches = {r.source.region.owner: r.source for r in records}
                    self.assertTrue(all(patches[s.owner].segments[s.segment].kind == 'curb'
                                        for s in loops[0].spans if s.owner in patches))
                    self.assertEqual(before, (a.source.points, b.source.points))

    def test_finite_scope_and_authored_radius_are_never_extended_or_reduced(self):
        a, b = fixture(length=1.)
        corner = CornerRadiusPolicy.solve(a, b, 10., 'STRICT')
        with self.assertRaises(ApproachExtentError):
            prepare(model.Key('patch', 1), corner, model.Settings())
        self.assertEqual(corner.radius, 10.)
        a, b = fixture()
        corner = CornerRadiusPolicy.solve(a, b, 10., 'STRICT')
        with self.assertRaises(ApproachExtentError):
            prepare(model.Key('patch', 1), corner, model.Settings(), (a.length, 0.))
        for invalid in (replace(corner, continuation=False), replace(corner, radius=0.)):
            with self.assertRaisesRegex(ValueError, 'request'):
                prepare(model.Key('patch', 1), invalid, model.Settings())
        b.offset = .2
        with self.assertRaisesRegex(ValueError, 'centred'):
            prepare(model.Key('patch', 1), corner, model.Settings())
        for change in ('half','width','absorbed_sites','mouth_piece','deferred_profile_reserve'):
            a,b=fixture();corner=CornerRadiusPolicy.solve(a,b,10.,'STRICT')
            if change=='half':a.profile=replace(a.profile,half=1)
            elif change=='width':a.width+=.1
            else:setattr(a.path,change,(model.Key('intersection',1),))
            with self.subTest(derived=change),self.assertRaisesRegex(ValueError,'unmodified original full-profile'):
                prepare(model.Key('patch',1),corner,model.Settings())


if __name__ == '__main__':
    unittest.main()
