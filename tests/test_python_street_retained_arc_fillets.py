"""Finite original-circle boundary primitives and explicit parent contacts."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.paths import Path, right
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.compounds import CompoundGroup, RimSpan, RimSegment,rim_label
from tested_python_streets.regions import RegionUnion, SourceSpan
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.generated_junction_curbs import GeneratedJunctionArea
from tested_python_streets.source_curb_conditions import SourceCurbConditions
from tested_python_streets.retained_arc_fillets import (RetainedArcFillet,
    OriginalAnalyticalCircleCurve, prepare_pair, reflex_switches, auxiliary_candidates,
    reconcile, _view, _topology, canonicalize_tangent_fragments,check_consumed_reservations,
    RetainedArcFitUnavailable, RetainedArcReconciliationError)
from test_python_street_source_cycles import prepare, triangle


def primitive_fixture(*, reverse=False, angle=0., translation=(0., 0.), mode='STYLE',
                      second_origin=(-3.,3.),radii=(4.,4.)):
    """Two actual local corner patches, not a complete junction surface fixture."""
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    sections, sites, sources, products = [], [], [], {}
    for node, origin in ((1, np.array((0., 0.))), (2, np.asarray(second_origin))):
        arms = []
        for j, direction in enumerate(((30., 0.), (0., 30.))):
            sid = 2*node-1+j
            xy = np.array((origin, origin+direction))@rotation.T+translation
            xyz = np.c_[xy, np.zeros(2)]
            if reverse:
                xyz = xyz[::-1]
            section = model.Section(sid, sid, tuple(map(tuple, xyz)), width=2., interpolation=0)
            path, profile = Path(section, .005), Crown(.08, 0.)
            sections.append(section); products[model.Key('path', sid)] = path
            products[model.Key('profile', sid)] = profile
            arms.append(Approach(path, int(reverse), profile))
        radius=float(radii[node-1])
        corner = CornerRadiusPolicy.solve(*arms, radius, mode)
        source = corner_patch(model.Key('cycle_corner', node), corner, .2, .00025)
        source.layout.original_junction = model.Key('intersection', node)
        sources.append(source)
        sites.append(model.Intersection(node, tuple(model.End(a.source.id, a.end) for a in arms),
                                        radius=radius, radius_modes=(mode, mode)))
    network = model.Network('retained-circle-primitives', tuple(sections), tuple(sites))
    context = BuildContext(network, products, PatternLibrary())
    group = CompoundGroup(tuple(site.key for site in sites), (), tuple(sources),
        RegionUnion().build(source.region for source in sources),
        SimpleNamespace(corners=tuple(source.layout.corners[0] for source in sources)))
    return context, group


def prepared_view(group):
    switch, = reflex_switches(group, tolerance=.00025)
    li, at, _, _, _ = switch
    fillet, wedge = prepare_pair(group, li, at, model.Key('retained_arc_fillet', 0), .5,
                                spacing=.2, tolerance=.00025)
    physical = copy(group.physical_layout); physical.retained_arc_fillets = (fillet,)
    result = replace(group, footprints=(*group.footprints, wedge), physical_layout=physical)
    return fillet, wedge, _view(result)


class RetainedArcFilletTests(unittest.TestCase):
    def test_retained_parent_contact_uses_its_real_crown_and_inward_normal(self):
        from tested_python_streets.bank_circle_contacts import OriginalCornerArc, ContactSampling
        from tested_python_streets.retained_circle_bank_join import CircleContact
        from tested_python_streets.sampling import AnalyticalContactError
        for reverse,angle in ((False,0.),(True,.71)):
            context,group=primitive_fixture(reverse=reverse,angle=angle,translation=(1200.,-1800.))
            source=group.footprints[0]
            corner,=source.layout.corners
            index=next(s.index for s in source.segments if s.kind=='curb')
            arc=OriginalCornerArc(source.region.owner,index,source.layout.boundaries[index],corner)
            contact=CircleContact(arc,(.2,.8),.4)
            contact.check_curve(ContactSampling(arc,.00025,.4))
            with self.assertRaises(AnalyticalContactError):
                contact.check_curve(ContactSampling(arc,.00025,.45))
            resolver=SimpleNamespace(context=context,group=group,sources={contact.owner:source})
            z,gradient=contact.profile_value(resolver)
            self.assertEqual(z,0.)
            slope=np.linalg.norm(corner.first.physical_gradient(corner.stations[0],-1))
            tangent=arc.derivative(.4)
            np.testing.assert_allclose(gradient,-slope*right(tangent)/np.linalg.norm(tangent),atol=1e-12)

    def test_explicit_circle_subinterval_keeps_auxiliary_leads_and_maps_true_circle_phase(self):
        from test_python_street_corner_patch import approaches
        from tested_python_streets.bank_circle_contacts import OriginalCornerArc
        from tested_python_streets.sampling import CertifiedCurbSampler
        for case in ('negative','mixed'):
            for routed in (False,True):
                first,last=approaches(case,.713,(1200.,-1800.),routed)
                corner=CornerRadiusPolicy.solve(first,last,4.,'STRICT')
                source=corner_patch(model.Key('cycle_corner',1),corner,.5,.00025)
                index,prepared=next((i,curve) for i,curve in enumerate(source.layout.boundaries)
                                    if hasattr(curve,'original_circle_interval'))
                curve=OriginalAnalyticalCircleCurve(source.region.owner,index,prepared,corner)
                declared=curve.circle_interval
                self.assertIsNotNone(declared)
                locus=OriginalCornerArc(source.region.owner,index,prepared,corner)
                knots=np.asarray(declared.prepared_parameters)
                t=.5*(knots[:-1]+knots[1:])
                step=.01*np.diff(knots)
                numerical=(locus(t+step)-locus(t-step))/(2*step[:,None])
                np.testing.assert_allclose(locus.derivative(t),numerical,atol=1e-5,rtol=1e-6)
                np.testing.assert_array_equal(locus.derivative(t[0]),locus.derivative(t)[0])
                phase=np.linspace(0.,1.,61)
                expected=corner.center+corner.radius*np.column_stack((
                    np.cos(corner.start_angle+corner.sweep*phase),np.sin(corner.start_angle+corner.sweep*phase)))
                np.testing.assert_allclose(curve(curve.source_parameter(phase)),expected,atol=2e-12,rtol=0.)
                np.testing.assert_allclose(curve.circle_parameter(curve.source_parameter(phase)),phase,atol=2e-15,rtol=0.)
                for lead in declared.auxiliary_leads:
                    at=np.linspace(*lead.prepared_bounds,17)
                    # The shared contact itself belongs to the analytical arc;
                    # every strictly auxiliary point keeps the original lead.
                    at=at[(at<declared.prepared_bounds[0])|(at>declared.prepared_bounds[1])]
                    np.testing.assert_array_equal(curve(at),prepared(at))
                values=curve.certified_parameters(CertifiedCurbSampler(.00025,.5,8192),0.,1.)
                self.assertTrue(all(value in values for value in declared.prepared_bounds))
                rule=SourceCurbConditions({},CompoundGroup((),(),(source,),())).derived_rule(curve,0.,1.)
                self.assertTrue(rule.corner_fragment);self.assertFalse(rule.source_anchor)
                heights,gradients=rule.sample(np.linspace(0.,1.,101))
                self.assertTrue(np.isfinite(gradients).all())
                self.assertGreaterEqual(heights.min(),0.);self.assertLessEqual(heights.max(),.08+1e-12)

    def test_new_pure_circle_metadata_keeps_the_legacy_wrapper_exact(self):
        from tested_python_streets.sampling import CertifiedCurbSampler
        _,group=primitive_fixture()
        source=group.footprints[0];prepared=source.layout.boundaries[0];corner=source.layout.corners[0]
        legacy=copy(prepared);del legacy.original_circle_interval
        old=OriginalAnalyticalCircleCurve(source.region.owner,0,legacy,corner)
        current=OriginalAnalyticalCircleCurve(source.region.owner,0,prepared,corner)
        self.assertIsNone(current.circle_interval)
        phase=np.linspace(0.,1.,301)
        np.testing.assert_array_equal(current(phase),old(phase))
        np.testing.assert_array_equal(current.lengths,old.lengths)
        sampler=CertifiedCurbSampler(.00025,.5,8192)
        np.testing.assert_array_equal(current.certified_parameters(sampler,.1,.9),
                                      old.certified_parameters(sampler,.1,.9))

    def test_overlapping_consumption_is_a_finite_fit_decline_with_exact_parent_identity(self):
        _,group=primitive_fixture(second_origin=(-1.,1.))
        li,at,*_=reflex_switches(group,tolerance=.00025)[0]
        first,_=prepare_pair(group,li,at,model.Key('retained_arc_fillet',0),1.,spacing=.2,tolerance=.00025)
        second,_=prepare_pair(group,li,at,model.Key('retained_arc_fillet',1),.5,spacing=.2,tolerance=.00025)
        self.assertIsNone(check_consumed_reservations(first,()))
        with self.assertRaisesRegex(RetainedArcFitUnavailable,'reserved finite original'):
            check_consumed_reservations(second,(first,))
        changed=replace(first.parent_arcs[0],original_curve=copy(first.parent_arcs[0].original_curve))
        wrong=replace(first,parent_arcs=(changed,first.parent_arcs[1]))
        with self.assertRaisesRegex(RetainedArcReconciliationError,'analytical ancestry'):
            check_consumed_reservations(second,(wrong,))

    def test_future_unavailable_strict_candidate_cannot_trigger_an_earlier_style_fit(self):
        from tested_python_streets import retained_arc_fillets as module
        context,group=primitive_fixture()
        switch=reflex_switches(group,tolerance=.00025)[0]
        parents=tuple(source.layout.corners[0] for source in group.footprints)
        strict=tuple(replace(parent,radius_mode='STRICT') for parent in parents)
        planner=SimpleNamespace(context=context,sites={node.key:node for node in context.network.intersections})
        # Two policy jobs expose a future hard reservation. The actual finite
        # circle fit still runs and must reject its unavailable 4 m candidate;
        # it may neither spend an earlier STYLE interval nor invent a smaller
        # hard radius. Source discovery/authority is tested independently.
        with patch.object(module,'reflex_switches',return_value=(switch,switch)), \
                patch.object(module,'_original_parents',side_effect=(parents,strict)), \
                patch.object(module,'prepare_pair',wraps=prepare_pair) as prepared:
            with self.assertRaisesRegex(RetainedArcReconciliationError,'STRICT'):
                reconcile(planner,group)
        self.assertEqual(prepared.call_count,1)
        self.assertEqual(prepared.call_args.args[3],model.Key('retained_arc_fillet',1))
        self.assertEqual(prepared.call_args.args[4],4.)
        self.assertEqual(tuple(auxiliary_candidates(*strict)),((4.,False),))
        self.assertTrue(all(parent.radius_mode=='STYLE' for parent in parents))

    def test_tangent_chord_aliases_use_exact_consumed_parents_without_repeating_contacts(self):
        for reverse,angle in ((False,0.),(True,.713)):
            _,group=primitive_fixture(reverse=reverse,angle=angle,translation=(1200.,-1800.),
                                      second_origin=(-.5,.5))
            li,at,*_=reflex_switches(group,tolerance=.00025)[0]
            fillet,wedge=prepare_pair(group,li,at,model.Key('retained_arc_fillet',0),3.2,
                                     spacing=.2,tolerance=.00025)
            layout=copy(group.physical_layout);layout.retained_arc_fillets=(fillet,)
            operands=(*group.footprints,wedge)
            view=_view(replace(group,footprints=operands,physical_layout=layout,
                loops=RegionUnion().build(source.region for source in operands)))
            sampler=AnalyticalBoundarySampler(tolerance=.005,spacing=.7)
            fixed,journal=canonicalize_tangent_fragments(view,sampler)
            self.assertEqual(len(journal),1)
            self.assertIs(fixed.footprints,view.footprints)
            self.assertIs(fixed.physical_layout,view.physical_layout)
            self.assertEqual(sum(loop.hole for loop in fixed.loops),sum(loop.hole for loop in view.loops))
            self.assertTrue(journal[0]['consumed_parent_fragments'])
            for loop in fixed.loops:
                ids=[i for i,span in enumerate(loop.spans) if span.owner==fillet.owner]
                if ids:
                    self.assertEqual(len(ids),1)
                    self.assertEqual(loop.spans[ids[0]].label,RimSegment('curb',0,0.,1.))
                    np.testing.assert_array_equal(loop.xy[ids[0]],fillet(0.))
                    np.testing.assert_array_equal(loop.xy[(ids[0]+1)%len(loop.xy)],fillet(1.))
            # The ordinary duplicate-contact guard remains strict. The helper
            # resolves the actual tangent-run topology rather than skipping a
            # second lookup or increasing its finite allocation.
            self.assertIs(canonicalize_tangent_fragments(fixed,sampler)[0],fixed)

    def test_tangent_resolution_declines_foreign_labels_and_retains_default_views(self):
        _,group=primitive_fixture(second_origin=(-.5,.5))
        li,at,*_=reflex_switches(group,tolerance=.00025)[0]
        fillet,wedge=prepare_pair(group,li,at,model.Key('retained_arc_fillet',0),3.2,
                                 spacing=.2,tolerance=.00025)
        layout=copy(group.physical_layout);layout.retained_arc_fillets=(fillet,)
        operands=(*group.footprints,wedge)
        view=_view(replace(group,footprints=operands,physical_layout=layout,
            loops=RegionUnion().build(source.region for source in operands)))
        sampler=AnalyticalBoundarySampler(tolerance=.005,spacing=.7)
        for loop_index,loop in enumerate(view.loops):
            ids=[i for i,span in enumerate(loop.spans) if span.owner==fillet.owner]
            if not ids:continue
            begin=min(ids,key=lambda i:rim_label(loop.spans[i],{p.region.owner:p for p in view.footprints}).first)
            end=max(ids,key=lambda i:rim_label(loop.spans[i],{p.region.owner:p for p in view.footprints}).last)
            order=np.roll(np.arange(len(loop.spans)),-begin)
            stop=int(np.flatnonzero(order==end)[0])+1
            alias=next(i for i in order[:stop] if loop.spans[i].owner!=fillet.owner)
            label=rim_label(loop.spans[alias],{p.region.owner:p for p in view.footprints})
            spans=list(loop.spans);spans[alias]=RimSpan(loop.spans[alias].owner,replace(label,kind='internal'))
            loops=list(view.loops);loops[loop_index]=replace(loop,spans=tuple(spans))
            with self.assertRaisesRegex(RetainedArcReconciliationError,'third source'):
                canonicalize_tangent_fragments(replace(view,loops=tuple(loops)),sampler)
        excessive=copy(view.physical_layout);excessive.retained_arc_fillets=(fillet,)*17
        with self.assertRaisesRegex(RetainedArcReconciliationError,'finite pair allocation'):
            canonicalize_tangent_fragments(replace(view,physical_layout=excessive),sampler)
        _,ordinary=primitive_fixture();fillet,wedge,view=prepared_view(ordinary)
        # Build the actual union here: prepared_view alone is a primitive
        # ancestry fixture, not a complete post-insertion physical boundary.
        view=replace(view,loops=RegionUnion().build(source.region for source in view.footprints))
        unchanged,journal=canonicalize_tangent_fragments(view,sampler)
        self.assertIs(unchanged,view);self.assertEqual(journal,())

    def test_exact_finite_tangent_primitive_rotates_and_reverses_with_real_sources(self):
        for reverse, angle in ((False, 0.), (True, .713)):
            context, group = primitive_fixture(reverse=reverse, angle=angle, translation=(1200., -1800.))
            original = context.network.to_dict()
            old_points = [source.region.xy.copy() for source in group.footprints]
            fillet, wedge, view = prepared_view(group)
            self.assertIsInstance(fillet, GeneratedJunctionArea)
            self.assertIsInstance(fillet, RetainedArcFillet)
            self.assertFalse(any(hasattr(fillet, name) for name in ('first', 'second', 'path', 'end')))
            self.assertEqual(wedge.layout.corners, ())
            self.assertIs(wedge.layout.generated_junction_area, fillet)
            self.assertIs(view.physical_layout.corners, group.physical_layout.corners)
            self.assertGreater(wedge.region.xy.shape[0], 4)
            for i, parent in enumerate(fillet.parent_arcs):
                actual = fillet(float(i))
                np.testing.assert_allclose(actual, parent.analytical_contact_xy, rtol=0., atol=1e-10)
                tangent = right((actual-parent.original_corner.center)/parent.original_corner.radius)
                np.testing.assert_allclose(tangent, fillet.tangent(float(i)), rtol=0., atol=1e-10)
                self.assertLessEqual(parent.analytical_to_original_chord_distance, .00025+1e-9)
                low, high = sorted(parent.exposed_circle_parameters)
                self.assertLess(low, parent.analytical_contact_parameter)
                self.assertLess(parent.analytical_contact_parameter, high)
            self.assertEqual(context.network.to_dict(), original)
            for source, old in zip(group.footprints, old_points):
                np.testing.assert_array_equal(source.region.xy, old)

    def test_finite_visible_interval_and_original_strict_minimum_are_not_overridden(self):
        _, group = primitive_fixture()
        li, at, _, _, _ = reflex_switches(group, tolerance=.00025)[0]
        with self.assertRaisesRegex(RetainedArcFitUnavailable, 'finite analytical arc'):
            prepare_pair(group, li, at, model.Key('retained_arc_fillet', 0), 1., spacing=.2, tolerance=.00025)
        context, strict = primitive_fixture(mode='STRICT')
        li, at, _, _, _ = reflex_switches(strict, tolerance=.00025)[0]
        with self.assertRaisesRegex(RetainedArcFitUnavailable, 'STRICT minimum'):
            prepare_pair(strict, li, at, model.Key('retained_arc_fillet', 0), .5, spacing=.2, tolerance=.00025)
        planner = SimpleNamespace(context=context, sites={site.key: site for site in context.network.intersections})
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'original radius constraints'):
            reconcile(planner, strict)

    def test_auxiliary_style_tail_reproduces_point_eight_after_unchanged_ordinary_sequence(self):
        _, group = primitive_fixture()
        parents = tuple(source.layout.corners[0] for source in group.footprints)
        ordinary = tuple(CornerRadiusPolicy.candidates(4., 'STYLE'))
        candidates = tuple(auxiliary_candidates(*parents))
        self.assertEqual(tuple(value for value, extra in candidates[:len(ordinary)]), ordinary)
        self.assertFalse(any(extra for value, extra in candidates[:len(ordinary)]))
        self.assertEqual(candidates[len(ordinary)], (.8, True))
        self.assertEqual(candidates[-1][0], .5)
        self.assertLessEqual(len(candidates), 17)
        self.assertTrue(all(value >= .5 for value, extra in candidates))

    def test_parent_parameters_are_explicit_and_wrong_ancestry_or_repeat_declines(self):
        _, group = primitive_fixture()
        _, _, view = prepared_view(group)
        fillet = view.physical_layout.retained_arc_fillets[0]
        parent = fillet.parent_arcs[0]
        layouts = {source.region.owner: source.layout for source in view.footprints}
        first = RimSpan(parent.owner, RimSegment('curb', parent.curve_index,
            max(0., parent.analytical_contact_parameter-.01), parent.analytical_contact_parameter))
        last = RimSpan(fillet.owner, RimSegment('curb', 0, 0., .01))
        sampler = AnalyticalBoundarySampler(tolerance=.005, spacing=.5)
        sampler.retained_arc_contact_tokens = set(); sampler.retained_arc_contact_handoffs = []
        values = sampler._contact(first, last, (first.label, last.label), layouts)
        np.testing.assert_array_equal(values, (parent.analytical_contact_parameter, 0.))
        self.assertEqual(len(sampler.retained_arc_contact_handoffs), 1)
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'pair allocation'):
            sampler._contact(first, last, (first.label, last.label), layouts)
        bad = copy(layouts[parent.owner]); curves = list(bad.boundaries)
        wrong = copy(curves[parent.curve_index]); wrong.original_prepared_curve = copy(parent.original_curve)
        curves[parent.curve_index] = wrong; bad.boundaries = tuple(curves)
        layouts[parent.owner] = bad
        sampler.retained_arc_contact_tokens = set()
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'ancestry'):
            sampler._contact(first, last, (first.label, last.label), layouts)

    def test_clipped_terminal_chord_uses_its_own_sampling_at_either_end(self):
        _, group = primitive_fixture()
        _, _, view = prepared_view(group)
        fillet, = view.physical_layout.retained_arc_fillets
        layouts = {source.region.owner: copy(source.layout) for source in view.footprints}
        knots = fillet.source_reference_parameters
        # Many closure samples make the wedge-wide row interval much shorter
        # than either terminal arc chord.
        layouts[fillet.owner].rows = 100 * len(knots)
        for which, parent in enumerate(fillet.parent_arcs):
            segment = len(knots)-2 if which else 0
            span = SourceSpan(fillet.owner, segment, .25, .75)
            a, b = knots[segment:segment+2]
            label = RimSegment('curb', 0, a+.25*(b-a), a+.75*(b-a))
            curve = layouts[parent.owner].boundaries[parent.curve_index]
            t = curve.source_parameter(parent.analytical_contact_parameter)
            neighbor = RimSpan(parent.owner, RimSegment('curb', parent.curve_index, t, t))
            def resolve(source, clipped):
                sampler = AnalyticalBoundarySampler(tolerance=.005, spacing=.5)
                sampler.retained_arc_contact_tokens = set()
                sampler.retained_arc_contact_handoffs = []
                spans = (source, neighbor) if which else (neighbor, source)
                labels = (clipped, neighbor.label) if which else (neighbor.label, clipped)
                return sampler._contact(*spans, labels, layouts)
            values = resolve(span, label)
            np.testing.assert_array_equal(values, (1., t) if which else (t, 0.))
            with self.assertRaisesRegex(RetainedArcReconciliationError, 'neighborhood'):
                resolve(replace(span, segment=segment+1), label)
            with self.assertRaisesRegex(RetainedArcReconciliationError, 'neighborhood'):
                resolve(span, replace(label, first=label.first+.001))

    def test_formal_generated_height_rule_uses_real_parent_levels_without_new_incidence(self):
        context, group = primitive_fixture()
        _, _, view = prepared_view(group)
        fillet = view.physical_layout.retained_arc_fillets[0]
        rule = SourceCurbConditions(context, view).source_rule(fillet.owner, 0, .1, .9)
        self.assertTrue(rule.corner_fragment)
        self.assertFalse(rule.source_anchor)
        height, gradient = rule.sample(np.linspace(0., 1., 17))
        np.testing.assert_array_equal(height, 0.)
        self.assertTrue(np.isfinite(gradient).all())

    def test_source_circle_equivalence_and_work_bounds_are_not_assumed(self):
        _, group = primitive_fixture()
        source = group.footprints[0]; curve = source.layout.boundaries[0]
        changed = copy(curve); changed.points = curve.points.copy(); changed.points[1] += (.001, 0.)
        with self.assertRaisesRegex(ValueError, 'uniform finite circular|source knots or vertices'):
            OriginalAnalyticalCircleCurve(source.region.owner, 0, changed, source.layout.corners[0])
        context, full = prepare(triangle(2.))
        planner = SimpleNamespace(context=context, sites={site.key: site for site in context.network.intersections})
        self.assertIs(reconcile(planner, full), full)
        for options in ({'maximum_pairs': 0}, {'maximum_pairs': 17}, {'maximum_batches': 0}, {'maximum_batches': 3}):
            with self.assertRaises(ValueError): reconcile(planner, full, **options)
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'hole topology'):
            _topology(full, replace(full, loops=tuple(loop for loop in full.loops if not loop.hole)), .00025)


if __name__ == '__main__': unittest.main()
