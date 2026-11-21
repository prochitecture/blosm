"""Ordinary prefix retries retain the original graph and finite ownership."""
from copy import copy
from dataclasses import replace
from functools import lru_cache
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package

from test_python_street_ordinary_singleton_sources import fixture
from test_python_street_source_cycle_atlas import grid
from tested_python_streets.model import End, Intersection, Key, Network, Section, fingerprint
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import Route, RouteTopology, route_key
from tested_python_streets.updates import BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.ordinary_mouth_clearance import correct
from tested_python_streets.ordinary_singleton_sources import build as source_union
from tested_python_streets.source_strip_containment import SourceStripContainmentError
from tested_python_streets.ordinary_source_prefix_atlas import (
    OrdinarySourcePrefixScope, OrdinarySourcePrefixCuts, OrdinarySourcePrefixUnavailable,
    eligible, signature)
from tested_python_streets.source_cycle_atlas import SourceCycleEnergyChart
from tested_python_streets.partition import CompoundQuadChart


@lru_cache(maxsize=1)
def admitted():
    from tested_python_streets.ordinary_mouth_clearance import _group
    from tested_python_streets.cycle_dimensional_retry import certify_layout
    planner, site, raw = fixture()
    try:
        # A fresh actual complete-width failure is sufficient for this scope
        # fixture. Exhaustive ordinary mouth search is covered by its own tests.
        certify_layout(planner.context, _group(site.key,raw), ())
    except SourceStripContainmentError as failure:
        return planner.context, source_union(planner.context, site.key, raw, planner.layouts, (site,), failure)
    raise AssertionError('Fixture must exercise a real failed original-root source union')


def changed_layout(group, **values):
    layout = copy(group.physical_layout)
    for key, value in values.items():
        setattr(layout, key, value)
    return replace(group, physical_layout=layout)


def straight_scope_fixture(reverse=False):
    """A known complete root union tests orientation independently of planning.

    The failed-only ordinary planner is exercised by admitted() and the normal
    native qualification. This analytical T has no root failure to repair.
    """
    from tested_python_streets.ordinary_mouth_clearance import _group
    from tested_python_streets.cycle_dimensional_retry import certify_layout
    controls=(((0.,0.,0.),(-15.,0.,0.)),((0.,0.,0.),(15.,0.,0.)),((0.,0.,0.),(0.,15.,0.)))
    sections=tuple(Section(sid,sid,points[::-1] if reverse else points,interpolation=0)
                   for sid,points in enumerate(controls,1))
    site=Intersection(1,tuple(End(sid,int(reverse)) for sid in range(1,4)))
    network=Network('ordinary-prefix-orientation',sections,(site,))
    products={Key('profile',s.id):Crown(s.rise,s.datum) for s in sections}
    context=BuildContext(network,products,PatternLibrary())
    for section in sections:products[Key('path',section.id)]=PathProvider.build(context,section,())
    topology=RouteTopology(network)
    for end in site.ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    layout=IntersectionProvider.footprint(context,site);group=_group(site.key,layout)
    layout.source_cycle_full_width_proof=certify_layout(context,group,())
    layout.ordinary_complete_mouth_clearance=dict(
        correction='FINITE_ORIGINAL_SINGLE_JUNCTION_SOURCE_UNION',
        source_semantics='one_original_junction_and_its_owned_prefixes',
        original_site=site.key.text(),original_members=[site.key.text()],consumed_sections=[],
        original_paths_profiles_incidence_preserved=True,complete_source_grid=True,
        strict_authored_radii_preserved=True)
    return context,group


class OrdinarySourcePrefixAtlasTests(unittest.TestCase):
    def test_finalized_singleton_keeps_only_its_own_prefix_authority(self):
        from tested_python_streets.ordinary_mouth_clearance import finish_groups

        context, group = straight_scope_fixture()
        site, = group.members
        finalized, = finish_groups((group,), {site: group.physical_layout})
        self.assertTrue(eligible(finalized))
        self.assertEqual(OrdinarySourcePrefixScope.build(context, group),
                         OrdinarySourcePrefixScope.build(context, finalized))

        wrapper = finalized.physical_layout.ordinary_complete_mouth_clearance
        row, = wrapper['corrections']
        for edit in (dict(original_members=['intersection:999']),
                     dict(corrections=[dict(row, original_site='intersection:999')]),
                     dict(corrections=[row, row])):
            with self.subTest(edit=edit):
                stale = changed_layout(finalized, ordinary_complete_mouth_clearance={**wrapper, **edit})
                self.assertFalse(eligible(stale))

    def test_reversed_original_ends_keep_complete_profiles_and_exact_owned_ranges(self):
        expected = None
        for reverse in (False, True):
            context, group = straight_scope_fixture(reverse); before = context.network.to_dict()
            scope = OrdinarySourcePrefixScope.build(context, group)
            scope.validate(context, group)
            self.assertEqual(set(scope.ends), set(context.network.intersections[0].ends))
            self.assertTrue({Key(kind, s.id) for s in context.network.sections
                             for kind in ('path', 'profile')} <= set(scope.dependencies))
            self.assertIn(context.network.intersections[0].key, scope.dependencies)
            policy = OrdinarySourcePrefixCuts(scope, scope.candidates()[0], source_corridors=True,
                                             allow_filled=True, existing_vertices=True)
            self.assertEqual(policy._source_inventory(group, context),
                             (set(scope.candidates()[0]), set(scope.candidates()[0])))
            self.assertEqual(policy.maximum_cuts, 128)
            self.assertEqual(policy.maximum_unpaired_span_widths, 8.)
            self.assertEqual(context.network.to_dict(), before)
            lengths = {sid:sum(b-a for a,b in ranges) for sid,ranges in scope.ranges}
            if expected is None:
                expected = lengths
            else:
                self.assertEqual(set(lengths), set(expected))
                np.testing.assert_allclose(list(lengths.values()), list(expected.values()), atol=1e-10, rtol=0.)

    def test_candidate_axes_rank_owned_lengths_never_entire_path_or_unowned_hull(self):
        context, group = admitted(); scope = OrdinarySourcePrefixScope.build(context, group)
        scope = replace(scope, ranges=((1, ((0.,2.), (100.,103.))), (2, ((0.,8.),)),
                                      (3, ((0.,7.),)), (4, ((0.,6.),))))
        self.assertEqual(scope.candidates(), ((2,), (3,), (4,), (2,3,4)))

    def test_stale_or_overlapping_range_proofs_and_changed_graph_are_rejected(self):
        context, group = admitted(); scope = OrdinarySourcePrefixScope.build(context, group)
        proof = dict(group.physical_layout.source_cycle_full_width_proof)
        proof['sections'] = [*proof['sections'], proof['sections'][0]]
        with self.assertRaisesRegex(ValueError, 'disjoint'):
            OrdinarySourcePrefixScope.build(context, changed_layout(group, source_cycle_full_width_proof=proof))
        proof = dict(group.physical_layout.source_cycle_full_width_proof)
        first = dict(proof['sections'][0]); first['interval'] = [first['interval'][0], .99*first['interval'][1]]
        proof['sections'] = [first, *proof['sections'][1:]]
        with self.assertRaisesRegex(ValueError, 'stale'):
            OrdinarySourcePrefixScope.build(context, changed_layout(group, source_cycle_full_width_proof=proof))
        original = context.network.intersections[0]
        other = replace(original, ends=tuple(End(end.section_id, 1-end.end) for end in original.ends))
        wrong = BuildContext(replace(context.network, intersections=(other,)), context.products, context.patterns)
        with self.assertRaisesRegex(ValueError, 'raw End'):
            OrdinarySourcePrefixScope.build(wrong, group)
        other_graph=BuildContext(replace(context.network,uuid='another-original-network'),context.products,context.patterns)
        with self.assertRaisesRegex(ValueError,'another original graph'):
            scope.validate(other_graph,group)
        outward = list(group.physical_layout.outward_stations); outward[0] += .001
        with self.assertRaisesRegex(ValueError, 'stale'):
            scope.validate(context, changed_layout(group, outward_stations=tuple(outward)))

    def test_full_profile_identity_and_original_union_admission_are_mandatory(self):
        context, group = admitted()
        with self.assertRaises(OrdinarySourcePrefixUnavailable):
            OrdinarySourcePrefixScope.build(context, changed_layout(group, ordinary_complete_mouth_clearance={}))
        changed = copy(group.physical_layout.approaches[0]); changed.offset = .1
        with self.assertRaisesRegex(ValueError, 'centred complete'):
            OrdinarySourcePrefixScope.build(context, changed_layout(group,
                approaches=(changed,*group.physical_layout.approaches[1:])))
        self.assertFalse(eligible(replace(group, claims=(object(),))))
        with self.assertRaisesRegex(ValueError, 'original Path identity'):
            products = dict(context.products)
            first = context.network.sections[0]
            products[Key('path',first.id)] = PathProvider.build(context,first,())
            OrdinarySourcePrefixScope.build(BuildContext(context.network,products,context.patterns),group)

    def test_signature_closes_every_actual_path_profile_range_and_original_site(self):
        context, group = admitted()
        original_digest = context.digest
        context.digest = lambda key: fingerprint(context.network.records()[key]) if key.kind=='intersection' else fingerprint(key.text())
        try:
            first = signature(context, group)
            fingerprint(first)
            context.digest = lambda key: 'changed-native-profile' if key==Key('profile',3) else 'old'
            second = signature(context, group)
            self.assertNotEqual(first, second)
            self.assertEqual({key for key,_ in first[-1]},
                {'intersection:1', *(Key(kind,s.id).text() for s in context.network.sections for kind in ('path','profile'))})
        finally:
            context.digest = original_digest

    def dispatch(self, accept, admission=None):
        context, group = admitted(); scope = OrdinarySourcePrefixScope.build(context, group)
        parent = CompoundQuadChart(context.patterns,source_height_admission=admission)
        calls = []
        def attempt(chart, value, actual, **kwargs):
            policy=chart.partitioner.policy; token=(policy.axes,policy.cadence_multiplier)
            calls.append(token)
            self.assertEqual(chart.max_cells,50000);self.assertEqual(policy.maximum_cuts,128)
            # A child cannot recursively re-enter this same new family.
            with self.assertRaisesRegex(ValueError,'recursively'):
                chart._ordinary_source_prefix_steps(value,actual,'previous')
            yield token
            if not accept(token):raise ValueError('finite candidate rejected')
            return SimpleNamespace(mesh=grid(),holes={},partition=SimpleNamespace(group=value),metadata={})
        with patch.object(SourceCycleEnergyChart,'build_steps',attempt):
            steps = parent._ordinary_source_prefix_steps(group,context,'old grading rejected')
            advances = []
            try:
                while True:
                    advances.append(next(steps))
            except StopIteration as finished:
                result = finished.value
            except ValueError as failure:
                if 'No qualified finite' not in str(failure):raise
                result=None
            finally:
                steps.close()
            self.assertEqual(advances, calls)
        return scope,calls,result

    def test_optional_height_admission_rejects_only_candidates_and_keeps_finite_order(self):
        attempts=[]
        def admission(group,chart):
            attempts.append(chart)
            if len(attempts)<3:raise ValueError('Exact immutable incoming floor exceeds configured limit')
            return dict(necessary_feasibility=True,solver_ran=False)
        scope,calls,result=self.dispatch(lambda token:True,admission)
        self.assertEqual(calls,[(axes,1) for axes in scope.candidates()[:3]])
        self.assertEqual(result.metadata['ordinary_source_prefix_atlas']['source_height_admission'],
                         dict(necessary_feasibility=True,solver_ran=False))
        for invalid in (True,False,0,'policy'):
            with self.assertRaisesRegex(ValueError,'optional callable'):
                CompoundQuadChart(PatternLibrary(),source_height_admission=invalid)

    def test_finite_cadence_order_and_success_preserve_old_owned_inventory(self):
        scope=OrdinarySourcePrefixScope.build(*admitted())
        accepted=(scope.candidates()[1],2)
        _,calls,result=self.dispatch(lambda token:token==accepted)
        wanted=[(axes,cadence) for cadence in (1,2) for axes in scope.candidates()]
        self.assertEqual(calls,wanted[:wanted.index(accepted)+1])
        details=result.metadata['ordinary_source_prefix_atlas']
        self.assertEqual(details['source_ranges'],dict(scope.ranges))
        self.assertFalse(details['physical_boundary_changed'])
        self.assertEqual(details['maximum_attempts'],8)
        _,calls,result=self.dispatch(lambda token:False)
        self.assertIsNone(result);self.assertEqual(calls,wanted);self.assertLessEqual(len(calls),8)

    def test_previously_successful_normal_chart_never_enters_new_retry(self):
        context,group=admitted()
        before=CompoundQuadChart(context.patterns).build(group,context)
        with patch.object(CompoundQuadChart,'_ordinary_source_prefix_steps',side_effect=AssertionError('Successful old chart retried')):
            after=CompoundQuadChart(context.patterns).build(group,context)
        np.testing.assert_array_equal(before.mesh.xyz,after.mesh.xyz)
        np.testing.assert_array_equal(before.mesh.quads,after.mesh.quads)
        self.assertEqual(before.metadata,after.metadata)


if __name__=='__main__':unittest.main()
