"""Owned-prefix computational axes never replace original physical ownership."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_paired_interval_atlas import fixture
from test_python_street_hole_fibres import unequal_lens
from tested_python_streets.model import Key,Section,Network,fingerprint
from tested_python_streets.paths import Path
from tested_python_streets.partition import BoundaryRef
from tested_python_streets.paired_interval_atlas import (
    OwnedPairedPrefixScope,OwnedPairedPrefixCuts,original_intervals,signature)
from tested_python_streets.physical_hole_cuts import PhysicalHoleCuts
from tested_python_streets.hole_fibre_chart import (
    RuledHoleBanks,_RuledHoleFrame,_ruled_paired_hole,build_fibre_hole)
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.bundle_axis import Boundary,Island
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import UpdateCoordinator,BuildContext


def finite_scope(context,group):
    rows=[dict(section=sid,interval=(0.,context[Key('path',sid)].length)) for sid in (3,6,9)]
    rows.extend(dict(section=sid,interval=(8.,10.)) for sid in (2,5))
    proof=dict(sections=rows,maximum_error=.002)
    group.physical_layout.source_cycle_full_width_proof=proof
    return OwnedPairedPrefixScope.build(context,group,original_intervals(context,group),proof)


class InventoryCuts(PhysicalHoleCuts):
    def __init__(self,*,all_intervals=False,**kwargs):
        super().__init__(**kwargs);self.all_intervals=all_intervals
    def _source_inventory(self,group,context):return {1},({1} if self.all_intervals else set())


def cut_region(*,hole=False):
    section=Section(1,1,((0.,0.,0.),(100.,0.,0.)),width=4.,interpolation=0)
    context={Key('path',1):Path(section)}
    loop=SimpleNamespace(xy=np.array(((0.,-2.),(100.,-2.),(100.,20.),(0.,20.))))
    loops=(loop,)
    if hole:loops+= (SimpleNamespace(xy=np.array(((20.,2.),(20.,10.),(80.,10.),(80.,2.)))),)
    group=SimpleNamespace(members=(Key('intersection',7),),loops=loops)
    labels=tuple(tuple(SimpleNamespace(kind='curb') for _ in value.xy) for value in loops)
    return context,group,labels


def point(group,ref):
    xy=group.loops[ref.loop].xy
    return xy[ref.edge]+ref.fraction*(xy[(ref.edge+1)%len(xy)]-xy[ref.edge])


class OwnedPairedPrefixChartTests(unittest.TestCase):
    def test_both_roles_include_only_actual_owned_member_prefixes_and_keep_claims(self):
        for reverse in (False,True):
            context,group=fixture(reverse);before=context.network.to_dict();claims=group.claims
            scope=finite_scope(context,group)
            for role,wanted,all_intervals in ((0,{2,3,9},{2,3}),(1,{5,6,9},{5,6})):
                policy=OwnedPairedPrefixCuts(scope,role,source_corridors=True,hole_bridges=True)
                self.assertEqual(policy._source_inventory(group,context),(wanted,all_intervals))
                self.assertEqual(policy.source_ranges[2],((8.,10.),))
                self.assertEqual(policy.cadence_multiplier,2)
                self.assertEqual(policy.maximum_unpaired_span_widths,8.)
                self.assertNotIn(1,policy.source_ranges)
                self.assertNotIn(4,policy.source_ranges)
            self.assertIs(group.claims,claims)
            self.assertEqual(context.network.to_dict(),before)
            self.assertTrue({Key(kind,sid) for sid in range(1,7) for kind in ('path','profile')}<=set(scope.dependencies))

    def test_original_disjoint_ranges_share_the_band_budget_and_do_not_sample_the_gap(self):
        context,group,labels=cut_region()
        policy=InventoryCuts(source_ranges={1:((0.,10.),(90.,100.))},cadence_multiplier=2)
        cuts=policy._source_cuts(group,context,labels)
        xs=[float(point(group,cut.first)[0]) for cut in cuts]
        self.assertEqual(xs,[5.,95.])
        self.assertTrue(all(cut.section_id is None for cut in cuts))
        fragmented={1:tuple((3.*i,3.*i+1.) for i in range(33))}
        with self.assertRaisesRegex(ValueError,'shared finite range band budget'):
            InventoryCuts(source_ranges=fragmented)._source_cuts(group,context,labels)
        for ranges in ({1:((0.,10.),(9.,11.))},{1:((0.,101.),)}):
            with self.assertRaises(ValueError):
                InventoryCuts(source_ranges=ranges)._source_cuts(group,context,labels)

    def test_all_pavement_intervals_and_metre_span_filter_change_only_computational_cuts(self):
        context,group,labels=cut_region(hole=True)
        before=[loop.xy.copy() for loop in group.loops]
        flags=dict(source_ranges={1:((44.,56.),)},cadence_multiplier=2)
        center=InventoryCuts(**flags)._source_cuts(group,context,labels)
        all_values=InventoryCuts(all_intervals=True,**flags)._source_cuts(group,context,labels)
        self.assertEqual(len(center),1);self.assertEqual(len(all_values),2)
        self.assertEqual(sorted(float(np.linalg.norm(point(group,c.first)-point(group,c.last))) for c in all_values),[4.,10.])
        plain=InventoryCuts(source_ranges={1:((0.,10.),)},maximum_unpaired_span_widths=4.)
        self.assertEqual(plain._source_cuts(group,context,labels),[])
        self.assertEqual(plain.source_cut_reports[0]['physical_span'],22.)
        for loop,xy in zip(group.loops,before):np.testing.assert_array_equal(loop.xy,xy)

    def test_surface_signature_includes_unowned_remote_member_paths_and_exact_ranges(self):
        context,group=fixture();finite_scope(context,group)
        patterns=PatternLibrary();updates=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns)
        first=updates.prepare(context.network);updates.commit(first)
        current=BuildContext(context.network,first.products,patterns,first.task_fingerprints,first.source_fingerprints)
        before=signature(current,group);fingerprint(before)
        section=context.network.section_index[1]
        points=(section.points[0],(.2,5.,0.),section.points[-1])
        network=replace(context.network,sections=tuple(replace(s,points=points) if s.id==1 else s for s in context.network.sections))
        second=updates.prepare(network)
        changed=BuildContext(network,second.products,patterns,second.task_fingerprints,second.source_fingerprints)
        self.assertNotEqual(before,signature(changed,group))
        self.assertNotIn(1,dict(OwnedPairedPrefixScope.build(current,group,original_intervals(current,group),
            group.physical_layout.source_cycle_full_width_proof).ranges))

    def test_ruled_bank_limits_are_continuous_and_keep_actual_adjacent_tip_vertices(self):
        boundary=unequal_lens(52,.713)
        pair=RuledHoleBanks.build(boundary);before=boundary.copy()
        plan=SimpleNamespace(origin=pair.origin,axis=pair.axis,left=pair.left,
            islands=(Island(Boundary(pair.first,[]),Boundary(pair.second,[]),(True,True)),))
        frame=_RuledHoleFrame(plan,plan.islands[0],0,.7)
        np.testing.assert_array_equal(frame.rows,frame.original)
        self.assertTrue(all(frame.columns(i)%4==2 and frame.columns(i)<=128 for i in range(len(frame.rows))))
        stations=frame.original[1:-1]
        for side,bank in ((0.,pair.first),(1.,pair.second)):
            actual=frame.point(stations,np.full(len(stations),side))
            near=frame.point(stations,np.full(len(stations),1e-9 if side==0. else 1.-1e-9))
            np.testing.assert_array_equal(actual,bank[1:-1])
            self.assertLess(float(np.linalg.norm(near-actual,axis=1).max()),1e-7)
        for end,row,near in ((0,0,1),(1,-1,-2)):
            tip=pair.tips[end]
            np.testing.assert_array_equal(tip['original_tip'],pair.first[row])
            np.testing.assert_array_equal(tip['adjacent_original_bank_vertices'],(pair.first[near],pair.second[near]))
        np.testing.assert_array_equal(boundary,before)
        self.assertTrue(boundary.flags.writeable)

    def test_existing_success_is_bit_exact_and_ruled_opt_in_cannot_increase_allocation(self):
        boundary=unequal_lens(51,0.);patterns=PatternLibrary()
        before=build_fibre_hole(patterns,boundary,.7,3000)
        with patch('tested_python_streets.hole_fibre_chart._ruled_paired_hole',side_effect=AssertionError('Unexpected fallback')):
            after=build_fibre_hole(patterns,boundary,.7,3000,ruled_paired=True)
        np.testing.assert_array_equal(after.mesh.xyz,before.mesh.xyz)
        np.testing.assert_array_equal(after.mesh.quads,before.mesh.quads)
        self.assertEqual(after.metadata,before.metadata)
        with self.assertRaises(ValueError):build_fibre_hole(patterns,boundary,.7,1,ruled_paired=True)
        with self.assertRaises(ValueError):build_fibre_hole(patterns,boundary,.7,3000,ruled_paired=1)

    def test_ruled_half_strip_chart_is_a_protected_positive_closed_disk(self):
        for count,angle in ((51,0.),(52,.713)):
            boundary=unequal_lens(count,angle);before=boundary.copy()
            result=_ruled_paired_hole(PatternLibrary(),boundary,.7,3000,ValueError('Original failed chart'))
            np.testing.assert_array_equal(boundary,before)
            self.assertTrue(result.metadata['exact_physical_polygon'])
            self.assertTrue(result.metadata['ruled_paired_hole']['half_strip_templates'])
            self.assertLessEqual(result.mesh.validate()['maximum_edge_ratio'],8.)
            self.assertLessEqual(result.metadata['stored_quality']['maximum_edge_ratio'],8.)
            self.assertLessEqual(len(result.mesh.quads),3000)
            rim=result.mesh.xyz[result.mesh.seams[0],:2]
            for point in boundary:self.assertEqual(int(np.all(rim==point,axis=1).sum()),1)
            self.assertGreater(result.metadata['ruled_paired_hole']['intrinsic_tip_sine'],0.)

    def test_finite_scope_rejects_overlapping_or_outside_original_ranges_and_nonpaired_claims(self):
        context,group=fixture();scope=finite_scope(context,group)
        proof=dict(scope.width_proof,sections=[*scope.width_proof['sections'],dict(section=2,interval=(9.,10.))])
        with self.assertRaisesRegex(ValueError,'disjoint'):
            OwnedPairedPrefixScope.build(context,group,scope.claims,proof)
        proof=dict(scope.width_proof,sections=[dict(section=2,interval=(9.,11.))])
        with self.assertRaisesRegex(ValueError,'finite original source'):
            OwnedPairedPrefixScope.build(context,group,scope.claims,proof)
        with self.assertRaisesRegex(ValueError,'typed original paired terminal'):
            OwnedPairedPrefixScope.build(context,group,(scope.claims[0],scope.claims[0]),scope.width_proof)


if __name__=='__main__':unittest.main()
