"""Original-bank wrapper coincidence preserves geometry authority and bounds."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Key
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_cycles import PolylineCurb
from tested_python_streets.compound_curbs import OffsetCurb,CompositeCurb
from tested_python_streets.compounds import CompoundGroup,SampledFootprint,RimSpan,RimSegment
from tested_python_streets.regions import Region,BoundaryLoop
from tested_python_streets.source_curb_coordinates import (
    CurbCoordinates,NativeBankCurb,original_bank_intervals,canonical_native_banks)
from tested_python_streets.coincident_source_curbs import contact
from tested_python_streets.sampling import CertifiedCurbSampler


def path(sid=1):
    return Path(Section(sid,sid,((0.,0.,0.),(8.,2.,0.),(16.,0.,0.)),
                        width=4.,interpolation=3),.005)


def approach(source,reverse=False):
    return Approach(Route((Piece(source,not reverse),)),int(reverse),Crown(.08,0.))


def prepared_bank(arm,first=0.,last=6.,offset=2.,tolerance=.005):
    source=arm.path.pieces[0]
    if source.forward:
        at=source.path.samples(first,last,.25,tolerance=tolerance,max_samples=20000)
    else:
        at=source.path.length-source.path.samples(source.path.length-last,
            source.path.length-first,.25,tolerance=tolerance,max_samples=20000)[::-1]
    curve=PolylineCurb(arm.point(at,offset))
    curve.source_coordinates=CurbCoordinates(arm.path,at,np.full(len(at),offset),
        np.full(len(at),1 if offset==arm.offset+.5*arm.width else -1),approach=arm,
        maximum_chord_error=tolerance)
    return curve


def composite(arm,first,last,offset):
    corner=SimpleNamespace(center=None)
    return CompositeCurb(corner,(OffsetCurb(arm,first,last,offset),))


def contact_group(first,last,*,first_contact=.5,last_contact=.5):
    # This test isolates eligibility metadata. Mesh validity is tested through
    # the normal singleton/ordinary-port fixtures and native import audits.
    owners=Key('prepared_bank',1),Key('retained_bank',2)
    triangle=np.array(((0.,0.),(1.,0.),(0.,1.)))
    footprints=tuple(SampledFootprint(Region(owner,triangle),(),
        SimpleNamespace(boundaries=(curve,),rows=100)) for owner,curve in zip(owners,(first,last)))
    spans=(RimSpan(owners[0],RimSegment('curb',0,0.,first_contact)),
           RimSpan(owners[1],RimSegment('curb',0,last_contact,1.)))
    return CompoundGroup((Key('intersection',1),),(),footprints,
        (BoundaryLoop(np.array(((0.,0.),(1.,0.))),spans),),SimpleNamespace())


class NativeCurbOwnershipTests(unittest.TestCase):
    def test_native_evaluation_retains_every_old_vertex_and_improves_certified_error(self):
        source=path();arm=approach(source);old=prepared_bank(arm)
        bank=NativeBankCurb(old);sampler=CertifiedCurbSampler(.005,.7,20000)
        t=bank.certified_parameters(sampler,0.,1.)
        np.testing.assert_array_equal(bank(bank.stations/bank.length),old.points)
        self.assertTrue(set(bank.stations/bank.length)<=set(t))
        bound=bank.curvature_bound*(bank.maximum_speed_bound*bank.length*np.max(np.diff(t)))**2/8.
        self.assertLessEqual(bound,bank.original_chord_error_bound*(1.+1e-12))
        self.assertLessEqual(bound,sampler.tolerance*(1.+1e-12))
        q=np.linspace(0.,1.,201)
        station=bank.outward_station(q)
        np.testing.assert_array_equal(bank(q),arm.point(station,2.))
        np.testing.assert_array_equal(old.points,old.contact_reference_points)

    def test_reversed_route_uses_the_identical_original_path_and_native_bank(self):
        source=path();a=approach(source);b=approach(source,True)
        old=prepared_bank(a,1.,5.);bank=NativeBankCurb(old)
        other=composite(b,source.length-5.,source.length-1.,-2.)
        lower=np.array((.1,.1));upper=np.array((.9,.9));initial=np.array((.45,.55))
        result=contact((bank,other),initial,lower,upper)
        self.assertIsNotNone(result)
        np.testing.assert_allclose(bank(result[0]),other(result[1]),atol=1e-10,rtol=0.)
        rows=original_bank_intervals(other,.1,.9,.5)
        self.assertTrue(rows)
        self.assertTrue(all(row['path'] is source and row['offset']==2. for row in rows))

    def test_equal_coordinates_other_sources_offsets_and_disjoint_domains_do_not_connect(self):
        source=path();a=approach(source);bank=NativeBankCurb(prepared_bank(a,1.,5.))
        initial=np.array((.5,.5));lower=np.array((.1,.1));upper=np.array((.9,.9))
        other_path=path(2)
        self.assertIsNone(contact((bank,composite(approach(other_path),1.,5.,2.)),initial,lower,upper))
        # An actual differently authored bank on the same source Path is still
        # distinct, even below the old geometric proximity tolerance.
        arm=approach(source)
        other=composite(arm,1.,5.,2.+1e-9)
        self.assertIsNone(contact((bank,other),initial,lower,upper))
        self.assertIsNone(contact((bank,composite(a,8.,10.,2.)),initial,lower,upper))

    def test_only_exposed_finite_same_bank_contacts_replace_a_prepared_curve(self):
        source=path();a=approach(source);old=prepared_bank(a,1.,5.)
        other=composite(a,1.,5.,2.)
        group=contact_group(old,other)
        result,records=canonical_native_banks(group)
        self.assertTrue(records)
        bank=result.footprints[0].layout.boundaries[0]
        self.assertIsInstance(bank,NativeBankCurb)
        self.assertIs(result.footprints[1],group.footprints[1])
        self.assertIs(group.footprints[0].layout.boundaries[0],old)
        self.assertEqual(result.members,group.members);self.assertEqual(result.claims,group.claims)
        np.testing.assert_array_equal(result.loops[0].xy,group.loops[0].xy)
        np.testing.assert_array_equal(result.footprints[0].region.xy,group.footprints[0].region.xy)
        disjoint=contact_group(old,composite(a,8.,10.,2.))
        self.assertIs(canonical_native_banks(disjoint)[0],disjoint)
        unrelated=contact_group(old,composite(approach(path(2)),1.,5.,2.))
        self.assertIs(canonical_native_banks(unrelated)[0],unrelated)
        no_ledger=PolylineCurb(old.points)
        ordinary=contact_group(no_ledger,other)
        self.assertIs(canonical_native_banks(ordinary)[0],ordinary)
        legacy=prepared_bank(a,1.,5.)
        legacy.source_coordinates=CurbCoordinates(a.path,legacy.source_coordinates.stations,
            legacy.source_coordinates.offsets,legacy.source_coordinates.sides,approach=a)
        legacy_group=contact_group(legacy,other)
        self.assertIs(canonical_native_banks(legacy_group)[0],legacy_group)
        with self.assertRaisesRegex(ValueError,'original declared source sampling bound'):
            NativeBankCurb(legacy)
        # A prepared chord's native map is not itself an exact analytical contact.
        self.assertFalse(original_bank_intervals(old,0.,1.,.5))

    def test_route_seams_and_exact_knots_stay_compulsory_with_native_reverse_offsets(self):
        first=Path(Section(1,1,((0.,0.,0.),(4.,0.,0.)),width=4.,interpolation=0),.005)
        second=Path(Section(2,1,((9.,0.,0.),(4.,0.,0.)),width=4.,interpolation=0),.005)
        arm=Approach(Route((Piece(first,True),Piece(second,False))),0,Crown(.08,0.))
        at=np.array((0.,1.5,6.3));old=PolylineCurb(arm.point(at,2.))
        old.source_coordinates=CurbCoordinates(arm.path,at,np.full(3,2.),np.ones(3),approach=arm,
                                               maximum_chord_error=.00025)
        bank=NativeBankCurb(old);sampler=CertifiedCurbSampler(.005,.7,20000)
        t=bank.certified_parameters(sampler,0.,1.)
        outward=bank.outward_station(t)
        self.assertTrue(np.any(outward==4.))
        rows=original_bank_intervals(bank,0.,1.,.5)
        self.assertEqual({(row['path'].section.id,row['offset']) for row in rows},{(1,2.),(2,-2.)})
        self.assertTrue(set(bank.stations/bank.length)<=set(t))
        with self.assertRaisesRegex(ValueError,'sample budget'):
            bank.certified_parameters(CertifiedCurbSampler(.005,.00001,4),0.,1.)

    def test_changed_or_nonbank_ledgers_are_not_accepted(self):
        arm=approach(path());old=prepared_bank(arm)
        changed=PolylineCurb(old.points.copy());changed.points[4,0]+=.001
        changed.source_coordinates=old.source_coordinates
        with self.assertRaisesRegex(ValueError,'physical vertices'):NativeBankCurb(changed)
        cap=PolylineCurb(arm.point(np.array((0.,0.)),np.array((-2.,2.))))
        cap.source_coordinates=CurbCoordinates(arm.path,(0.,0.),(-2.,2.),(-1,1),approach=arm)
        self.assertFalse(original_bank_intervals(cap,0.,1.,.5,prepared=True))

    def test_original_evaluation_batch_is_exact_through_arbitrary_subsets_without_accepting_moved_vertices(self):
        source=path();grid=source.samples(0.,source.length,.17,tolerance=.00025,max_samples=20000)
        points=np.vstack((source.point(grid,2.),source.point(grid,-2.)[::-1]))
        coordinates=CurbCoordinates(source,np.r_[grid,grid[::-1]],
            np.r_[np.full(len(grid),2.),np.full(len(grid),-2.)],np.r_[np.ones(len(grid)),-np.ones(len(grid))],
            maximum_chord_error=.00025,evaluation_grid=grid)
        self.assertTrue(coordinates.validate(points))
        indices=np.r_[np.arange(1,len(grid),3),len(grid)+np.arange(0,len(grid),5)]
        subset=coordinates.subset(indices);self.assertTrue(subset.validate(points[indices]))
        moved=points[indices].copy();moved[0,0]=np.nextafter(moved[0,0],np.inf)
        with self.assertRaisesRegex(ValueError,'physical vertices'):subset.validate(moved)
        with self.assertRaises(ValueError):subset._vertex_evaluation.points[0,0]+=.001
        with self.assertRaises(AttributeError):subset._vertex_evaluation.points=moved
        subset.stations=subset.stations.copy();subset.stations[0]=np.nextafter(subset.stations[0],np.inf)
        with self.assertRaisesRegex(ValueError,'declared source evaluation'):subset.validate(points[indices])
        with self.assertRaisesRegex(ValueError,'evaluation batch'):
            CurbCoordinates(source,(grid[0],grid[1]),(2.,2.),(1,1),evaluation_grid=grid[1:])
        with self.assertRaises(TypeError):
            CurbCoordinates(source,(0.,1.),(2.,2.),(1,1),evaluated_points=points[:2])
        class WrongOwner:
            def __getitem__(self,key):return path()
        with self.assertRaisesRegex(ValueError,'Path ownership'):coordinates.validate(points,WrongOwner())

    def test_declared_source_tolerances_survive_subsets_and_control_all_new_chords(self):
        arm=approach(path())
        for tolerance in (.001,.00025,.00001):
            old=prepared_bank(arm,tolerance=tolerance)
            coordinates=old.source_coordinates
            self.assertEqual(coordinates.maximum_chord_error,tolerance)
            self.assertEqual(coordinates.subset(np.arange(len(old.points))).maximum_chord_error,tolerance)
            with self.assertRaises(AttributeError):coordinates.maximum_chord_error=.1
            bank=NativeBankCurb(old)
            # Deliberately much looser caller: its tolerance cannot replace the
            # original source construction's finer mandatory certificate.
            sampler=CertifiedCurbSampler(.1,1.,20000)
            t=bank.certified_parameters(sampler,0.,1.)
            bound=bank.curvature_bound*(bank.maximum_speed_bound*bank.length*np.max(np.diff(t)))**2/8.
            self.assertLessEqual(bound,tolerance*(1.+1e-12))
            np.testing.assert_array_equal(bank(bank.stations/bank.length),old.points)
            self.assertTrue(set(bank.stations/bank.length)<=set(t))
        for bound in (True,0.,-1.,float('inf'),float('nan')):
            with self.assertRaisesRegex(ValueError,'chord-error bound'):
                CurbCoordinates(arm.path,(0.,1.),(2.,2.),(1,1),approach=arm,
                                maximum_chord_error=bound)

    def test_offset_lookup_breaks_are_compulsory_and_each_actual_jet_piece_fits(self):
        from tested_python_streets.conditions import BoundaryData
        from tested_python_streets.source_curb_conditions import _OffsetRule
        from tested_python_streets.surface_fields import BoundaryJetField
        arm=approach(path());curve=OffsetCurb(arm,2.5,1.75,2.)
        sampler=CertifiedCurbSampler(.00025,.7,20000)
        before=curve(np.linspace(0.,1.,101)).copy()
        t=curve.sample(sampler,.03,.97)
        cuts=curve.cumulative/curve.length
        self.assertTrue(set(cuts[(cuts>.03)&(cuts<.97)])<=set(t))
        k,speed,_=sampler.lead_metrics(arm,curve.parameters,curve.cumulative,curve.offset,0.,curve.length)
        self.assertLessEqual(k*(speed*curve.length*np.max(np.diff(t)))**2/8.,sampler.tolerance*(1.+1e-12))
        for a,b in zip(t,t[1:]):
            self.assertFalse(np.any((cuts>a)&(cuts<b)))
            # Exercise the real fitter on each finite analytical jet piece.
            # This is a two-sided local fitting probe, not a pavement mesh.
            rules=(_OffsetRule(curve,a,b),_OffsetRule(curve,b,a))
            xy=curve(np.array((a,b)))
            data=BoundaryData(xy,np.zeros(2),np.zeros((2,2)),rules,(),0.,.08,False)
            field=BoundaryJetField((data,))
            self.assertLessEqual(field.metadata['maximum_sampled_fit_error'],1e-8)
        np.testing.assert_array_equal(curve(np.linspace(0.,1.,101)),before)


if __name__=='__main__':unittest.main()
