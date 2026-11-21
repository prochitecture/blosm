"""Finite parent-chain contacts, sampling references and whole-arc ownership."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.continuations import WidthContinuation
from tested_python_streets.sampling import CertifiedCurbSampler
from tested_python_streets.source_owned_connector_fillets import NativeCurbContact
from tested_python_streets.source_parent_chain_fillets import (offset_solution,ParentChainUnavailable,
    PreparedParentChainFillet,CertifiedCoreContact,OriginalCircleChain,_regularity,validate,discover)
from tested_python_streets.compounds import CompoundGroup,RimSpan,RimSegment
from tested_python_streets.regions import BoundaryLoop
from tested_python_streets.source_parent_chain_preparation import ordered_curb_closure


def parents(reverse=False,angle=0.,translation=(0.,0.)):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    def transform(points):return np.asarray(points)@rotation.T+translation
    core=WidthContinuation(*transform(((0.,0.),(10.,0.))),rotation@np.array((1.,0.)),rotation@np.array((1.,0.)),4.)
    xy=transform(((9.,0.),(9.,-10.)))
    if reverse:xy=xy[::-1]
    section=model.Section(1,1,tuple((x,y,0.) for x,y in xy),width=2.,interpolation=0)
    path=Path(section);arm=Approach(path,int(reverse),Crown(.08,0.))
    return core,arm,transform


def primitive(core,arm,solved):
    owner=model.Key('source_parent_chain_fillet',0)
    contact=CertifiedCoreContact(model.Key('core',1),0,core,.6,(0.,1.),None,model.Key('intersection',1))
    native=NativeCurbContact(arm.path,-1 if arm.end else 1,arm.length-2. if arm.end else 2.,((0.,arm.length),))
    parent=OriginalCircleChain(model.Key('circle',1),0,None,None,None,model.Key('intersection',1))
    knots=tuple(np.linspace(0.,1.,51))
    return PreparedParentChainFillet(owner,solved['center'],2.,solved['start_angle'],solved['sweep'],
        solved['first'],solved['last'],contact,native,parent,(),2.,'STYLE',(2.,),knots,.00025,(0.,0.))


class ParentChainTests(unittest.TestCase):
    def test_closure_traces_declared_curbs_across_loop_origin(self):
        points = np.array([(0.,0.), (1.,0.), (2.,1.), (1.,2.), (0.,1.)])
        owners = [model.Key('source', i) for i in range(5)]
        spans = [RimSpan(owner, RimSegment('curb', 0, 0., 1.)) for owner in owners]
        start = (owners[4], 0, 0., points[4])
        middle = ((owners[0], 0), (owners[1], 0))
        # The final contact lies halfway along its own source curb.
        spans[2:4] = [RimSpan(owners[3], RimSegment('curb', 0, a, b)) for a,b in ((0.,.5),(.5,1.))]
        loop = BoundaryLoop(points, tuple(spans))
        finish = (owners[3], 0, .5, points[3])
        xy, ancestry = ordered_curb_closure((loop,), {}, start, finish, middle)
        np.testing.assert_array_equal(xy, points[[3,2,1,0,4]])
        self.assertEqual([owner for owner,_ in ancestry], [owners[i] for i in (4,0,1,3)])
        for wrong in (middle[::-1], middle[:1]):
            with self.assertRaises(ParentChainUnavailable):
                ordered_curb_closure((loop,), {}, start, finish, wrong)

    def test_finite_normal_offset_contacts_preserve_geometry_under_rigid_motion_and_reverse_storage(self):
        for reverse,angle,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            core,arm,transform=parents(reverse,angle,translation)
            before=(core.controls.copy(),arm.source.points,arm.profile)
            solved=offset_solution(core,(0.,1.),arm,(0.,arm.length),2.)
            self.assertLessEqual(solved['iterations'],32)
            self.assertLess(solved['residual'],1e-10)
            self.assertAlmostEqual(solved['poly_contact'],.6,places=8)
            self.assertAlmostEqual(solved['bank_outward_contact'],2.,places=8)
            np.testing.assert_allclose(solved['first'],transform((6.,0.)),atol=2e-9,rtol=0.)
            np.testing.assert_allclose(solved['last'],transform((8.,-2.)),atol=2e-9,rtol=0.)
            np.testing.assert_array_equal(before[0],core.controls)
            self.assertEqual(before[1:],(arm.source.points,arm.profile))
            for poly_range,bank_range in (((0.,.4),(0.,arm.length)),((0.,1.),(0.,1.5))):
                with self.assertRaises(ParentChainUnavailable):
                    offset_solution(core,poly_range,arm,bank_range,2.)
            for radius in (0.,-1.,float('nan')):
                with self.assertRaises(ValueError):offset_solution(core,(0.,1.),arm,(0.,10.),radius)

    def test_original_source_knots_and_chord_budget_survive_a_coarser_final_sampler(self):
        core,arm,_=parents();record=primitive(core,arm,offset_solution(core,(0.,1.),arm,(0.,10.),2.))
        sampler=CertifiedCurbSampler(.005,2.,8192)
        at=record.certified_parameters(sampler,0.,1.)
        self.assertTrue(all(t in at for t in record.source_reference_parameters))
        error=record.radius*(1.-np.cos(.5*abs(record.sweep)*np.diff(at)))
        self.assertLessEqual(float(error.max()),record.source_chord_tolerance)
        partial=record.certified_parameters(sampler,.17,.83)
        self.assertEqual(partial[0],.17);self.assertEqual(partial[-1],.83)
        self.assertTrue(all(t in partial for t in record.source_reference_parameters if .17<t<.83))
        np.testing.assert_array_equal(record(0.),record.first);np.testing.assert_array_equal(record(1.),record.last)
        self.assertFalse(hasattr(record.corner,'first'))

    def test_regular_offsets_refuse_folded_core_and_exact_native_miter_crossings(self):
        core,arm,_=parents()
        candidate=SimpleNamespace(polynomial=core,bank=SimpleNamespace(path=arm.path,side=1),bank_arm=arm,bank_outward_interval=(0.,10.))
        self.assertEqual(_regularity(candidate,2.),(0.,0.))
        with self.assertRaisesRegex(ParentChainUnavailable,'regularity'):
            _regularity(replace_namespace(candidate,polynomial=SimpleNamespace(curvature_bound=.5)),2.)
        section=model.Section(2,2,((9.,0.,0.),(9.,-5.,0.),(14.,-5.,0.)),width=2.,interpolation=0)
        bank=Path(section);b=Approach(bank,0,Crown(.08,0.))
        crossing=replace_namespace(candidate,bank=SimpleNamespace(path=bank,side=1),bank_arm=b,bank_outward_interval=(0.,bank.length))
        with self.assertRaisesRegex(ParentChainUnavailable,'miter'):_regularity(crossing,2.)

    def test_actual_endpoint_jets_cannot_be_clamped_to_a_generated_arc_normal(self):
        from tested_python_streets.source_parent_chain_conditions import endpoint_slopes
        from tested_python_streets.source_curb_conditions import SourceCurbConditionUnavailable
        core,arm,_=parents();record=primitive(core,arm,offset_solution(core,(0.,1.),arm,(0.,10.),2.))
        tangent=record.tangent(np.array((0.,1.)));inward=np.column_stack((-tangent[:,1],tangent[:,0]))
        jets=np.array((.03,.08))[:,None]*inward
        np.testing.assert_allclose(endpoint_slopes(record,jets),(.03,.08),atol=1e-15)
        for invalid in (jets+np.array((.001,0.))[:,None]*tangent,-jets):
            with self.assertRaisesRegex(SourceCurbConditionUnavailable,'tangent/inward'):
                endpoint_slopes(record,invalid)

    def test_whole_arc_guard_rejects_clipped_curbs_and_exposed_internal_closure(self):
        core,arm,_=parents();record=primitive(core,arm,offset_solution(core,(0.,1.),arm,(0.,10.),2.))
        area=SimpleNamespace(region=SimpleNamespace(owner=record.owner),layout=SimpleNamespace(source_parent_chain_fillet=record))
        group=CompoundGroup((),(),(area,),())
        def view(label):
            loop=BoundaryLoop(np.array(((0.,0.),(1.,0.),(0.,1.))),tuple(RimSpan(record.owner,item) for item in label))
            return replace(group,loops=(loop,))
        good=view((RimSegment('curb',0,0.,.5),RimSegment('curb',0,.5,1.)))
        self.assertIs(validate(good),good)
        for labels in ((RimSegment('curb',0,.01,1.),),(RimSegment('curb',0,0.,.99),),
                       (RimSegment('curb',0,0.,.4),RimSegment('curb',0,.5,1.)),
                       (RimSegment('curb',0,0.,1.),RimSegment('internal',1,0.,1.))):
            with self.assertRaisesRegex(ParentChainUnavailable,'clips|closure'):validate(view(labels))


def replace_namespace(value,**changes):
    return SimpleNamespace(**dict(vars(value),**changes))


if __name__=='__main__':unittest.main()
