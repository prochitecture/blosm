"""A tiny exposed circular auxiliary lead permits sampling, not radius changes."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Key,Section,Settings
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.compounds import CompoundGroup,SharedClaim
from tested_python_streets.regions import RegionUnion
from tested_python_streets.source_exterior_sampling import tiny_auxiliary_edges,recover_sampling
from tested_python_streets.quad_charts import QuadChartLayerBudgetError
from tested_python_streets.partition import QuadChartGradingError


def fixture(half=False,angle=0.,shift=(0.,0.),turn=.004,radius=3.4):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    def path(sid,direction):
        points=np.array(((0.,0.),direction))@rotation.T+shift
        return Path(Section(sid,sid,tuple((float(x),float(y),0.) for x,y in points),width=7.))
    a=Approach(path(1,(40.,0.)),0,Crown(.08,0.,-1 if half else 0))
    gap=math.pi+turn
    b=Approach(path(2,(40.*math.cos(gap),40.*math.sin(gap))),0,Crown(.08,0.,1 if half else 0))
    # R=3.4 leaves a genuine -0.2 mm original offset lead at both contacts.
    corner=CornerRadiusPolicy.solve(a,b,radius,'STRICT')
    source=corner_patch(Key('test_corner',1),corner,.7,.00025)
    loops=RegionUnion().build((source.region,))
    layout=SimpleNamespace(source_cycle_full_width_proof={'maximum_error':.002})
    claim=SharedClaim(1,Key('intersection',1),Key('intersection',2),1.,0.)
    return CompoundGroup((claim.first_owner,claim.last_owner),(claim,),(source,),loops,layout),corner


def positive_fixture(angle=0.,shift=(0.,0.),authored_knot=False):
    from tested_python_streets.routes import Route,Piece
    from tested_python_streets.compound_cycles import _port_footprint
    base,old=fixture(angle=angle,shift=shift,turn=-.0004)
    def approach(old):
        path=old.path
        if authored_knot and old.source.id==1:
            points=list(path.section.points)
            middle=np.asarray(points[0])+.0001/40.*(np.asarray(points[-1])-points[0])
            path=Path(replace(path.section,points=(points[0],tuple(middle),points[-1])))
        return Approach(Route((Piece(path,True),)),0,old.profile)
    a,b=approach(old.first),approach(old.second)
    corner=CornerRadiusPolicy.solve(a,b,3.4,'STRICT')
    source=corner_patch(Key('test_corner',1),corner,.7,.00025)
    ports=[]
    for i,app in enumerate((a,b)):
        s=corner.stations[i]
        port=_port_footprint(Key('cycle_port',app.source.id),app,5.,4,.7,.00025,
            {app.source.id:np.array((0.,s/2.,s,5.))})
        if i==1:
            # Retain its genuine strip but isolate this test to the first
            # original source's positive contact-phase provenance.
            layout=copy(port.layout);layout.original_root_caps={}
            port=replace(port,layout=layout)
        ports.append(port)
    # These are original source-prefix operands for the incidence test only;
    # this fixture does not claim or qualify a combined physical exterior.
    return tuple(ports),corner


class SourceExteriorSamplingTests(unittest.TestCase):
    def test_real_negative_circle_leads_and_equal_plateau_survive_rotation(self):
        for half in (False,True):
            for angle in (0.,.713):
                group,corner=fixture(half,angle,(1200.,-1800.))
                before=corner.center.copy();stations=corner.stations
                rows=tiny_auxiliary_edges(group,Settings())
                self.assertTrue(rows)
                self.assertTrue(all(row['original_radius']==3.4 for row in rows))
                self.assertTrue(all(row['facing_height_policy']==('ORIGINAL_EQUAL_PLATEAU' if half
                    else 'ORIGINAL_ZERO_OUTER') for row in rows))
                np.testing.assert_array_equal(corner.center,before)
                self.assertEqual(corner.stations,stations)

    def test_up_to_five_mm_requires_generated_negative_root_and_contact_knots(self):
        group,corner=fixture(turn=.01,radius=2.8)
        rows=tiny_auxiliary_edges(group,Settings())
        self.assertTrue(rows)
        self.assertTrue(all(.003<abs(row['phase_knots']['source_contact_station'])<.005 for row in rows))
        self.assertTrue(all(row['generated_contact_phase']=='NEGATIVE_AUXILIARY_CONTACT' for row in rows))
        beyond,_=fixture(turn=.02,radius=2.8)
        self.assertFalse(tiny_auxiliary_edges(beyond,Settings()))

    def test_positive_phase_needs_same_original_source_profile_side_and_stored_knots(self):
        from tested_python_streets.source_exterior_sampling import _positive_source_phase
        from tested_python_streets.compounds import RimSegment
        for angle in (0.,.713):
            ports,corner=positive_fixture(angle,(1200.,-1800.))
            port=ports[0]
            curve=port.layout.boundaries[1]
            label=RimSegment('curb',1,float(curve.stations[-2]/curve.length),1.)
            phase=lambda source,span=label:_positive_source_phase(source,span,corner.first,corner.stations[0],-1)
            self.assertIsNotNone(phase(port))
            self.assertIsNone(phase(port,replace(label,index=0)))
            layout=copy(port.layout);layout.original_root_caps={}
            self.assertIsNone(phase(replace(port,layout=layout)))
            layout=copy(port.layout);wrong=copy(layout.approaches[0]);wrong.path=copy(wrong.path)
            layout.approaches=(wrong,)
            self.assertIsNone(phase(replace(port,layout=layout)))
            from tested_python_streets.compound_cycles import PolylineCurb
            layout=copy(port.layout)
            # The generated contact, not a nearby arbitrary point, is compulsory.
            contact=corner.first.point(corner.stations[0],-.5*corner.first.width)
            at=int(np.linalg.norm(curve.points-contact,axis=1).argmin())
            changed=PolylineCurb(np.delete(curve.points,at,axis=0))
            layout.boundaries=(layout.boundaries[0],changed,layout.boundaries[2])
            self.assertIsNone(phase(replace(port,layout=layout)))
        authored,c=positive_fixture(authored_knot=True)
        curve=authored[0].layout.boundaries[1]
        label=RimSegment('curb',1,float(curve.stations[-2]/curve.length),1.)
        self.assertIsNone(_positive_source_phase(authored[0],label,c.first,c.stations[0],-1))

    def test_reversed_curved_authored_control_uses_metric_station(self):
        from tested_python_streets.source_exterior_sampling import _generated_original_phase
        from tested_python_streets.routes import Route,Piece
        path=Path(Section(8,8,((0.,0.,0.),(10.,10.,0.),(20.,0.,0.),(20.0001,0.,0.)),width=2.))
        near=float(path.length-path.station(path.knots[-2]))
        self.assertLess(near,.005)
        self.assertGreater(path.length-float(path.knots[-2]),.005)
        for a in (Approach(path,1,Crown(.08,0.)),
                  Approach(Route((Piece(path,False),)),1,Crown(.08,0.))):
            self.assertFalse(_generated_original_phase(a,near*1.1))

    def test_large_negative_leads_hidden_sources_wrong_node_and_unequal_plateau_reject(self):
        group,corner=fixture()
        self.assertFalse(tiny_auxiliary_edges(replace(group,claims=()),Settings()))
        source=group.footprints[0]
        unexposed=replace(source,region=replace(source.region,owner=Key('hidden_corner',1)))
        plain=copy(source.layout);plain.corners=()
        self.assertFalse(tiny_auxiliary_edges(replace(group,footprints=(
            replace(source,layout=plain),unexposed)),Settings()))
        for field,value in (('center',None),('stations',(-16.,-28.))):
            changed=replace(corner,**{field:value})
            layout=copy(group.footprints[0].layout);layout.corners=(changed,)
            altered=replace(group,footprints=(replace(group.footprints[0],layout=layout),))
            self.assertFalse(tiny_auxiliary_edges(altered,Settings()))
        shifted=Approach(Path(Section(2,2,((1.,0.,0.),(-40.,0.,0.)),width=7.)),0,Crown(.08,0.))
        layout=copy(group.footprints[0].layout);layout.corners=(replace(corner,second=shifted),)
        self.assertFalse(tiny_auxiliary_edges(replace(group,footprints=(replace(group.footprints[0],layout=layout),)),Settings()))
        plateau,c=fixture(True)
        b=copy(c.second);b.profile=Crown(.09,0.,1)
        layout=copy(plateau.footprints[0].layout);layout.corners=(replace(c,second=b),)
        self.assertFalse(tiny_auxiliary_edges(replace(plateau,footprints=(replace(plateau.footprints[0],layout=layout),)),Settings()))

    def test_unrelated_errors_and_repeated_sampling_do_not_enter_rebuild(self):
        group,_=fixture();context=SimpleNamespace(network=SimpleNamespace(settings=Settings()))
        failure=ValueError('actual original source escaped')
        with self.assertRaises(ValueError) as caught:recover_sampling(None,group,context,failure)
        self.assertIs(caught.exception,failure)
        group.physical_layout.source_rim_normalization=({'additional_error_allowance':1e-6},)
        self.assertFalse(tiny_auxiliary_edges(group,Settings()))
        failure=QuadChartLayerBudgetError('Direct quad chart exceeds its layer budget')
        with self.assertRaises(QuadChartLayerBudgetError):recover_sampling(None,group,context,failure)

    def test_source_coverage_failure_propagates_before_any_sampling_or_rebuild(self):
        group,_=fixture();context=SimpleNamespace(network=SimpleNamespace(settings=Settings()))
        from tested_python_streets.source_strip_containment import SourceStripContainmentError
        failure=SourceStripContainmentError('original finite source width escapes')
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',side_effect=failure):
            with patch('tested_python_streets.source_rim_normalization.normalize_exterior') as normalize:
                with self.assertRaises(SourceStripContainmentError) as caught:
                    recover_sampling(None,group,context,QuadChartGradingError('finite chart budget'))
        self.assertIs(caught.exception,failure)
        normalize.assert_not_called()

    def test_one_retry_requires_fresh_consumed_and_prefix_proof_before_after_and_final(self):
        group,_=fixture();context=SimpleNamespace(network=SimpleNamespace(settings=Settings()))
        normalized=replace(group,physical_layout=copy(group.physical_layout))
        normalized.physical_layout.source_rim_normalization=({'maximum_projection_error':.5e-6},)
        normalized.physical_layout.source_exterior_sampling={'original_holes_exact':True}
        output=SimpleNamespace(metadata={},partition=SimpleNamespace(group=normalized))
        parent=SimpleNamespace(max_cells=50000,build=lambda g,c:output)
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',return_value={'maximum_error':.002}) as proof:
            with patch('tested_python_streets.source_rim_normalization.normalize_exterior',return_value=normalized):
                result=recover_sampling(parent,group,context,QuadChartGradingError('finite grading budget'))
        self.assertEqual(proof.call_count,3)
        report=result.metadata['derived_exterior_sampling_recovery']
        self.assertEqual(report['maximum_rebuilds'],1)
        self.assertLessEqual(report['composed_curve_error_bound'],report['authored_curve_budget'])


if __name__=='__main__':unittest.main()
