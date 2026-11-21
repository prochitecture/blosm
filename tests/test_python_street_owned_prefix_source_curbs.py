"""Typed source-prefix ports retain their actual native half-profile banks."""
from copy import copy
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_owned_prefix_corner_authority import case
from tested_python_streets.model import Key,Network
from tested_python_streets.updates import BuildContext
from tested_python_streets.compound_cycles import _port_footprint
from tested_python_streets.source_curb_coordinates import CurbCoordinates,CurbCoordinatesUnavailable
from tested_python_streets.source_curb_conditions import SourceCurbConditions,_native_gradient
from tested_python_streets.owned_prefix_source_curbs import validate_curve


def prepared(*,reverse=False,angle=0.):
    context,group,node,corner=case(reverse=reverse,angle=angle,translation=(1300.,-1800.))
    arm=next(a for a in (corner.first,corner.second) if hasattr(a.path,'claim'))
    mouth=next(at for a,at in zip(group.physical_layout.approaches,group.physical_layout.outward_stations)
        if getattr(a.path,'claim',None) is arm.path.claim)
    partitions={piece.section.id:piece.path.samples(0.,piece.path.length,.5,tolerance=.00025,max_samples=20000)
        for piece in arm.path.pieces}
    args=(Key('cycle_port',100),arm,mouth,4,.5,.00025,partitions)
    port=_port_footprint(*args)
    return context,replace(group,footprints=(*group.footprints,port)),port,args


class OwnedPrefixSourceCurbTests(unittest.TestCase):
    def test_explicit_native_ledger_adds_no_vertex_or_geometry_change(self):
        context,group,port,args=prepared()
        with patch('tested_python_streets.source_curb_coordinates.CurbCoordinates',
                side_effect=CurbCoordinatesUnavailable('old unsupported source type')):
            old=_port_footprint(*args)
        np.testing.assert_array_equal(port.region.xy,old.region.xy)
        self.assertEqual(port.segments,old.segments)
        for mapped,original in zip(port.layout.boundaries,old.layout.boundaries):
            np.testing.assert_array_equal(mapped.points,original.points)
            self.assertIs(type(mapped.source_coordinates),CurbCoordinates)
            self.assertTrue(mapped.source_coordinates.validate(mapped.points,context))
            self.assertEqual(mapped.source_coordinates.maximum_chord_error,.00025)

    def test_normal_dispatch_reads_actual_native_half_profile_in_both_storage_directions(self):
        for reverse,angle in ((False,.713),(True,0.)):
            context,group,port,_=prepared(reverse=reverse,angle=angle)
            resolver=SourceCurbConditions(context,group,maximum_slope=.2);levels=[]
            for index in (0,1):
                curve=port.layout.boundaries[index]
                proof=validate_curve(resolver,port,index,curve,0.,1.)
                self.assertTrue(proof.final_mouth_width_proof['maximum_error']<=.002)
                rule=resolver.source_rule(port.region.owner,index,0.,1.)
                self.assertTrue(rule.source_anchor)
                for first,last in zip((0.,*rule.polynomial_breaks()),(*rule.polynomial_breaks(),1.)):
                    height,gradient=rule.sample_piece(np.array((0.,.5,1.)),first,last)
                    metric=(rule.first+.5*(first+last)*(rule.last-rule.first))*curve.length
                    record=rule.records[int(np.clip(np.searchsorted(rule.ends,metric,side='right'),0,len(rule.records)-1))]
                    p0,p1,path,profile,s0,s1,side=record
                    at=(rule.first+(first+np.array((0.,.5,1.))*(last-first))*(rule.last-rule.first))*curve.length
                    native=s0+np.clip((at-p0)/(p1-p0),0.,1.)*(s1-s0)
                    native_interval=s0+np.clip(((rule.first+np.array((first,last))*(rule.last-rule.first))*curve.length-p0)/(p1-p0),0.,1.)*(s1-s0)
                    np.testing.assert_array_equal(height,np.full(3,profile.datum+profile.value(side)))
                    np.testing.assert_allclose(gradient,_native_gradient(path,profile,native,side,interval=native_interval),atol=1e-12,rtol=0.)
                    self.assertNotEqual(path.section.id,curve.source_coordinates.approach.source.id)
                    levels.extend(height)
            self.assertGreater(max(levels)-min(levels),.07)

    def test_equal_reconstructed_network_retains_canonical_native_coordinates(self):
        context,group,port,_=prepared();network=Network.from_dict(context.network.to_dict())
        context=BuildContext(network,context.products,context.patterns,context.task_digests,context.source_digests)
        resolver=SourceCurbConditions(context,group,maximum_slope=.2)
        self.assertIsNotNone(validate_curve(resolver,port,0,port.layout.boundaries[0],0.,1.))
        self.assertTrue(np.isfinite(resolver.source_rule(port.region.owner,0,0.,1.).sample(np.array((0.,1.)))[1]).all())

    def test_changed_vertices_claims_and_borrowed_native_paths_are_refused(self):
        for failure in ('vertices','claim','path'):
            context,group,port,_=prepared();curve=port.layout.boundaries[0]
            if failure=='vertices':curve.points[1,0]+=.001
            elif failure=='claim':curve.source_coordinates.path.claim=replace(curve.source_coordinates.path.claim)
            else:
                route=curve.source_coordinates.path;piece=route.pieces[0]
                route.pieces=(replace(piece,path=copy(piece.path)),*route.pieces[1:])
            with self.assertRaises(ValueError):
                validate_curve(SourceCurbConditions(context,group,maximum_slope=.2),port,0,curve,0.,1.)


if __name__=='__main__':unittest.main()
