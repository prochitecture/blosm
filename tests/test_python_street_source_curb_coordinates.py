"""Canonical source identity survives prepared strip and port polylines."""
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Key
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_cycles import _footprint,_port_footprint
from tested_python_streets.source_curb_coordinates import CurbCoordinates,CurbCoordinatesUnavailable


def bent(mode=0,reverse=False):
    points=((0.,0.,0.),(8.,0.,0.),(8.,8.,0.))
    if reverse:points=points[::-1]
    section=Section(1,1,points,width=4.,interpolation=mode)
    return Path(section,.005)


class SourceCurbCoordinateTests(unittest.TestCase):
    def test_full_strip_keeps_every_source_station_side_and_geometry_bit(self):
        for mode in (0,3):
            for reverse in (False,True):
                path=bent(mode,reverse)
                at=path.samples(0.,path.length,.4,tolerance=.00025,max_samples=20000)
                right,left=path.point(at,2.),path.point(at,-2.)[::-1]
                points=np.vstack((right,left));count=len(at)
                ledger=CurbCoordinates(path,np.r_[at,at[::-1]],
                    np.r_[np.full(count,2.),np.full(count,-2.)],
                    np.r_[np.ones(count),-np.ones(count)])
                legacy=_footprint(Key('strip',1),points,internal_edges=(count-1,2*count-1))
                mapped=_footprint(Key('strip',1),points,internal_edges=(count-1,2*count-1),
                                  source_coordinates=ledger)
                np.testing.assert_array_equal(mapped.region.xy,legacy.region.xy)
                self.assertEqual(mapped.segments,legacy.segments)
                for old,new in zip(legacy.layout.boundaries,mapped.layout.boundaries):
                    np.testing.assert_array_equal(old.points,new.points)
                    np.testing.assert_array_equal(old(np.linspace(0.,1.,101)),new(np.linspace(0.,1.,101)))
                    coordinates=new.source_coordinates
                    self.assertTrue(coordinates.validate(new.points,{Key('path',1):path}))
                    self.assertFalse(coordinates.stations.flags.writeable)
                    for edge in range(len(new.points)-1):
                        piece,=coordinates.edge_sources(edge)
                        self.assertIs(piece.path,path)
                        self.assertEqual(piece.first_station,coordinates.stations[edge])
                        self.assertEqual(piece.last_station,coordinates.stations[edge+1])
                        if coordinates.sides[edge]==coordinates.sides[edge+1]:
                            self.assertEqual(piece.authored_side,int(coordinates.sides[edge]))
                        else:self.assertIsNone(piece.authored_side)

    def test_route_chord_is_split_into_actual_oriented_sources_without_new_vertices(self):
        first=Path(Section(1,1,((0.,0.,0.),(4.,0.,0.)),width=4.,interpolation=0),.005)
        second=Path(Section(2,1,((9.,0.,0.),(4.,0.,0.)),width=4.,interpolation=0),.005)
        route=Route((Piece(first,True),Piece(second,False)))
        approach=Approach(route,0,Crown(.08,0.))
        partitions={1:np.array((0.,1.5,4.)),2:np.array((0.,2.,5.))}
        mapped=_port_footprint(Key('port',1),approach,6.3,4,.7,.00025,partitions)
        with patch('tested_python_streets.source_curb_coordinates.CurbCoordinates',
                   side_effect=CurbCoordinatesUnavailable('explicit unsupported test')):
            legacy=_port_footprint(Key('port',1),approach,6.3,4,.7,.00025,partitions)
        np.testing.assert_array_equal(mapped.region.xy,legacy.region.xy)
        self.assertEqual(mapped.segments,legacy.segments)
        for new,old in zip(mapped.layout.boundaries,legacy.layout.boundaries):
            np.testing.assert_array_equal(new.points,old.points)
            new.source_coordinates.validate(new.points,{Key('path',1):first,Key('path',2):second})
        curve=mapped.layout.boundaries[0]
        crossing=[curve.source_coordinates.edge_sources(i) for i in range(len(curve.points)-1)
                  if len(curve.source_coordinates.edge_sources(i))>1]
        self.assertEqual(len(crossing),1)
        left,right=crossing[0]
        self.assertIs(left.path,first);self.assertIs(right.path,second)
        self.assertEqual(left.authored_side,1);self.assertEqual(right.authored_side,-1)
        self.assertEqual(left.last_station,4.);self.assertEqual(right.first_station,5.)
        self.assertEqual(left.last_fraction,right.first_fraction)
        self.assertEqual(left.first_fraction,0.);self.assertEqual(right.last_fraction,1.)
        self.assertEqual(left.first_offset,2.);self.assertEqual(right.first_offset,-2.)

    def test_false_source_vertex_maps_and_nonfinite_ranges_are_rejected(self):
        path=bent();at=np.array((0.,4.,8.));offset=np.full(3,2.)
        ledger=CurbCoordinates(path,at,offset,np.ones(3))
        points=path.point(at,offset)
        self.assertTrue(ledger.validate(points))
        changed=points.copy();changed[1,0]+=.001
        with self.assertRaisesRegex(ValueError,'physical vertices'):ledger.validate(changed)
        with self.assertRaisesRegex(ValueError,'ownership'):
            ledger.validate(points,{Key('path',1):bent()})
        with self.assertRaisesRegex(ValueError,'bank side'):
            CurbCoordinates(path,at,offset,-np.ones(3))
        for stations in ((0.,np.nan,8.),(-.001,4.,8.),(0.,4.,path.length+.001)):
            with self.assertRaises(ValueError):CurbCoordinates(path,stations,offset,np.ones(3))
        for indices in ((0,-1),(0,3),(.0,1.)):
            with self.assertRaises(ValueError):ledger.subset(indices)
        for edge in (-1,2,True):
            with self.assertRaises(ValueError):ledger.edge_sources(edge)
        with self.assertRaises(CurbCoordinatesUnavailable):
            CurbCoordinates(SimpleNamespace(length=10.),at,offset,np.ones(3))
        with self.assertRaises(CurbCoordinatesUnavailable):
            CurbCoordinates(Route((Piece(path,True),)),at,offset,np.ones(3))

    def test_reversed_native_approach_retains_authored_station_and_side(self):
        path=bent();approach=Approach(path,1,Crown(.08,0.))
        at=np.array((0.,4.,8.));offset=np.full(3,2.)
        ledger=CurbCoordinates(path,at,offset,np.ones(3),approach=approach)
        self.assertTrue(ledger.validate(approach.point(at,offset),{Key('path',1):path}))
        first,=ledger.edge_sources(0)
        self.assertEqual(first.first_station,path.length)
        self.assertEqual(first.last_station,path.length-4.)
        self.assertEqual(first.first_offset,-2.)
        self.assertEqual(first.authored_side,-1)

    def test_policy_change_rebuilds_real_plan_and_projected_surfaces_with_exact_geometry(self):
        from _python_street_fixtures import fixture_networks
        from tested_python_streets import model
        from tested_python_streets.backend import PythonMeshBackend
        from tested_python_streets.ownership import PLAN
        network=fixture_networks(model)[0]
        backend=PythonMeshBackend(network.settings)
        initial=backend.prepare(network);backend.commit(initial)
        surfaces={key:value for key,value in initial.products.items() if key.kind=='junction_surface'}
        path=initial.products[Key('path',network.sections[0].id)]
        self.assertTrue(surfaces)
        with patch('tested_python_streets.ownership.SOURCE_CURB_COORDINATE_POLICY','next-source-parameter-policy'):
            changed=backend.prepare(network)
        self.assertIn(PLAN,changed.rebuilt)
        self.assertIs(changed.products[Key('path',network.sections[0].id)],path)
        for key,product in surfaces.items():
            self.assertIn(key,changed.rebuilt)
            if product is not None:
                np.testing.assert_array_equal(changed.products[key].pavement.mesh.xyz,product.pavement.mesh.xyz)
                np.testing.assert_array_equal(changed.products[key].pavement.mesh.quads,product.pavement.mesh.quads)


if __name__=='__main__':unittest.main()
