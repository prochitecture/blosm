"""Run with standalone Python. Test geometry conventions independently of bpy."""
import json
import math
import sys
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mapillary_test_support import fixture, mvt, field, BOUNDS, BASE_ID
from mapillary.tiles import image_points, identifier, tile_coordinates, fields, MAX_TILE_BYTES
from mapillary.pose import interpret, rotation_vector, matrix_product, transpose
import importlib.util
spec = importlib.util.spec_from_file_location('blosm_projection', Path(__file__).resolve().parents[1] / 'util/transverse_mercator.py')
projection_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(projection_module)
TransverseMercator = projection_module.TransverseMercator
from progressive.core import Job, Download, Drain


class PoseChecks(unittest.TestCase):
    def setUp(self):
        self.projection = TransverseMercator(lat=60.17, lon=24.94)
        self.coverage = {'coverage_geometry': [24.94, 60.17], 'compass_angle': 0.}

    def pose(self, **metadata):
        metadata.setdefault('exif_orientation',1)
        return interpret(metadata, self.coverage, self.projection, 20.)

    def vector(self, rotation, local):
        return [sum(a*b for a,b in zip(row,local)) for row in rotation]

    def assertVector(self, a, b):
        for x, y in zip(a, b):
            self.assertAlmostEqual(x, y, places=6)

    def test_cardinal_heading(self):
        for heading, direction in ((0,(0,1,0)), (90,(1,0,0)), (180,(0,-1,0)), (270,(-1,0,0))):
            p = self.pose(compass_angle=heading)
            self.assertVector(self.vector(p.rotation, (0,0,-1)), direction)
            self.assertTrue(p.properties['pose_approximate'])

    def test_full_rotation_inverts_world_to_camera(self):
        # Independent inverse Rodrigues expectation: spatial.viewingDirection
        # rotates +Z by the negative rotation vector. Never Euler angles.
        for vector in ((0,0,0), (.7,-.3,.4), (math.pi/2,0,0)):
            p = self.pose(computed_rotation=vector, computed_altitude=32.)
            inverse = rotation_vector(tuple(-v for v in vector))
            self.assertVector(self.vector(p.rotation, (0,0,-1)), self.vector(inverse,(0,0,1)))
            self.assertVector(self.vector(p.rotation, (0,1,0)), self.vector(inverse,(0,-1,0)))
            self.assertEqual(p.position, (0.,0.,12.))
            self.assertTrue(p.properties['pose_complete'])

    def test_exif_applied_once(self):
        # Transform.upVector() gives these explicit raw-camera up directions.
        up = {1:(0,-1,0), 3:(0,1,0), 6:(-1,0,0), 8:(1,0,0)}
        for exif, direction in up.items():
            p = self.pose(computed_rotation=[0,0,0], exif_orientation=exif)
            self.assertVector(self.vector(p.rotation,(0,1,0)), direction)
            self.assertVector(self.vector(p.rotation,(0,0,-1)), (0,0,1))
        for exif in up:
            self.assertVector(self.vector(self.pose(compass_angle=90,exif_orientation=exif).rotation,(0,0,-1)), (1,0,0))

    def test_intrinsics_and_missing_data(self):
        p = self.pose(camera_type='perspective', camera_parameters=[.85,0,0], width=4000,height=3000,exif_orientation=6)
        self.assertAlmostEqual(p.lens, 40.8)
        self.assertAlmostEqual(p.sensor_height, 48.)
        self.assertEqual(p.position[2], 0.)
        self.assertEqual(p.properties['altitude_source'], 'missing_scene_plane')
        self.assertEqual(self.pose(altitude=0).position[2], -20.)
        self.assertEqual(self.pose(camera_type='spherical').camera_type,'PANO')
        self.assertIn('proxy', self.pose(camera_type='fisheye',camera_parameters=[1,0,0],width=10,height=10).properties['projection_status'])
        self.assertTrue(self.pose(computed_rotation=[0,0,0],exif_orientation=2).properties['exif_mirrored'])
        self.assertTrue(self.pose(computed_rotation=[0,0,0],altitude=32,exif_orientation=None).properties['pose_approximate'])

    def test_projection_and_position_preference(self):
        p = self.pose(geometry={'coordinates':[24.9401,60.1701]}, computed_geometry={'coordinates':[24.9402,60.1702]}, altitude=30)
        expected = self.projection.fromGeographic(60.1702,24.9402)
        self.assertVector(p.position,(expected[0],expected[1],10.))
        self.assertEqual(p.properties['position_source'],'computed')


class TileChecks(unittest.TestCase):
    def test_dense_live_style_dictionaries_and_unused_layers(self):
        tile,data,_=fixture(1)
        layer=bytes(next(fields(memoryview(data)))[2])
        # Real coverage exceeded both the original 4 MiB resource cap and the
        # original small dictionary cap. Keep this regression payload in memory.
        dense=layer+field(4,field(5,1))*280000
        irrelevant=field(1,b'sequence')+field(4,field(5,1))*1050000
        payload=field(3,irrelevant)+field(3,dense)
        self.assertGreater(len(payload),4*1024*1024)
        self.assertLess(len(payload),MAX_TILE_BYTES)
        records=[r for r in image_points(payload,tile) if r is not None]
        self.assertEqual(records[0]['id'],str(BASE_ID))

    def test_decode_extent_ids_duplicates_and_checkpoints(self):
        tile, data, entries = fixture(96)
        decoded = list(image_points(data,tile))
        self.assertIn(None, decoded)
        records = [r for r in decoded if r is not None]
        self.assertEqual(len(records),97)
        self.assertEqual(records[0]['id'], str(BASE_ID))
        self.assertEqual(records[-1]['id'], str(BASE_ID))
        for extent in (512,8192):
            record = next(r for r in image_points(mvt([(BASE_ID,24.94,60.17)],tile,extent),tile) if r)
            self.assertAlmostEqual(record['coverage_geometry'][0],24.94,places=4)
            self.assertAlmostEqual(record['coverage_geometry'][1],60.17,places=4)

    def test_invalid_and_bounds(self):
        for value in (float(BASE_ID), True, '1e12', '-1',None):
            with self.assertRaises(ValueError): identifier(value)
        tile = next(tile_coordinates(BOUNDS))
        with self.assertRaises(ValueError): list(image_points(b'\x1a\xff',tile))
        with self.assertRaises(ValueError): list(tile_coordinates((-180,-80,180,80)))


class ConsumerChecks(unittest.TestCase):
    def test_empty_authorization_redaction_is_bounded(self):
        from downloads.diagnostics import redact_text
        self.assertEqual(redact_text('failure',headers={'Authorization':'OAuth '}),'failure')

    def test_cooperative_consumer_pause_drain_and_close(self):
        class Transport:
            max_requests = 1
            pending = None
            def submit(self, request): self.pending = (1,b'data',None); return 1
            def poll(self): result,self.pending = self.pending,None; return result
            def close(self): self.closed=True
        events=[]
        def consume(data, error):
            try:
                for i in range(100):
                    events.append(i)
                    yield None
            finally: events.append('consumer_closed')
        def steps():
            try:
                yield Download('request',consume)
                yield Drain()
                events.append('finished')
            finally: events.append('producer_closed')
        job=Job(steps(),Transport())
        for _ in range(3): job.step(max_steps=2)
        before=list(events)
        job.poll()
        self.assertEqual(events,before)
        self.assertNotIn('finished',events)
        job.cancel()
        self.assertIn('consumer_closed',events)
        self.assertIn('producer_closed',events)
        self.assertTrue(job.transport.closed)


if __name__=='__main__':
    unittest.main()
