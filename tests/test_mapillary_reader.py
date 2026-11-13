"""Standalone regressions for the bounded count/replay tile reader."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from mapillary_test_support import fixture, mvt, field, BOUNDS, BASE_ID
from mapillary import tiles


def records(iterator):
    return [value for value in iterator if value is not None]


class ReaderChecks(unittest.TestCase):
    def test_count_replay_matches_uncached_records_and_releases(self):
        tile, data, _ = fixture(96)
        expected = records(tiles.image_points(data,tile,BOUNDS))
        reader = tiles.ImageTile(data,tile,BOUNDS)
        try:
            identities = records(reader.image_ids())
            self.assertTrue(reader.index_ready)
            self.assertEqual(identities,[row['id'] for row in expected])
            self.assertEqual(records(reader.image_points()),expected)
            self.assertEqual(records(reader.image_ids()),identities)
            self.assertGreater(reader.index_bytes,0)
            self.assertLess(reader.index_bytes,len(data))
            self.assertEqual(reader.layers[0].values.offsets.typecode,'I')
            self.assertEqual(reader.layers[0].selected.typecode,'I')
        finally:
            reader.close()
        self.assertIsNone(reader.data)
        self.assertEqual(reader.index_bytes,0)
        with self.assertRaisesRegex(ValueError,'closed'):
            records(reader.image_points())

    def test_replay_does_not_rescan_layer_or_rebuild_dictionary(self):
        tile, data, _ = fixture(16)
        reader = tiles.ImageTile(data,tile,BOUNDS)
        try:
            identities = records(reader.image_ids())
            with patch.object(tiles,'_spans',side_effect=AssertionError('Unexpected full scan')):
                actual = records(reader.image_points())
            self.assertEqual([row['id'] for row in actual],identities)
        finally:
            reader.close()

    def test_bounds_extents_duplicates_and_multiple_image_layers(self):
        for extent in (512,4096,8192):
            tile = next(tiles.tile_coordinates(BOUNDS))
            points = [(BASE_ID+i,24.938+i*.0001,60.168+i*.0001) for i in range(40)]
            data = mvt(points[:20],tile,extent,duplicate=True)+mvt(points[20:],tile,extent)
            reference = records(tiles.image_points(data,tile))
            for bounds in (BOUNDS,(24.938,60.168,24.94,60.170),(24.941,60.171,24.944,60.174)):
                selected = [row for row in reference if tiles.contains(bounds,*row['coverage_geometry'])]
                reader = tiles.ImageTile(data,tile,bounds)
                try:
                    self.assertEqual(records(reader.image_ids()),[row['id'] for row in selected])
                    self.assertEqual(records(reader.image_points()),selected)
                finally:
                    reader.close()

    def test_selected_value_validation_and_unused_scalar_skipping(self):
        tile, data, _ = fixture(1)
        layer = bytes(next(tiles.fields(memoryview(data)))[2])
        # A dictionary entry's outer framing is valid, its nested scalar is not.
        malformed = field(4,b'\x28\x80')
        unused = field(3,layer+malformed)
        self.assertEqual(records(tiles.image_ids(unused,tile,BOUNDS)),[str(BASE_ID),str(BASE_ID)])
        raw_feature = field(2,b'\x00\x00')+field(3,1)+field(4,b'\x09\x00\x00')
        selected = field(3,field(1,b'image')+field(2,raw_feature)+field(3,b'id')+malformed)
        with self.assertRaisesRegex(ValueError,'Truncated'):
            records(tiles.image_ids(selected,tile,(-180,-85,180,85)))
        # A malformed non-ID attribute referenced by a selected feature also
        # fails during counting, before the tile becomes a usable sample.
        raw_feature = field(2,b'\x00\x00\x01\x01')+field(3,1)+field(4,b'\x09\x00\x00')
        selected = field(3,field(1,b'image')+field(2,raw_feature)+field(3,b'id')+
                         field(3,b'captured_at')+field(4,field(5,BASE_ID))+malformed)
        with self.assertRaisesRegex(ValueError,'Truncated'):
            records(tiles.image_ids(selected,tile,(-180,-85,180,85)))

    def test_large_integer_ids_and_invalid_float_ids(self):
        import struct
        tile = next(tiles.tile_coordinates(BOUNDS))
        raw_feature = field(2,b'\x00\x00')+field(3,1)+field(4,b'\x09\x00\x00')
        for value, valid in ((field(5,BASE_ID),True),(field(3,struct.pack('<d',float(BASE_ID))),False)):
            data = field(3,field(1,b'image')+field(2,raw_feature)+field(3,b'id')+field(4,value))
            if valid:
                self.assertEqual(records(tiles.image_ids(data,tile,(-180,-85,180,85))),[str(BASE_ID)])
            else:
                with self.assertRaisesRegex(ValueError,'identifier'):
                    records(tiles.image_ids(data,tile,(-180,-85,180,85)))

    def test_dictionary_selection_layer_and_payload_caps(self):
        tile, data, _ = fixture(2)
        with patch.object(tiles,'MAX_TABLE',7):
            with self.assertRaisesRegex(ValueError,'dictionary'):
                records(tiles.image_ids(data,tile,BOUNDS))
        _, one, _ = fixture(1)
        layer = bytes(next(tiles.fields(memoryview(one)))[2])
        raw_feature = next(value for number, wire, value in tiles.fields(memoryview(layer)) if number==2)
        many = field(3,layer+field(2,bytes(raw_feature))*8)
        with patch.object(tiles,'MAX_TABLE',5):
            with self.assertRaisesRegex(ValueError,'selected feature'):
                records(tiles.image_ids(many,tile,BOUNDS))
        with patch.object(tiles,'MAX_IMAGE_LAYERS',1):
            with self.assertRaisesRegex(ValueError,'layer count'):
                records(tiles.image_ids(one+one,tile,BOUNDS))
        with patch.object(tiles,'MAX_TILE_BYTES',len(data)-1):
            with self.assertRaisesRegex(ValueError,'32 MiB'):
                tiles.ImageTile(data,tile,BOUNDS)

    def test_cancelled_count_does_not_publish_partial_replay_index(self):
        tile, data, _ = fixture(96)
        reader = tiles.ImageTile(data,tile,BOUNDS)
        iterator = reader.image_ids()
        next(value for value in iterator if value is not None)
        iterator.close()
        self.assertFalse(reader.index_ready)
        self.assertEqual(sum(len(layer.selected) for layer in reader.layers),0)
        try:
            self.assertEqual(len(records(reader.image_points())),97)
            self.assertFalse(reader.index_ready)
        finally:
            reader.close()

    def test_cancelled_dictionary_preparation_releases_partial_local_layer(self):
        import weakref
        tile, data, _ = fixture(96)
        reader = tiles.ImageTile(data,tile,BOUNDS)
        observed = []
        original = tiles._Layer
        def capture(payload):
            layer = original(payload)
            observed.append(weakref.ref(layer))
            return layer
        iterator = reader.image_ids()
        try:
            with patch.object(tiles,'_Layer',capture):
                for _ in range(1000):
                    next(iterator)
                    if observed and len(observed[0]().values.offsets):
                        break
                else:
                    self.fail('Did not pause within dictionary preparation')
                self.assertFalse(reader.prepared)
                self.assertEqual(reader.layers,[])
                iterator.close()
                self.assertIsNone(observed[0]())
            reader.close()
            self.assertIsNone(reader.data)
            self.assertEqual(reader.index_bytes,0)
        finally:
            iterator.close()
            reader.close()

    def test_coarse_rejection_preserves_geographic_edges_and_tile_buffers(self):
        import math
        from mapillary_test_support import vint
        for tile in ((14,0,0),(14,8192,8192),(14,16383,16383)):
            for extent in (512,4096,1 << 24):
                positions = (-extent//10,-1,0,1,extent//2,extent-1,extent,extent+1,extent*11//10)
                features, values = [], []
                for i, (x,y) in enumerate((x,y) for x in positions for y in positions):
                    zig = lambda value: value*2 if value >= 0 else -value*2-1
                    geometry = b''.join(map(vint,(9,zig(x),zig(y))))
                    features.append(field(2,field(2,vint(0)+vint(i))+field(3,1)+field(4,geometry)))
                    values.append(field(4,field(5,BASE_ID+i)))
                data = field(3,field(1,b'image')+b''.join(features)+field(3,b'id')+
                             b''.join(values)+field(5,extent))
                reference = records(tiles.image_points(data,tile))
                for point in reference[::10]:
                    lon, lat = point['coverage_geometry']
                    # Exact decoded edges and adjacent floating-point values
                    # retain the geographic comparison as final authority.
                    for w, s in ((lon,lat),(math.nextafter(lon,math.inf),lat),
                                 (lon,math.nextafter(lat,math.inf))):
                        bounds = (w,s,lon+1e-7,lat+1e-7)
                        selected = [row for row in reference if tiles.contains(bounds,*row['coverage_geometry'])]
                        reader = tiles.ImageTile(data,tile,bounds)
                        try:
                            self.assertEqual(records(reader.image_ids()),[row['id'] for row in selected])
                            self.assertEqual(records(reader.image_points()),selected)
                        finally:
                            reader.close()

    def test_dictionary_tail_is_not_traversed_twice_and_late_features_work(self):
        tile, data, _ = fixture(2)
        reader = tiles.ImageTile(data,tile,BOUNDS)
        span_lengths = []
        original = tiles._spans
        def observed(payload):
            span_lengths.append(len(payload))
            yield from original(payload)
        try:
            with patch.object(tiles,'_spans',observed):
                identities = records(reader.image_ids())
            self.assertEqual(len(span_lengths),2)
            self.assertLess(span_lengths[1],span_lengths[0])
            layer = bytes(next(tiles.fields(memoryview(data)))[2])
            last_feature = next(value for number,wire,value in tiles.fields(memoryview(layer)) if number==2)
            late = field(3,layer+field(2,bytes(last_feature)))
            reader.close()
            reader = tiles.ImageTile(late,tile,BOUNDS)
            self.assertEqual(records(reader.image_ids()),identities+[identities[0]])
            self.assertEqual([row['id'] for row in records(reader.image_points())],identities+[identities[0]])
        finally:
            reader.close()

    def test_fixed_point_decoder_rejects_truncated_and_multiple_geometry(self):
        tile = next(tiles.tile_coordinates(BOUNDS))
        for geometry in (b'',b'\x09',b'\x09\x00',b'\x09\x00\x00\x00',b'\x0a\x00\x00'):
            data = field(3,field(1,b'image')+field(2,field(2,b'\x00\x00')+field(3,1)+field(4,geometry))+
                         field(3,b'id')+field(4,field(5,BASE_ID)))
            with self.assertRaises(ValueError):
                records(tiles.image_ids(data,tile,BOUNDS))


if __name__=='__main__':
    unittest.main()
