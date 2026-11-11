"""Synthetic MVT/Graph fixtures built in memory, never saved in the repository."""
from pathlib import Path
import json
import math
import struct
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mapillary.tiles import ZOOM, tile_coordinates

BOUNDS = (24.939, 60.169, 24.941, 60.171)
BASE_ID = 9007199254740993


def vint(v):
    data = bytearray()
    while v > 127:
        data.append((v & 127) | 128)
        v >>= 7
    data.append(v)
    return bytes(data)


def field(index, data):
    if isinstance(data, int):
        return vint(index << 3) + vint(data)
    return vint(index << 3 | 2) + vint(len(data)) + data


def mvt(points, tile, extent=4096, duplicate=False):
    keys = ('id', 'sequence_id', 'compass_angle', 'captured_at')
    values, features = [], []
    z, tx, ty = tile
    for identity, lon, lat in points:
        properties = (identity, 'fixture-sequence', 0, 1700000000000)
        tags = []
        for index, value in enumerate(properties):
            encoded = field(1, value.encode()) if isinstance(value, str) else field(5, value)
            tags.extend((index, len(values)))
            values.append(field(4, encoded))
        x = round(((lon+180)/360*(1 << z)-tx)*extent)
        y = round(((1-math.log(math.tan(math.pi/4+math.radians(lat)/2))/math.pi)/2*(1 << z)-ty)*extent)
        zig = lambda v: v*2 if v >= 0 else -v*2-1
        geometry = b''.join(map(vint, (9, zig(x), zig(y))))
        feature = field(2, b''.join(map(vint, tags))) + field(3, 1) + field(4, geometry)
        features.append(field(2, feature))
    if duplicate and features:
        features.append(features[0])
    # Deliberately put dictionaries AFTER features, as permitted by protobuf.
    layer = field(15, 2) + field(1, b'image') + b''.join(features)
    layer += b''.join(field(3, k.encode()) for k in keys) + b''.join(values) + field(5, extent)
    return field(3, layer)


def fixture(count=96):
    tile = next(tile_coordinates(BOUNDS))
    points, entries = [], {}
    for i in range(count):
        identity = str(BASE_ID+i)
        lon, lat = 24.94 + (i % 11)*.00001, 60.17 + (i % 13)*.00001
        points.append((int(identity), lon, lat))
        entry = dict(id=identity, geometry={'type': 'Point', 'coordinates': [lon, lat]},
                     computed_geometry={'type': 'Point', 'coordinates': [lon, lat]},
                     altitude=30., computed_altitude=32., computed_rotation=[math.pi/2, 0., 0.],
                     compass_angle=0., camera_type='perspective', camera_parameters=[.85, 0., 0.],
                     width=4000, height=3000, exif_orientation=1, sequence='fixture-sequence', captured_at=1700000000000)
        if i % 7 == 1:
            entry.pop('computed_rotation')
            entry['compass_angle'] = 90.
        if i % 7 == 2:
            entry.pop('altitude')
            entry.pop('computed_altitude')
        if i % 7 == 3:
            entry['exif_orientation'] = 6
        if i % 7 == 4:
            entry['camera_type'] = 'spherical'
            entry['camera_parameters'] = []
        if i % 7 == 5:
            entry['computed_rotation'] = [0., 0., 0.]
        entries[identity] = entry
    return tile, mvt(points, tile, duplicate=True), entries


class MemoryTransport:
    max_requests = 4

    def __init__(self, count=96, missing=False, reverse=True):
        self.tile, self.content, self.entries = fixture(count)
        self.pending = {}
        self.reverse, self.missing = reverse, missing
        self.requests = self.next_id = self.peak = 0
        self.closed = False

    def submit(self, request):
        from urllib.parse import urlparse, parse_qs
        assert len(self.pending) < self.max_requests
        self.next_id += 1
        self.requests += 1
        query = parse_qs(urlparse(request.url).query)
        if 'image_ids' in query:
            data = [self.entries[i] for i in query['image_ids'][0].split(',') if i in self.entries]
            content = json.dumps({'data': data[1:] if self.missing else data}).encode()
        else:
            parts=urlparse(request.url).path.rstrip('/').split('/')
            tile=tuple(map(int,parts[-3:]))
            if tile==self.tile:
                content=self.content
            else:
                points=[(int(i),*entry['geometry']['coordinates']) for i,entry in self.entries.items()]
                content=mvt(points,tile,duplicate=True)
        self.pending[self.next_id] = (self.next_id, content, None)
        self.peak = max(self.peak, len(self.pending))
        return self.next_id

    def poll(self):
        if self.pending:
            identity = (next(reversed(self.pending)) if self.reverse else next(iter(self.pending)))
            return self.pending.pop(identity)

    def close(self):
        self.closed = True
        self.pending.clear()
