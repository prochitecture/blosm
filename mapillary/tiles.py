"""Small, cooperative MVT point reader; integer identifiers remain exact.

Only the documented image point layer is consumed. Unknown protobuf fields are
skipped. No protobuf/geometry dependency or intermediate full GeoJSON is needed.
"""
from math import atan, ceil, degrees, floor, isfinite, log, pi, radians, sinh, tan
from array import array
import struct

ZOOM = 14
MAX_TILE_BYTES = 32 * 1024 * 1024
MAX_TABLE = 1048576
MAX_IMAGE_LAYERS = 256


def identifier(value):
    # Never accept even an apparently integral float: its low bits may be lost.
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise ValueError('Image identifier must be an integer or decimal string')
    value = str(value)
    if not value.isascii() or not value.isdecimal() or len(value) > 24:
        raise ValueError('Invalid image identifier')
    return value


def validate_bounds(bounds):
    w, s, e, n = bounds
    if not all(isfinite(v) for v in bounds) or not (-180 <= w < e <= 180 and -85.05112878 <= s < n <= 85.05112878):
        raise ValueError('Select finite west < east / south < north bounds within Web Mercator coverage')


def tile_coordinates(bounds):
    validate_bounds(bounds)
    w, s, e, n = bounds
    size = 1 << ZOOM
    def y(lat):
        return max(0, min(size - 1, floor((1 - log(tan(pi/4 + radians(lat)/2))/pi)/2 * size)))
    left, right = max(0, floor((w+180)/360*size)), min(size-1, floor((e+180)/360*size))
    top, bottom = y(n), y(s)
    if (right-left+1)*(bottom-top+1) > 4096:
        raise ValueError('Mapillary area exceeds 4096 coverage tiles; select a smaller area')
    for x in range(left, right+1):
        for row in range(top, bottom+1):
            yield ZOOM, x, row


def contains(bounds, lon, lat):
    return bounds[0] <= lon <= bounds[2] and bounds[1] <= lat <= bounds[3]


def varint(data, position):
    limit = len(data)
    if position >= limit:
        raise ValueError('Truncated MVT varint')
    byte = data[position]
    position += 1
    if byte < 128:
        return byte, position
    value = byte & 127
    for shift in range(7, 70, 7):
        if position >= limit:
            raise ValueError('Truncated MVT varint')
        byte = data[position]
        position += 1
        value |= (byte & 127) << shift
        if byte < 128:
            if value > 0xffffffffffffffff:
                raise ValueError('MVT integer exceeds uint64')
            return value, position
    raise ValueError('Invalid MVT varint')


def fields(data):
    position = 0
    while position < len(data):
        key, position = varint(data, position)
        number, wire = key >> 3, key & 7
        if not number:
            raise ValueError('Invalid protobuf field')
        if wire == 0:
            value, position = varint(data, position)
        elif wire in (1, 2, 5):
            length = 8 if wire == 1 else 4
            if wire == 2:
                length, position = varint(data, position)
            if position + length > len(data):
                raise ValueError('Truncated MVT field')
            value = data[position:position+length]
            position += length
        else:
            raise ValueError('Unsupported protobuf wire type')
        yield number, wire, value


def text(value):
    if len(value) > 65536:
        raise ValueError('MVT string exceeds limit')
    return bytes(value).decode('utf-8')


def scalar(data):
    for number, wire, value in fields(data):
        if number == 1 and wire == 2:
            return text(value)
        if number == 2 and wire == 5:
            return struct.unpack('<f', value)[0]
        if number == 3 and wire == 1:
            return struct.unpack('<d', value)[0]
        if wire == 0:
            if number == 4:
                return value if value < 1 << 63 else value - (1 << 64)
            if number == 5:
                return value
            if number == 6:
                return (value >> 1) ^ -(value & 1)
            if number == 7:
                return bool(value)
    return None


def packed(data):
    position = 0
    while position < len(data):
        value, position = varint(data, position)
        yield value


def feature(data, keys, values, extent, tile, bounds=None, ids_only=False, local_bounds=None):
    if len(data) > 65536:
        raise ValueError('MVT feature exceeds limit')
    tag_blocks, kind, geometry = [], None, None
    for number, wire, value in fields(data):
        if number == 2 and wire == 2:
            tag_blocks.append(value)
        elif number == 3 and wire == 0:
            kind = value
        elif number == 4 and wire == 2:
            geometry = value
    if kind != 1 or geometry is None:
        return None
    command, position = varint(geometry,0)
    encoded_x, position = varint(geometry,position)
    encoded_y, position = varint(geometry,position)
    # Mapillary documents one point per image. Coordinates are tile-local,
    # origin at top-left, y down; buffered coordinates may lie outside extent.
    if command != 9 or position != len(geometry):
        raise ValueError('Expected one MVT image point')
    x = (encoded_x >> 1) ^ -(encoded_x & 1)
    y = (encoded_y >> 1) ^ -(encoded_y & 1)
    if local_bounds is not None and not contains(local_bounds,x,y):
        return None
    z, tx, ty = tile
    size = 1 << z
    lon = (tx + x/extent)/size*360 - 180
    lat = degrees(atan(sinh(pi*(1-2*(ty+y/extent)/size))))
    # Dense city tiles mostly lie outside a small import extent. Decode point
    # geometry first; irrelevant feature attributes need never be materialized.
    if bounds is not None and not contains(bounds,lon,lat):
        return None
    props = {}
    for block in tag_blocks:
        tags = iter(packed(block))
        for key in tags:
            index = next(tags,None)
            if key >= len(keys) or index is None or index >= len(values):
                raise ValueError('Invalid MVT tags')
            name = keys[key]
            selected_value = values[index]
            if not ids_only or name == 'id':
                props[name] = selected_value
    if 'id' not in props:
        return None
    props['id'] = identifier(props['id'])
    if ids_only:
        return props['id']
    props['coverage_geometry'] = [lon, lat]
    return props


def image_points(content, tile, bounds=None):
    """Yield records or None checkpoints; retain bytes, not a decoded tile list."""
    reader = ImageTile(content,tile,bounds)
    try:
        yield from reader.image_points()
    finally:
        reader.close()


def image_ids(content, tile, bounds):
    """Count-path IDs/checkpoints, sharing the full reader's geometry and checks."""
    reader = ImageTile(content,tile,bounds)
    try:
        yield from reader.image_ids()
    finally:
        reader.close()


def _spans(data):
    """Yield payload offsets without allocating a memoryview for every value."""
    position, limit = 0, len(data)
    while position < limit:
        key, position = varint(data,position)
        number, wire = key >> 3, key & 7
        if not number:
            raise ValueError('Invalid protobuf field')
        if wire == 0:
            start = position
            _, position = varint(data,position)
        elif wire in (1,2,5):
            length = 8 if wire == 1 else 4
            if wire == 2:
                length, position = varint(data,position)
            start = position
            position += length
            if position > limit:
                raise ValueError('Truncated MVT field')
        else:
            raise ValueError('Unsupported protobuf wire type')
        yield number, wire, start, position


class _Values:
    """A compact dictionary index; only referenced values are decoded."""
    def __init__(self,data):
        self.data = data
        self.offsets = array('I')

    def __len__(self):
        return len(self.offsets)//2

    def __getitem__(self,index):
        start, end = self.offsets[index*2:index*2+2]
        return scalar(self.data[start:end])


class _Layer:
    def __init__(self,data):
        self.data = data
        self.keys = []
        self.values = _Values(data)
        self.extent = 4096
        self.selected = array('I')
        self.feature_end = 0
        self.local_bounds = None


def _local_bounds(bounds,extent,tile):
    """Conservative integer rejection; survivors still use geographic bounds.

    At documented zoom 14 and extent <= 2**24, expand floor/ceil by another
    complete tile unit to cover inverse-projection floating point roundoff.
    Other zooms or out-of-Mercator bounds use only the original exact check.
    Buffered feature coordinates can be negative or exceed the tile extent.
    """
    if bounds is None or tile[0] != ZOOM:
        return None
    w, s, e, n = bounds
    if not all(isfinite(value) for value in bounds) or not (
            -180 <= w <= e <= 180 and -85.05112878 <= s <= n <= 85.05112878):
        return None
    size = 1 << tile[0]
    if not (0 <= tile[1] < size and 0 <= tile[2] < size):
        return None
    left = ((w+180)/360*size-tile[1])*extent
    right = ((e+180)/360*size-tile[1])*extent
    def y(lat):
        return ((1-log(tan(pi/4+radians(lat)/2))/pi)/2*size-tile[2])*extent
    return floor(left)-1, floor(y(n))-1, ceil(right)+1, ceil(y(s))+1


class ImageTile:
    """Bounded per-tile dictionary/selection index shared by count and replay.

    Each value and selected feature uses two uint32 payload offsets, not a
    memoryview or decoded record. Whole-tile caps are MAX_TABLE values/keys and
    MAX_TABLE selected features (8 MiB each offset index), plus 256 image layers.
    The raw tile retains its existing 32 MiB cap. Scalar values are validated
    when selected features reference them; unused scalar contents are skipped.
    Close iterators before close() when cancelling an active traversal.
    """
    def __init__(self,content,tile,bounds=None):
        if len(content) > MAX_TILE_BYTES:
            raise ValueError('Mapillary coverage tile exceeds 32 MiB')
        self.data = memoryview(content)
        self.tile = tuple(tile)
        self.bounds = tuple(bounds) if bounds is not None else None
        self.layers = []
        self.prepared = self.index_ready = self.closed = False

    @property
    def index_bytes(self):
        return sum(layer.values.offsets.itemsize*len(layer.values.offsets) +
                   layer.selected.itemsize*len(layer.selected) for layer in self.layers)

    def close(self):
        self.closed = True
        self.layers.clear()
        self.data = None
        self.prepared = self.index_ready = False

    def _check(self):
        if self.closed:
            raise ValueError('Mapillary tile reader is closed')

    def _prepare(self):
        self._check()
        if self.prepared:
            return
        key_count = value_count = 0
        try:
            for number, wire, data in fields(self.data):
                if number != 3 or wire != 2:
                    continue
                name = None
                for index, (number, wire, value) in enumerate(fields(data)):
                    if number == 1 and wire == 2:
                        name = text(value)
                        break
                    if index % 64 == 0:
                        yield None
                if name != 'image':
                    continue
                if len(self.layers) >= MAX_IMAGE_LAYERS:
                    raise ValueError('MVT image layer count exceeds limit')
                layer = _Layer(data)
                for index, (number, wire, start, end) in enumerate(_spans(data)):
                    if number == 1 and wire == 2:
                        name = text(data[start:end])
                    elif number == 3 and wire == 2:
                        layer.keys.append(text(data[start:end]))
                        key_count += 1
                    elif number == 4 and wire == 2:
                        layer.values.offsets.extend((start,end))
                        value_count += 1
                    elif number == 2 and wire == 2:
                        layer.feature_end = end
                    elif number == 5 and wire == 0:
                        layer.extent, _ = varint(data,start)
                    if key_count > MAX_TABLE or value_count > MAX_TABLE:
                        raise ValueError('MVT dictionary exceeds limit')
                    if index % 64 == 0:
                        yield None
                if name != 'image':
                    continue
                if not 1 <= layer.extent <= 1 << 24:
                    raise ValueError('Invalid MVT extent')
                layer.local_bounds = _local_bounds(self.bounds,layer.extent,self.tile)
                self.layers.append(layer)
            self.prepared = True
        finally:
            if not self.prepared:
                self.layers.clear()

    def image_ids(self):
        if self.index_ready:
            yield from self._replay(True)
        else:
            yield from self._scan(True)

    def image_points(self):
        if self.index_ready:
            yield from self._replay(False)
        else:
            yield from self._scan(False)

    def _scan(self,ids_only):
        self._check()
        yield from self._prepare()
        selected_count = 0
        if ids_only:
            for layer in self.layers:
                del layer.selected[:]
        completed = False
        try:
            for layer in self.layers:
                # Preparation checked all framing already. Dictionaries often
                # follow all features, so do not traverse that large tail twice.
                for index, (number, wire, start, end) in enumerate(_spans(layer.data[:layer.feature_end])):
                    if number == 2 and wire == 2:
                        record = feature(layer.data[start:end],layer.keys,layer.values,
                                         layer.extent,self.tile,self.bounds,ids_only,layer.local_bounds)
                        if record is not None:
                            if ids_only:
                                if selected_count >= MAX_TABLE:
                                    raise ValueError('MVT selected feature index exceeds limit')
                                layer.selected.extend((start,end))
                                selected_count += 1
                            yield record
                    if index % 32 == 0:
                        yield None
            completed = True
            self.index_ready = ids_only
        finally:
            if ids_only and not completed:
                for layer in self.layers:
                    del layer.selected[:]

    def _replay(self,ids_only):
        self._check()
        for layer in self.layers:
            for index in range(0,len(layer.selected),2):
                start, end = layer.selected[index:index+2]
                record = feature(layer.data[start:end],layer.keys,layer.values,
                                 layer.extent,self.tile,self.bounds,ids_only,layer.local_bounds)
                if record is not None:
                    yield record
                if index % 64 == 0:
                    yield None
