"""Small, cooperative MVT point reader; integer identifiers remain exact.

Only the documented image point layer is consumed. Unknown protobuf fields are
skipped. No protobuf/geometry dependency or intermediate full GeoJSON is needed.
"""
from math import atan, degrees, floor, isfinite, log, pi, radians, sinh, tan
import struct

ZOOM = 14
MAX_TILE_BYTES = 32 * 1024 * 1024
MAX_TABLE = 1048576


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


def feature(data, keys, values, extent, tile, bounds=None, ids_only=False):
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
    coords = list(packed(geometry))
    # Mapillary documents one point per image. Coordinates are tile-local,
    # origin at top-left, y down; buffered coordinates may lie outside extent.
    if len(coords) != 3 or coords[0] != 9:
        raise ValueError('Expected one MVT image point')
    x, y = ((v >> 1) ^ -(v & 1) for v in coords[1:])
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
            if not ids_only or name == 'id':
                props[name] = values[index]
    if 'id' not in props:
        return None
    props['id'] = identifier(props['id'])
    if ids_only:
        return props['id']
    props['coverage_geometry'] = [lon, lat]
    return props


def image_points(content, tile, bounds=None):
    """Yield records or None checkpoints; retain bytes, not a decoded tile list."""
    yield from _image_entries(content,tile,bounds,False)


def image_ids(content, tile, bounds):
    """Count-path IDs/checkpoints, sharing the full reader's geometry and checks."""
    yield from _image_entries(content,tile,bounds,True)


def _image_entries(content, tile, bounds, ids_only):
    if len(content) > MAX_TILE_BYTES:
        raise ValueError('Mapillary coverage tile exceeds 32 MiB')
    data = memoryview(content)
    for number, wire, layer in fields(data):
        if number != 3 or wire != 2:
            continue
        name, extent, keys, values = None, 4096, [], []
        # Other coverage layers (e.g. sequences) can also be very large. Read
        # their name without allocating dictionaries we will never consume.
        for index, (number, wire, value) in enumerate(fields(layer)):
            if number == 1 and wire == 2:
                name = text(value)
                break
            if index % 64 == 0:
                yield None
        if name != 'image':
            continue
        # MVT permits feature fields before dictionaries. First pass gathers
        # only dictionaries; the second pass streams feature records.
        for index, (number, wire, value) in enumerate(fields(layer)):
            if number == 1 and wire == 2:
                name = text(value)
            elif number == 3 and wire == 2:
                keys.append(text(value))
            elif number == 4 and wire == 2:
                values.append(scalar(value))
            elif number == 5 and wire == 0:
                extent = value
            if len(keys) > MAX_TABLE or len(values) > MAX_TABLE:
                raise ValueError('MVT dictionary exceeds limit')
            if index % 64 == 0:
                yield None
        if name != 'image':
            continue
        if not 1 <= extent <= 1 << 24:
            raise ValueError('Invalid MVT extent')
        for index, (number, wire, value) in enumerate(fields(layer)):
            if number == 2 and wire == 2:
                record = feature(value, keys, values, extent, tile, bounds, ids_only)
                if record is not None:
                    yield record
            if index % 32 == 0:
                yield None
