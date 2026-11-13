"""Bounded, Blender-independent validation of photo metadata and image headers.

Graph orientation describes the original capture, not necessarily the downloaded
thumbnail. Only the downloaded file's own EXIF orientation controls its pixels.
These checks bound decoded size before Blender's native decoder is invoked;
the native decoder still validates image compression data.
"""
import json
import struct
from urllib.parse import urlsplit
import zlib

MAX_METADATA_BYTES = 256 * 1024
MAX_PHOTO_BYTES = 32 * 1024 * 1024
MAX_PIXELS = 8 * 1024 * 1024
MAX_SIDE = 8192
MAX_HEADER_BYTES = 1024 * 1024
MAX_SEGMENTS = 4096
_HOSTS = ('fbcdn.net', 'fbsbx.com', 'mapillary.com')
_BAD_METADATA = 'Mapillary did not return a valid photo download.'
_BAD_IMAGE = 'Mapillary returned an invalid or unsupported photo.'


def photo_url(content, image_id):
    """Return the selected image's HTTPS download URL and bounded attribution.

    Do not include provider payloads, signed URLs or JSON errors in diagnostics.
    The caller should keep the returned URL transient and omit API credentials
    when downloading the CDN resource.
    """
    try:
        if len(content) > MAX_METADATA_BYTES:
            raise ValueError()
        document = json.loads(content)
        if not isinstance(document, dict) or 'error' in document:
            raise ValueError()
        rows = document.get('data')
        if not isinstance(rows, list) or len(rows) != 1:
            raise ValueError()
        row = rows[0]
        if not isinstance(row, dict) or not isinstance(row.get('id'), (str, int)) or isinstance(row.get('id'), bool):
            raise ValueError()
        if str(row['id']) != str(image_id):
            raise ValueError()
        url = row.get('thumb_2048_url')
        if not isinstance(url, str) or not 1 <= len(url) <= 16384 or any(ord(c) < 33 or ord(c) == 127 for c in url):
            raise ValueError()
        parsed = urlsplit(url)
        host = (parsed.hostname or '').lower()
        if (parsed.scheme != 'https' or parsed.username is not None or parsed.password is not None
                or parsed.port not in (None, 443) or parsed.fragment
                or not any(host == allowed or host.endswith('.' + allowed) for allowed in _HOSTS)):
            raise ValueError()
        creator = row.get('creator')
        if isinstance(creator, dict):
            creator = creator.get('username') or creator.get('id') or ''
        if not isinstance(creator, (str, int)) or isinstance(creator, bool):
            creator = ''
        creator = ''.join(c for c in str(creator) if c.isprintable())[:256]
        return url, creator
    except (ValueError, TypeError, UnicodeError, RecursionError, OverflowError):
        raise ValueError(_BAD_METADATA) from None


def _dimensions(width, height):
    if not width or not height or max(width,height) > MAX_SIDE or width * height > MAX_PIXELS:
        raise ValueError(_BAD_IMAGE)
    return width, height


def _orientation(tiff):
    """Read only IFD0's inline SHORT Orientation, without following TIFF links."""
    if len(tiff) < 8 or bytes(tiff[:2]) not in (b'II', b'MM'):
        raise ValueError(_BAD_IMAGE)
    order = '<' if bytes(tiff[:2]) == b'II' else '>'
    if struct.unpack_from(order + 'H', tiff, 2)[0] != 42:
        raise ValueError(_BAD_IMAGE)
    offset = struct.unpack_from(order + 'I', tiff, 4)[0]
    if offset < 8 or offset + 2 > len(tiff):
        raise ValueError(_BAD_IMAGE)
    count = struct.unpack_from(order + 'H', tiff, offset)[0]
    if count > MAX_SEGMENTS or offset + 2 + count * 12 + 4 > len(tiff):
        raise ValueError(_BAD_IMAGE)
    orientation = None
    for index in range(count):
        entry = offset + 2 + index * 12
        tag, kind, amount = struct.unpack_from(order + 'HHI', tiff, entry)
        if tag != 0x112:
            continue
        if orientation is not None or kind != 3 or amount != 1:
            raise ValueError(_BAD_IMAGE)
        orientation = struct.unpack_from(order + 'H', tiff, entry + 8)[0]
        if not 1 <= orientation <= 8:
            raise ValueError(_BAD_IMAGE)
    return orientation or 1


def _jpeg(data):
    if bytes(data[-2:]) != b'\xff\xd9':
        raise ValueError(_BAD_IMAGE)
    position, segments = 2, 0
    dimensions = None
    orientation, had_exif = 1, False
    while position < len(data):
        segments += 1
        if segments > MAX_SEGMENTS or position > MAX_HEADER_BYTES:
            raise ValueError(_BAD_IMAGE)
        if data[position] != 255:
            raise ValueError(_BAD_IMAGE)
        while position < len(data) and data[position] == 255:
            position += 1
        if position >= len(data) or position > MAX_HEADER_BYTES:
            raise ValueError(_BAD_IMAGE)
        marker = data[position]
        position += 1
        if marker in (0, 0xd8, 0xd9, 1) or 0xd0 <= marker <= 0xd7:
            raise ValueError(_BAD_IMAGE)
        if position + 2 > len(data):
            raise ValueError(_BAD_IMAGE)
        length = struct.unpack_from('>H', data, position)[0]
        end = position + length
        if length < 2 or end > len(data) or end > MAX_HEADER_BYTES:
            raise ValueError(_BAD_IMAGE)
        payload = data[position+2:end]
        if marker in (0xc0, 0xc1, 0xc2):
            if dimensions is not None or len(payload) < 6:
                raise ValueError(_BAD_IMAGE)
            depth, height, width, channels = struct.unpack_from('>BHHB', payload)
            if depth != 8 or channels not in (1, 3, 4) or len(payload) != 6 + 3 * channels:
                raise ValueError(_BAD_IMAGE)
            dimensions = _dimensions(width, height)
        elif 0xc0 <= marker <= 0xcf and marker not in (0xc4, 0xc8, 0xcc):
            raise ValueError(_BAD_IMAGE)
        elif marker == 0xe1 and bytes(payload[:6]) == b'Exif\0\0':
            if had_exif:
                raise ValueError(_BAD_IMAGE)
            orientation, had_exif = _orientation(payload[6:]), True
        elif marker == 0xda:
            if dimensions is None or len(payload) < 6:
                raise ValueError(_BAD_IMAGE)
            channels = payload[0]
            if channels not in (1, 2, 3, 4) or len(payload) != 4 + 2 * channels:
                raise ValueError(_BAD_IMAGE)
            return *dimensions, orientation
        position = end
    raise ValueError(_BAD_IMAGE)


def _png(data):
    position, chunks = 8, 0
    dimensions = None
    had_pixels, had_exif, orientation = False, False, 1
    while position < len(data):
        chunks += 1
        if chunks > MAX_SEGMENTS or position + 12 > len(data):
            raise ValueError(_BAD_IMAGE)
        length = struct.unpack_from('>I', data, position)[0]
        kind = bytes(data[position+4:position+8])
        end = position + 12 + length
        if end > len(data) or (kind != b'IDAT' and length > 65536):
            raise ValueError(_BAD_IMAGE)
        payload = data[position+8:end-4]
        crc = zlib.crc32(data[position+4:end-4]) & 0xffffffff
        if crc != struct.unpack_from('>I', data, end-4)[0]:
            raise ValueError(_BAD_IMAGE)
        if kind == b'IHDR':
            if dimensions is not None or position != 8 or length != 13:
                raise ValueError(_BAD_IMAGE)
            width, height, depth, color, compression, filtering, interlace = struct.unpack('>IIBBBBB', payload)
            depths = {0:(1,2,4,8,16), 2:(8,16), 3:(1,2,4,8), 4:(8,16), 6:(8,16)}
            if depth not in depths.get(color, ()) or compression or filtering or interlace not in (0,1):
                raise ValueError(_BAD_IMAGE)
            dimensions = _dimensions(width, height)
        elif dimensions is None:
            raise ValueError(_BAD_IMAGE)
        elif kind == b'IDAT':
            had_pixels = True
        elif kind == b'eXIf':
            if had_exif:
                raise ValueError(_BAD_IMAGE)
            orientation, had_exif = _orientation(payload), True
        elif kind == b'IEND':
            if length or not had_pixels or end != len(data):
                raise ValueError(_BAD_IMAGE)
            return *dimensions, orientation
        position = end
    raise ValueError(_BAD_IMAGE)


def image_info(content):
    """Return encoded width, height and file EXIF orientation before decoding.

    JPEG entropy data and PNG compressed pixels are left to Blender. Header
    framing and dimensions are checked here to bound the native allocation.
    """
    try:
        if not 16 <= len(content) <= MAX_PHOTO_BYTES:
            raise ValueError()
        data = memoryview(content)
        if bytes(data[:2]) == b'\xff\xd8':
            return _jpeg(data)
        if bytes(data[:8]) == b'\x89PNG\r\n\x1a\n':
            return _png(data)
    except (ValueError, TypeError, IndexError, struct.error, OverflowError):
        pass
    raise ValueError(_BAD_IMAGE) from None
