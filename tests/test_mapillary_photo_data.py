"""Pure-Python regressions for bounded photo metadata/header handling."""
import json
from pathlib import Path
import struct
import sys
import unittest
from unittest.mock import patch
import zlib

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mapillary import photo_data as photos


def tiff(orientation=1, order='<'):
    return (b'II' if order == '<' else b'MM') + struct.pack(order+'HIHHHIH',42,8,1,0x112,3,1,orientation) + b'\0'*6


def segment(marker, payload):
    return bytes((255,marker))+struct.pack('>H',len(payload)+2)+payload


def jpeg(width=640,height=480,orientation=None,order='<'):
    exif = segment(0xe1,b'Exif\0\0'+tiff(orientation,order)) if orientation is not None else b''
    frame = struct.pack('>BHHB',8,height,width,3)+b'\x01\x11\x00\x02\x11\x01\x03\x11\x01'
    scan = b'\x03\x01\x00\x02\x11\x03\x11\x00\x3f\x00'
    return b'\xff\xd8'+exif+segment(0xc0,frame)+segment(0xda,scan)+b'\x00\xff\xd9'


def chunk(kind,payload):
    return struct.pack('>I',len(payload))+kind+payload+struct.pack('>I',zlib.crc32(kind+payload)&0xffffffff)


def png(width=3,height=2,extra=b''):
    header=struct.pack('>IIBBBBB',width,height,8,2,0,0,0)
    return b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',header)+extra+chunk(b'IDAT',zlib.compress((b'\0'+b'\x10\x20\x30'*width)*height))+chunk(b'IEND',b'')


class PhotoMetadataChecks(unittest.TestCase):
    def metadata(self, **kw):
        row={'id':'9007199254740993','thumb_2048_url':'https://scontent.example.fbcdn.net/photo.jpg'}
        row.update(kw)
        return json.dumps({'data':[row]}).encode()

    def test_url_and_attribution_preserve_exact_ids(self):
        url,creator=photos.photo_url(self.metadata(creator={'username':'Alice\nExample'}),'9007199254740993')
        self.assertEqual(url,'https://scontent.example.fbcdn.net/photo.jpg')
        self.assertEqual(creator,'AliceExample')
        self.assertEqual(photos.photo_url(self.metadata(creator={'id':123}),'9007199254740993')[1],'123')
        self.assertEqual(len(photos.photo_url(self.metadata(creator={'username':'a'*1000}),'9007199254740993')[1]),256)

    def test_rejects_untrusted_downloads_without_echoing_payload(self):
        for url in ('http://scontent.example.fbcdn.net/a','https://fbcdn.net.evil.test/a',
                    'https://evil.test/a','https://user:SECRET@fbcdn.net/a',
                    'https://fbcdn.net:8443/a','https://fbcdn.net/a#fragment','https://fbcdn.net/a\n'):
            with self.assertRaises(ValueError) as error:
                photos.photo_url(self.metadata(thumb_2048_url=url),'9007199254740993')
            self.assertEqual(str(error.exception),photos._BAD_METADATA)

    def test_invalid_metadata_is_bounded_and_does_not_echo_secrets(self):
        for content in (b'not JSON SECRET',b'{"error":{"message":"SECRET"}}',b'{}',
                        b'{"data":[]}',self.metadata(id=9007199254740992.),self.metadata(id='1'),
                        self.metadata(thumb_2048_url=None)):
            with self.assertRaisesRegex(ValueError,'valid photo download') as error:
                photos.photo_url(content,'9007199254740993')
            self.assertNotIn('SECRET',str(error.exception))
        with patch.object(photos,'MAX_METADATA_BYTES',3):
            with self.assertRaises(ValueError):photos.photo_url(self.metadata(),'9007199254740993')


class PhotoHeaderChecks(unittest.TestCase):
    def test_jpeg_dimensions_and_every_embedded_exif_orientation(self):
        self.assertEqual(photos.image_info(jpeg()),(640,480,1))
        for orientation in range(1,9):
            for order in ('<','>'):
                self.assertEqual(photos.image_info(jpeg(orientation=orientation,order=order)),(640,480,orientation))

    def test_png_dimensions_crc_and_embedded_exif(self):
        self.assertEqual(photos.image_info(png()),(3,2,1))
        self.assertEqual(photos.image_info(png(extra=chunk(b'eXIf',tiff(6)))),(3,2,6))
        damaged=bytearray(png()); damaged[20]^=1
        with self.assertRaises(ValueError):photos.image_info(damaged)

    def test_pixel_and_download_caps_apply_before_native_decode(self):
        for content in (jpeg(65535,65535),jpeg(0,12),png(0,1),jpeg(photos.MAX_SIDE+1,1),jpeg(1,photos.MAX_SIDE+1)):
            with self.assertRaises(ValueError):photos.image_info(content)
        with patch.object(photos,'MAX_PIXELS',5):
            with self.assertRaises(ValueError):photos.image_info(png())
        with patch.object(photos,'MAX_PHOTO_BYTES',15):
            with self.assertRaises(ValueError):photos.image_info(jpeg())

    def test_malformed_and_truncated_containers(self):
        for content in (b'not image SECRET',jpeg()[:-2],jpeg()[:15],png()[:-1],png()+b'trailing',
                        b'\xff\xd8'+segment(0xe1,b'Exif\0\0broken')+jpeg()[2:]):
            with self.assertRaises(ValueError) as error:photos.image_info(content)
            self.assertEqual(str(error.exception),photos._BAD_IMAGE)

    def test_invalid_or_ambiguous_exif_and_header_limits(self):
        for orientation in (0,9):
            with self.assertRaises(ValueError):photos.image_info(jpeg(orientation=orientation))
        duplicate=b'\xff\xd8'+segment(0xe1,b'Exif\0\0'+tiff(1))+jpeg(orientation=1)[2:]
        with self.assertRaises(ValueError):photos.image_info(duplicate)
        with patch.object(photos,'MAX_HEADER_BYTES',16):
            with self.assertRaises(ValueError):photos.image_info(jpeg())
        with patch.object(photos,'MAX_SEGMENTS',1):
            with self.assertRaises(ValueError):photos.image_info(jpeg())


if __name__=='__main__':
    unittest.main()
