"""Offline imagery stitching, cache decoding, URL expansion and 404 fallback."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
from io import BytesIO
import unittest
from urllib.error import HTTPError
from unittest.mock import patch
import bpy
import numpy as np
from release_support import arguments, enable, reset, REPO, PACKAGE, run

args = arguments()
enable(args)
module = importlib.import_module(PACKAGE + '.overlay')
Overlay = module.Overlay


class OverlayChecks(unittest.TestCase):
    def setUp(self):
        reset(args)
        module.app.maxNumTiles = 16
        module.app.area = None
        module.app.terrain = None
        module.app.saveOverlayToFile = False
        module.app.setOverlayMaterial = False
        self.overlay = Overlay('https://fixture.invalid/{z}/{x}/{y}.png', 2, PACKAGE)
        self.overlay.tileWidth = self.overlay.tileHeight = 4
        self.overlay.overlayDir = str(args.artifacts / self._testMethodName)
        self.overlay.prepareImport(-10, -10, 10, 10)

    def test_quadrant_orientation_and_packing(self):
        overlay = self.overlay
        self.assertEqual((overlay.numTilesX, overlay.numTilesY), (2, 2))
        colors = {(0, 0): (1, 0, 0, 1), (1, 0): (0, 1, 0, 1),
                  (0, 1): (0, 0, 1, 1), (1, 1): (1, 1, 0, 1)}
        def tile(zoom, x, y):
            return np.tile(colors[x-overlay.l, y-overlay.t], 16)
        with patch.object(overlay, 'getTileData', side_effect=tile):
            while overlay.importNextTile():
                pass
        self.assertTrue(overlay.finalizeImport())
        image = bpy.data.images.get(Overlay.blenderImageName)
        self.assertEqual(tuple(image.size), (8, 8))
        self.assertIsNotNone(image.packed_file)
        pixels = np.array(image.pixels[:]).reshape(8, 8, 4)
        for (x, y), color in colors.items():
            np.testing.assert_allclose(pixels[(1-y)*4:(2-y)*4, x*4:(x+1)*4],
                                       np.broadcast_to(color, (4, 4, 4)), atol=1/255)

    def test_download_then_cache(self):
        image = bpy.data.images.new('fixture', width=4, height=4)
        image.pixels[:] = [1, 0, 0, 1] * 16
        image.filepath_raw = str(args.artifacts / 'red.png')
        image.file_format = 'PNG'
        image.save()
        payload = Path(image.filepath_raw).read_bytes()
        bpy.data.images.remove(image)
        before = len(bpy.data.images)
        with patch.object(module.request, 'urlopen', return_value=BytesIO(payload)) as request:
            first = self.overlay.getTileData(2, 1, 1)
            request.assert_called_once()
        # Outer offline guard now rejects network: the second request must hit disk.
        second = self.overlay.getTileData(2, 1, 1)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_allclose(first.reshape(-1, 4), np.tile((1, 0, 0, 1), (16, 1)), atol=1/255)
        self.assertEqual(len(bpy.data.images), before, 'Temporary tile image leaked')

    def test_404_reduces_zoom_then_stops(self):
        with patch.object(module.request, 'urlopen', side_effect=HTTPError('fixture', 404, '', {}, None)):
            self.assertIs(self.overlay.getTileData(2, 1, 1), True)
            self.assertEqual(self.overlay.zoom, 1)
            self.assertIs(self.overlay.getTileData(1, 0, 0), False)
        self.assertIsNone(self.overlay.imageData)

    def test_url_and_geographic_controls(self):
        overlay = Overlay('https://[a,b].fixture.invalid', 19, PACKAGE)
        self.assertEqual(overlay.getTileUrl(4, 3, 2), 'https://a.fixture.invalid/4/3/2.png')
        overlay.tileCounter = 1
        self.assertEqual(overlay.getTileUrl(4, 3, 2), 'https://b.fixture.invalid/4/3/2.png')
        y, x = Overlay.toSphericalMercator(0, 0, True)
        self.assertAlmostEqual(x, module.halfEquator)
        self.assertAlmostEqual(y, module.halfEquator)
        self.assertEqual(Overlay.toTileCoord(x, 2), 2)
        self.assertAlmostEqual(Overlay.fromTileCoord(2, 2), x)


run(args, OverlayChecks)
