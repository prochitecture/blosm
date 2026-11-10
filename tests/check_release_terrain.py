"""Big-endian HGT elevation, voids, reduction and two-axis seam topology."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from collections import Counter
import gzip
import importlib
import unittest
import numpy as np
from release_support import arguments, enable, run

args = arguments()
addon = enable(args)
Terrain = importlib.import_module(addon.__name__ + '.terrain').Terrain


class IdentityProjection:
    def fromGeographic(self, lat, lon):
        return lon, lat, 0


class TerrainChecks(unittest.TestCase):
    def check_grid(self, reduction, primitive, void):
        # A miniature grid exercises the same byte addressing and stitching as SRTM.
        app = addon.blenderApp.BlenderApp()
        app.terrainSize = 4
        app.terrainReductionRatio = str(reduction)
        app.terrainPrimitiveType = primitive
        app.projection = IdentityProjection()
        app.latIntervals = ((1, 2), (0, 1))
        app.lonIntervals = ((0, 1), (1, 2))
        app.terrainDir = str(args.artifacts / 'hgt')
        Path(app.terrainDir).mkdir(exist_ok=True)
        for lat in (0, 1):
            for lon in (0, 1):
                grid = np.array([[100 + (lat*4 + 4-y)*10 + lon*4 + x
                                  for x in range(5)] for y in range(5)], dtype='>i2')
                if void and (lat, lon) == (0, 0):
                    grid[4, 0] = -32768
                with gzip.open(Path(app.terrainDir) / Terrain.getHgtFileName(lat, lon), 'wb') as stream:
                    stream.write(grid.tobytes())
        vertices, faces = [], []
        minimum = app.buildTerrain(vertices, faces, None)
        width = 8 // reduction + 1
        self.assertEqual(len(vertices), width**2)
        self.assertEqual(len(faces), (width-1)**2 * (1 if primitive == 'quad' else 2))
        self.assertEqual(len({tuple(v[:2]) for v in vertices}), len(vertices), 'Duplicated seam vertices')
        for x, y, z in vertices:
            self.assertEqual(z, 0 if void and x == y == 0 else 100 + y*40 + x*4)
        self.assertEqual(minimum, 0 if void else 100)
        edges = Counter(tuple(sorted((face[i], face[(i+1) % len(face)])))
                        for face in faces for i in range(len(face)))
        self.assertEqual(sum(count == 1 for count in edges.values()), 4*(width-1))
        self.assertTrue(all(count in (1, 2) for count in edges.values()), 'Nonmanifold seam')
        self.assertTrue(all(0 <= i < len(vertices) for face in faces for i in face))
        v2, f2 = [], []
        self.assertEqual(app.buildTerrain(v2, f2, 42), 42)
        self.assertEqual(vertices, v2)
        self.assertEqual(faces, f2)


for reduction in (1, 2):
    for primitive in ('quad', 'triangle'):
        for void in (False, True):
            def check(self, reduction=reduction, primitive=primitive, void=void):
                self.check_grid(reduction, primitive, void)
            setattr(TerrainChecks, f'test_{primitive}_reduction_{reduction}_void_{void}', check)

run(args, TerrainChecks)
