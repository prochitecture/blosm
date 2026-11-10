"""Selected real OSM extracts: grouping parity, layers, repeat import and persistence."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import unittest
import random
import bpy
from release_support import (arguments, enable, catalog, configure_osm, meshes,
                             geometry, assert_geometry, run)

args = arguments()
addon = enable(args)


class OsmChecks(unittest.TestCase):
    def compare_grouping(self, entry, mode):
        snapshots, counts, areas = [], [], []
        for single in (False, True):
            settings = configure_osm(args, entry, mode, single)
            roads = entry['id'] in ('frankfurter_allee', 'kuzmice_railways')
            if roads:
                settings.buildings = False
                settings.highways = entry['id'] == 'frankfurter_allee'
                settings.railways = True
            self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
            self.assertIn('lat', bpy.context.scene)
            if roads:
                curves = [o for o in bpy.context.scene.objects if o.type == 'CURVE']
                self.assertTrue(curves, 'No road/rail curves imported')
                snapshots.append(sorted(tuple(round(c, 3) for c in o.matrix_world @ p.co.xyz)
                                        for o in curves for s in o.data.splines for p in s.points))
                counts.append(len(curves))
            else:
                self.assertTrue(meshes(), 'No buildings imported')
                snapshots.append(geometry(meshes()))
                counts.append(len(meshes()))
                areas.append(sum(p.area for o in meshes() for p in o.data.polygons))
                if mode == '2D':
                    z = [p[2] for f in snapshots[-1] for p in f]
                    self.assertLess(max(z) - min(z), 0.001, 'Footprints must be planar')
                else:
                    z = [p[2] for f in snapshots[-1] for p in f]
                    self.assertGreater(max(z) - min(z), 1.0)
        if roads:
            self.assertEqual(snapshots[0], snapshots[1])
        elif mode == '2D':
            # Courtyard face decomposition can differ when joining footprints.
            # Compare surface area and all boundary/control vertices instead.
            self.assertAlmostEqual(areas[0], areas[1], delta=max(0.05, areas[0]*1e-6))
            from mathutils.kdtree import KDTree
            for source, target in (snapshots, list(reversed(snapshots))):
                points = [p for face in target for p in face]
                tree = KDTree(len(points))
                for i, p in enumerate(points):
                    tree.insert(p, i)
                tree.balance()
                self.assertTrue(all(tree.find(p)[2] < 0.01 for face in source for p in face))
        else:
            assert_geometry(self, *snapshots)
        self.assertLessEqual(counts[1], counts[0])
        if entry['id'] == 'shared_edges':
            self.assertEqual(counts, [3, 1])

    def test_reimport_and_save_reopen(self):
        entry = next(e for e in catalog() if e['id'] == 'rectangle')
        configure_osm(args, entry)
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        original = geometry(meshes())
        origin = (bpy.context.scene['lat'], bpy.context.scene['lon'])
        random.seed(1729)
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertEqual((bpy.context.scene['lat'], bpy.context.scene['lon']), origin)
        assert_geometry(self, sorted(original * 2), geometry(meshes()))
        path = args.artifacts / 'repeated_import.blend'
        bpy.ops.wm.save_as_mainfile(filepath=str(path))
        bpy.ops.wm.open_mainfile(filepath=str(path), load_ui=False)
        assert_geometry(self, sorted(original * 2), geometry(meshes()))

    def test_disabled_layers(self):
        entry = next(e for e in catalog() if e['id'] == 'karl_marx_allee')
        settings = configure_osm(args, entry)
        settings.buildings = False
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertFalse(list(bpy.context.scene.objects))


for entry in catalog(args.extended):
    for mode in ('2D', '3Dsimple'):
        def check(self, entry=entry, mode=mode):
            self.compare_grouping(entry, mode)
        setattr(OsmChecks, f'test_{entry["id"]}_{mode}', check)

run(args, OsmChecks)
