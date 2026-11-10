"""Real GPX: independent XML point oracle, curve/mesh parity and shared origin."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import math
import unittest
import xml.etree.ElementTree as ET
import bpy
from release_support import arguments, enable, reset, real_fixture, run

args = arguments()
addon = enable(args)


class GpxChecks(unittest.TestCase):
    def check_import(self, kind, project, existing_origin):
        path = real_fixture(args, 'gpx')
        root = ET.parse(path).getroot()
        segments = [[(float(p.attrib['lat']), float(p.attrib['lon']),
                      float(p.findtext('{*}ele', '0'))) for p in segment.findall('{*}trkpt')]
                    for segment in root.findall('.//{*}trkseg')]
        points = [p for segment in segments for p in segment]
        self.assertTrue(points)
        settings = reset(args)
        settings.dataType = 'gpx'
        settings.gpxFilepath = str(path)
        settings.gpxImportType = kind
        settings.gpxProjectOnTerrain = project  # Without terrain: must retain elevations.
        if existing_origin:
            bpy.context.scene['lat'], bpy.context.scene['lon'] = points[0][:2]
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        obj = bpy.context.view_layer.objects.active
        if kind == 'curve':
            self.assertEqual([len(s.points) for s in obj.data.splines], [len(s) for s in segments])
            actual = [p.co.xyz for s in obj.data.splines for p in s.points]
            self.assertIsNotNone(obj.data.bevel_object)
            self.assertEqual(obj.data.bevel_mode, 'OBJECT', 'GPX width profile is not active')
        else:
            self.assertEqual(len(obj.data.edges), sum(max(0, len(s) - 1) for s in segments))
            self.assertEqual(len(obj.data.polygons), 0)
            actual = [v.co for v in obj.data.vertices]
            expected_edges, offset = set(), 0
            for segment in segments:
                expected_edges.update((offset + i, offset + i + 1) for i in range(len(segment) - 1))
                offset += len(segment)
            self.assertEqual({tuple(sorted(e.vertices)) for e in obj.data.edges}, expected_edges)
        self.assertEqual(len(actual), len(points))
        app = addon.blenderApp.app
        for source, target in zip(points, actual):
            self.assertTrue(all(math.isfinite(v) for v in target))
            self.assertAlmostEqual(source[2], target.z, delta=0.001)
            # Inverse geographic control points catch swapped axes and incorrect origin.
            lat, lon = app.projection.toGeographic(target.x, target.y)
            self.assertAlmostEqual(lat, source[0], delta=1e-6)
            self.assertAlmostEqual(lon, source[1], delta=1e-6)
        expected_origin = points[0][:2] if existing_origin else (
            (min(p[0] for p in points) + max(p[0] for p in points))/2,
            (min(p[1] for p in points) + max(p[1] for p in points))/2)
        for key, value in zip(('lat', 'lon'), expected_origin):
            self.assertAlmostEqual(bpy.context.scene[key], value)

    def check_terrain(self, kind, project):
        path = real_fixture(args, 'gpx')
        points = [(float(p.attrib['lat']), float(p.attrib['lon']), float(p.findtext('{*}ele', '0')))
                  for p in ET.parse(path).getroot().findall('.//{*}trkpt')]
        settings = reset(args)
        settings.dataType = 'gpx'
        settings.gpxFilepath = str(path)
        settings.gpxImportType = kind
        settings.gpxProjectOnTerrain = project
        bpy.context.scene['lat'], bpy.context.scene['lon'] = points[0][:2]
        bpy.ops.mesh.primitive_plane_add(size=100000, location=(0, 0, 10))
        terrain = bpy.context.object
        terrain['height_offset'] = 100.0
        settings.terrainObject = terrain.name
        bpy.context.view_layer.update()
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        obj = bpy.context.object
        if kind == 'curve' and project:
            shrinkwrap = [m for m in obj.modifiers if m.type == 'SHRINKWRAP']
            self.assertEqual(len(shrinkwrap), 1)
            self.assertEqual(shrinkwrap[0].target, terrain)
            self.assertAlmostEqual(shrinkwrap[0].offset, 0.5)
            return  # Curve deformation is evaluated after tessellation; width is checked separately.
        evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get()) if project else obj
        if kind == 'curve':
            actual = [p.co.z for s in evaluated.data.splines for p in s.points]
        else:
            actual = [v.co.z for v in evaluated.data.vertices]
        self.assertEqual(len(actual), len(points))
        expected = [10.5] * len(points) if project else [p[2] - 100 for p in points]
        for value, z in zip(actual, expected):
            self.assertAlmostEqual(value, z, delta=0.01)


for kind in ('curve', 'mesh'):
    for project in (False, True):
        for origin in (False, True):
            def check(self, kind=kind, project=project, origin=origin):
                self.check_import(kind, project, origin)
            setattr(GpxChecks, f'test_{kind}_project_{project}_existing_origin_{origin}', check)
        def terrain_check(self, kind=kind, project=project):
            self.check_terrain(kind, project)
        setattr(GpxChecks, f'test_{kind}_terrain_project_{project}', terrain_check)

run(args, GpxChecks)
