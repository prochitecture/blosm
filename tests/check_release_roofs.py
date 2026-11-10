"""Small generated OSM controls supplement the selected real extracts."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import unittest
import xml.etree.ElementTree as ET
import bpy
from release_support import arguments, enable, reset, geometry, meshes, run

args = arguments()
enable(args)


class RoofChecks(unittest.TestCase):
    def check_roof(self, shape, single):
        root = ET.Element('osm', version='0.6')
        for i, (lat, lon) in enumerate(((52, 13), (52, 13.0003), (52.0002, 13.0003), (52.0002, 13)), 1):
            ET.SubElement(root, 'node', id=str(i), lat=str(lat), lon=str(lon))
        way = ET.SubElement(root, 'way', id='100')
        for ref in (1, 2, 3, 4, 1):
            ET.SubElement(way, 'nd', ref=str(ref))
        for key, value in {'building': 'yes', 'height': '12', 'roof:height': '3',
                           'roof:shape': shape, 'building:min_level': '0'}.items():
            ET.SubElement(way, 'tag', k=key, v=value)
        path = args.artifacts / f'{shape}.osm'
        ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)
        settings = reset(args)
        settings.dataType = 'osm'
        settings.osmSource = 'file'
        settings.osmFilepath = str(path)
        settings.mode = '3Dsimple'
        settings.singleObject = single
        settings.loadMissingMembers = False
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertEqual(len(meshes()), 1)
        faces = geometry(meshes())
        z = [p[2] for face in faces for p in face]
        self.assertAlmostEqual(min(z), 0, delta=0.001)
        self.assertAlmostEqual(max(z), 12, delta=0.001)
        self.assertGreaterEqual(len(faces), 5)


for shape in ('flat', 'gabled', 'pyramidal', 'skillion', 'hipped', 'dome', 'onion',
              'round', 'half-hipped', 'gambrel', 'saltbox', 'mansard', 'unknown-shape'):
    for single in (False, True):
        def check(self, shape=shape, single=single):
            self.check_roof(shape, single)
        setattr(RoofChecks, f'test_{shape.replace("-", "_")}_single_{single}', check)

run(args, RoofChecks)
