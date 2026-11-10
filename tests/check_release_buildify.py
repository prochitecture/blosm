"""Use the supplied Buildify binary on real OSM footprints, without executing scripts."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import unittest
import xml.etree.ElementTree as ET
import bpy
from release_support import (arguments, enable, catalog, configure_osm, real_fixture,
                             meshes, geometry, assert_geometry, run)

args = arguments()
enable(args)


class BuildifyChecks(unittest.TestCase):
    def check_buildify(self, fixture_id):
        entry = next(e for e in catalog() if e['id'] == fixture_id)
        configure_osm(args, entry, '2D', True)
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        footprints = geometry(meshes())
        settings = configure_osm(args, entry, '2D', True)
        settings.gnBlendFile2d = str(real_fixture(args, 'gn2d'))
        settings.gnSetup2d = 'building'
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        objects = meshes()
        if fixture_id == 'rectangle':
            assert_geometry(self, footprints, geometry(objects))
        else:
            # The current GN renderer supports closed building ways; ordinary 2D
            # also triangulates courtyard relations, so face counts need not match.
            root = ET.parse(settings.osmFilepath).getroot()
            ways = [w for w in root.findall('way') if w.find("tag[@k='building']") is not None]
            closed = [w for w in ways if w.find('nd').get('ref') == w.findall('nd')[-1].get('ref')]
            # Every closed building way must survive. Relations may contribute
            # additional outer polygons; courtyard tessellation is not identical.
            self.assertGreaterEqual(sum(len(o.data.polygons) for o in objects), len(closed))
            self.assertLessEqual(sum(len(o.data.polygons) for o in objects), len(footprints))
        self.assertEqual(len(objects), 1)
        obj = objects[0]
        attribute = obj.data.attributes.get('building:levels')
        self.assertIsNotNone(attribute)
        self.assertEqual(attribute.domain, 'FACE')
        self.assertEqual(attribute.data_type, 'INT')
        self.assertEqual(len(attribute.data), len(obj.data.polygons))
        self.assertTrue(all(v.value > 0 for v in attribute.data))
        modifiers = [m for m in obj.modifiers if m.type == 'NODES']
        self.assertEqual(len(modifiers), 1)
        self.assertTrue(modifiers[0].node_group.name.startswith('building'))
        bpy.context.view_layer.update()
        # Buildify returns instances. Realize them in a temporary copy to measure
        # geometry; object.dimensions alone does not include instances.
        group = modifiers[0].node_group.copy()
        modifiers[0].node_group = group
        output = next(n for n in group.nodes if n.type == 'GROUP_OUTPUT' and n.is_active_output)
        socket = next(s for s in output.inputs if s.type == 'GEOMETRY')
        original = socket.links[0].from_socket
        realize = group.nodes.new('GeometryNodeRealizeInstances')
        group.links.new(original, realize.inputs['Geometry'])
        group.links.new(realize.outputs['Geometry'], socket)
        obj.update_tag()
        bpy.context.view_layer.update()
        evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
        self.assertGreater(evaluated.dimensions.z, 0.1, 'Buildify produced no building height')
        path = args.artifacts / f'{fixture_id}_buildify.blend'
        bpy.ops.wm.save_as_mainfile(filepath=str(path))
        bpy.ops.wm.open_mainfile(filepath=str(path), load_ui=False)
        self.assertIsNotNone(meshes()[0].data.attributes.get('building:levels'))
        self.assertIsNotNone(meshes()[0].modifiers[0].node_group)

    def test_rectangle(self):
        self.check_buildify('rectangle')

    def test_karl_marx_allee(self):
        self.check_buildify('karl_marx_allee')


run(args, BuildifyChecks)
