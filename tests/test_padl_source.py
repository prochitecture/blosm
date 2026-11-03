"""Multipolygon association, stable identity and source-edge preservation."""
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
PREFIX = '_padl_source_tests'
for suffix, path in (('', ROOT), ('.building', ROOT / 'building'),
                     ('.building.padl_integration', ROOT / 'building/padl_integration'),
                     ('.padl', ROOT / 'padl'), ('.parse', ROOT / 'parse')):
    package = ModuleType(PREFIX + suffix)
    package.__path__ = [str(path)]
    sys.modules[package.__name__] = package
osm_module = ModuleType(PREFIX + '.parse.osm')
osm_module.Osm = SimpleNamespace(outer=4, inner=5, way=2, relation=3)
sys.modules[osm_module.__name__] = osm_module
source = importlib.import_module(PREFIX + '.building.padl_integration.source')
mapping = importlib.import_module(PREFIX + '.building.padl_integration.mapping')


def ring(x=0):
    return ((x, 0), (x+10, 0), (x+10, 10), (x, 10))


def line(way, points, role=4):
    return SimpleNamespace(wayId=str(way), parts=None, role=role, points=points)


def relation(lines):
    return SimpleNamespace(ls=lines, tags={'building': 'apartments'},
                           getLinestringData=lambda item, osm: item.points)


class SourceTests(unittest.TestCase):
    def ways(self, items):
        report = source.ImportReport()
        ways = {f'osm:way:{key}': SimpleNamespace(tags=tags, getData=lambda osm, p=points: p,
                                                **({'o': parent} if parent else {}))
                for key, points, tags, parent in items}
        return source.OsmFootprintSource(None, mapping.OsmFaceAttributeMapper(), report).collect(ways, {}), report

    def test_parts_shared_edges_classification_and_grouping(self):
        result, report = self.ways([
            (10, ring(), {'building': 'apartments', 'height': '30'}, None),
            (11, ((0,0),(5,0),(5,10),(0,10)), {'building:part': 'yes', 'height': '12'}, None),
            (12, ((5,0),(10,0),(10,10),(5,10)), {'building:part': 'office', 'height': '20'}, None)])
        self.assertFalse(report.issues)
        self.assertEqual([f.face_type for f in result], [5, 3, 3])
        self.assertEqual([f.parent_id for f in result], [0, result[0].footprint_id, result[0].footprint_id])
        self.assertEqual([f.attributes['building_use'] for f in result], ['apartments', 'apartments', 'office'])
        self.assertEqual([f.attributes['max_elevation'] for f in result], [30, 12, 20])
        self.assertEqual(source.group_buildings(result, False), (tuple(result),))

    def test_dual_outline_and_explicit_overhanging_part(self):
        result, report = self.ways([
            (10, ring(), {'building': 'yes', 'building:part': 'yes'}, None),
            (11, ring(5), {'building:part': 'yes'}, ('10', 2))])
        self.assertFalse(report.issues)
        self.assertEqual([f.face_type for f in result], [6, 3])
        self.assertEqual(result[1].parent_id, result[0].footprint_id)
        standalone, _ = self.ways([(10, ring(), {'building': 'yes', 'building:part': 'yes'}, None)])
        self.assertEqual(standalone[0].face_type, 6)

    def test_orphan_ambiguous_and_rejected_part_do_not_create_outline_volume(self):
        result, report = self.ways([
            (10, ring(), {'building': 'yes'}, None),
            (11, ring(), {'building:part': 'yes', 'roof:shape': 'dome'}, None),
            (12, ring(30), {'building:part': 'yes'}, None)])
        self.assertFalse(result)
        self.assertEqual({r['code'] for r in report.issues}, {'unsupported_roof', 'part_ownership'})
        result, report = self.ways([
            (10, ring(), {'building': 'yes'}, None),
            (11, ring(), {'building': 'yes'}, None),
            (12, ring(), {'building:part': 'yes'}, None)])
        self.assertEqual(len(result), 2)
        self.assertEqual(report.issues[0]['code'], 'part_ownership')

    def test_invalid_part_geometry_does_not_enable_outline_fallback(self):
        result, report = self.ways([
            (10, ring(), {'building':'yes'}, None),
            (11, ((1,1),(2,2),(3,3)), {'building:part':'yes'}, ('10',2))])
        self.assertFalse(result)
        self.assertEqual(report.issues[0]['code'], 'invalid_geometry')

    def test_multipolygon_part_owns_its_courtyard(self):
        report = source.ImportReport()
        part = relation([line(20, ((1,1),(9,1),(9,9),(1,9))),
                         line(21, ((3,3),(3,5),(5,5),(5,3)),role=5)])
        part.tags = {'building:part':'yes'}
        way = SimpleNamespace(tags={'building':'apartments'},getData=lambda osm:ring())
        result = source.OsmFootprintSource(None,mapping.OsmFaceAttributeMapper(),report).collect(
            {'osm:way:10':way},{'osm:relation:7':part})
        building = next(f for f in result if f.face_type == 5)
        child = next(f for f in result if f.face_type == 3)
        self.assertEqual(child.parent_id,building.footprint_id)
        self.assertEqual(len(child.holes),1)
        self.assertEqual(child.attributes['building_use'],'apartments')
        self.assertFalse(report.issues)

    def test_part_can_share_a_building_courtyard_boundary(self):
        report = source.ImportReport()
        parent = relation([line(10,ring()),line(30,((4,4),(4,6),(6,6),(6,4)),role=5)])
        way = SimpleNamespace(tags={'building:part':'yes'},
            getData=lambda osm:((0,0),(4,0),(4,10),(0,10)))
        result = source.OsmFootprintSource(None,mapping.OsmFaceAttributeMapper(),report).collect(
            {'osm:way:50':way},{'osm:relation:7':parent})
        self.assertEqual(sorted(f.face_type for f in result),[3,5])
        self.assertFalse(report.issues)

    def collect(self, lines, ways=None):
        report = source.ImportReport()
        result = source.OsmFootprintSource(None, mapping.OsmFaceAttributeMapper(), report).collect(
            ways or {}, {'osm:relation:7': relation(lines)})
        return result, report

    def test_multiple_outers_and_owned_courtyard(self):
        outer, other = line(10, ring()), line(20, ring(20))
        hole = line(30, ((2, 2), (2, 4), (4, 4), (4, 2)), role=5)
        result, report = self.collect([outer, other, hole])
        self.assertEqual(len(result), 2)
        self.assertEqual([len(f.holes) for f in result], [1, 0])
        self.assertEqual([f.identity for f in result], ['osm:relation:7:outer:10', 'osm:relation:7:outer:20'])
        reordered, _ = self.collect([other, hole, outer])
        self.assertEqual(result, reordered)
        self.assertFalse(report.issues)

    def test_tagged_members_are_not_rendered_twice(self):
        result, report = self.collect([line(10, ring())], {'osm:way:10': object()})
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].identity, 'osm:relation:7')

    def test_invalid_hole_does_not_fall_back_to_member_way(self):
        result, report = self.collect([line(10, ring()), line(30, ring(30), role=5)],
                                      {'osm:way:10': object()})
        self.assertFalse(result)
        self.assertEqual(report.issues[0]['code'], 'invalid_geometry')

    def test_only_closing_coordinate_is_normalized(self):
        points = ((0,0), (5,0), (10,0), (10,10), (0,10), (0,0))
        self.assertEqual(source.clean_ring(points), points[:-1])


if __name__ == '__main__':
    unittest.main()
