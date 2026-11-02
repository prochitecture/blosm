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
osm_module.Osm = SimpleNamespace(outer=4, inner=5)
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
