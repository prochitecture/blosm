"""OSM source identity versus geometric coincidence, in isolated Blender."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_street_progressive as support


def way(identifier, points, nodes, tags=None):
    from mathutils import Vector
    segments = [SimpleNamespace(v1=Vector(a),v2=Vector(b),id1=na,id2=nb)
                for a,b,na,nb in zip(points,points[1:],nodes,nodes[1:])]
    return SimpleNamespace(source_id=str(identifier),category='residential',
        element=SimpleNamespace(tags=dict(highway='residential',**(tags or {}))),segments=segments)


def prepare(ways, bounds=(-100.,100.,-100.,100.)):
    from blosm.action.generate_streets import StreetGenerator
    generator=object.__new__(StreetGenerator)
    generator.app=SimpleNamespace(minX=bounds[0],maxX=bounds[1],minY=bounds[2],maxY=bounds[3])
    generator.wayManager=SimpleNamespace(getAllWays=lambda:ways)
    generator.allWays=True;generator.leftHandTraffic=True
    generator.findSelfIntersections()
    generator.createWaySectionNetwork()
    return generator


def sections(generator):
    return [s for s in generator.sectionNetwork.iterAllForwardSegments() if s.category!='scene_border']


class SourceTopologyTests(unittest.TestCase):
    def test_projected_crossings_never_create_unmapped_osm_junctions(self):
        for tags in ({},{'bridge':'yes','layer':'1'}):
            with self.subTest(tags=tags):
                generator=prepare([way(1,[(-20.,0.),(20.,0.)],[1,2]),
                                   way(2,[(0.,-20.),(0.,20.)],[3,4],tags)])
                roads=sections(generator)
                self.assertEqual(len(roads),2)
                self.assertEqual([len(s.path) for s in roads],[2,2])
                self.assertNotIn((0.,0.),[tuple(p) for s in roads for p in (s.s,s.t)])
                audit=generator.wayManager.streetPreparationAudit
                self.assertEqual(audit['policy'],'OSM_NODE_IDENTITY')
                self.assertTrue(audit['events'])
                self.assertTrue(all(e['connectivity']=='SEPARATE_SOURCE_INCIDENCES' for e in audit['events']))

    def test_shared_osm_node_remains_a_junction_across_grade_tag_transitions(self):
        generator=prepare([way(1,[(-20.,0.),(0.,0.),(20.,0.)],[1,5,2]),
            way(2,[(0.,-20.),(0.,0.)],[3,5],{'bridge':'yes','layer':'1'}),
            way(3,[(0.,0.),(0.,20.)],[5,4])])
        roads=sections(generator)
        self.assertEqual(len(roads),4)
        center=[p for p in generator.sectionNetwork if tuple(p)==(0.,0.)]
        self.assertEqual(len(center),1)
        self.assertEqual(generator.sectionNetwork.order(center[0]),4)
        self.assertFalse(generator.wayManager.streetPreparationAudit['events'])

    def test_distinct_coincident_nodes_keep_two_junctions_without_displacing_geometry(self):
        from mathutils import Vector
        from blosm.way.source_topology import SourcePoint
        from blosm.way.waymap.waymap import WayMap
        from blosm.way.item import Intersection, Street
        ways=[way(i,[(0.,0.),point],[center,i]) for i,point,center in (
            (1,(-20.,0.),100),(2,(20.,0.),100),(3,(0.,20.),100),
            (4,(-20.,-20.),200),(5,(20.,-20.),200),(6,(0.,-20.),200))]
        generator=prepare(ways)
        roads=sections(generator)
        self.assertEqual(len(roads),6)
        centers=[p for p in generator.sectionNetwork if tuple(p)==(0.,0.)]
        self.assertEqual(len(centers),2)
        self.assertTrue(all(isinstance(p,SourcePoint) for p in centers))
        self.assertNotEqual(centers[0],centers[1])
        self.assertNotEqual(centers[0],Vector((0.,0.)))
        self.assertNotEqual(Vector((0.,0.)),centers[0])
        self.assertEqual({generator.sectionNetwork.order(p) for p in centers},{3})
        waymap=WayMap()
        for road in roads:
            waymap.addNode(Intersection(road.s));waymap.addNode(Intersection(road.t))
            waymap.addEdge(Street(road.s,road.t))
        self.assertEqual(waymap.number_of_nodes(),8)
        self.assertEqual([sum(map(len,waymap.getInOutEdges(p))) for p in centers],[3,3])
        for point in centers:
            self.assertEqual(tuple(point),(0.,0.))
            self.assertEqual(type(point + Vector((1.,2.))),Vector)
            self.assertEqual(type(point - centers[0]),Vector)
            self.assertEqual(type(point * 2.),Vector)
            self.assertEqual(type(-point),Vector)
            self.assertEqual(point.copy(),point)

    def test_clipped_source_edges_with_one_coincident_cut_stay_separate(self):
        generator=prepare([way(1,[(-10.,-10.),(10.,10.)],[1,2]),
                           way(2,[(-10.,10.),(10.,-10.)],[3,4])],(0.,20.,-20.,20.))
        roads=sections(generator)
        self.assertEqual(len(roads),2)
        cuts=[p for p in generator.sectionNetwork if tuple(p)==(0.,0.)
              and generator.sectionNetwork.borderlessOrder(p)]
        self.assertEqual(len(cuts),2)
        self.assertEqual([generator.sectionNetwork.order(p) for p in cuts],[1,1])
        self.assertEqual({tuple(s.s) for s in roads},{(0.,0.)})
        self.assertEqual({s.sourceSegments[0][3] for s in roads},{.5})

    def test_nonconnections_are_not_later_rejected_as_unverified_junctions(self):
        from blosm.way.street_network.adapter import NetworkSnapshot
        from blosm.way.street_network.connectivity import SurfaceConnectivityReview
        generator=prepare([way(1,[(-20.,0.),(20.,0.)],[1,2]),
                           way(2,[(0.,-20.),(0.,20.)],[3,4])])
        audit=generator.wayManager.streetPreparationAudit
        def snapshot():
            return NetworkSnapshot(preparation_audit=audit,streets={1:dict(id=1,status='SUPPORTED',
                items=[dict(kind='Section',category='residential',source_way_ids=('1',))])})
        current=snapshot();SurfaceConnectivityReview().apply(current,{'residential'})
        self.assertEqual(current.streets[1]['status'],'SUPPORTED')
        self.assertEqual(current.preparation_audit['surface_review_counts']['separate'],1)
        self.assertFalse(current.issues)
        legacy=snapshot();legacy.preparation_audit=dict(audit,policy='legacy')
        SurfaceConnectivityReview().apply(legacy,{'residential'})
        self.assertEqual(legacy.streets[1]['status'],'UNSUPPORTED')


if __name__=='__main__':
    support.enable()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SourceTopologyTests))
    if not result.wasSuccessful():raise RuntimeError('OSM source topology regression failed')
