"""Directed transition-chain ownership regressions in background Blender."""
from itertools import permutations, product
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import blosm
from mathutils import Vector
from blosm.lib.CompGeom.PolyLine import PolyLine
from blosm.way.item import Intersection, Section, Street
from blosm.way.way_network import NetSection
from blosm.way.waymap.waymap import WayMap
from blosm.way.transition_chains import TransitionStreetAssembler
from blosm.way.minor_cycles import MinorCycleReview
from blosm.way.manager import WayManager


def source(a, b, pattern='NN'):
    a,b=Vector(a).freeze(),Vector(b).freeze()
    section=Section(NetSection(a,b,'residential',{}),PolyLine([a,b]),None)
    section.setSectionAttributes(True,pattern,'',0,dict(laneWidth=2.5))
    street=Street(a,b);street.insertEnd(section)
    return street


def graph(streets):
    result=WayMap()
    for street in streets:
        for p in (street.src,street.dst):
            if p not in result:result.addNode(Intersection(p))
        result.addEdge(street)
    for p,node in result.iterNodes(Intersection):
        arriving,leaving=result.getInOutEdges(p)
        node.update(arriving,leaving)
    return result


def minor_manager(streets):
    waymap=graph(streets);minor={};major={}
    for point,node in waymap.iterNodes(Intersection):
        node.processIntersection()
        if len(node.leaveWays)>=3 and node.isMinorIntersection():node.transformToMinor()
        if node.isMinor:minor[point]=node
        elif node.order>1:major[point]=node
    style=SimpleNamespace(styleBlocks={})
    return SimpleNamespace(waymap=waymap,minorIntersections=minor,majorIntersections=major,
                          styleStore=SimpleNamespace(get=lambda key:style),getStyle=lambda street:'test')


def minor_loop():
    points=((0.,0.),(60.,0.),(30.,50.))
    through=[source(points[i],points[(i+1)%3]) for i in range(3)]
    sides=[source(p,(p[0]*1.5-15.,p[1]*1.5-10.)) for p in points]
    for street in sides:street.head.category='service'
    return through,sides


class StreetPreparationTests(unittest.TestCase):
    def test_bundle_terminals_exclude_internal_directed_continuation(self):
        from blosm.way.item.bundle import StreetGroup, orderHeadTail
        for reverse in (False, True):
            for order in permutations(range(3)):
                with self.subTest(reverse=reverse, order=order):
                    points = ((0., 10.), (50., 11.), (100., 10.))
                    if reverse:
                        points = points[::-1]
                    upper = [source(a, b) for a, b in zip(points, points[1:])]
                    lower = source((0., 0.), (100., 0.))
                    streets = upper + [lower]
                    branch = source(points[1], (50., 40.))
                    waymap = graph(streets + [branch])
                    for _, node in waymap.iterNodes(Intersection):
                        node.processIntersection()
                    links = [(s.pred, s.succ) for s in streets + [branch]]
                    controls = [tuple(s.head.centerline) for s in streets + [branch]]
                    group = StreetGroup([streets[i] for i in order])
                    head, tail = orderHeadTail(group)
                    first, last = upper[::-1] if reverse else upper
                    self.assertEqual([row['street'] for row in head], [first, lower])
                    self.assertEqual([row['street'] for row in tail], [last, lower])
                    self.assertNotIn(points[1], [tuple(row['firstVert']) for row in head + tail])
                    self.assertEqual(set(group), set(streets))
                    self.assertEqual(links, [(s.pred, s.succ) for s in streets + [branch]])
                    self.assertEqual(controls, [tuple(s.head.centerline) for s in streets + [branch]])

    def test_bundle_continuations_do_not_infer_joins_or_cut_cycles(self):
        from blosm.way.item.bundle import _bundleContinuationLocations
        examples = (
            [source((0., 0.), (10., 0.)), source((10., .001), (20., 0.))],
            [source((0., 0.), (10., 0.)), source((20., 0.), (10., 0.))],
            [source((0., 0.), (10., 0.)), source((10., 0.), (0., 1.))],
            [source((0., 0.), (10., 0.)), source((10., 0.), (20., 0.)),
             source((10., 0.), (20., 1.))],
        )
        for streets in examples:
            self.assertFalse(_bundleContinuationLocations(streets))
        points = ((0., 0.), (10., 0.), (10., 10.), (0., 10.))
        cycle = [source(points[i], points[(i + 1) % 4]) for i in range(4)]
        self.assertFalse(_bundleContinuationLocations(cycle))
        self.assertFalse(_bundleContinuationLocations(cycle + [source(points[0], (-10., 0.))]))
        # A broad U-turn is not a hairpin at any one join, but the established
        # forwardOrder frame reverses. Keep that join visible to preparation.
        returning = [source((0., 0.), (10., 9.)), source((10., 9.), (0., 18.))]
        self.assertFalse(_bundleContinuationLocations(returning))

    def test_clipping_point_contact_does_not_create_a_zero_length_section(self):
        from blosm.action.generate_streets import StreetGenerator
        def way(identifier, a, b):
            segment=SimpleNamespace(v1=Vector(a),v2=Vector(b),id1=identifier*2,id2=identifier*2+1)
            return SimpleNamespace(source_id=str(identifier),category='residential',
                element=SimpleNamespace(tags={'highway':'residential'}),segments=[segment])
        ways=(way(1,(-1.,-1.),(0.,0.)),way(2,(1.,1.),(2.,1.)),
              way(3,(3.,3.),(3.000001,3.)))
        generator=object.__new__(StreetGenerator)
        generator.app=SimpleNamespace(minX=0.,maxX=10.,minY=0.,maxY=10.)
        generator.wayManager=SimpleNamespace(getAllWays=lambda:ways)
        generator.allWays=True;generator.leftHandTraffic=False;generator.intersectingSegments={}
        generator.createWaySectionNetwork()
        retained=[s for s in generator.networkGraph.iterAllForwardSegments() if s.category!='scene_border']
        self.assertEqual({sid for s in retained for sid in s.sourceWayIds},{'2','3'})
        self.assertTrue(all(s.s!=s.t and s.length>0 for s in retained))

    def test_parallel_chain_keeps_its_outer_connectors_and_free_ends(self):
        from blosm.way.item.bundle import StreetGroup, mergePseudoMinors
        for count in (2, 3):
            for order in permutations(range(count)):
                for connected_start, connected_end in ((False,False),(False,True),(True,False),(True,True)):
                    with self.subTest(count=count,order=order,start=connected_start,end=connected_end):
                        points=[(float(i*40),0.) for i in range(count+1)]
                        through=[source(a,b) for a,b in zip(points,points[1:])]
                        branches=[source(p,(p[0],25.)) for p in points[1:-1]]
                        for connected,p in ((connected_start,points[0]),(connected_end,points[-1])):
                            if connected:
                                branches.extend(source(p,(p[0],y)) for y in (-25.,25.))
                        waymap=graph(through+branches)
                        nodes=[node for _,node in waymap.iterNodes(Intersection)]
                        for node in nodes:node.processIntersection()
                        manager=SimpleNamespace(streets={s.id:s for s in through+branches},
                            majorIntersections={n.location:n for n in nodes if n.startConnector},
                            minorIntersections={})
                        outer=(through[0].pred,through[-1].succ)
                        sections=[s.head for s in through]
                        controls=[tuple(s.centerline) for s in sections]
                        branch_links=[(s.pred,s.succ) for s in branches]
                        group=StreetGroup([through[i] for i in order])
                        membership={p:[group.id] for s in through for p in (s.src,s.dst)}
                        mergePseudoMinors(manager,group,membership)
                        self.assertEqual(len(group),1)
                        merged=group[0]
                        self.assertIs(merged.pred,outer[0])
                        self.assertIs(merged.succ,outer[1])
                        self.assertEqual((merged.src,merged.dst),(through[0].src,through[-1].dst))
                        items=list(merged.iterItems())
                        self.assertEqual(items[::2],sections)
                        self.assertEqual([tuple(s.centerline) for s in sections],controls)
                        self.assertTrue(all(item.street is merged for item in items))
                        self.assertTrue(all(j.isMinor for j in items[1::2]))
                        self.assertEqual(set(manager.streets),{merged.id}|{s.id for s in branches})
                        for connector in outer:
                            if connector:self.assertIs(connector.item,merged)
                        for branch,links in zip(branches,branch_links):
                            self.assertEqual((branch.pred,branch.succ),links)
                            for connector in links:
                                if connector:self.assertIs(connector.item,branch)
                        for junction in items[1::2]:
                            self.assertNotIn(junction.location,manager.majorIntersections)
                            self.assertIs(manager.minorIntersections[junction.location],junction)

    def test_bundle_terminal_copies_each_original_incidence_once(self):
        from blosm.way.item.bundle import endBundleIntersection
        for head in (True, False):
            with self.subTest(head=head):
                p=Vector((0.,0.)).freeze()
                members=[source(p,(x,20.)) if head else source((x,-20.),p) for x in (-2.,2.)]
                branch=source(p,(0.,-20. if head else 20.))
                j=Intersection(p);j.update([] if head else members,members+[branch] if head else [branch])
                j.processIntersection()
                bundle=SimpleNamespace(streetsHead=members,streetsTail=members,
                    headLocs=[p,p] if head else [s.src for s in members],
                    tailLocs=[s.dst for s in members] if head else [p,p],
                    pred=None if head else object(),succ=object() if head else None)
                for member in members:member.bundle=bundle
                before=[tuple(s.head.centerline) for s in members+[branch]]
                manager=SimpleNamespace(majorIntersections={p:j})
                endBundleIntersection(manager,bundle)
                junction=next(iter(manager.majorIntersections.values()))
                connectors=list(junction)
                self.assertEqual(len(connectors),2)
                outgoing=next(c for c in connectors if c.item is branch)
                self.assertIs(branch.pred,outgoing)
                self.assertTrue(outgoing.leaving)
                self.assertEqual(before,[tuple(s.head.centerline) for s in members+[branch]])
                for member in members:
                    self.assertIs(member.pred if head else member.succ,connectors[0])

    def test_bundle_terminal_preserves_both_incident_ends_of_one_street(self):
        from blosm.way.item.bundle import endBundleIntersection
        from blosm.way.item.connectors import IntConnector
        ends=[Vector(p).freeze() for p in ((-4.,0.),(4.,0.))]
        members=[source(p,(p.x,20.)) for p in ends]
        branch=source(*ends);major={}
        for i,p in enumerate(ends):
            j=Intersection(p);major[p]=j
            for item,leaving in ((members[i],True),(branch,i==0)):
                c=IntConnector(j);c.item=item;c.leaving=leaving;j.insertConnector(c)
                if leaving:item.pred=c
                else:item.succ=c
        bundle=SimpleNamespace(streetsHead=members,streetsTail=members,headLocs=ends,
            tailLocs=[s.dst for s in members],pred=None,succ=object())
        for member in members:member.bundle=bundle
        manager=SimpleNamespace(majorIntersections=major)
        endBundleIntersection(manager,bundle)
        junction=next(iter(major.values()));refs=[c for c in junction if c.item is branch]
        self.assertEqual(len(refs),2)
        self.assertEqual({c.leaving for c in refs},{False,True})
        self.assertIsNot(branch.pred,branch.succ)
        self.assertIs(branch.pred.intersection,junction)
        self.assertIs(branch.succ.intersection,junction)

    def test_bundle_terminal_order_is_independent_of_physical_ring_starts(self):
        """Retain the exterior walk across a connecting Street at either end."""
        from blosm.way.item.bundle import endBundleIntersection
        from blosm.way.item.connectors import IntConnector
        for head, left_start, right_start in product((False, True), range(4), range(4)):
            with self.subTest(head=head, left=left_start, right=right_start):
                left, right = (Vector(p).freeze() for p in ((-4., 0.), (4., 0.)))
                members = [source(p, (p.x, -20.)) if head else source((p.x, -20.), p)
                           for p in (left, right)]
                inner = source(left, right)
                a, b = (source(left, (-20., y)) for y in (12., -12.))
                c, d = (source(right, (20., y)) for y in (12., -12.))
                streets = members + [inner, a, b, c, d]
                controls = [tuple(s.head.centerline) for s in streets]
                rings = (((a, True), (b, True), (members[0], head), (inner, True)),
                         ((c, True), (inner, False), (members[1], head), (d, True)))
                major = {}
                for p, ring, shift in zip((left, right), rings, (left_start, right_start)):
                    junction = Intersection(p); major[p] = junction
                    for item, leaving in ring:
                        connector = IntConnector(junction)
                        connector.item = item; connector.leaving = leaving
                        junction.insertConnector(connector)
                        if leaving: item.pred = connector
                        else: item.succ = connector
                    for _ in range(shift): junction.startConnector = junction.startConnector.succ
                original = [(j, tuple(j)) for j in major.values()]
                ordered = members[::-1] if head else members
                bundle = SimpleNamespace(streetsHead=ordered, streetsTail=ordered,
                    headLocs=[s.src for s in ordered], tailLocs=[s.dst for s in ordered],
                    pred=None if head else object(), succ=object() if head else None)
                for member in members: member.bundle = bundle
                endBundleIntersection(SimpleNamespace(majorIntersections=major), bundle)
                logical, = major.values()
                result = list(logical)
                self.assertIs(result[0].item, bundle)
                self.assertEqual([(r.item, r.leaving) for r in result[1:]],
                    [(d, True), (c, True), (inner, False), (inner, True), (a, True), (b, True)])
                self.assertIsNot(inner.pred, inner.succ)
                for original_junction, ring in original:
                    self.assertEqual(tuple(original_junction), ring)
                    self.assertTrue(all(r.succ.pred is r and r.pred.succ is r for r in ring))
                self.assertEqual(controls, [tuple(s.head.centerline) for s in streets])
                for r in result[1:]:
                    self.assertIs(r.item.pred if r.leaving else r.item.succ, r)

    def test_bundle_tail_order_uses_tail_locations_after_a_bend(self):
        """A quarter-turn must not sort the tail by the remote head points."""
        from itertools import product
        from math import cos, sin, pi
        from blosm.way.item.bundle import orderHeadTail
        paths=(((0.,8.),(30.,8.),(38.,16.),(38.,60.)),
               ((-5.,0.),(35.,0.),(46.,11.),(46.,60.)))
        cases=0
        for angle,reverse0,reverse1,swap in product(range(8),(False,True),(False,True),(False,True)):
            theta=angle*pi/4
            streets=[]
            for points,backward in zip(paths,(reverse0,reverse1)):
                points=tuple(Vector((cos(theta)*x-sin(theta)*y,sin(theta)*x+cos(theta)*y)).freeze()
                             for x,y in (reversed(points) if backward else points))
                section=Section(NetSection(points[0],points[-1],'residential',{}),PolyLine(points),None)
                street=Street(points[0],points[-1]);street.insertEnd(section);streets.append(street)
            if swap:streets.reverse()
            with self.subTest(angle=angle,reverse0=reverse0,reverse1=reverse1,swap=swap):
                head,tail=orderHeadTail(streets)
                self.assertEqual([r['street'] for r in head],[r['street'] for r in tail])
                self.assertEqual({r['street'] for r in head},set(streets))
                for h,t in zip(head,tail):
                    self.assertEqual(h['firstVert'],t['lastVert'])
                    self.assertEqual(h['lastVert'],t['firstVert'])
            cases+=1
        self.assertEqual(cases,64)

    def test_bundle_ordering_does_not_hide_a_real_order_change(self):
        from blosm.way.item.bundle import orderHeadTail
        streets=[]
        for coordinates in (((0.,8.),(30.,8.),(46.,60.)),((-5.,0.),(35.,0.),(38.,60.))):
            points=tuple(Vector(p).freeze() for p in coordinates)
            section=Section(NetSection(points[0],points[-1],'residential',{}),PolyLine(points),None)
            street=Street(points[0],points[-1]);street.insertEnd(section);streets.append(street)
        head,tail=orderHeadTail(streets)
        self.assertEqual([r['street'] for r in head],streets)
        self.assertEqual([r['street'] for r in tail],list(reversed(streets)))

    def test_minor_cycle_retains_all_original_roads_and_complete_junctions(self):
        for order in permutations(range(3)):
            with self.subTest(order=order):
                through,sides=minor_loop()
                original=through+sides
                manager=minor_manager([through[i] for i in order]+sides)
                self.assertEqual(len(manager.minorIntersections),3)
                nodes=list(manager.minorIntersections.values())
                rings={n:tuple(n.iterConnectors()) for n in nodes}
                sections={s:s.head for s in original}
                coordinates={s:tuple(tuple(p) for p in s.head.centerline) for s in original}
                prepared=list(WayManager.iterStreetsFromWaymap(manager))
                self.assertEqual(set(prepared),set(original))
                self.assertEqual(len(prepared),6)
                self.assertFalse(manager.minorIntersections)
                self.assertEqual(set(manager.majorIntersections.values()),set(nodes))
                self.assertEqual(len(manager.minorCycleAudit),1)
                self.assertEqual(set(manager.minorCycleAudit[0]['street_ids']),{s.id for s in through})
                for street in prepared:
                    self.assertEqual(list(street.iterItems()),[sections[street]])
                    self.assertIs(street.head.street,street)
                    self.assertEqual(tuple(tuple(p) for p in street.head.centerline),coordinates[street])
                for node in nodes:
                    self.assertEqual(tuple(node),rings[node])
                    for conn in node:
                        self.assertIs(conn.item.pred if conn.leaving else conn.item.succ,conn)
                # No mutation or duplicate generation on another read.
                self.assertEqual(set(WayManager.iterStreetsFromWaymap(manager)),set(original))
                self.assertEqual(len(manager.minorCycleAudit),1)

    def test_open_minor_chain_still_joins_and_preserves_side_connectors(self):
        through=[source((30.*i,0.),(30.*(i+1),0.)) for i in range(3)]
        sides=[source((30.*i,0.),(30.*i,25.)) for i in (1,2)]
        for street in sides:street.head.category='service'
        manager=minor_manager([through[1],*sides,through[2],through[0]])
        original=[s.head for s in through];nodes=list(manager.minorIntersections.values())
        self.assertFalse(MinorCycleReview(manager).find())
        prepared=list(WayManager.iterStreetsFromWaymap(manager))
        self.assertEqual(len(prepared),3)
        combined=next(s for s in prepared if s not in sides)
        self.assertEqual([i for i in combined.iterItems() if isinstance(i,Section)],original)
        self.assertEqual(len(list(combined.iterItems())),5)
        self.assertTrue(all(n.isMinor and n.street is combined for n in nodes))
        self.assertTrue(all(s.pred.item is s for s in sides))

    def test_return_to_major_junction_is_retained_instead_of_silently_dropped(self):
        # Distinct paths share their two junction locations; the top path has
        # an extra control so this is a real enclosed loop, not a doubled edge.
        forward=source((0.,0.),(60.,0.));back=source((60.,0.),(0.,0.))
        back.head.polyline=PolyLine([back.src,Vector((30.,40.)).freeze(),back.dst])
        at_minor=source((60.,0.),(90.,0.));at_minor.head.category='service'
        at_major=[source((0.,0.),(-30.,0.)),source((0.,0.),(0.,-30.))]
        manager=minor_manager([forward,back,at_minor,*at_major])
        self.assertEqual(len(manager.minorIntersections),1)
        self.assertFalse(MinorCycleReview(manager).find())
        original={s.head for s in (forward,back,at_minor,*at_major)}
        prepared=list(WayManager.iterStreetsFromWaymap(manager))
        self.assertEqual(len(prepared),4)
        retained=[i for s in prepared for i in s.iterItems() if isinstance(i,Section)]
        self.assertEqual(len(retained),len(original));self.assertEqual(set(retained),original)
        loop=next(s for s in prepared if s.src==s.dst)
        self.assertEqual(len(list(loop.iterItems())),3)
        self.assertIs(loop.pred.item,loop);self.assertIs(loop.succ.item,loop)

    def test_minor_cycle_validation_does_not_partially_reclassify_on_failure(self):
        through,sides=minor_loop();manager=minor_manager(through+sides)
        nodes=list(manager.minorIntersections.values())
        nodes[-1].startConnector=None
        with self.assertRaisesRegex(ValueError,'complete original connector ring'):
            MinorCycleReview(manager).apply()
        self.assertEqual(len(manager.minorIntersections),3)
        self.assertFalse(manager.majorIntersections)
        self.assertTrue(all(n.isMinor for n in nodes))

    def test_minor_cycle_validation_is_bounded_for_a_damaged_connector_ring(self):
        through,sides=minor_loop();manager=minor_manager(through+sides)
        node=next(iter(manager.minorIntersections.values()))
        node.startConnector.succ.succ=node.startConnector.succ
        with self.assertRaisesRegex(ValueError,'complete original connector ring'):
            MinorCycleReview(manager).apply()
        self.assertTrue(all(n.isMinor for n in manager.minorIntersections.values()))

    def test_long_chain_is_owned_once_in_every_input_order(self):
        for order in permutations(range(4)):
            with self.subTest(order=order):
                streets=[source((40.*i,0.),(40.*(i+1),0.),'N'*(2+i%2)) for i in range(4)]
                original=[s.head for s in streets]
                controls=[tuple(tuple(p) for p in s.centerline) for s in original]
                waymap=graph([streets[i] for i in order]);symmetric=[];sided=[]
                combined,=TransitionStreetAssembler(waymap,False).assemble(symmetric,sided)
                self.assertEqual(waymap.number_of_edges(),1)
                self.assertEqual(waymap.number_of_nodes(),2)
                self.assertEqual(len(symmetric),3);self.assertFalse(sided)
                items=list(combined.iterItems())
                self.assertEqual(items[::2],original)
                self.assertEqual(len(items),7)
                self.assertTrue(all(item.street is combined for item in items))
                self.assertEqual([tuple(tuple(p) for p in s.centerline) for s in original],controls)
                for before,transition,after in zip(original,symmetric,original[1:]):
                    self.assertIs(transition.pred,before);self.assertIs(transition.succ,after)
                    self.assertIs(before.succ,transition);self.assertIs(after.pred,transition)
                # Calling the stage again cannot re-collect stale old owners.
                self.assertEqual(TransitionStreetAssembler(waymap,False).assemble([],[]),())

    def test_repeated_side_lane_joins_keep_the_pixel_correction_setting(self):
        streets=[source((40.*i,0.),(40.*(i+1),0.),pattern)
                 for i,pattern in enumerate(('NN','NN','NNR','NNRR'))]
        original=[s.head for s in streets]
        waymap=graph([streets[i] for i in (1,3,0,2)]);symmetric=[];sided=[]
        combined,=TransitionStreetAssembler(waymap,False).assemble(symmetric,sided)
        self.assertEqual(len(symmetric),1);self.assertEqual(len(sided),2)
        self.assertEqual([s.width for s in original],[5.,5.,7.5,10.])
        self.assertEqual([s.doPixelCorrection for s in sided],[False,False])
        self.assertEqual([s for s in combined.iterItems() if isinstance(s,Section)],original)

    def test_outer_junction_remaps_each_chain_to_its_own_street(self):
        left=[source((-80.,0.),(-40.,0.)),source((-40.,0.),(0.,0.))]
        right=[source((0.,0.),(40.,0.)),source((40.,0.),(80.,0.))]
        top=[source((0.,0.),(0.,40.)),source((0.,40.),(0.,80.))]
        waymap=graph(left+right+top);original={s.head:s for s in left+right+top}
        result=TransitionStreetAssembler(waymap,False).assemble([],[])
        self.assertEqual(len(result),3)
        node=waymap.getNode(Vector((0.,0.)))['object']
        self.assertEqual(len({w.street for w in node.leaveWays}),3)
        for way in node.leaveWays:
            self.assertIn(way.street,result)
            self.assertIs(way.section.street,way.street)
            self.assertIs(way.street.head if way.leaving else way.street.tail,way.section)
        node.processIntersection()
        self.assertEqual({c.item for c in node},set(result))
        self.assertEqual({s for street in result for s in street.iterItems() if isinstance(s,Section)},set(original))

    def test_cycle_is_retained_without_mutating_its_streets_or_nodes(self):
        a,b,c=(0.,0.),(40.,0.),(20.,35.)
        streets=[source(a,b),source(b,c),source(c,a)];waymap=graph(streets)
        original=[s.head for s in streets]
        symmetric=[];sided=[]
        self.assertEqual(TransitionStreetAssembler(waymap,False).assemble(symmetric,sided),())
        self.assertFalse(symmetric or sided)
        self.assertEqual(waymap.number_of_nodes(),3);self.assertEqual(waymap.number_of_edges(),3)
        for street,section in zip(streets,original):
            self.assertEqual(list(street.iterItems()),[section])
            self.assertIs(section.street,street)
            self.assertIsNone(section.pred);self.assertIsNone(section.succ)

    def test_incompatible_directions_keep_distinct_incidents(self):
        streets=[source((0.,0.),(40.,0.)),source((0.,0.),(-40.,0.))]
        waymap=graph(streets)
        self.assertFalse(TransitionStreetAssembler(waymap,False).assemble([],[]))
        self.assertEqual(waymap.number_of_edges(),2)
        self.assertEqual({w.street for w in waymap.getNode(Vector((0.,0.)))['object'].leaveWays},set(streets))


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(StreetPreparationTests))
    if not result.wasSuccessful():raise RuntimeError('Street preparation regression failed')
