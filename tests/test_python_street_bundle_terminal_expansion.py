"""Whole short-member ownership acquires exact opposite terminal sites."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import _python_street_package
from tested_python_streets.model import Section, End, Intersection, Bundle, BundleMember, Traversal, Network, Key
from tested_python_streets.bundle_terminal_expansion import missing_terminals
from tested_python_streets.compound_preparation import DeferredCompoundPlanner, CompoundExpansionRequired


def fixture(reverse=False):
    points = ((0., 0.), (0., 10.), (10., 0.), (10., 10.))
    def section(sid, a, b):
        return Section(sid, sid, ((*a, 0.), (*b, 0.)), width=3.5)
    sections = [section(1, points[0], points[1]), section(2, points[2], points[3]),
                section(3, points[0], points[2]), section(4, points[1], points[3])]
    for i, point in enumerate(points):
        sections.append(section(5+i, point, (point[0]+(-20. if i<2 else 20.), point[1])))
    sites = (Intersection(10, (End(1, 0), End(3, 0), End(5, 0))),
             Intersection(11, (End(1, 1), End(4, 0), End(6, 0))),
             Intersection(12, (End(2, 0), End(3, 1), End(7, 0))),
             Intersection(13, (End(2, 1), End(4, 1), End(8, 0))))
    if reverse:
        sections = [replace(s, points=s.points[::-1]) for s in sections]
        sites = tuple(replace(j, ends=tuple(End(e.section_id, 1-e.end) for e in j.ends)) for j in sites)
    bundle = Bundle(7, tuple(BundleMember((Traversal(sid, not reverse),)) for sid in (1, 2)))
    network = Network('whole_member_terminal_closure', tuple(sections), sites, (bundle,))
    return SimpleNamespace(network=network), bundle, {j.key: j for j in sites}


class BundleTerminalExpansionTests(unittest.TestCase):
    def test_exact_source_terminal_expansion_precedes_unchanged_whole_pair_claim(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                context, bundle, layouts = fixture(reverse)
                before = context.network.to_dict()
                members = {Key('intersection', i) for i in (10, 11, 13)}
                self.assertEqual(missing_terminals(context, bundle, {1, 4}, members, layouts),
                                 (Key('intersection', 12),))
                planner = DeferredCompoundPlanner(context, context.network.intersections, layouts)
                with self.assertRaises(CompoundExpansionRequired) as caught:
                    planner.check_bundle_terminals({1, 4}, (), members)
                self.assertEqual(caught.exception.neighbors, (Key('intersection', 12),))
                _, claims = planner.check_bundle_terminals({1, 2, 3, 4}, (), set(layouts))
                self.assertEqual([c.bundle_end for c in claims], ['WHOLE', 'WHOLE'])
                self.assertEqual(before, context.network.to_dict())

    def test_missing_scope_selection_and_unowned_source_attachments_reject(self):
        context, bundle, layouts = fixture()
        members = {Key('intersection', i) for i in (10, 11, 13)}
        with self.assertRaisesRegex(ValueError, 'available original physical'):
            missing_terminals(context, bundle, {1}, members, {k:v for k,v in layouts.items() if k.id != 12})
        for attribute in ('scope_ends', 'selection_ends'):
            network = replace(context.network, sections=tuple(replace(s, **{attribute:(0,)}) if s.id == 2 else s
                                                             for s in context.network.sections),
                              intersections=tuple(j for j in context.network.intersections if j.id != 12))
            with self.assertRaisesRegex(ValueError, 'available original physical|scope or selection'):
                missing_terminals(SimpleNamespace(network=network), bundle, {1}, members, layouts)
        self.assertEqual(missing_terminals(context, bundle, {1}, {Key('intersection', 10)}, layouts), ())
        sections=tuple(replace(s,points=((10.,0.,0.),(10.,5.,0.))) if s.id==2 else s
                       for s in context.network.sections)+(Section(9,9,((10.,5.,0.),(10.,10.,0.)),width=3.5),
                        Section(10,10,((10.,5.,0.),(30.,5.,0.)),width=3.5))
        sites=tuple(replace(j,ends=(End(9,1),End(4,1),End(8,0))) if j.id==13 else j
                    for j in context.network.intersections)+(Intersection(14,(End(2,1),End(9,0),End(10,0))),)
        branched=replace(bundle,members=(bundle.members[0],BundleMember((Traversal(2,True),Traversal(9,True)))))
        network=replace(context.network,sections=sections,intersections=sites,bundles=(branched,))
        with self.assertRaisesRegex(ValueError,'skip an unowned original attachment'):
            missing_terminals(SimpleNamespace(network=network),branched,{1},members,{j.key:j for j in sites})


if __name__ == '__main__':
    unittest.main()
