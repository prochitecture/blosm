"""A partial paired head keeps the real opposite terminal and its incidence."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import _python_street_package
from test_python_street_bundle_terminal_expansion import fixture as base_fixture
from tested_python_streets.model import Section,End,Intersection,BundleMember,Traversal,Key
from tested_python_streets.bundle_partial_terminal_expansion import missing_terminal
from tested_python_streets.compound_preparation import DeferredCompoundPlanner,CompoundExpansionRequired


def fixture(reverse=False,tail=False):
    context,bundle,_=base_fixture(reverse)
    points=((0.,10.,0.),(0.,20.,0.))
    added=Section(9,9,points[::-1] if reverse else points,width=3.5)
    sites=tuple(replace(j,ends=(*j.ends,End(9,int(reverse)))) if j.id==11 else j
                for j in context.network.intersections)
    members=(BundleMember((*bundle.members[0].path,Traversal(9,not reverse))),bundle.members[1])
    if tail:
        members=tuple(BundleMember(tuple(Traversal(v.section_id,not v.forward) for v in m.path[::-1])) for m in members)
    bundle=replace(bundle,members=members)
    n=replace(context.network,sections=(*context.network.sections,added),intersections=sites,bundles=(bundle,))
    return SimpleNamespace(network=n),bundle,{j.key:j for j in sites}


class PartialTerminalExpansionTests(unittest.TestCase):
    def test_head_tail_and_reversed_storage_keep_exact_original_neighbor(self):
        for reverse in (False,True):
            for tail in (False,True):
                with self.subTest(reverse=reverse,tail=tail):
                    c,bundle,sites=fixture(reverse,tail);before=c.network.to_dict()
                    members={Key('intersection',10),Key('intersection',11)}
                    self.assertEqual(missing_terminal(c,bundle,{1},members,sites),(Key('intersection',12),))
                    planner=DeferredCompoundPlanner(c,c.network.intersections,sites)
                    with self.assertRaises(CompoundExpansionRequired) as caught:
                        planner.check_bundle_terminals({1},(),members)
                    self.assertEqual(caught.exception.neighbors,(Key('intersection',12),))
                    self.assertEqual(missing_terminal(c,bundle,{1},members|{Key('intersection',12)},sites),())
                    self.assertEqual(c.network.to_dict(),before)

    def test_missing_owner_and_coincident_unconnected_end_do_not_expand(self):
        c,bundle,sites=fixture();members={Key('intersection',10),Key('intersection',11)}
        with self.assertRaisesRegex(ValueError,'available original physical'):
            missing_terminal(c,bundle,{1},members,{k:v for k,v in sites.items() if k.id!=12})
        added=Section(10,10,((10.,0.,0.),(20.,-20.,0.)),width=3.5)
        disconnected=tuple(replace(j,ends=(End(2,0),End(7,0),End(10,0))) if j.id==12 else j
                           for j in c.network.intersections)
        n=replace(c.network,sections=(*c.network.sections,added),intersections=disconnected)
        with self.assertRaisesRegex(ValueError,'original neighboring source route'):
            missing_terminal(SimpleNamespace(network=n),bundle,{1},members,{j.key:j for j in disconnected})

    def test_short_real_route_closure_acquires_every_intermediate_owner_and_rejects_distant_routes(self):
        c,bundle,sites=fixture();members={Key('intersection',10),Key('intersection',11)}
        source=c.network.section_index[3]
        sections=tuple(replace(s,points=((0.,0.,0.),(5.,0.,0.))) if s.id==3 else s
                       for s in c.network.sections)
        sections+= (Section(10,source.street_id,((5.,0.,0.),(10.,0.,0.)),width=source.width),
                    Section(11,11,((5.,0.,0.),(5.,-20.,0.)),width=source.width))
        changed=tuple(replace(j,ends=(End(2,0),End(10,1),End(7,0))) if j.id==12 else j
                      for j in c.network.intersections)
        changed+= (Intersection(99,(End(3,1),End(10,0),End(11,0))),)
        network=replace(c.network,sections=sections,intersections=changed)
        before=network.to_dict()
        self.assertEqual(missing_terminal(SimpleNamespace(network=network),bundle,{1},members,
            {j.key:j for j in changed}),(Key('intersection',12),Key('intersection',99)))
        self.assertEqual(network.to_dict(),before)
        distant=replace(network,sections=tuple(replace(s,points=tuple((10*x,y,z) for x,y,z in s.points))
                                              for s in network.sections))
        with self.assertRaisesRegex(ValueError,'original neighboring source route'):
            missing_terminal(SimpleNamespace(network=distant),bundle,{1},members,{j.key:j for j in changed})


if __name__=='__main__':unittest.main()
