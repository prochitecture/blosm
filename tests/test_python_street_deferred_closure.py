"""Finite original ownership closure has a per-pool, not per-network, budget."""
from types import SimpleNamespace, MethodType
from collections import Counter
import unittest
import _python_street_package
from tested_python_streets.model import Section, End, Intersection, Network, Key
from tested_python_streets.compound_preparation import (
    DeferredCompoundPlanner, DeferredCompoundFootprint, CompoundExpansionRequired)


def fixture(count):
    sections=[];sites=[];layouts={};triples=[]
    for i in range(count):
        a,b,c=(Key('intersection',3*i+j) for j in (1,2,3))
        triples.append((a,b,c));s=7*i
        origins=((20.*i,0.,0.),(20.*i,5.,0.),(20.*i,10.,0.))
        sections.extend((Section(s+1,s+1,origins[:2]),Section(s+2,s+2,origins[1:])))
        ends=[[End(s+1,0)],[End(s+1,1),End(s+2,0)],[End(s+2,1)]]
        sid=s+3
        for row,origin in zip(ends,origins):
            while len(row)<3:
                destination=(origin[0]+10.,origin[1]+2.*len(row),0.)
                sections.append(Section(sid,sid,(origin,destination)));row.append(End(sid,0));sid+=1
        local=[Intersection(k.id,tuple(row)) for k,row in zip((a,b,c),ends)]
        sites.extend(local)
        layouts[a]=DeferredCompoundFootprint(local[0],(),(b,),'actual original short-route extent')
        for key in (b,c):layouts[key]=SimpleNamespace(mouths={},absorbed_transitions=())
    network=Network('closure-graph',tuple(sections),tuple(sites))
    return DeferredCompoundPlanner(SimpleNamespace(network=network),tuple(sites),layouts),tuple(triples)


def accepted(members):
    # This test isolates closure, not geometry: real source topology is above,
    # while the native exterior callback stands in for an already certified rim.
    return SimpleNamespace(members=tuple(members),physical_layout=SimpleNamespace(mouths={}),footprints=(),loops=())


class DeferredClosureTests(unittest.TestCase):
    def test_nine_independent_real_expansions_keep_local_limit_and_reuse_unchanged_rims(self):
        planner,triples=fixture(9);calls=Counter()
        target={a:c for a,b,c in triples}
        def exterior(self,members):
            yield
            calls[members]+=1
            if len(members)==2:raise CompoundExpansionRequired((target[members[0]],))
            return accepted(members)
        planner.exterior_steps=MethodType(exterior,planner)
        groups=planner.build(())
        self.assertEqual(tuple(g.members for g in groups),triples)
        self.assertEqual(len(calls),18)
        self.assertEqual(set(calls.values()),{1})
        self.assertEqual(planner.maximum_members,8)

    def test_non_growing_empty_or_unavailable_requests_reject_without_restarting(self):
        for mode in ('same','empty','missing'):
            planner,triples=fixture(1);calls=[]
            def exterior(self,members):
                yield
                calls.append(members)
                request={'same':members,'empty':(),'missing':(Key('intersection',999),)}[mode]
                raise CompoundExpansionRequired(request)
            planner.exterior_steps=MethodType(exterior,planner)
            with self.subTest(mode=mode),self.assertRaisesRegex(ValueError,
                    'new original neighbor' if mode!='missing' else 'available original physical'):
                planner.build(())
            self.assertEqual(len(calls),1)

    def test_existing_two_site_success_and_no_pending_fast_path_keep_order(self):
        planner,triples=fixture(2);calls=[]
        def exterior(self,members):
            yield
            calls.append(members)
            return accepted(members)
        planner.exterior_steps=MethodType(exterior,planner)
        groups=planner.build(())
        self.assertEqual(tuple(g.members for g in groups),tuple(t[:2] for t in triples))
        self.assertEqual(calls,[t[:2] for t in triples])
        for key in planner.layouts:planner.layouts[key]=SimpleNamespace(mouths={})
        old=(accepted((triples[1][2],)),accepted((triples[0][2],)))
        self.assertEqual(planner.build(old),old)

    def test_real_growth_still_rejects_a_physical_pool_over_eight_sites(self):
        planner,triples=fixture(3)
        def exterior(self,members):
            yield
            if len(members)==2:raise CompoundExpansionRequired(tuple(k for t in triples for k in t))
            return accepted(members)
        planner.exterior_steps=MethodType(exterior,planner)
        with self.assertRaisesRegex(ValueError,'local member budget'):planner.build(())


if __name__=='__main__':unittest.main()
