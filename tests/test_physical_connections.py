"""Joining actual endpoint owners must retain source identities and road loops."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

from test_physical_junctions import contract as k,module as physical

spec=importlib.util.spec_from_file_location('tested_physical_connections',
    Path(__file__).resolve().parents[1]/'way/street_network/physical_connections.py')
module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)


def fixture():
    sites={1:(0.,0.),2:(40.,0.),3:(40.,40.),4:(0.,40.)}
    snapshot=NS(streets={},bundles={},junction_origins={'physical':[]})
    sections={};rings={jid:[] for jid in sites}
    for sid,(a,b) in enumerate(((1,2),(2,3),(3,4),(4,1)),1):
        sections[sid]=(k.SectionInput(100+sid,(sites[a],sites[b]),5.+sid,category='residential'),)
        snapshot.streets[sid]=dict(id=sid,status='SUPPORTED',bundle_id=None,endpoints=(a,b))
        for jid,end in ((a,0),(b,1)):
            rings[jid].append(dict(connector_id=10*jid+len(rings[jid]),item_kind='Street',
                endpoints=[dict(street_id=sid,end=end,terminal_kind='Section',terminal_id=100+sid)]))
    for jid,location in sites.items():
        snapshot.junction_origins['physical'].append(dict(junction_id=jid,location=location,ring=rings[jid]))
    reader=physical.PreparedPhysicalJunctions(k,lambda snapshot,sid:sections[sid],{'residential':{'radius':3.}})
    return snapshot,reader


class PhysicalConnectionTests(unittest.TestCase):
    def test_two_distinct_links_close_a_hole_without_absorbing_either_street(self):
        snapshot,reader=fixture();before=deepcopy(snapshot)
        seed=reader.prepare_sites(snapshot,(1,2,3)).junctions
        result=module.PreparedPhysicalConnections(reader).extend(snapshot,seed,(3,4))
        self.assertEqual(result.added_sites,(4,));self.assertEqual(result.cycle_rank,1)
        self.assertEqual(result.internal_links,((1,1,2),(2,2,3),(3,3,4),(4,4,1)))
        self.assertEqual(result.junctions[:3],seed);self.assertEqual(snapshot,before)

    def test_following_one_link_adds_the_complete_remote_ring(self):
        snapshot,reader=fixture();seed=reader.prepare_sites(snapshot,(1,2)).junctions
        result=module.PreparedPhysicalConnections(reader).extend(snapshot,seed,(2,))
        self.assertEqual(result.added_sites,(3,));self.assertEqual(result.cycle_rank,0)
        self.assertEqual([s.source_id for s in result.junctions[-1].streets],[2,3])

    def test_absent_and_duplicate_original_far_owners_reject(self):
        for change in ('absent','duplicate'):
            snapshot,reader=fixture();seed=reader.prepare_sites(snapshot,(1,2,3)).junctions
            sites=snapshot.junction_origins['physical']
            if change=='absent':sites.pop()
            else:sites.append(deepcopy(sites[-1]))
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'one exact'):
                module.PreparedPhysicalConnections(reader).extend(snapshot,seed,(3,))

    def test_existing_sources_and_profiles_cannot_be_silently_changed(self):
        snapshot,reader=fixture();seed=list(reader.prepare_sites(snapshot,(1,2,3)).junctions)
        j=seed[0];seed[0]=replace(j,profiles=tuple(replace(p,rise=.1) for p in j.profiles))
        with self.assertRaisesRegex(ValueError,'profile/height transition'):
            module.PreparedPhysicalConnections(reader).extend(snapshot,seed,(3,))

    def test_invalid_selection_and_missing_site_reject(self):
        snapshot,reader=fixture();seed=reader.prepare_sites(snapshot,(1,2)).junctions
        for ids in ((),(2,2),(1,),(99,)):
            with self.subTest(ids=ids),self.assertRaises(ValueError):
                module.PreparedPhysicalConnections(reader).extend(snapshot,seed,ids)
        for ids in ((),(1,1),(99,)):
            with self.subTest(ids=ids),self.assertRaises(ValueError):reader.prepare_sites(snapshot,ids)


if __name__=='__main__':unittest.main()
