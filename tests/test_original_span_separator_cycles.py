"""Internal regions retain original incidence; topology is not surface coverage."""
from copy import deepcopy
from dataclasses import asdict, replace
import importlib
from types import SimpleNamespace as S
import unittest

from test_original_span_separators import package

m = importlib.import_module(package.__name__+'.original_span_separator_cycles')
p = importlib.import_module(package.__name__+'.physical_ports')
c = importlib.import_module(package.__name__+'.port_clusters')


def fixture():
    sites = (40, 60, 50, 70)
    streets, bundles, spans, rings = {}, {}, [], {jid: [] for jid in sites}
    for i, sid in enumerate((101, 102, 103, 104)):
        section, bid = sid+200, 900+i % 2
        jids = (sites[i], sites[(i+1) % 4])
        if i % 2:
            jids = jids[::-1]
        ends = tuple(p.SpanEnd(p.SectionEnd(sid, section, at), 'PHYSICAL', jid)
                     for at, jid in enumerate(jids))
        span = p.PhysicalSpan(sid, (section,), *ends, 0, 1, 2)
        spans.append(span)
        for end in ends:
            rings[end.junction_id].append(p.PhysicalPort(end.junction_id, end.endpoint, 'PHYSICAL', sid))
        streets[sid] = dict(id=sid, status='SUPPORTED', bundle_id=bid,
            endpoint_space='PHYSICAL_JUNCTION', endpoints=jids,
            items=[dict(kind='Section', id=section, points=[(0, 0, 0), (1, 1, 0)],
                        category='primary', offset=0., width=8.)])
        bundle = bundles.setdefault(bid, dict(id=bid, status='SUPPORTED', owned_streets=[],
            topology=dict(bundle_id=bid, status='MAPPED', members=[], interior_streets=[],
                          junctions=[], required_features=[])))
        bundle['owned_streets'].append(sid)
        members = bundle['topology']['members']; order = len(members)
        forward = bool(i % 2); a, b = (0., 1.) if forward else (1., 0.)
        members.append(dict(street_id=sid, forward=forward, head_order=order,
            tail_order=order, complete=True, spans=[dict(ordinal=0,
                start=dict(bundle_end='HEAD', junction_id=None),
                end=dict(bundle_end='TAIL', junction_id=None), sections=[dict(
                    start=dict(section_id=section, segment_index=0, factor=a),
                    end=dict(section_id=section, segment_index=0, factor=b))])]))
    cluster = c.LogicalPortCluster(800, sites, tuple(spans), (), (), (), 1)
    component = S(ports=S(source_streets=tuple(streets), spans=tuple(spans),
        sites=tuple(p.PhysicalPortSite(jid, (0., 0.), tuple(ports)) for jid, ports in rings.items())),
        network=S(logical_clusters=(cluster,), ordinary_site_ids=(), compounds=(), routes=()))
    return S(streets=streets, bundles=bundles), component


class InternalCyclePreparationTests(unittest.TestCase):
    def setUp(self):
        self.snapshot, self.component = fixture()
        self.planner = m.InternalCycleRegionPlanner(('primary',))

    def plan(self):
        return self.planner.prepare(self.snapshot, self.component)

    def test_original_identity_only_and_no_input_mutation(self):
        before = deepcopy((self.snapshot, self.component))
        result = self.plan()
        self.assertFalse(result.pending)
        region, = result.regions
        self.assertEqual(region.bundle_ids, (900, 901))
        self.assertEqual(region.physical_junction_ids, (40, 50, 60, 70))
        self.assertEqual([m.span_id[0] for m in region.members], [101, 102, 103, 104])
        self.assertEqual([m.from_end.end for m in region.members], [0, 1, 0, 1])
        self.assertEqual((self.snapshot, self.component), before)
        self.assertNotIn('points', repr(asdict(region)))
        self.assertNotIn('width', repr(asdict(region)))
        self.assertFalse(hasattr(region, 'geometry_accepted'))

    def test_table_order_positions_widths_and_names_do_not_define_identity(self):
        expected = self.plan()
        self.component.ports.spans = self.component.ports.spans[::-1]
        self.component.ports.sites = self.component.ports.sites[::-1]
        self.component.network.logical_clusters = tuple(replace(row,
            internal_spans=row.internal_spans[::-1], physical_junction_ids=row.physical_junction_ids[::-1])
            for row in self.component.network.logical_clusters)
        for source in self.snapshot.streets.values():
            source['items'][0].update(points=[(4., 7., 0.), (-30., 20., 0.)], width=15., name='edited')
        for bundle in self.snapshot.bundles.values():
            bundle['topology']['members'].reverse()
            # Live BundleTopology.report() retains tuples; JSON uses lists.
            # Both are the same contract and must yield the same declaration.
            for member in bundle['topology']['members']:
                for extent in member['spans']:
                    extent['sections'] = tuple(extent['sections'])
                member['spans'] = tuple(member['spans'])
            bundle['topology']['members'] = tuple(bundle['topology']['members'])
        self.assertEqual(self.plan(), expected)

    def test_stale_source_and_bundle_claims_are_explicitly_pending(self):
        for case in ('unsupported', 'inactive', 'wrong_end', 'offset', 'curved', 'category',
                     'partial_section', 'partial_bundle', 'direction', 'competing', 'missing_member'):
            with self.subTest(case=case):
                self.setUp(); source = self.snapshot.streets[101]
                member = self.snapshot.bundles[900]['topology']['members'][0]
                if case == 'unsupported': source['status'] = 'UNSUPPORTED'
                if case == 'inactive': self.component.ports.source_streets = (102, 103, 104)
                if case == 'wrong_end': source['endpoints'] = source['endpoints'][::-1]
                if case == 'offset': source['items'][0]['offset'] = 0.1
                if case == 'curved': source['items'][0]['points'].append((2, 2, 0))
                if case == 'category': source['items'][0]['category'] = 'footway'
                if case == 'partial_section': member['spans'][0]['sections'][0]['start']['factor'] = .1
                if case == 'partial_bundle': member['complete'] = False
                if case == 'direction': member['forward'] = not member['forward']
                if case == 'competing': self.snapshot.bundles[999] = dict(owned_streets=[101])
                if case == 'missing_member': self.snapshot.bundles[900]['owned_streets'].pop()
                result = self.plan()
                self.assertFalse(result.regions)
                self.assertEqual(len(result.pending), 1)
                self.assertEqual(result.pending[0]['logical_junction_ids'], (800,))

    def test_incidence_and_original_span_changes_are_not_hidden(self):
        for case in ('missing_port', 'duplicate_port', 'wrong_port_site', 'missing_span', 'missing_site',
                     'repeated_span', 'multiple_cycles', 'extra_footprint'):
            with self.subTest(case=case):
                self.setUp(); plan = self.component.ports; cluster, = self.component.network.logical_clusters
                site, *rest = plan.sites
                if case == 'missing_port': site = replace(site, ports=site.ports[1:])
                if case == 'duplicate_port': site = replace(site, ports=site.ports+site.ports[:1])
                if case == 'wrong_port_site': site = replace(site, ports=(replace(site.ports[0], junction_id=999),)+site.ports[1:])
                plan.sites = (site, *rest)
                if case == 'missing_span': plan.spans = plan.spans[1:]
                if case == 'missing_site': plan.sites = plan.sites[1:]
                if case == 'repeated_span': cluster = replace(cluster, internal_spans=cluster.internal_spans+cluster.internal_spans[:1])
                if case == 'multiple_cycles': cluster = replace(cluster, source_cycle_rank=2)
                if case == 'extra_footprint': cluster = replace(cluster, footprint_spans=cluster.internal_spans[:1])
                self.component.network.logical_clusters = (cluster,)
                result = self.plan()
                self.assertFalse(result.regions); self.assertTrue(result.pending)

    def test_whole_shared_logical_owner_retained_and_acyclic_clusters_skipped(self):
        old, = self.component.network.logical_clusters
        shared = c.SharedLogicalPortCluster(**old.__dict__, logical_junction_ids=(800, 801))
        self.component.network.logical_clusters = (shared,)
        self.assertEqual(self.plan().regions[0].logical_junction_ids, (800, 801))
        self.component.network.logical_clusters = (replace(old, source_cycle_rank=0),)
        self.assertEqual(self.plan(), m.InternalCyclePreparation((), ()))
        self.component.network.logical_clusters = (old,)
        self.component.network.ordinary_site_ids = (40,)
        with self.assertRaises(ValueError): self.plan()

    def ordinary_cycle(self):
        self.snapshot.bundles.clear()
        for source in self.snapshot.streets.values():
            source['bundle_id'] = None
        cluster, = self.component.network.logical_clusters
        self.component.network.logical_clusters = (replace(cluster,
            footprint_spans=((101, 301, 301), (103, 303, 303))),)

    def test_ordinary_intersection_cycle_does_not_infer_a_bundle_island(self):
        self.ordinary_cycle()
        original = deepcopy((self.snapshot, self.component))
        result = self.plan()
        self.assertFalse(result.regions)
        self.assertFalse(result.pending)
        row, = result.not_applicable
        self.assertEqual(row['code'], 'NO_INTERNAL_BUNDLE_MEMBERS')
        self.assertEqual(row['street_ids'], (101, 102, 103, 104))
        self.assertEqual(row['physical_junction_ids'], (40, 50, 60, 70))
        self.assertEqual((self.snapshot, self.component), original)
        self.assertNotIn('geometry_accepted', row)
        self.assertNotIn('bundle_ids', row)
        self.component.network.logical_clusters = tuple(replace(c,
            internal_spans=c.internal_spans[::-1]) for c in self.component.network.logical_clusters)
        for source in self.snapshot.streets.values():
            source['items'][0].update(points=[(30, 40, 0), (-3, 6, 0)], width=13)
        self.assertEqual(self.plan(), result)

    def test_missing_reverse_bundle_membership_cannot_hide_a_pending_region(self):
        for role in ('owned_streets', 'head', 'tail', 'members', 'topology_members',
                     'chains', 'interior_streets', 'street_bundle_id'):
            with self.subTest(role=role):
                self.setUp(); self.ordinary_cycle()
                bundle = dict(id=900, topology={})
                self.snapshot.bundles[900] = bundle
                if role in ('owned_streets', 'head', 'tail'): bundle[role] = [101]
                elif role == 'members': bundle['members'] = [dict(street_id=101)]
                elif role == 'topology_members': bundle['topology']['members'] = [dict(street_id=101)]
                elif role == 'chains': bundle['topology']['chains'] = [[101]]
                elif role == 'interior_streets': bundle['topology']['interior_streets'] = [101]
                else: self.snapshot.streets[101]['bundle_id'] = 900
                result = self.plan()
                self.assertFalse(result.not_applicable)
                self.assertFalse(result.regions)
                self.assertEqual(len(result.pending), 1)

    def test_unknown_cycle_sources_remain_pending(self):
        for case in ('missing', 'unsupported', 'inactive', 'duplicate', 'wrong_id'):
            with self.subTest(case=case):
                self.setUp(); self.ordinary_cycle()
                if case == 'missing': del self.snapshot.streets[103]
                elif case == 'unsupported': self.snapshot.streets[103]['status'] = 'UNSUPPORTED'
                elif case == 'inactive': self.component.ports.source_streets = (101, 102, 104)
                elif case == 'wrong_id': self.snapshot.streets[103]['id'] = 999
                else:
                    cluster, = self.component.network.logical_clusters
                    self.component.network.logical_clusters = (replace(cluster,
                        internal_spans=cluster.internal_spans+cluster.internal_spans[:1]),)
                result = self.plan()
                self.assertFalse(result.not_applicable)
                self.assertFalse(result.regions)
                self.assertEqual(len(result.pending), 1)


if __name__ == '__main__':
    unittest.main()
